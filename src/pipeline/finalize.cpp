//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "finalize.h"
#include "seqio.h"
#include "../common/logger.h"

namespace finalize
{

namespace
{
	const int SCAFFOLD_GAP = 100;

	std::string rc(const std::string& sign) {return sign == "-" ? "+" : "-";}
	std::string unsignedName(const std::string& ctg) {return ctg.substr(1);}

	//Python's str.strip(chars)
	std::string stripChars(const std::string& s, const std::string& chars)
	{
		size_t b = 0, e = s.size();
		while (b < e && chars.find(s[b]) != std::string::npos) ++b;
		while (e > b && chars.find(s[e - 1]) != std::string::npos) --e;
		return s.substr(b, e - b);
	}

	struct SeqStats
	{
		std::string name;
		std::string length;
		std::string coverage;
		std::string circular = "N";
		std::string repeat = "N";
		std::string mult = "1";
		std::string telomere = "none";
		std::string altGroup = "*";
		std::string graphPath;

		void printOut(std::ostream& out) const
		{
			out << name << "\t" << length << "\t" << coverage << "\t" << circular << "\t"
				<< repeat << "\t" << mult << "\t" << altGroup << "\t" << graphPath << "\n";
		}
	};

	int64_t calcN50(std::vector<int64_t> lengths, int64_t assemblyLen)
	{
		std::stable_sort(lengths.begin(), lengths.end(), std::greater<int64_t>());
		int64_t sum = 0;
		for (int64_t l : lengths)
		{
			sum += l;
			if (sum > assemblyLen / 2) return l;
		}
		return 0;
	}

	//Python's int() of a decimal string (graph path elements)
	int64_t pyInt(const std::string& s)
	{
		size_t pos = 0;
		int64_t value = std::stoll(pycompat::pyStrip(s), &pos);
		if (pos != pycompat::pyStrip(s).size()) throw std::runtime_error("invalid literal for int(): " + s);
		return value;
	}
}

Scaffolds generateScaffolds(const std::string& contigsFile, const std::string& linksFile,
							const std::string& outScaffolds)
{
	pycompat::OrderedDict<std::string> contigsFasta = seqio::readSequenceDict(contigsFile);
	std::unordered_set<std::string> usedContigs;
	std::unordered_map<std::string, std::string> connections;

	if (!linksFile.empty())
	{
		std::ifstream in(linksFile);
		std::string line;
		while (std::getline(in, line))
		{
			line = pycompat::pyStrip(line);
			if (line.empty()) continue;
			std::vector<std::string> tokens = pycompat::pySplitSep(line, '\t');
			if (tokens.size() != 4) throw std::runtime_error("Bad scaffold link: " + line);
			const std::string& ctg1 = tokens[0];
			const std::string& sign1 = tokens[1];
			const std::string& ctg2 = tokens[2];
			const std::string& sign2 = tokens[3];
			if (contigsFasta.contains(ctg1) && contigsFasta.contains(ctg2))
			{
				connections[sign1 + ctg1] = sign2 + ctg2;
				connections[rc(sign2) + ctg2] = rc(sign1) + ctg1;
			}
		}
	}

	pycompat::OrderedDict<std::string> scaffoldsFasta;
	Scaffolds scaffoldsSeq;
	for (auto& ctgKv : contigsFasta)
	{
		std::string ctg = ctgKv.first;
		if (usedContigs.count(ctg)) continue;

		usedContigs.insert(ctg);
		std::vector<std::string> scf = {"-" + ctg};
		//extending right
		while (connections.count(scf.back()) &&
			   !usedContigs.count(unsignedName(connections.at(scf.back()))))
		{
			scf.push_back(connections.at(scf.back()));
			usedContigs.insert(unsignedName(scf.back()));
		}

		for (size_t i = 0; i < scf.size(); ++i)
		{
			ctg = scf[i];	//the loop variable is reused afterwards
			scf[i] = rc(std::string(1, ctg[0])) + unsignedName(ctg);
		}
		std::reverse(scf.begin(), scf.end());

		//extending left
		while (connections.count(scf.back()) &&
			   !usedContigs.count(unsignedName(connections.at(scf.back()))))
		{
			scf.push_back(connections.at(scf.back()));
			usedContigs.insert(unsignedName(scf.back()));
		}

		//generating sequence interleaved by Ns
		if (scf.size() == 1)
		{
			scaffoldsFasta[unsignedName(ctg)] = contigsFasta.at(unsignedName(ctg));
			scaffoldsSeq[unsignedName(ctg)] = scf;
		}
		else
		{
			std::string scfName = "scaffold_" + stripChars(unsignedName(scf[0]), "contig_");
			scaffoldsSeq[scfName] = scf;
			std::string gap(SCAFFOLD_GAP, 'N');
			std::string scfSeq;
			for (size_t i = 0; i < scf.size(); ++i)
			{
				if (i > 0) scfSeq += gap;
				const std::string& seq = contigsFasta.at(unsignedName(scf[i]));
				scfSeq += scf[i][0] == '+' ? seq : seqio::reverseComplement(seq);
			}
			scaffoldsFasta[scfName] = std::move(scfSeq);
		}
	}

	seqio::writeFastaDict(scaffoldsFasta, outScaffolds);
	return scaffoldsSeq;
}

void generateStats(const std::string& repeatFile, const std::string& polishedFile,
				   const Scaffolds& scaffolds, const std::string& outStats)
{
	std::unordered_map<std::string, SeqStats> contigsStats;
	const std::string headerLine = "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path";
	{
		std::ifstream in(repeatFile);
		std::string line;
		while (std::getline(in, line))
		{
			if (!line.empty() && line[0] == '#') continue;
			std::vector<std::string> tokens = pycompat::pySplitSep(pycompat::pyStrip(line), '\t');
			if (tokens.size() < 1 || tokens[0].empty()) continue;
			SeqStats stats;
			std::string* fields[] = {&stats.name, &stats.length, &stats.coverage,
									 &stats.circular, &stats.repeat, &stats.mult,
									 &stats.telomere, &stats.altGroup, &stats.graphPath};
			for (size_t i = 0; i < tokens.size() && i < 9; ++i) *fields[i] = tokens[i];
			contigsStats[tokens[0]] = stats;
		}
	}

	if (!polishedFile.empty())
	{
		std::ifstream in(polishedFile);
		std::string line;
		while (std::getline(in, line))
		{
			if (!line.empty() && line[0] == '#') continue;
			std::vector<std::string> tokens = pycompat::pySplitSep(pycompat::pyStrip(line), '\t');
			if (tokens.size() < 3) continue;
			SeqStats& stats = contigsStats.at(tokens[0]);

			//update multiplicity proportionally
			double covRate = std::stod(tokens[2]) / (std::stod(stats.coverage) + 1);
			stats.mult = std::to_string(std::max<int64_t>(1,
								(int64_t)(std::stod(stats.mult) * covRate)));

			//update length and coverage
			stats.length = tokens[1];
			stats.coverage = tokens[2];
		}
	}

	pycompat::OrderedDict<SeqStats> scaffoldsStats;
	for (const auto& kv : scaffolds)
	{
		const std::string& scf = kv.first;
		const std::vector<std::string>& scfSeq = kv.second;
		SeqStats& st = scaffoldsStats[scf];
		st = SeqStats();
		st.name = scf;

		int64_t scfLength = 0;
		std::vector<int64_t> covs;
		int64_t scfMult = INT64_MAX;
		for (const std::string& c : scfSeq)
		{
			const SeqStats& cs = contigsStats.at(unsignedName(c));
			scfLength += pyInt(cs.length);
			covs.push_back(pyInt(cs.coverage));
			scfMult = std::min(scfMult, pyInt(cs.mult));
		}
		scfLength += ((int64_t)scfSeq.size() - 1) * SCAFFOLD_GAP;
		st.length = std::to_string(scfLength);

		int64_t covSum = 0;
		for (int64_t c : covs) covSum += c;
		st.coverage = std::to_string(covs.empty() ? 0 : covSum / (int64_t)covs.size());

		const SeqStats& first = contigsStats.at(unsignedName(scfSeq.front()));
		const SeqStats& last = contigsStats.at(unsignedName(scfSeq.back()));
		st.repeat = first.repeat;
		st.circular = first.circular;
		st.altGroup = first.altGroup;
		st.mult = std::to_string(scfMult);

		//telomere information
		bool scfLeft = false;
		bool scfRight = false;
		if (scfSeq.front()[0] == '+') scfLeft = first.telomere == "left" || first.telomere == "both";
		else scfLeft = first.telomere == "right" || first.telomere == "both";
		if (scfSeq.back()[0] == '+') scfRight = last.telomere == "right" || last.telomere == "both";
		else scfRight = last.telomere == "left" || last.telomere == "both";

		//graph path
		std::vector<std::string> path;
		for (const std::string& ctg : scfSeq)
		{
			std::string ctgPath = contigsStats.at(unsignedName(ctg)).graphPath;
			if (ctg[0] == '-')
			{
				std::vector<std::string> parts = pycompat::pySplitSep(ctgPath, ',');
				std::reverse(parts.begin(), parts.end());
				std::string joined;
				for (size_t i = 0; i < parts.size(); ++i)
				{
					if (i > 0) joined += ",";
					joined += std::to_string(-pyInt(parts[i]));
				}
				ctgPath = joined;
			}
			path.push_back(ctgPath);
		}
		std::string joinedPath;
		for (size_t i = 0; i < path.size(); ++i)
		{
			if (i > 0) joinedPath += ",??,";
			joinedPath += path[i];
		}
		st.graphPath = std::string(scfLeft ? "*," : "") + joinedPath + (scfRight ? ",*" : "");
	}

	std::vector<std::string> order;
	for (auto& kv : scaffoldsStats) order.push_back(kv.first);
	std::stable_sort(order.begin(), order.end(),
		[&scaffoldsStats](const std::string& a, const std::string& b)
		{return pyInt(scaffoldsStats.at(a).length) > pyInt(scaffoldsStats.at(b).length);});
	{
		std::ofstream out(outStats);
		out << headerLine << "\n";
		for (const std::string& scf : order) scaffoldsStats.at(scf).printOut(out);
	}

	int64_t totalLength = 0;
	std::vector<int64_t> scaffoldLengths;
	for (auto& kv : scaffoldsStats)
	{
		totalLength += pyInt(kv.second.length);
		scaffoldLengths.push_back(pyInt(kv.second.length));
	}
	if (totalLength == 0) return;

	int64_t numScaffolds = (int64_t)scaffoldsStats.size();
	int64_t numContigs = 0;
	for (auto& kv : scaffolds) numContigs += (int64_t)kv.second.size();
	int64_t largestScf = *std::max_element(scaffoldLengths.begin(), scaffoldLengths.end());
	int64_t scfN50 = calcN50(scaffoldLengths, totalLength);
	int64_t meanReadCov = 0;
	for (auto& kv : scaffoldsStats)
	{
		meanReadCov += pyInt(kv.second.length) * pyInt(kv.second.coverage);
	}
	meanReadCov /= totalLength;

	Logger::get().info() << "Assembly statistics:\n\n"
		<< "\tTotal length:\t" << totalLength << "\n"
		<< "\tFragments:\t" << numScaffolds << "\n"
		<< "\tFragments N50:\t" << scfN50 << "\n"
		<< "\tLargest frg:\t" << largestScf << "\n"
		<< "\tScaffolds:\t" << numContigs - numScaffolds << "\n"
		<< "\tMean coverage:\t" << meanReadCov << "\n";
}

std::pair<int64_t, int64_t> shortStatistics(const std::string& fastaFile)
{
	std::vector<int64_t> lengths;
	for (auto& kv : seqio::readSequenceLengths(fastaFile)) lengths.push_back((int64_t)kv.second);
	if (lengths.empty()) return {0, 0};
	int64_t total = 0;
	for (int64_t l : lengths) total += l;
	return {total, calcN50(lengths, total)};
}

}	//namespace finalize
