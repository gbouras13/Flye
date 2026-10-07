//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <zlib.h>
#include <sys/stat.h>

#include "polish.h"
#include "aligner.h"
#include "samalign.h"
#include "bubbles.h"
#include "seqio.h"
#include "py_compat.h"
#include "../polishing/bubble_processor.h"
#include "../common/logger.h"

namespace polish
{

namespace
{
	const int MAX_READ_COVERAGE = 1000;
	const int HARD_MINIMUM_COVERAGE = 3;
	const int RELATIVE_MINIMUM_COVERAGE = 5;

	std::string joinPath(const std::string& dir, const std::string& name)
	{
		if (dir.empty() || dir.back() == '/') return dir + name;
		return dir + "/" + name;
	}

	size_t fileSize(const std::string& path)
	{
		struct stat st;
		if (stat(path.c_str(), &st) != 0) return 0;
		return (size_t)st.st_size;
	}

	std::string formatFloat6(double value)
	{
		char buf[64];
		snprintf(buf, sizeof(buf), "%f", value);	//Python's "%f"
		return buf;
	}

	void writeGzipText(const std::string& path, const std::string& text)
	{
		gzFile out = gzopen(path.c_str(), "wb");
		if (!out) throw std::runtime_error("Can't open " + path);
		if (!text.empty() && gzwrite(out, text.data(), (unsigned)text.size()) == 0)
		{
			gzclose(out);
			throw std::runtime_error("Error writing " + path);
		}
		gzclose(out);
	}

	struct ComposedBubble
	{
		int64_t pos;
		int64_t subPos;
		int64_t coverage;
		std::string seq;
	};

	//polish._compose_sequence
	void composeSequence(const std::string& consensusFile,
						 pycompat::OrderedDict<std::string>& polishedFasta,
						 pycompat::OrderedDict<int64_t>& polishedLengths,
						 pycompat::OrderedDict<std::vector<std::pair<int64_t, int64_t>>>& coverages)
	{
		pycompat::OrderedDict<std::vector<ComposedBubble>> consensuses;
		std::ifstream in(consensusFile);
		if (!in.is_open()) throw std::runtime_error("Can't open " + consensusFile);
		std::string line;
		bool header = true;
		std::string ctgId;
		int64_t ctgPos = 0, coverage = 0, ctgSubPos = 0;
		while (std::getline(in, line))
		{
			if (header)
			{
				std::vector<std::string> tokens = pycompat::pySplitSep(pycompat::pyStrip(line), ' ');
				if (tokens.size() != 4) throw std::runtime_error("Bubble format error");
				ctgId = tokens[0].substr(1);
				ctgPos = std::stoll(tokens[1]);
				coverage = std::stoll(tokens[2]);
				ctgSubPos = std::stoll(tokens[3]);
			}
			else
			{
				consensuses[ctgId].push_back({ctgPos, ctgSubPos, coverage, pycompat::pyStrip(line)});
			}
			header = !header;
		}

		for (auto& kv : consensuses)
		{
			std::vector<ComposedBubble>& seqs = kv.second;
			std::stable_sort(seqs.begin(), seqs.end(),
				[](const ComposedBubble& a, const ComposedBubble& b)
				{
					if (a.pos != b.pos) return a.pos < b.pos;
					return a.subPos < b.subPos;
				});
			std::string concat;
			std::vector<std::pair<int64_t, int64_t>> bubbleCoverages;
			for (const ComposedBubble& b : seqs)
			{
				concat += b.seq;
				bubbleCoverages.emplace_back((int64_t)b.seq.size(), b.coverage);
			}
			polishedLengths[kv.first] = (int64_t)concat.size();
			polishedFasta[kv.first] = std::move(concat);
			coverages[kv.first] = std::move(bubbleCoverages);
		}
	}

	//order of contigs as they first appear in the bubbles file (the order of the
	//polisher output with a single thread)
	std::vector<std::string> bubbleContigOrder(const std::string& bubblesFile)
	{
		std::vector<std::string> order;
		std::unordered_map<std::string, bool> seen;
		std::ifstream in(bubblesFile);
		std::string line;
		while (std::getline(in, line))
		{
			if (line.empty() || line[0] != '>') continue;
			size_t space = line.find(' ');
			if (space == std::string::npos) continue;	//branch header
			std::string ctg = line.substr(1, space - 1);
			if (!seen.count(ctg)) {seen[ctg] = true; order.push_back(ctg);}
		}
		return order;
	}
}

std::string polish(const std::string& contigs, const std::vector<std::string>& reads,
				   const std::string& workDir, int numIters,
				   const PolishSettings& settings, std::string& statsFile)
{
	const bool nano = settings.platform == "nano";
	const std::string subsMatrix = joinPath(settings.configDir, nano ?
		"nano_r94_substitutions.mat" : "pacbio_chm13_substitutions.mat");
	const std::string hopoMatrix = joinPath(settings.configDir, nano ?
		"nano_r94_g36_homopolymers.mat" : "pacbio_chm13_homopolymers.mat");
	const bool useHopo = false;	//hopo_enabled is False for both platforms
	statsFile = joinPath(workDir, "contigs_stats.txt");
	const std::string bedCoverage = joinPath(workDir, "base_coverage.bed.gz");

	std::string prevAssembly = contigs;
	pycompat::OrderedDict<int64_t> contigLengths;
	pycompat::OrderedDict<int64_t> coverageStats;
	pycompat::OrderedDict<std::vector<std::pair<int64_t, int64_t>>> bubbleCoverages;
	for (int i = 0; i < numIters; ++i)
	{
		Logger::get().info() << "Polishing genome (" << i + 1 << "/" << numIters << ")";
		Logger::get().info() << "Running minimap2";
		aligner::AlignmentSet alignment =
			aligner::makeAlignment(prevAssembly, reads, settings.numThreads,
								   settings.platform, settings.readType);

		Logger::get().info() << "Separating alignment into bubbles";
		samalign::SeqDict contigsFasta = seqio::readSequenceDict(prevAssembly);
		const std::string bubblesFile = joinPath(workDir, "bubbles_" + std::to_string(i + 1) + ".fasta");
		double meanAlnError = 0;
		coverageStats = bubbles::makeBubbles(alignment, contigsFasta, settings.platform,
											 settings.numThreads, bubblesFile, meanAlnError);
		alignment = aligner::AlignmentSet();	//free memory
		Logger::get().info() << "Alignment error rate: " << formatFloat6(meanAlnError);

		const std::string consensusOut = joinPath(workDir, "consensus_" + std::to_string(i + 1) + ".fasta");
		const std::string polishedFile = joinPath(workDir, "polished_" + std::to_string(i + 1) + ".fasta");
		if (fileSize(bubblesFile) == 0)
		{
			Logger::get().info() << "No reads were aligned during polishing";
			std::ofstream(statsFile) << "#seq_name\tlength\tcoverage\n";
			std::ofstream(polishedFile).close();
			writeGzipText(bedCoverage, "");
			std::remove(bubblesFile.c_str());
			return polishedFile;
		}

		Logger::get().info() << "Correcting bubbles";
		{
			BubbleProcessor bp(subsMatrix, hopoMatrix, /*show progress*/ true, useHopo);
			bp.polishAll(bubblesFile, consensusOut, settings.numThreads);
		}

		pycompat::OrderedDict<std::string> polishedFasta;
		pycompat::OrderedDict<int64_t> polishedLengths;
		pycompat::OrderedDict<std::vector<std::pair<int64_t, int64_t>>> composedCoverages;
		composeSequence(consensusOut, polishedFasta, polishedLengths, composedCoverages);
		seqio::writeFastaDict(polishedFasta, polishedFile);

		//the polisher writes contigs in the order threads finish them; use the
		//order of the bubbles file instead (= Flye's order with one thread)
		contigLengths = pycompat::OrderedDict<int64_t>();
		bubbleCoverages = pycompat::OrderedDict<std::vector<std::pair<int64_t, int64_t>>>();
		for (const std::string& ctg : bubbleContigOrder(bubblesFile))
		{
			if (!polishedLengths.contains(ctg)) continue;
			contigLengths[ctg] = polishedLengths.at(ctg);
			bubbleCoverages[ctg] = composedCoverages.at(ctg);
		}

		std::remove(bubblesFile.c_str());
		std::remove(consensusOut.c_str());
		prevAssembly = polishedFile;
	}

	{
		std::ofstream out(statsFile);
		out << "#seq_name\tlength\tcoverage\n";
		for (auto& kv : contigLengths)
		{
			if (!coverageStats.contains(kv.first)) coverageStats[kv.first] = 0;
			out << kv.first << "\t" << kv.second << "\t" << coverageStats.at(kv.first) << "\n";
		}
	}

	std::ostringstream bed;
	bed << "#seq_name\tstart\tend\tcoverage\n";
	for (auto& kv : bubbleCoverages)
	{
		int64_t curPos = 0;
		for (auto& bub : kv.second)
		{
			bed << kv.first << "\t" << curPos << "\t" << curPos + bub.first << "\t"
				<< bub.second << "\n";
			curPos += bub.first;
		}
	}
	writeGzipText(bedCoverage, bed.str());
	return prevAssembly;
}

void filterByCoverage(bool meta, const std::string& readType,
					  const std::string& statsIn, const std::string& contigsIn,
					  const std::string& statsOut, const std::string& contigsOut)
{
	const int SUBASM_MIN_COVERAGE = 1;

	struct CtgStat {int64_t length; int64_t coverage;};
	std::unordered_map<std::string, CtgStat> ctgStats;
	int64_t sumCov = 0;
	int64_t sumLength = 0;
	{
		std::ifstream in(statsIn);
		std::string line;
		while (std::getline(in, line))
		{
			if (!line.empty() && line[0] == '#') continue;
			if (line.empty()) continue;
			std::vector<std::string> tokens = pycompat::pySplitSep(line, '\t');
			CtgStat stat {std::stoll(tokens.at(1)), std::stoll(tokens.at(2))};
			ctgStats[tokens[0]] = stat;
			sumCov += stat.coverage * stat.length;
			sumLength += stat.length;
		}
	}

	if (sumLength == 0) throw std::runtime_error("float division by zero");	//as in Flye
	int64_t meanCoverage = (int64_t)((double)sumCov / (double)sumLength);
	int64_t coverageThreshold = 0;
	if (readType == "subasm")
	{
		coverageThreshold = SUBASM_MIN_COVERAGE;
	}
	else if (meta)
	{
		coverageThreshold = HARD_MINIMUM_COVERAGE;
	}
	else
	{
		coverageThreshold = (int64_t)pycompat::pyRound((double)meanCoverage /
													   RELATIVE_MINIMUM_COVERAGE);
		coverageThreshold = std::max<int64_t>(HARD_MINIMUM_COVERAGE, coverageThreshold);
	}
	Logger::get().debug() << "Mean contig coverage: " << meanCoverage
		<< ", selected threshold: " << coverageThreshold;

	int64_t filteredNum = 0;
	int64_t filteredSeq = 0;
	pycompat::OrderedDict<std::string> goodFasta;
	seqio::streamSequence(contigsIn, [&](std::string& hdr, std::string& seq)
	{
		const CtgStat& stat = ctgStats.at(hdr);
		if (stat.coverage >= coverageThreshold)
		{
			goodFasta[hdr] = seq;
		}
		else
		{
			++filteredNum;
			filteredSeq += stat.length;
		}
	});
	Logger::get().debug() << "Filtered " << filteredNum << " contigs of total length "
		<< filteredSeq;

	seqio::writeFastaDict(goodFasta, contigsOut);
	std::ofstream out(statsOut);
	out << "#seq_name\tlength\tcoverage\n";
	for (auto& kv : goodFasta)
	{
		const CtgStat& stat = ctgStats.at(kv.first);
		out << kv.first << "\t" << stat.length << "\t" << stat.coverage << "\n";
	}
}

void generatePolishedEdges(const std::string& edgesFile, const std::string& gfaFile,
						   const std::string& polishedContigs, const std::string& workDir,
						   const PolishSettings& settings, const std::string& polishedStats)
{
	Logger::get().debug() << "Generating polished GFA";

	std::unordered_map<std::string, int64_t> edgesNewCoverage;
	{
		std::ifstream in(polishedStats);
		std::string line;
		while (std::getline(in, line))
		{
			if (!line.empty() && line[0] == '#') continue;
			std::vector<std::string> tokens = pycompat::pySplit(line);
			if (tokens.size() != 3) continue;
			std::string ctgId = pycompat::pySplitSep(tokens[0], '_').at(1);
			edgesNewCoverage[ctgId] = std::stoll(tokens[2]);
		}
	}

	samalign::SeqDict polishedDict = seqio::readSequenceDict(polishedContigs);
	aligner::AlignmentSet alignment =
		aligner::makeAlignment(polishedContigs, {edgesFile}, settings.numThreads,
							   settings.platform, settings.readType);
	samalign::SamReader reader(alignment, polishedDict, MAX_READ_COVERAGE,
							   /*use secondary*/ false);
	pycompat::OrderedDict<std::vector<samalign::Alignment>> alnByEdge;
	for (samalign::Alignment& aln : reader.getAllAlignments())
	{
		alnByEdge[aln.qryId].push_back(std::move(aln));
	}

	const double MIN_CONTAINMENT = 0.9;
	int64_t updatedSeqs = 0;
	samalign::SeqDict edgesDict = seqio::readSequenceDict(edgesFile);
	for (auto& edgeKv : edgesDict)
	{
		const std::string& edge = edgeKv.first;
		if (!alnByEdge.contains(edge)) continue;
		std::vector<samalign::Alignment>& alns = alnByEdge.at(edge);
		//sort(key=qry span, reverse=True): stable descending
		std::stable_sort(alns.begin(), alns.end(),
			[](const samalign::Alignment& a, const samalign::Alignment& b)
			{return (a.qryEnd - a.qryStart) > (b.qryEnd - b.qryStart);});
		const samalign::Alignment& mainAln = alns[0];
		int64_t mapStart = mainAln.trgStart;
		int64_t mapEnd = mainAln.trgEnd;
		for (const samalign::Alignment& aln : alns)
		{
			if (aln.trgId == mainAln.trgId && aln.trgSign == mainAln.trgSign)
			{
				mapStart = std::min(mapStart, aln.trgStart);
				mapEnd = std::max(mapEnd, aln.trgEnd);
			}
		}

		const std::string& trgSeq = polishedDict.at(mainAln.trgId);
		int64_t from = std::min<int64_t>(mapStart, (int64_t)trgSeq.size());
		int64_t to = std::min<int64_t>(mapEnd, (int64_t)trgSeq.size());
		std::string newSeq = to > from ? trgSeq.substr(from, to - from) : "";
		if (mainAln.qrySign == '-') newSeq = seqio::reverseComplement(newSeq);

		//(Flye compares with the length of the last alignment in the loop)
		if ((double)newSeq.size() / (double)alns.back().qryLen > MIN_CONTAINMENT)
		{
			edgeKv.second = newSeq;
			++updatedSeqs;
		}
	}

	{
		std::ifstream gfaIn(gfaFile);
		std::ofstream gfaOut(joinPath(workDir, "polished_edges.gfa"));
		std::string line;
		while (std::getline(gfaIn, line))
		{
			if (!line.empty() && line[0] == 'S')
			{
				std::vector<std::string> tokens = pycompat::pySplit(line);
				std::string seqId = tokens.at(1);
				std::string coverageTag = tokens.at(3);
				std::string seqNum = pycompat::pySplitSep(seqId, '_').at(1);
				auto covIt = edgesNewCoverage.find(seqNum);
				if (covIt != edgesNewCoverage.end())
				{
					coverageTag = "dp:i:" + std::to_string(covIt->second);
				}
				gfaOut << "S\t" << seqId << "\t" << edgesDict.at(seqId) << "\t"
					   << coverageTag << "\n";
			}
			else
			{
				gfaOut << line << "\n";
			}
		}
	}
	Logger::get().debug() << (int64_t)edgesDict.size() - updatedSeqs
		<< " sequences remained unpolished";
}

}	//namespace polish
