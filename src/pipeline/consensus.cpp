//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <array>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "consensus.h"
#include "parallel_util.h"

namespace consensus
{

namespace
{
	//symbols that can appear in the aligned strings (SAM SEQ after the BAM
	//round trip, the reference, and '-'), in Python's sorted() order
	const std::string SYMBOLS = "-=ABCDGHKMNRSTVWY";

	struct SymbolIndex
	{
		int8_t index[256];
		SymbolIndex()
		{
			for (int i = 0; i < 256; ++i) index[i] = -1;
			for (size_t i = 0; i < SYMBOLS.size(); ++i) index[(uint8_t)SYMBOLS[i]] = (int8_t)i;
		}
	};

	int symbolIndex(char c)
	{
		static SymbolIndex table;
		int idx = table.index[(uint8_t)c];
		if (idx < 0) throw std::runtime_error(std::string("Unexpected symbol in alignment: ") + c);
		return idx;
	}

	//consensus.Profile
	struct Profile
	{
		std::vector<std::pair<uint32_t, std::string>> insertions;	//read -> inserted string
		std::array<int32_t, 17> matches;
		char nucl = '-';

		Profile() {matches.fill(0);}
	};

	//consensus._contig_profile
	std::vector<Profile> contigProfile(const std::vector<samalign::Alignment>& alignment,
									   std::vector<double>& alnErrors)
	{
		if (alignment.empty()) return {};
		int64_t genomeLen = alignment[0].trgLen;
		std::vector<Profile> profile(genomeLen);
		std::unordered_map<std::string, uint32_t> readIds;

		for (const samalign::Alignment& aln : alignment)
		{
			alnErrors.push_back(aln.errRate);
			auto idIt = readIds.find(aln.qryId);
			uint32_t readId = 0;
			if (idIt == readIds.end())
			{
				readId = (uint32_t)readIds.size();
				readIds[aln.qryId] = readId;
			}
			else
			{
				readId = idIt->second;
			}

			//after gap shifting it is possible that two gaps are aligned against each other
			std::string qrySeq = samalign::shiftGaps(aln.trgSeq, aln.qrySeq);
			std::string trgSeq = samalign::shiftGaps(qrySeq, aln.trgSeq);

			int64_t trgPos = aln.trgStart;
			size_t len = std::min(trgSeq.size(), qrySeq.size());
			for (size_t i = 0; i < len; ++i)
			{
				char trgNuc = trgSeq[i];
				char qryNuc = qrySeq[i];
				if (trgNuc == '-') --trgPos;
				if (trgPos >= genomeLen) trgPos -= genomeLen;

				//Python indexing: profile[-1] is the last element
				Profile& elem = profile[trgPos < 0 ? trgPos + genomeLen : trgPos];
				if (trgNuc == '-' && qryNuc != '-')
				{
					bool found = false;
					for (auto& ins : elem.insertions)
					{
						if (ins.first == readId) {ins.second.push_back(qryNuc); found = true; break;}
					}
					if (!found) elem.insertions.emplace_back(readId, std::string(1, qryNuc));
				}
				else
				{
					elem.nucl = trgNuc;
					++elem.matches[symbolIndex(qryNuc)];
				}
				++trgPos;
			}
		}
		return profile;
	}

	//consensus._flatten_profile
	std::string flattenProfile(const std::vector<Profile>& profile)
	{
		std::string growingSeq;
		growingSeq.reserve(profile.size());
		std::map<std::string, int64_t> insGroup;	//sorted, like sorted(ins_group)

		for (const Profile& elem : profile)
		{
			insGroup.clear();
			for (const auto& ins : elem.insertions) ++insGroup[ins.second];

			int64_t matchAndDelNum = 0;
			for (int32_t c : elem.matches) matchAndDelNum += c;
			int64_t delNum = elem.matches[0];		//'-' (accessing it adds the key)
			int64_t numIns = (int64_t)elem.insertions.size();

			//max(sorted(pos_matches), key=pos_matches.get): '-' is always a key now,
			//other keys exist only with positive counts
			char maxMatch = '-';
			int64_t maxCount = elem.matches[0];
			for (size_t s = 1; s < SYMBOLS.size(); ++s)
			{
				if (elem.matches[s] > maxCount)
				{
					maxCount = elem.matches[s];
					maxMatch = SYMBOLS[s];
				}
			}

			const std::string* maxInsert = nullptr;
			int64_t maxInsCount = -1;
			for (const auto& kv : insGroup)
			{
				if (kv.second > maxInsCount) {maxInsCount = kv.second; maxInsert = &kv.first;}
			}

			bool isDeletion = maxMatch == '-' || delNum > matchAndDelNum / 3;
			bool isInsertion = maxInsert && !maxInsert->empty() && *maxInsert != "-" &&
							   numIns > matchAndDelNum / 3;

			if (!isDeletion) growingSeq.push_back(maxMatch);
			if (isInsertion) growingSeq += *maxInsert;
		}
		return growingSeq;
	}
}

samalign::SeqDict getConsensus(const aligner::AlignmentSet& alignment,
							   const samalign::SeqDict& contigs,
							   int numThreads, double& meanAlnError)
{
	const int64_t CHUNK_SIZE = 1000000;
	const int MAX_READ_COVERAGE = 1000;

	samalign::SamReader reader(alignment, contigs, MAX_READ_COVERAGE, /*use secondary*/ true);
	std::vector<samalign::ContigRegion> regions = samalign::chunkRegions(contigs, CHUNK_SIZE);

	struct ChunkResult
	{
		bool done = false;
		std::string sequence;
		std::vector<double> alnErrors;
	};
	std::vector<ChunkResult> results(regions.size());
	std::vector<size_t> tasks(regions.size());
	for (size_t i = 0; i < tasks.size(); ++i) tasks[i] = i;

	std::function<void(const size_t&)> worker =
	[&reader, &regions, &results] (const size_t& idx)
	{
		const samalign::ContigRegion& region = regions[idx];
		std::vector<samalign::Alignment> ctgAln =
			reader.getAlignments(region.ctgId, region.start, region.end);
		if (ctgAln.empty()) return;

		ctgAln = samalign::trimAndTranspose(ctgAln, region.start, region.end);
		double medianCov = 0;
		if (ctgAln.empty())
		{
			//Flye unpacks the empty list returned by get_uniform_alignments and fails
			throw std::runtime_error("No alignments left after trimming region " +
									 region.ctgId);
		}
		ctgAln = samalign::getUniformAlignments(ctgAln, medianCov);

		ChunkResult& res = results[idx];
		std::vector<Profile> profile = contigProfile(ctgAln, res.alnErrors);
		res.sequence = flattenProfile(profile);
		res.done = true;
	};
	processInParallelRethrow(tasks, worker, (size_t)numThreads, false);

	//chunks are concatenated in the order of their start positions
	pycompat::OrderedDict<std::vector<std::pair<int64_t, std::string>>> chunkConsensus;
	std::vector<double> totalAlnErrors;
	for (size_t i = 0; i < regions.size(); ++i)
	{
		if (!results[i].done) continue;
		totalAlnErrors.insert(totalAlnErrors.end(), results[i].alnErrors.begin(),
							  results[i].alnErrors.end());
		if (!results[i].sequence.empty())
		{
			chunkConsensus[regions[i].ctgId].emplace_back(regions[i].start,
														  std::move(results[i].sequence));
		}
	}

	samalign::SeqDict outFasta;
	for (auto& kv : chunkConsensus)
	{
		std::stable_sort(kv.second.begin(), kv.second.end(),
			[](const std::pair<int64_t, std::string>& a,
			   const std::pair<int64_t, std::string>& b) {return a.first < b.first;});
		std::string seq;
		for (auto& chunk : kv.second) seq += chunk.second;
		outFasta[kv.first] = std::move(seq);
	}

	double sumErrors = 0;
	for (double e : totalAlnErrors) sumErrors += e;
	meanAlnError = sumErrors / (double)(totalAlnErrors.size() + 1);
	return outFasta;
}

}	//namespace consensus
