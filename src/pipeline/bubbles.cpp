//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cstdlib>

#include "bubbles.h"
#include "seqio.h"
#include "parallel_util.h"
#include "../common/logger.h"

namespace bubbles
{

namespace
{
	//py_cfg.vals
	const int SIMPLE_KMER_LENGTH = 4;
	const int SOLID_KMER_LENGTH = 10;
	const int MAX_BUBBLE_LENGTH = 500;
	const int MAX_BUBBLE_BRANCHES = 50;
	const int MAX_READ_COVERAGE = 1000;
	const int MIN_POLISH_ALN_LEN = 500;
	const double SOLID_MISSMATCH = 0.3;	//same for "pacbio" and "nano"
	const double SOLID_INDEL = 0.3;

	struct ProfileInfo
	{
		char nucl = 0;
		std::vector<std::pair<uint32_t, std::string>> insertions;	//read -> inserted string
		int64_t propagatedIns = 0;
		int64_t numDeletions = 0;
		int64_t numMissmatch = 0;
		int64_t coverage = 0;
	};

	struct Bubble
	{
		std::string contigId;
		int64_t position;
		int64_t subPosition = 0;
		std::vector<std::string> branches;
		std::string consensus;

		Bubble(const std::string& ctg, int64_t pos): contigId(ctg), position(pos) {}
	};

	//bubbles._compute_profile
	std::vector<ProfileInfo> computeProfile(const std::vector<samalign::Alignment>& alignment,
											const std::string& refSequence,
											std::vector<double>& alnErrors)
	{
		if (alignment.empty()) throw std::runtime_error("No alignmemnts!");
		int64_t genomeLen = alignment[0].trgLen;
		int64_t minAlnLen = std::min<int64_t>(MIN_POLISH_ALN_LEN, genomeLen / 2);
		std::vector<ProfileInfo> profile(genomeLen);
		for (int64_t i = 0; i < genomeLen; ++i) profile[i].nucl = refSequence.at(i);

		std::unordered_map<std::string, uint32_t> readIds;
		for (const samalign::Alignment& aln : alignment)
		{
			if ((int64_t)aln.qrySeq.size() < minAlnLen) continue;
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

			std::string qrySeq = samalign::shiftGaps(aln.trgSeq, aln.qrySeq);
			std::string trgSeq = samalign::shiftGaps(qrySeq, aln.trgSeq);

			int64_t trgPos = aln.trgStart;
			size_t len = std::min(trgSeq.size(), qrySeq.size());
			for (size_t i = 0; i < len; ++i)
			{
				char trgNuc = trgSeq[i];
				char qryNuc = qrySeq[i];
				if (trgNuc == '-') --trgPos;

				//Python indexing: profile[-1] is the last element
				ProfileInfo& elem = profile[trgPos < 0 ? trgPos + genomeLen : trgPos];
				if (trgNuc == '-')
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
					++elem.coverage;
					if (qryNuc == '-') ++elem.numDeletions;
					else if (trgNuc != qryNuc) ++elem.numMissmatch;
				}
				++trgPos;
			}
		}

		for (int64_t i = 0; i < genomeLen; ++i)
		{
			for (const auto& ins : profile[i].insertions)
			{
				++profile[i].propagatedIns;
				int64_t span = (int64_t)ins.second.size();
				for (int64_t j = std::max<int64_t>(0, i - span); j < i; ++j)
				{
					++profile[j].propagatedIns;
				}
				for (int64_t j = i + 1; j < std::min(i + span + 1, genomeLen); ++j)
				{
					++profile[j].propagatedIns;
				}
			}
		}
		return profile;
	}

	//bubbles._is_solid_kmer
	bool isSolidKmer(const std::vector<ProfileInfo>& profile, int64_t position)
	{
		for (int64_t i = position; i < position + SOLID_KMER_LENGTH; ++i)
		{
			if (profile[i].coverage == 0) return false;
			double localMissmatch = (double)(profile[i].numMissmatch + profile[i].numDeletions) /
									(double)profile[i].coverage;
			double localIns = (double)profile[i].propagatedIns / (double)profile[i].coverage;
			if (localMissmatch > SOLID_MISSMATCH || localIns > SOLID_INDEL) return false;
		}
		return true;
	}

	//bubbles._is_simple_kmer
	bool isSimpleKmer(const std::vector<ProfileInfo>& profile, int64_t position)
	{
		const int extendedLen = SIMPLE_KMER_LENGTH * 2;
		std::vector<char> nuclStr;
		int64_t from = position - extendedLen / 2;
		int64_t to = std::min<int64_t>(position + extendedLen / 2, (int64_t)profile.size());
		for (int64_t i = from; i < to; ++i) nuclStr.push_back(profile[i].nucl);
		auto at = [&nuclStr](int64_t i) {return nuclStr.at(i);};

		//single nucleotide homopolymers
		for (int i = extendedLen / 2 - SIMPLE_KMER_LENGTH / 2;
			 i < extendedLen / 2 + SIMPLE_KMER_LENGTH / 2 - 1; ++i)
		{
			if (at(i) == at(i + 1)) return false;
		}

		//dinucleotide homopolymers
		for (int shift = 0; shift <= 1; ++shift)
		{
			for (int i = 0; i < SIMPLE_KMER_LENGTH - shift - 1; ++i)
			{
				int pos = extendedLen / 2 - SIMPLE_KMER_LENGTH + shift + i * 2;
				if (at(pos) == at(pos + 2) && at(pos + 1) == at(pos + 3)) return false;
			}
		}
		return true;
	}

	//bubbles._get_partition
	std::vector<int64_t> getPartition(const std::vector<ProfileInfo>& profile,
									  int64_t& longBubbles)
	{
		int64_t profLen = (int64_t)profile.size();
		std::vector<bool> solidFlags(profLen, false);
		int64_t profPos = 0;
		while (profPos < profLen - SOLID_KMER_LENGTH)
		{
			if (isSolidKmer(profile, profPos))
			{
				for (int64_t i = profPos; i < profPos + SOLID_KMER_LENGTH; ++i) solidFlags[i] = true;
				profPos += SOLID_KMER_LENGTH;
			}
			else
			{
				++profPos;
			}
		}

		std::vector<int64_t> partition;
		int64_t prevPartition = SOLID_KMER_LENGTH;
		longBubbles = 0;
		profPos = SOLID_KMER_LENGTH;
		while (profPos < profLen - SOLID_KMER_LENGTH)
		{
			int64_t curPartition = profPos + SIMPLE_KMER_LENGTH / 2;
			bool allSolid = true;
			for (int64_t i = profPos; i < std::min<int64_t>(profPos + SIMPLE_KMER_LENGTH, profLen); ++i)
			{
				if (!solidFlags[i]) {allSolid = false; break;}
			}
			bool landmark = allSolid && isSimpleKmer(profile, curPartition);

			if (profPos - prevPartition > MAX_BUBBLE_LENGTH) ++longBubbles;

			if (landmark || profPos - prevPartition > MAX_BUBBLE_LENGTH)
			{
				partition.push_back(curPartition);
				prevPartition = curPartition;
				profPos += SOLID_KMER_LENGTH;
			}
			else
			{
				++profPos;
			}
		}
		return partition;
	}

	//Python's str.replace("-", "")
	std::string removeGaps(const std::string& s, size_t from, size_t to)
	{
		std::string out;
		to = std::min(to, s.size());
		for (size_t i = from; i < to; ++i)
		{
			if (s[i] != '-') out.push_back(s[i]);
		}
		return out;
	}

	//bubbles._get_bubble_seqs
	std::vector<Bubble> getBubbleSeqs(const std::vector<samalign::Alignment>& alignment,
									  const std::vector<ProfileInfo>& profile,
									  const std::vector<int64_t>& partition,
									  const std::string& contigId)
	{
		if (partition.empty() || alignment.empty()) return {};
		int64_t ctgLen = alignment[0].trgLen;

		std::vector<int64_t> extPartition;
		extPartition.push_back(0);
		extPartition.insert(extPartition.end(), partition.begin(), partition.end());
		extPartition.push_back(ctgLen);

		std::vector<Bubble> bubbles;
		for (size_t i = 0; i + 1 < extPartition.size(); ++i)
		{
			bubbles.emplace_back(contigId, extPartition[i]);
			std::string consensus;
			for (int64_t p = extPartition[i];
				 p < std::min<int64_t>(extPartition[i + 1], (int64_t)profile.size()); ++p)
			{
				consensus.push_back(profile[p].nucl);
			}
			bubbles.back().consensus = consensus;
		}

		//bisect.bisect (right)
		auto bisect = [&extPartition](int64_t x)
		{
			return (int64_t)(std::upper_bound(extPartition.begin(), extPartition.end(), x) -
							 extPartition.begin());
		};

		for (const samalign::Alignment& aln : alignment)
		{
			int64_t bubbleId = bisect(aln.trgStart) - 1;
			int64_t nextBubbleStart = extPartition.at(bubbleId + 1);
			bool chromosomeEnd = aln.trgEnd >= extPartition.back();

			bool incompleteSegment = aln.trgStart > extPartition[bubbleId];
			int64_t trgPos = aln.trgStart;
			size_t branchStart = 0;
			for (size_t i = 0; i < aln.trgSeq.size(); ++i)
			{
				if (aln.trgSeq[i] == '-') continue;

				if (trgPos >= nextBubbleStart)
				{
					if (!incompleteSegment)
					{
						std::string branchSeq = seqio::toAcgt(removeGaps(aln.qrySeq,
																		 branchStart, i));
						bubbles.at(bubbleId).branches.push_back(std::move(branchSeq));
					}
					incompleteSegment = false;
					bubbleId = bisect(trgPos) - 1;
					nextBubbleStart = extPartition.at(bubbleId + 1);
					branchStart = i;
				}
				++trgPos;
			}

			if (chromosomeEnd)
			{
				std::string branchSeq = seqio::toAcgt(removeGaps(aln.qrySeq, branchStart,
																 aln.qrySeq.size()));
				bubbles.back().branches.push_back(std::move(branchSeq));
			}
		}
		return bubbles;
	}

	//sorted(branches, key=len)[len // 2]
	const std::string& medianBranch(const std::vector<std::string>& branches)
	{
		std::vector<size_t> order(branches.size());
		for (size_t i = 0; i < order.size(); ++i) order[i] = i;
		std::stable_sort(order.begin(), order.end(),
			[&branches](size_t a, size_t b) {return branches[a].size() < branches[b].size();});
		return branches[order[branches.size() / 2]];
	}

	//bubbles._postprocess_bubbles
	std::vector<Bubble> postprocessBubbles(std::vector<Bubble>& bubbles, int64_t& emptyBubbles)
	{
		std::vector<Bubble> newBubbles;
		emptyBubbles = 0;
		for (Bubble& bubble : bubbles)
		{
			if (bubble.branches.empty())
			{
				++emptyBubbles;
				continue;
			}

			std::string median = medianBranch(bubble.branches);
			if (median.empty())
			{
				++emptyBubbles;
				continue;
			}

			//only take branches that are not significantly differ in length from the median
			std::vector<std::string> newBranches;
			for (const std::string& branch : bubble.branches)
			{
				double inconsRate = std::abs((double)branch.size() - (double)median.size()) /
									(double)median.size();
				if (inconsRate < 0.5 && !branch.empty()) newBranches.push_back(branch);
			}

			if (newBranches.empty())
			{
				++emptyBubbles;
				continue;
			}

			//if bubble consensus has very different length from all the branches, replace
			//consensus with the median branch instead
			if (std::abs((int64_t)median.size() - (int64_t)bubble.consensus.size()) >
				(int64_t)median.size() / 2)
			{
				bubble.consensus = median;
			}

			//finally, keep only MAX_BRANCHES
			if ((int64_t)newBranches.size() > MAX_BUBBLE_BRANCHES)
			{
				int64_t left = (int64_t)newBranches.size() / 2 - MAX_BUBBLE_BRANCHES / 2;
				int64_t right = left + MAX_BUBBLE_BRANCHES;
				std::stable_sort(newBranches.begin(), newBranches.end(),
					[](const std::string& a, const std::string& b) {return a.size() < b.size();});
				newBranches = std::vector<std::string>(newBranches.begin() + left,
													   newBranches.begin() + right);
			}

			newBubbles.emplace_back(bubble.contigId, bubble.position);
			newBubbles.back().consensus = bubble.consensus;
			newBubbles.back().branches = std::move(newBranches);
		}
		return newBubbles;
	}

	//bubbles._split_long_bubbles
	std::vector<Bubble> splitLongBubbles(std::vector<Bubble>& bubbles, int64_t& longBranches)
	{
		std::vector<Bubble> newBubbles;
		longBranches = 0;
		for (Bubble& bubble : bubbles)
		{
			const std::string& median = medianBranch(bubble.branches);
			int64_t numChunks = (int64_t)median.size() / MAX_BUBBLE_LENGTH;
			if (numChunks > 1)
			{
				++longBranches;
				for (int64_t partNum = 0; partNum < numChunks; ++partNum)
				{
					std::vector<std::string> newBranches;
					for (const std::string& b : bubble.branches)
					{
						int64_t chunkLen = (int64_t)b.size() / numChunks;
						int64_t start = partNum * chunkLen;
						int64_t end = (partNum != numChunks - 1) ? (partNum + 1) * chunkLen
																 : (int64_t)b.size();
						newBranches.push_back(b.substr(start, end - start));
					}
					newBubbles.emplace_back(bubble.contigId, bubble.position);
					newBubbles.back().consensus = newBranches[0];
					newBubbles.back().branches = std::move(newBranches);
					newBubbles.back().subPosition = partNum;
				}
			}
			else
			{
				newBubbles.push_back(std::move(bubble));
			}
		}
		return newBubbles;
	}

	//bubbles._output_bubbles
	void outputBubbles(const std::vector<Bubble>& bubbles, std::ostream& out)
	{
		for (const Bubble& bubble : bubbles)
		{
			if (bubble.branches.empty()) throw std::runtime_error("No branches in a bubble");
			out << ">" << bubble.contigId << " " << bubble.position << " "
				<< bubble.branches.size() << " " << bubble.subPosition << "\n";
			out << bubble.consensus << "\n";
			for (size_t branchId = 0; branchId < bubble.branches.size(); ++branchId)
			{
				out << ">" << branchId << "\n";
				out << bubble.branches[branchId] << "\n";
			}
		}
	}
}

pycompat::OrderedDict<int64_t> makeBubbles(const aligner::AlignmentSet& alignment,
										   const samalign::SeqDict& contigs,
										   const std::string& errMode, int numThreads,
										   const std::string& bubblesOut,
										   double& meanAlnError)
{
	(void)errMode;	//all error modes use the same solid k-mer rates
	const int64_t CHUNK_SIZE = 1000000;

	samalign::SamReader reader(alignment, contigs, MAX_READ_COVERAGE, /*use secondary*/ true);
	std::vector<samalign::ContigRegion> regions = samalign::chunkRegions(contigs, CHUNK_SIZE);

	struct ChunkResult
	{
		bool done = false;
		std::string bubblesText;
		int64_t numBubbles = 0;
		int64_t numLongBubbles = 0;
		int64_t numEmpty = 0;
		int64_t numLongBranch = 0;
		std::vector<double> alnErrors;
		double meanCov = 0;
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
		std::string refSeq = reader.getRegionSequence(region.ctgId, region.start, region.end);

		//since we are working with contig chunks, tranform alignment coorinates
		ctgAln = samalign::trimAndTranspose(ctgAln, region.start, region.end);
		if (ctgAln.empty())
		{
			//Flye unpacks the empty list returned by get_uniform_alignments and fails
			throw std::runtime_error("No alignments left after trimming region " +
									 region.ctgId);
		}
		ChunkResult& res = results[idx];
		double meanCov = 0;
		ctgAln = samalign::getUniformAlignments(ctgAln, meanCov);

		std::vector<ProfileInfo> profile = computeProfile(ctgAln, refSeq, res.alnErrors);
		std::vector<int64_t> partition = getPartition(profile, res.numLongBubbles);
		std::vector<Bubble> ctgBubbles = getBubbleSeqs(ctgAln, profile, partition,
													   region.ctgId);

		const double coverageCap = 0.9 * MAX_READ_COVERAGE;
		if (meanCov > coverageCap)
		{
			meanCov = reader.getMedianDepth(region.ctgId, region.start, region.end);
		}

		ctgBubbles = postprocessBubbles(ctgBubbles, res.numEmpty);
		ctgBubbles = splitLongBubbles(ctgBubbles, res.numLongBranch);

		//transform coordinates back
		for (Bubble& b : ctgBubbles) b.position += region.start;

		std::ostringstream text;
		outputBubbles(ctgBubbles, text);
		res.bubblesText = text.str();
		res.numBubbles = (int64_t)ctgBubbles.size();
		res.meanCov = meanCov;
		res.done = true;
	};
	processInParallelRethrow(tasks, worker, (size_t)numThreads, false);

	//chunks are written in region order (Flye's order with a single thread)
	std::ofstream fout(bubblesOut, std::ios::out | std::ios::trunc);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + bubblesOut);

	int64_t totalBubbles = 0, totalLongBubbles = 0, totalLongBranches = 0, totalEmpty = 0;
	std::vector<double> totalAlnErrors;
	pycompat::OrderedDict<std::vector<double>> coverageLists;
	for (size_t i = 0; i < regions.size(); ++i)
	{
		ChunkResult& res = results[i];
		if (!res.done) continue;
		fout << res.bubblesText;
		totalLongBubbles += res.numLongBubbles;
		totalLongBranches += res.numLongBranch;
		totalEmpty += res.numEmpty;
		totalAlnErrors.insert(totalAlnErrors.end(), res.alnErrors.begin(), res.alnErrors.end());
		totalBubbles += res.numBubbles;
		coverageLists[regions[i].ctgId].push_back(res.meanCov);
	}
	fout.close();

	pycompat::OrderedDict<int64_t> coverageStats;
	for (auto& kv : coverageLists)
	{
		double sum = 0;
		for (double c : kv.second) sum += c;
		coverageStats[kv.first] = (int64_t)(sum / (double)kv.second.size());
	}

	double sumErrors = 0;
	for (double e : totalAlnErrors) sumErrors += e;
	meanAlnError = sumErrors / (double)(totalAlnErrors.size() + 1);
	Logger::get().debug() << "Generated " << totalBubbles << " bubbles";
	Logger::get().debug() << "Split " << totalLongBubbles << " long bubbles";
	Logger::get().debug() << "Skipped " << totalEmpty << " empty bubbles";
	Logger::get().debug() << "Skipped " << totalLongBranches << " bubbles with long branches";
	return coverageStats;
}

}	//namespace bubbles
