//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <algorithm>
#include <unordered_map>
#include <tuple>

#include "samalign.h"

namespace samalign
{

namespace
{
	enum {OP_M = 0, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_EQ, OP_X};

	char upper(char c) {return (c >= 'a' && c <= 'z') ? (char)(c - 32) : c;}

	void appendUpper(std::string& out, const std::string& src, int64_t pos, int64_t size)
	{
		//Python slicing clamps to the string bounds
		int64_t begin = std::min<int64_t>(std::max<int64_t>(pos, 0), (int64_t)src.size());
		int64_t end = std::min<int64_t>(std::max<int64_t>(pos + size, 0), (int64_t)src.size());
		for (int64_t i = begin; i < end; ++i) out.push_back(upper(src[i]));
	}
}

std::string SamReader::getRegionSequence(const std::string& ctg,
										 int64_t start, int64_t end) const
{
	const std::string& seq = _refFasta.at(ctg);
	start = std::min<int64_t>(start, (int64_t)seq.size());
	end = std::min<int64_t>(end, (int64_t)seq.size());
	if (end <= start) return "";
	return seq.substr(start, end - start);
}

//record -> Alignment, as in get_alignments / get_all_alignments
//(returns false for records these functions skip)
bool SamReader::makeAlignment(const aligner::SamRecord& rec, const std::string& ctgId,
							  const std::string& ctgStr, Alignment& aln) const
{
	bool isUnmapped = rec.flag & 0x4;
	bool isSecondary = rec.flag & 0x100;
	bool isSupplementary = rec.flag & 0x800;
	bool isReversed = rec.flag & 0x16;

	if (isUnmapped) return false;
	if (isSecondary && !_useSecondary) return false;
	if (rec.seq == "*") return false;

	//_parse_cigar
	const std::string& readStr = rec.seq;
	std::string& trgSeq = aln.trgSeq;
	std::string& qrySeq = aln.qrySeq;
	trgSeq.clear();
	qrySeq.clear();
	int64_t trgStart = rec.pos - 1;
	int64_t trgPos = rec.pos - 1;
	int64_t qryStart = 0;
	int64_t qryPos = 0;
	bool leftHard = true;
	bool leftSoft = true;
	int64_t hardLeft = 0, hardRight = 0, softLeft = 0, softRight = 0;

	for (uint32_t c : rec.cigar)
	{
		int64_t size = c >> 4;
		int op = c & 0xf;
		if (op == OP_H)
		{
			if (leftHard) {qryStart += size; hardLeft += size;}
			else hardRight += size;
		}
		else if (op == OP_S)
		{
			qryPos += size;
			if (leftSoft) softLeft += size;
			else softRight += size;
		}
		else if (op == OP_M || op == OP_X || op == OP_EQ)
		{
			appendUpper(qrySeq, readStr, qryPos, size);
			appendUpper(trgSeq, ctgStr, trgPos, size);
			qryPos += size;
			trgPos += size;
		}
		else if (op == OP_I)
		{
			appendUpper(qrySeq, readStr, qryPos, size);
			trgSeq.append(size, '-');
			qryPos += size;
		}
		else if (op == OP_D)
		{
			qrySeq.append(size, '-');
			appendUpper(trgSeq, ctgStr, trgPos, size);
			trgPos += size;
		}
		else
		{
			const char OPS[] = "MIDNSHP=XB";
			throw AlignmentException(std::string("Unsupported CIGAR operation: ") + OPS[op]);
		}
		leftHard = false;
		if (op != OP_H) leftSoft = false;
	}

	int64_t matches = 0;
	size_t alnLen = std::min(trgSeq.size(), qrySeq.size());
	for (size_t i = 0; i < alnLen; ++i)
	{
		if (trgSeq[i] == qrySeq[i]) ++matches;
	}
	if (trgSeq.empty()) throw AlignmentException("Empty alignment in SAM record");
	double errRate = 1.0 - (double)matches / (double)trgSeq.size();

	int64_t trgEnd = trgPos;
	int64_t qryEnd = qryPos + hardLeft;
	int64_t qryLen = qryEnd + hardRight;
	qryStart += softLeft;
	qryEnd -= softRight;

	aln.qryId = rec.qname;
	aln.trgId = ctgId;
	aln.qryStart = qryStart;
	aln.qryEnd = qryEnd;
	aln.qrySign = isReversed ? '-' : '+';
	aln.qryLen = qryLen;
	aln.trgStart = trgStart;
	aln.trgEnd = trgEnd;
	aln.trgSign = '+';
	aln.trgLen = (int64_t)ctgStr.size();
	aln.errRate = errRate;
	aln.isSecondary = isSecondary;
	aln.isSupplementary = isSupplementary;
	aln.mapQv = rec.mapq;
	return true;
}

std::vector<Alignment> SamReader::getAlignments(const std::string& ctg,
												int64_t start, int64_t end) const
{
	const std::string& ctgStr = _refFasta.at(ctg);
	std::vector<const aligner::SamRecord*> chunkBuffer = _alignment.region(ctg, start, end);

	//shuffle alignments so that they uniformly distributed (same seed as Flye)
	pycompat::PyRandom(42).shuffle(chunkBuffer);

	int64_t sequenceLength = 0;
	std::vector<Alignment> alignments;
	for (const aligner::SamRecord* rec : chunkBuffer)
	{
		Alignment aln;
		if (!this->makeAlignment(*rec, ctg, ctgStr, aln)) continue;
		sequenceLength += aln.qryEnd - aln.qryStart;
		alignments.push_back(std::move(aln));
		//Python floor division (both sides non-negative here)
		if (sequenceLength / (int64_t)ctgStr.size() > _maxCoverage) break;
	}

	//sort by read and by aligned length (stable, like Python's sort)
	std::stable_sort(alignments.begin(), alignments.end(),
		[](const Alignment& a, const Alignment& b)
		{
			if (a.qryId != b.qryId) return a.qryId < b.qryId;
			return -(a.qryEnd - a.qryStart) < -(b.qryEnd - b.qryStart);
		});
	return alignments;
}

std::vector<Alignment> SamReader::getAllAlignments() const
{
	std::vector<Alignment> alignments;
	for (const aligner::SamRecord* rec : _alignment.all())
	{
		if (rec->tid < 0) continue;
		const std::string& ctgId = _alignment.refNames[rec->tid];
		Alignment aln;
		if (!this->makeAlignment(*rec, ctgId, _refFasta.at(ctgId), aln)) continue;
		alignments.push_back(std::move(aln));
	}
	return alignments;
}

//`samtools depth <bam> -r ctg:start-end -a -m 0 -Q 10 -l 100` -> median
double SamReader::getMedianDepth(const std::string& ctg, int64_t start, int64_t end) const
{
	const std::string& ctgStr = _refFasta.at(ctg);
	int64_t beg = std::max<int64_t>(start - 1, 0);
	int64_t stop = std::min<int64_t>(end, (int64_t)ctgStr.size());
	if (beg >= stop) return 0;

	std::vector<int64_t> diff(stop - beg + 1, 0);
	bool anyCovered = false;
	for (const aligner::SamRecord* rec : _alignment.region(ctg, start, end))
	{
		if (rec->flag & (0x4 | 0x100 | 0x200 | 0x400)) continue;
		if (rec->mapq < 10) continue;
		int64_t qlen = 0;
		for (uint32_t c : rec->cigar)
		{
			int op = c & 0xf;
			if (op == OP_M || op == OP_I || op == OP_S || op == OP_EQ || op == OP_X)
			{
				qlen += c >> 4;
			}
		}
		if (qlen < 100) continue;

		int64_t refPos = rec->pos - 1;
		for (uint32_t c : rec->cigar)
		{
			int64_t len = c >> 4;
			int op = c & 0xf;
			if (op == OP_M || op == OP_EQ || op == OP_X || op == OP_D || op == OP_N)
			{
				int64_t segBeg = std::max(refPos, beg);
				int64_t segEnd = std::min(refPos + len, stop);
				if (segBeg < segEnd)
				{
					anyCovered = true;	//position reported by the pileup
					if (op == OP_M || op == OP_EQ || op == OP_X)
					{
						diff[segBeg - beg] += 1;
						diff[segEnd - beg] -= 1;
					}
				}
				refPos += len;
			}
		}
	}
	//with -a, all region positions are reported once any position is covered
	if (!anyCovered) return 0;
	std::vector<double> coverage(stop - beg);
	int64_t running = 0;
	for (int64_t i = 0; i < stop - beg; ++i)
	{
		running += diff[i];
		coverage[i] = (double)running;
	}
	return pycompat::getMedian(coverage);
}

std::vector<Alignment> trimAndTranspose(const std::vector<Alignment>& alignments,
										int64_t regionStart, int64_t regionEnd)
{
	const int64_t MIN_ALN = 100;
	std::vector<Alignment> trimmed;
	for (const Alignment& aln : alignments)
	{
		if (aln.trgStart >= regionStart && aln.trgEnd <= regionEnd)
		{
			trimmed.push_back(aln);
			continue;
		}

		//trimming from left
		int64_t newQryStart = aln.qryStart;
		int64_t newTrgStart = aln.trgStart;
		int64_t leftOffset = 0;
		for (leftOffset = 0; leftOffset < (int64_t)aln.trgSeq.size(); ++leftOffset)
		{
			if (newTrgStart >= regionStart) break;
			if (aln.trgSeq[leftOffset] != '-') ++newTrgStart;
			if (aln.qrySeq[leftOffset] != '-') ++newQryStart;
		}
		//Python keeps the last loop value when the loop does not break
		if (!aln.trgSeq.empty() && leftOffset == (int64_t)aln.trgSeq.size()) --leftOffset;

		//trimming from right
		int64_t newQryEnd = aln.qryEnd;
		int64_t newTrgEnd = aln.trgEnd;
		int64_t rightOffset = 0;
		int64_t len = (int64_t)aln.trgSeq.size();
		for (rightOffset = 0; rightOffset < len; ++rightOffset)
		{
			if (newTrgEnd <= regionEnd) break;
			if (aln.trgSeq[len - 1 - rightOffset] != '-') --newTrgEnd;
			if (aln.qrySeq[len - 1 - rightOffset] != '-') --newQryEnd;
		}
		if (len > 0 && rightOffset == len) --rightOffset;

		if (newTrgEnd - newQryEnd > MIN_ALN)
		{
			Alignment newAln = aln;
			int64_t qLen = (int64_t)aln.qrySeq.size();
			int64_t qEnd = std::max(qLen - rightOffset, leftOffset);
			int64_t tEnd = std::max(len - rightOffset, leftOffset);
			newAln.qrySeq = aln.qrySeq.substr(leftOffset, qEnd - leftOffset);
			newAln.trgSeq = aln.trgSeq.substr(leftOffset, tEnd - leftOffset);
			newAln.qryStart = newQryStart;
			newAln.qryEnd = newQryEnd;
			newAln.trgStart = newTrgStart;
			newAln.trgEnd = newTrgEnd;
			trimmed.push_back(std::move(newAln));
		}
	}

	for (Alignment& aln : trimmed)
	{
		aln.trgStart -= regionStart;
		aln.trgEnd -= regionStart;
		aln.trgLen = regionEnd - regionStart;
	}
	return trimmed;
}

std::vector<Alignment> getUniformAlignments(const std::vector<Alignment>& alignments,
											double& medianCov)
{
	medianCov = 0;
	if (alignments.empty()) return {};

	const int64_t WINDOW = 100;
	const int MIN_COV = 20;
	const double GOOD_RATE = 0.66;
	const int MIN_QV = 20;

	auto isReliable = [](const Alignment& aln)
		{return !aln.isSecondary && aln.mapQv >= MIN_QV;};

	int64_t seqLen = alignments[0].trgLen;
	std::vector<int64_t> wndPrimaryCov(seqLen / WINDOW + 1, 0);

	for (const Alignment& aln : alignments)
	{
		for (int64_t i = aln.trgStart / WINDOW; i < aln.trgEnd / WINDOW + 1; ++i)
		{
			if (isReliable(aln)) ++wndPrimaryCov[i];
		}
	}

	int64_t covThreshold = std::max<int64_t>((int64_t)pycompat::getMedian(wndPrimaryCov),
											 MIN_COV);

	auto alnScore = [&wndPrimaryCov, covThreshold](const Alignment& aln,
														   int64_t& good, int64_t& bad)
	{
		good = 0;
		bad = 0;
		for (int64_t i = aln.trgStart / WINDOW; i < aln.trgEnd / WINDOW + 1; ++i)
		{
			if (wndPrimaryCov[i] < covThreshold) ++good;
			else ++bad;
		}
	};

	std::vector<Alignment> selected;
	//sec_aln_scores: dict keyed by read id (insertion order of first key,
	//value of the last assignment)
	struct SecScore {int64_t good; int64_t bad; size_t alnIdx;};
	std::vector<std::string> secOrder;
	std::unordered_map<std::string, SecScore> secScores;

	for (size_t idx = 0; idx < alignments.size(); ++idx)
	{
		const Alignment& aln = alignments[idx];
		if (isReliable(aln))
		{
			selected.push_back(aln);
		}
		else
		{
			SecScore score;
			alnScore(aln, score.good, score.bad);
			score.alnIdx = idx;
			if (!secScores.count(aln.qryId)) secOrder.push_back(aln.qryId);
			secScores[aln.qryId] = score;
		}
	}

	//sorted(..., reverse=True, key=(good - 2 * bad, trg span)): stable descending
	std::stable_sort(secOrder.begin(), secOrder.end(),
		[&secScores, &alignments](const std::string& a, const std::string& b)
		{
			const SecScore& sa = secScores.at(a);
			const SecScore& sb = secScores.at(b);
			const Alignment& aa = alignments[sa.alnIdx];
			const Alignment& ab = alignments[sb.alnIdx];
			auto ka = std::make_tuple(sa.good - 2 * sa.bad, aa.trgEnd - aa.trgStart);
			auto kb = std::make_tuple(sb.good - 2 * sb.bad, ab.trgEnd - ab.trgStart);
			return ka > kb;
		});

	for (const std::string& alnId : secOrder)
	{
		const Alignment& aln = alignments[secScores.at(alnId).alnIdx];
		int64_t good, bad;
		alnScore(aln, good, bad);
		bool toTake = (double)good / (double)(good + bad) > GOOD_RATE;
		if (toTake)
		{
			selected.push_back(aln);
			for (int64_t i = aln.trgStart / WINDOW; i < aln.trgEnd / WINDOW + 1; ++i)
			{
				++wndPrimaryCov[i];
			}
		}
	}

	medianCov = pycompat::getMedian(wndPrimaryCov);
	return selected;
}

std::string shiftGaps(const std::string& seqTrg, const std::string& seqQry)
{
	std::string lstTrg = "$" + seqTrg + "$";
	std::string lstQry = "$" + seqQry + "$";
	bool isGap = false;
	int64_t gapStart = 0;
	for (int64_t i = 0; i < (int64_t)lstTrg.size(); ++i)
	{
		if (isGap && lstQry[i] != '-')
		{
			isGap = false;
			int64_t swapLeft = gapStart - 1;
			int64_t swapRight = i - 1;

			while (swapLeft > 0 && swapRight >= gapStart &&
				   lstQry[swapLeft] == lstTrg[swapRight])
			{
				std::swap(lstQry[swapLeft], lstQry[swapRight]);
				--swapLeft;
				--swapRight;
			}
		}

		if (!isGap && lstQry[i] == '-')
		{
			isGap = true;
			gapStart = i;
		}
	}
	return lstQry.substr(1, lstQry.size() - 2);
}

std::vector<ContigRegion> chunkRegions(const SeqDict& reference, int64_t chunkSize)
{
	std::vector<ContigRegion> regions;
	for (const auto& kv : reference)
	{
		int64_t ctgLen = (int64_t)kv.second.size();
		int64_t chunk = chunkSize > 0 ? chunkSize : ctgLen;
		int64_t numChunks = std::max<int64_t>(chunk > 0 ? ctgLen / chunk : 1, 1);
		for (int64_t i = 0; i < numChunks; ++i)
		{
			int64_t regStart = i * chunk;
			int64_t regEnd = (i + 1) * chunk;
			if (ctgLen - regEnd < chunk) regEnd = ctgLen;
			regions.push_back({kv.first, regStart, regEnd});
		}
	}
	return regions;
}

}	//namespace samalign
