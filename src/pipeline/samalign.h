//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//Port of flye/utils/sam_parser.py (SynchronizedSamReader) and the alignment
//helpers of flye/polishing/alignment.py

#pragma once

#include <string>
#include <vector>
#include <stdexcept>

#include "aligner.h"
#include "py_compat.h"

namespace samalign
{

struct AlignmentException : public std::runtime_error
{
	explicit AlignmentException(const std::string& msg): std::runtime_error(msg) {}
};

//sam_parser.Alignment
struct Alignment
{
	std::string qryId;
	std::string trgId;
	int64_t qryStart;
	int64_t qryEnd;
	char qrySign;
	int64_t qryLen;
	int64_t trgStart;
	int64_t trgEnd;
	char trgSign;
	int64_t trgLen;
	std::string qrySeq;
	std::string trgSeq;
	double errRate;
	bool isSecondary;
	bool isSupplementary;
	int mapQv;
};

typedef pycompat::OrderedDict<std::string> SeqDict;

//SynchronizedSamReader
class SamReader
{
public:
	SamReader(const aligner::AlignmentSet& alignment, const SeqDict& refFasta,
			  int maxCoverage, bool useSecondary):
		_alignment(alignment), _refFasta(refFasta),
		_maxCoverage(maxCoverage), _useSecondary(useSecondary) {}

	std::string getRegionSequence(const std::string& ctg, int64_t start, int64_t end) const;
	std::vector<Alignment> getAlignments(const std::string& ctg, int64_t start,
										 int64_t end) const;
	std::vector<Alignment> getAllAlignments() const;
	double getMedianDepth(const std::string& ctg, int64_t start, int64_t end) const;

private:
	bool makeAlignment(const aligner::SamRecord& rec, const std::string& ctgId,
					   const std::string& ctgStr, Alignment& out) const;

	const aligner::AlignmentSet& _alignment;
	const SeqDict& _refFasta;
	int _maxCoverage;
	bool _useSecondary;
};

//SynchronizedSamReader.trim_and_transpose
std::vector<Alignment> trimAndTranspose(const std::vector<Alignment>& alignments,
										int64_t regionStart, int64_t regionEnd);

//alignment.get_uniform_alignments: returns selected alignments and median coverage
std::vector<Alignment> getUniformAlignments(const std::vector<Alignment>& alignments,
											double& medianCov);

//alignment.shift_gaps
std::string shiftGaps(const std::string& seqTrg, const std::string& seqQry);

//SynchonizedChunkManager regions
struct ContigRegion
{
	std::string ctgId;
	int64_t start;
	int64_t end;
};
std::vector<ContigRegion> chunkRegions(const SeqDict& reference, int64_t chunkSize);

}	//namespace samalign
