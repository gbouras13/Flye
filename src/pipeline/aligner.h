//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//In-process replacement for `flye-minimap2 ... | samtools view -u | samtools sort`
//as run by flye/polishing/alignment.py, plus region queries as done by
//`samtools view <bam> ctg:start-end`.

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace aligner
{

struct SamRecord
{
	std::string qname;
	uint16_t flag;
	int32_t tid;
	int64_t pos;					//1-based, as in SAM
	uint8_t mapq;
	std::vector<uint32_t> cigar;	//BAM encoding: length << 4 | op
	std::string seq;				//after the BAM round trip; "*" if absent
	int64_t refEnd;					//0-based exclusive end (bam_endpos)

	bool isReverse() const {return flag & 0x10;}
};

class AlignmentSet
{
public:
	std::vector<std::string> refNames;
	std::unordered_map<std::string, int32_t> refIds;
	std::vector<SamRecord> records;		//coordinate-sorted, like samtools sort
	std::vector<size_t> refOffsets;		//records of tid t: [refOffsets[t], refOffsets[t + 1])

	//records overlapping region `ctg:start-end` (Python passes 0-based start
	//into the 1-based region string), in file order
	std::vector<const SamRecord*> region(const std::string& ctg,
										 int64_t start, int64_t end) const;
	std::vector<const SamRecord*> all() const;
};

//alignment.make_alignment(reference_file, reads_files, num_proc, platform, read_type, ...)
AlignmentSet makeAlignment(const std::string& referenceFile,
						   const std::vector<std::string>& readFiles,
						   int numThreads, const std::string& platform,
						   const std::string& readType);

}	//namespace aligner
