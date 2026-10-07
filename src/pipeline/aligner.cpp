//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <sys/stat.h>

#include "aligner.h"
#include "minimap.h"

namespace aligner
{

namespace
{
	//htslib's seq_nt16_table / seq_nt16_str: SEQ goes through 4-bit encoding
	//when samtools converts SAM to BAM
	const unsigned char SEQ_NT16_TABLE[256] = {
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
		15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
		15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
		15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
		15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
		15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
	};
	const char SEQ_NT16_STR[] = "=ACMGRSVTWYHKDBN";
	const char CIGAR_OPS[] = "MIDNSHP=XB";

	enum {OP_M = 0, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_EQ, OP_X};

	int cigarOpCode(char c)
	{
		const char* p = strchr(CIGAR_OPS, c);
		if (!p || !c) throw std::runtime_error(std::string("Bad CIGAR op: ") + c);
		return (int)(p - CIGAR_OPS);
	}

	bool consumesRef(int op)
	{
		return op == OP_M || op == OP_D || op == OP_N || op == OP_EQ || op == OP_X;
	}

	struct Collector
	{
		const std::unordered_map<std::string, int32_t>* refIds;
		std::vector<SamRecord>* records;
		std::string error;
	};

	void splitTabs(const char* line, std::vector<std::pair<const char*, size_t>>& fields)
	{
		fields.clear();
		const char* start = line;
		for (const char* p = line; ; ++p)
		{
			if (*p == '\t' || *p == '\0')
			{
				fields.emplace_back(start, (size_t)(p - start));
				if (*p == '\0') break;
				start = p + 1;
			}
		}
	}

	//parses a SAM line as written by minimap2 into a record as samtools
	//would give it back after the BAM round trip
	void onSamLine(const char* line, void* data)
	{
		Collector& col = *static_cast<Collector*>(data);
		if (!col.error.empty() || line[0] == '@') return;

		thread_local std::vector<std::pair<const char*, size_t>> f;
		splitTabs(line, f);
		if (f.size() < 11)
		{
			col.error = "Malformed SAM line from minimap2";
			return;
		}
		SamRecord rec;
		rec.qname.assign(f[0].first, f[0].second);
		rec.flag = (uint16_t)strtol(std::string(f[1].first, f[1].second).c_str(), nullptr, 10);
		std::string rname(f[2].first, f[2].second);
		if (rname == "*")
		{
			rec.tid = -1;
		}
		else
		{
			auto it = col.refIds->find(rname);
			if (it == col.refIds->end())
			{
				col.error = "Unknown reference in SAM: " + rname;
				return;
			}
			rec.tid = it->second;
		}
		rec.pos = strtoll(std::string(f[3].first, f[3].second).c_str(), nullptr, 10);
		rec.mapq = (uint8_t)strtol(std::string(f[4].first, f[4].second).c_str(), nullptr, 10);

		//CIGAR
		std::string cigarStr(f[5].first, f[5].second);
		if (cigarStr != "*")
		{
			uint32_t num = 0;
			for (char c : cigarStr)
			{
				if (c >= '0' && c <= '9') num = num * 10 + (c - '0');
				else
				{
					rec.cigar.push_back(num << 4 | (uint32_t)cigarOpCode(c));
					num = 0;
				}
			}
		}

		//SEQ through 4-bit encoding
		std::string seqStr(f[9].first, f[9].second);
		if (seqStr == "*")
		{
			rec.seq = "*";
		}
		else
		{
			rec.seq.resize(seqStr.size());
			for (size_t i = 0; i < seqStr.size(); ++i)
			{
				rec.seq[i] = SEQ_NT16_STR[SEQ_NT16_TABLE[(uint8_t)seqStr[i]]];
			}
		}
		size_t lQseq = rec.seq == "*" ? 0 : rec.seq.size();

		//minimap2 -L moves CIGARs with >65535 operations into the CG tag;
		//htslib (bam_tag2cigar) moves them back when reading the BAM
		if (!rec.cigar.empty() && rec.tid >= 0 && rec.pos > 0 &&
			(rec.cigar[0] & 0xf) == (uint32_t)OP_S && (rec.cigar[0] >> 4) == lQseq)
		{
			for (size_t i = 11; i < f.size(); ++i)
			{
				if (f[i].second > 7 && strncmp(f[i].first, "CG:B:I,", 7) == 0)
				{
					std::vector<uint32_t> realCigar;
					std::string vals(f[i].first + 7, f[i].second - 7);
					size_t start = 0;
					while (start < vals.size())
					{
						size_t comma = vals.find(',', start);
						if (comma == std::string::npos) comma = vals.size();
						realCigar.push_back((uint32_t)strtoul(vals.substr(start,
														comma - start).c_str(), nullptr, 10));
						start = comma + 1;
					}
					if (realCigar.size() >= rec.cigar.size()) rec.cigar.swap(realCigar);
					break;
				}
			}
		}

		//bam_endpos
		int64_t refLen = 0;
		for (uint32_t c : rec.cigar)
		{
			if (consumesRef(c & 0xf)) refLen += c >> 4;
		}
		if (refLen == 0) refLen = 1;
		rec.refEnd = rec.pos - 1 + refLen;
		col.records->push_back(std::move(rec));
	}

	size_t fileSize(const std::string& path)
	{
		struct stat st;
		if (stat(path.c_str(), &st) != 0) return 0;
		return (size_t)st.st_size;
	}
}

AlignmentSet makeAlignment(const std::string& referenceFile,
						   const std::vector<std::string>& readFiles,
						   int numThreads, const std::string& platform,
						   const std::string& readType)
{
	//options as in alignment._run_minimap
	const std::string mode = platform + "-" + readType;
	const char* preset = nullptr;
	bool kmer17 = false;
	if (mode == "pacbio-raw" || mode == "pacbio-corrected" || mode == "pacbio-subasm")
	{
		preset = "map-pb";
	}
	if (mode == "pacbio-hifi")
	{
		preset = "map-hifi";
	}
	else if (mode == "nano-raw" || mode == "nano-corrected")
	{
		preset = "map-ont";
	}
	else if (mode == "nano-nano_hq")
	{
		preset = "map-ont";
		kmer17 = true;
	}
	if (!preset) throw std::runtime_error("Unknown read type for alignment: " + mode);

	const bool bigReference = fileSize(referenceFile) > 100 * 1024 * 1024;
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_set_opt(0, &iopt, &mopt);
	if (mm_set_opt(preset, &iopt, &mopt) < 0) throw std::runtime_error("minimap2 preset");
	if (kmer17) iopt.k = 17;										//-k 17
	mopt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;							//-a
	mopt.pri_ratio = 0.5f;											//-p 0.5
	mopt.best_n = 10;												//-N 10
	mopt.flag |= MM_F_SAM_HIT_ONLY;									//--sam-hit-only
	mopt.flag |= MM_F_LONG_CIGAR;									//-L
	mopt.mini_batch_size = bigReference ? 5000000000LL : 1000000000LL;	//-K
	mopt.zdrop = mopt.zdrop_inv = 1000;								//-z 1000
	mopt.flag |= MM_F_NO_QUAL;										//-Q
	//--secondary-seq: Flye's minimap2 registers it under the code of --rmq,
	//so in Flye 2.9.6 it has no effect (secondary records keep SEQ = '*')
	iopt.batch_size = 64000000000ULL;								//-I 64G
	if (mm_check_opt(&iopt, &mopt) < 0) throw std::runtime_error("Bad minimap2 options");

	AlignmentSet result;
	int oldVerbose = mm_verbose;
	mm_verbose = 1;

	mm_idx_reader_t* reader = mm_idx_reader_open(referenceFile.c_str(), &iopt, 0);
	if (!reader) throw std::runtime_error("Can't open reference: " + referenceFile);
	std::vector<SamRecord> unsorted;
	Collector collector;
	collector.refIds = &result.refIds;
	collector.records = &unsorted;

	mm_idx_t* index = nullptr;
	int parts = 0;
	while ((index = mm_idx_reader_read(reader, numThreads)) != 0)
	{
		if (++parts > 1) throw std::runtime_error("Unexpected multi-part minimap2 index");
		for (uint32_t i = 0; i < index->n_seq; ++i)
		{
			result.refIds[index->seq[i].name] = (int32_t)result.refNames.size();
			result.refNames.push_back(index->seq[i].name);
		}
		mm_mapopt_update(&mopt, index);
		mm_output_hook = onSamLine;
		mm_output_hook_data = &collector;
		for (const std::string& readFile : readFiles)
		{
			if (mm_map_file(index, readFile.c_str(), &mopt, numThreads) < 0)
			{
				mm_output_hook = nullptr;
				mm_idx_destroy(index);
				mm_idx_reader_close(reader);
				mm_verbose = oldVerbose;
				throw std::runtime_error("minimap2 failed to map " + readFile);
			}
		}
		mm_output_hook = nullptr;
		mm_output_hook_data = nullptr;
		mm_idx_destroy(index);
	}
	mm_idx_reader_close(reader);
	mm_verbose = oldVerbose;
	if (!collector.error.empty()) throw std::runtime_error(collector.error);

	//samtools sort: stable by (tid, pos, strand)
	std::vector<size_t> order(unsorted.size());
	for (size_t i = 0; i < order.size(); ++i) order[i] = i;
	std::stable_sort(order.begin(), order.end(),
		[&unsorted](size_t a, size_t b)
		{
			const SamRecord& ra = unsorted[a];
			const SamRecord& rb = unsorted[b];
			uint64_t pa = (uint64_t)(uint32_t)ra.tid << 32 | (uint64_t)(uint32_t)ra.pos;
			uint64_t pb = (uint64_t)(uint32_t)rb.tid << 32 | (uint64_t)(uint32_t)rb.pos;
			if (pa != pb) return pa < pb;
			return (int)ra.isReverse() < (int)rb.isReverse();
		});
	result.records.reserve(unsorted.size());
	for (size_t i : order) result.records.push_back(std::move(unsorted[i]));

	result.refOffsets.assign(result.refNames.size() + 1, result.records.size());
	for (size_t t = 0, i = 0; t <= result.refNames.size(); ++t)
	{
		while (i < result.records.size() && result.records[i].tid >= 0 &&
			   (size_t)result.records[i].tid < t) ++i;
		result.refOffsets[t] = i;
	}
	return result;
}

std::vector<const SamRecord*> AlignmentSet::region(const std::string& ctg,
												   int64_t start, int64_t end) const
{
	std::vector<const SamRecord*> out;
	auto it = refIds.find(ctg);
	if (it == refIds.end()) return out;
	//htslib hts_parse_reg: beg = max(start - 1, 0), end as given
	int64_t beg = std::max<int64_t>(start - 1, 0);
	if (beg >= end) return out;
	int32_t tid = it->second;
	for (size_t i = refOffsets[tid]; i < refOffsets[tid + 1]; ++i)
	{
		const SamRecord& rec = records[i];
		int64_t pos0 = rec.pos - 1;
		if (pos0 >= end) break;
		if (rec.refEnd > beg) out.push_back(&rec);
	}
	return out;
}

std::vector<const SamRecord*> AlignmentSet::all() const
{
	std::vector<const SamRecord*> out;
	out.reserve(records.size());
	for (const SamRecord& rec : records) out.push_back(&rec);
	return out;
}

}	//namespace aligner
