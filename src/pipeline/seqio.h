//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//FASTA/FASTQ I/O with the semantics of flye/utils/fasta_parser.py

#pragma once

#include <string>
#include <functional>
#include <stdexcept>

#include "py_compat.h"

namespace seqio
{

struct FastaError : public std::runtime_error
{
	explicit FastaError(const std::string& msg): std::runtime_error(msg) {}
};

//fasta_parser.stream_sequence: calls fn(header, sequence) per record
//(header = first word; non-ACGT symbols converted with to_acgt)
void streamSequence(const std::string& filename,
					const std::function<void(std::string& header,
											 std::string& seq)>& fn);

//fasta_parser.read_sequence_dict / read_sequence_lengths
pycompat::OrderedDict<std::string> readSequenceDict(const std::string& filename);
pycompat::OrderedDict<size_t> readSequenceLengths(const std::string& filename);

//fasta_parser.write_fasta_dict: sorted by header, 60 symbols per line
void writeFastaDict(const pycompat::OrderedDict<std::string>& fasta,
					const std::string& filename);

std::string reverseComplement(const std::string& seq);

//fasta_parser.to_acgt
std::string toAcgt(const std::string& seq);

}	//namespace seqio
