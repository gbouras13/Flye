//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

#include <fstream>
#include <atomic>
#include <cstring>
#include <zlib.h>

#include "seqio.h"
#include "../common/logger.h"

namespace seqio
{

namespace
{
	//reads a (possibly gzipped) file line by line
	class LineReader
	{
	public:
		explicit LineReader(const std::string& filename):
			_file(gzopen(filename.c_str(), "rb")), _buffer(1 << 20)
		{
			if (!_file) throw FastaError("Can't open file: " + filename);
			gzbuffer(_file, 1 << 20);
		}
		~LineReader() {if (_file) gzclose(_file);}

		bool getLine(std::string& line)
		{
			line.clear();
			bool gotAny = false;
			while (true)
			{
				if (!gzgets(_file, _buffer.data(), (int)_buffer.size()))
				{
					int err = 0;
					gzerror(_file, &err);
					if (err != Z_OK && err != Z_STREAM_END)
					{
						throw FastaError("Error reading compressed input");
					}
					return gotAny;
				}
				gotAny = true;
				size_t len = strlen(_buffer.data());
				line.append(_buffer.data(), len);
				if (len > 0 && _buffer[len - 1] == '\n') return true;
			}
		}

	private:
		gzFile _file;
		std::vector<char> _buffer;
	};

	bool isFastq(const std::string& filename)
	{
		//fasta_parser._is_fastq, including its rstrip(".gz") quirk
		auto lastSuffix = [](const std::string& s)
		{
			size_t pos = s.rfind('.');
			return pos == std::string::npos ? s : s.substr(pos + 1);
		};
		std::string withoutGz = filename;
		if (lastSuffix(filename) == "gz")
		{
			while (!withoutGz.empty() && (withoutGz.back() == '.' ||
				   withoutGz.back() == 'g' || withoutGz.back() == 'z'))
			{
				withoutGz.pop_back();
			}
		}
		std::string suffix = lastSuffix(withoutGz);
		if (suffix == "fasta" || suffix == "fa") return false;
		if (suffix == "fastq" || suffix == "fq") return true;
		throw FastaError("Unknown file extension: " + filename);
	}

	struct CharTables
	{
		bool valid[256];
		bool acgt[256];
		char toAcgt[256];
		char compl_[256];

		CharTables()
		{
			const char* VALID = "ACGTURYKMSWBDHVNXatgcurykmswbvdhnx";
			const char* ACGT = "ACGTacgt";
			const char* FROM = "URYKMSWBVDHNXurykmswbvdhnx";
			const char* TO   = "ACGTACGTACGTAacgtacgtacgta";
			const char* CFROM = "ATGCURYKMSWBVDHNXatgcurykmswbvdhnx";
			const char* CTO   = "TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx";
			for (int i = 0; i < 256; ++i)
			{
				valid[i] = false;
				acgt[i] = false;
				toAcgt[i] = (char)i;
				compl_[i] = (char)i;
			}
			for (const char* c = VALID; *c; ++c) valid[(uint8_t)*c] = true;
			for (const char* c = ACGT; *c; ++c) acgt[(uint8_t)*c] = true;
			for (size_t i = 0; FROM[i]; ++i) toAcgt[(uint8_t)FROM[i]] = TO[i];
			for (size_t i = 0; CFROM[i]; ++i) compl_[(uint8_t)CFROM[i]] = CTO[i];
		}
	};

	const CharTables& tables()
	{
		static CharTables t;
		return t;
	}

	std::atomic<bool> g_acgtWarned(false);

	//fasta_parser._to_acgt_bytes (in place)
	void toAcgtInPlace(std::string& seq)
	{
		const CharTables& t = tables();
		bool allAcgt = true;
		for (char c : seq)
		{
			if (!t.acgt[(uint8_t)c]) {allAcgt = false; break;}
		}
		if (allAcgt) return;
		if (!g_acgtWarned.exchange(true))
		{
			Logger::get().warning() << "Input contain non-ACGT characters - "
				"they will be converted to arbitrary ACGTs";
		}
		for (char& c : seq) c = t.toAcgt[(uint8_t)c];
	}

	bool validateSeq(const std::string& seq)
	{
		const CharTables& t = tables();
		for (char c : seq)
		{
			if (!t.valid[(uint8_t)c]) return false;
		}
		return true;
	}

	std::string firstWord(const std::string& s)
	{
		std::vector<std::string> words = pycompat::pySplit(s);
		if (words.empty()) throw FastaError("Empty sequence header");
		return words[0];
	}
}

void streamSequence(const std::string& filename,
					const std::function<void(std::string& header,
											 std::string& seq)>& fn)
{
	bool fastq = isFastq(filename);
	LineReader reader(filename);
	std::string rawLine;
	size_t entry = 0;

	auto emit = [&fn, &filename, &entry](std::string& hdr, std::string& seq)
	{
		if (hdr.find('@') != std::string::npos)
		{
			throw FastaError("Fasta/q sequence header has '@' symbol in file: " +
							 filename + ", entry " + std::to_string(entry));
		}
		if (!validateSeq(seq))
		{
			throw FastaError("Invalid sequence symbol in file: " + filename +
							 ", entry " + std::to_string(entry));
		}
		toAcgtInPlace(seq);
		fn(hdr, seq);
		++entry;
	};

	if (fastq)
	{
		int state = 0;
		size_t lineNo = 0;
		std::string header;
		std::string seq;
		for (; reader.getLine(rawLine); ++lineNo)
		{
			std::string line = pycompat::pyStrip(rawLine);
			if (line.empty()) continue;
			if (state == 0)
			{
				if (line[0] != '@')
				{
					throw FastaError("Fastq format error: " + filename +
									 " at line " + std::to_string(lineNo));
				}
				header = firstWord(line.substr(1));
			}
			else if (state == 1)
			{
				seq = line;
			}
			else if (state == 2)
			{
				if (line[0] != '+')
				{
					throw FastaError("Fastq format error: " + filename +
									 " at line " + std::to_string(lineNo));
				}
			}
			else
			{
				emit(header, seq);
			}
			state = (state + 1) % 4;
		}
	}
	else
	{
		std::string header;
		std::string seq;
		bool haveHeader = false;
		while (reader.getLine(rawLine))
		{
			std::string line = pycompat::pyStrip(rawLine);
			if (line.empty()) continue;
			if (line[0] == '>')
			{
				if (haveHeader) emit(header, seq);
				seq.clear();
				header = firstWord(line.substr(1));
				haveHeader = true;
			}
			else
			{
				seq += line;
			}
		}
		if (haveHeader && !seq.empty()) emit(header, seq);
	}
}

pycompat::OrderedDict<std::string> readSequenceDict(const std::string& filename)
{
	pycompat::OrderedDict<std::string> dict;
	streamSequence(filename, [&dict](std::string& hdr, std::string& seq)
		{dict[hdr] = std::move(seq);});
	return dict;
}

pycompat::OrderedDict<size_t> readSequenceLengths(const std::string& filename)
{
	pycompat::OrderedDict<size_t> dict;
	streamSequence(filename, [&dict](std::string& hdr, std::string& seq)
		{dict[hdr] = seq.length();});
	return dict;
}

void writeFastaDict(const pycompat::OrderedDict<std::string>& fasta,
					const std::string& filename)
{
	std::ofstream out(filename);
	if (!out.is_open()) throw std::runtime_error("Can't open " + filename);
	for (const std::string& header : fasta.sortedKeys())
	{
		const std::string& seq = fasta.at(header);
		out << ">" << header << "\n";
		for (size_t i = 0; i < seq.length(); i += 60)
		{
			out.write(seq.data() + i, std::min((size_t)60, seq.length() - i));
			out << "\n";
		}
	}
}

std::string reverseComplement(const std::string& seq)
{
	const CharTables& t = tables();
	std::string rc(seq.rbegin(), seq.rend());
	for (char& c : rc) c = t.compl_[(uint8_t)c];
	return rc;
}

std::string toAcgt(const std::string& seq)
{
	std::string out = seq;
	toAcgtInPlace(out);
	return out;
}

}	//namespace seqio
