//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//`flye-modules pipeline`: the whole Flye pipeline (flye/main.py) in a single
//process, without Python, minimap2/samtools subprocesses or SAM/BAM files.
//The Python logic (configure, consensus, polishing, finalize) is ported to
//C++ with identical results; the C++ modules are called in-process.

#include <string>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>
#include <limits.h>

#include "py_compat.h"
#include "seqio.h"
#include "aligner.h"
#include "consensus.h"
#include "polish.h"
#include "finalize.h"
#include "../common/logger.h"
#include "../common/portable_rand.h"
#include "../sequence/sequence_container.h"

int assemble_main(int argc, char** argv);
int repeat_main(int argc, char** argv);
int contigger_main(int argc, char** argv);

namespace
{
	const char* FLYE_VERSION = "2.9.6-b1802";
	const int PIPELINE_VERSION = 3;

	struct PipelineError : public std::runtime_error
	{
		explicit PipelineError(const std::string& msg): std::runtime_error(msg) {}
	};

	struct Args
	{
		std::vector<std::string> reads;
		std::string platform;
		std::string readType;
		int64_t genomeSize = 0;
		std::string outDir;
		int threads = 1;
		int numIters = 1;
		int minOverlap = 0;		//0 = auto
		int asmCoverage = 0;	//0 = not set
		double readError = 0;
		std::string extraParams;
		bool meta = false;
		bool keepHaplotypes = false;
		bool noAltContigs = false;
		bool scaffold = false;
		bool debug = false;
		bool deterministic = false;
		std::string stopAfter;
		bool resume = false;
		std::string resumeFrom;
		std::string configDir;
		std::string logFile;
		std::string asmConfig;
	};

	bool isDir(const std::string& path)
	{
		struct stat st;
		return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
	}

	bool fileExists(const std::string& path)
	{
		struct stat st;
		return stat(path.c_str(), &st) == 0;
	}

	int64_t fileSize(const std::string& path)
	{
		struct stat st;
		if (stat(path.c_str(), &st) != 0) return 0;
		return (int64_t)st.st_size;
	}

	void makeDir(const std::string& path)
	{
		if (!isDir(path) && mkdir(path.c_str(), 0755) != 0)
		{
			throw PipelineError("Can't create directory " + path);
		}
	}

	std::string joinPath(const std::string& dir, const std::string& name)
	{
		if (dir.empty() || dir.back() == '/') return dir + name;
		return dir + "/" + name;
	}

	std::string absPath(const std::string& path)
	{
		if (!path.empty() && path[0] == '/') return path;
		char buf[PATH_MAX];
		if (!getcwd(buf, sizeof(buf))) return path;
		std::string p = path;
		while (p.compare(0, 2, "./") == 0) p = p.substr(2);
		return joinPath(buf, p);
	}

	void copyFile(const std::string& from, const std::string& to)
	{
		std::ifstream in(from, std::ios::binary);
		if (!in.is_open()) throw PipelineError("Can't open " + from);
		std::ofstream out(to, std::ios::binary);
		out << in.rdbuf();
	}

	//flye.utils.bytes2human.human2bytes (after upper())
	int64_t human2bytes(const std::string& s)
	{
		std::string num;
		size_t i = 0;
		while (i < s.size() && ((s[i] >= '0' && s[i] <= '9') || s[i] == '.')) num += s[i++];
		std::string letter = pycompat::pyStrip(s.substr(i));
		const std::string SYMBOLS = "BKMGTPEZY";
		if (letter.size() != 1 || SYMBOLS.find(letter[0]) == std::string::npos || num.empty())
		{
			throw PipelineError("can't interpret '" + s + "'");
		}
		double value = std::stod(num);
		double prefix = std::pow(10.0, 3 * (double)SYMBOLS.find(letter[0]));
		return (int64_t)(value * prefix);
	}

	//configurator._calc_nx
	int64_t calcNx(std::vector<int64_t> lengths, int64_t total, double rate)
	{
		std::stable_sort(lengths.begin(), lengths.end(), std::greater<int64_t>());
		int64_t sum = 0;
		for (int64_t l : lengths)
		{
			sum += l;
			if ((double)sum > rate * (double)total) return l;
		}
		return 0;
	}

	//configurator._get_downsample_threshold
	int64_t downsampleThreshold(std::vector<int64_t> lengths, int64_t targetLen)
	{
		std::stable_sort(lengths.begin(), lengths.end(), std::greater<int64_t>());
		int64_t sum = 0;
		for (int64_t l : lengths)
		{
			sum += l;
			if (sum > targetLen) return l;
		}
		return 0;
	}

	struct RunParams
	{
		int64_t minOverlap = 0;
		int64_t minReadLength = 0;
		bool configured = false;
	};

	//configurator.setup_params
	RunParams setupParams(const Args& args)
	{
		Logger::get().info() << "Configuring run";
		RunParams params;
		int64_t totalLength = 0;
		std::vector<int64_t> readLengths;
		const int64_t MAX_READ_LEN = (1LL << 31) - 1;
		const int64_t MIN_OVLP = args.readType == "subasm" ? 1000 : 1000;
		const int64_t MAX_OVLP = args.readType == "subasm" ? 1000 : 10000;

		int64_t lowestReadLen = args.minOverlap ? args.minOverlap : MIN_OVLP;
		int64_t passingReads = 0;
		for (const std::string& readFile : args.reads)
		{
			for (auto& kv : seqio::readSequenceLengths(readFile))
			{
				int64_t seqLen = (int64_t)kv.second;
				if (seqLen > MAX_READ_LEN)
				{
					throw PipelineError("Length of single read in '" + readFile +
										"' exceeded maximum (" + std::to_string(MAX_READ_LEN) + ")");
				}
				if (seqLen > lowestReadLen) ++passingReads;
				totalLength += seqLen;
				readLengths.push_back(seqLen);
			}
		}
		if (!passingReads)
		{
			throw PipelineError("No reads above minimum length threshold (" +
								std::to_string(lowestReadLen) + ")");
		}

		int64_t readsN50 = calcNx(readLengths, totalLength, 0.50);
		int64_t readsN90 = calcNx(readLengths, totalLength, 0.90);

		Logger::get().info() << "Total read length: " << totalLength;
		int64_t coverage = 0;
		if (args.genomeSize)
		{
			coverage = totalLength / args.genomeSize;
			Logger::get().info() << "Input genome size: " << args.genomeSize;
			Logger::get().info() << "Estimated coverage: " << coverage;
			if (coverage < 5 || coverage > 1000)
			{
				Logger::get().warning() << "Expected read coverage is " << coverage
					<< ", the assembly is not guaranteed to be optimal in this setting."
					<< " Are you sure that the genome size was entered correctly?";
			}
		}

		Logger::get().info() << "Reads N50/N90: " << readsN50 << " / " << readsN90;
		if (!args.minOverlap)
		{
			const int64_t GRADE = 1000;
			int64_t intMinOvlp = (int64_t)pycompat::pyRound((double)readsN90 / GRADE) * GRADE;
			int64_t maxOvlp = MAX_OVLP;
			if (args.meta) maxOvlp = std::min<int64_t>(maxOvlp, 10000);
			params.minOverlap = std::max(MIN_OVLP, std::min(maxOvlp, intMinOvlp));
			Logger::get().info() << "Minimum overlap set to " << params.minOverlap;
		}
		else
		{
			params.minOverlap = args.minOverlap;
			Logger::get().info() << "Selected minimum overlap: " << params.minOverlap;
		}

		//downsampling reads for the first assembly stage to save memory
		if (args.asmCoverage && args.asmCoverage < coverage)
		{
			Logger::get().info() << "Using longest " << args.asmCoverage << "x reads for contig assembly";
			params.minReadLength = downsampleThreshold(readLengths, args.genomeSize * args.asmCoverage);
			Logger::get().debug() << "Min read length cutoff: " << params.minReadLength;
		}
		params.configured = true;
		return params;
	}

	//runs a flye-modules module in this process
	int runModule(int (*moduleMain)(int, char**), std::vector<std::string> argv)
	{
		//state that a fresh process would start with
		portableRandReset();
		SequenceContainer::g_nextSeqId = 0;
#if defined(__APPLE__) || defined(__FreeBSD__)
		optreset = 1;
#endif
		optind = 0;

		std::string cmd;
		for (auto& a : argv) cmd += " " + a;
		Logger::get().debug() << "Running: flye-modules" << cmd;

		std::vector<char*> cargv;
		for (auto& a : argv) cargv.push_back(&a[0]);
		cargv.push_back(nullptr);
		int ret = moduleMain((int)argv.size(), cargv.data());
		optind = 0;
		return ret;
	}

	std::string joinReads(const std::vector<std::string>& reads)
	{
		std::string out;
		for (size_t i = 0; i < reads.size(); ++i) out += (i ? "," : "") + reads[i];
		return out;
	}

	void saveParams(const Args& args, const RunParams& params, const std::string& stage)
	{
		std::ofstream out(joinPath(args.outDir, "params.json"));
		if (params.configured)
		{
			out << "{\"pipeline_version\": " << PIPELINE_VERSION << ", \"min_overlap\": "
				<< params.minOverlap << ", \"min_read_length\": " << params.minReadLength
				<< ", \"stage_name\": \"" << stage << "\"}";
		}
		else
		{
			out << "{\"stage_name\": \"" << stage << "\"}";
		}
	}

	std::string findConfigDir(const std::string& requested, const char* argv0)
	{
		std::vector<std::string> candidates;
		if (!requested.empty()) candidates.push_back(requested);
		if (const char* env = getenv("FLYE_CONFIG_DIR")) candidates.push_back(env);
		std::string self = argv0 ? argv0 : "";
		size_t slash = self.rfind('/');
		std::string selfDir = slash == std::string::npos ? "." : self.substr(0, slash);
		candidates.push_back(joinPath(selfDir, "../flye/config/bin_cfg"));
		candidates.push_back("/flye/config/bin_cfg");		//embedded in WebAssembly builds
		for (const std::string& dir : candidates)
		{
			if (fileExists(joinPath(dir, "asm_nano_hq.cfg"))) return dir;
		}
		throw PipelineError("Flye configuration directory not found (use --config-dir)");
	}

	void printUsage()
	{
		std::cerr << "Usage: flye-modules pipeline (--pacbio-raw | --pacbio-corr | --pacbio-hifi |\n"
			"\t--nano-raw | --nano-corr | --nano-hq ) file1 [file_2 ...] --out-dir PATH\n"
			"\t[--genome-size SIZE] [--threads int] [--iterations int] [--meta]\n"
			"\t[--min-overlap SIZE] [--asm-coverage int] [--read-error float]\n"
			"\t[--keep-haplotypes] [--no-alt-contigs] [--scaffold] [--extra-params str]\n"
			"\t[--resume] [--resume-from stage] [--stop-after stage] [--deterministic]\n"
			"\t[--debug] [--config-dir PATH]\n";
	}

	void addExtraParam(Args& args, const std::string& param)
	{
		args.extraParams = args.extraParams.empty() ? param : args.extraParams + "," + param;
	}

	//argparse-like parsing of flye's command line
	bool parseArgs(int argc, char** argv, Args& args)
	{
		struct ReadMode {const char* flag; const char* platform; const char* type;};
		const ReadMode MODES[] = {
			{"--pacbio-raw", "pacbio", "raw"}, {"--pacbio-corr", "pacbio", "corrected"},
			{"--pacbio-hifi", "pacbio", "hifi"}, {"--nano-raw", "nano", "raw"},
			{"--nano-corr", "nano", "corrected"}, {"--nano-hq", "nano", "nano_hq"},
			{"--subassemblies", "pacbio", "subasm"}};

		std::string genomeSize;
		double hifiError = 0;
		auto needValue = [argc, argv](int& i) -> std::string
		{
			if (i + 1 >= argc) throw PipelineError(std::string("argument ") + argv[i] + ": expected one argument");
			return argv[++i];
		};
		auto checkRange = [](const std::string& name, int value, int lo, int hi)
		{
			if (value < lo || value > hi)
			{
				throw PipelineError("argument " + name + ": value should be in the range [" +
									std::to_string(lo) + ", " + std::to_string(hi) + "]");
			}
			return value;
		};

		for (int i = 1; i < argc; ++i)
		{
			std::string arg = argv[i];
			bool isMode = false;
			for (const ReadMode& m : MODES)
			{
				if (arg != m.flag) continue;
				if (!args.readType.empty()) throw PipelineError("only one read type option is allowed");
				isMode = true;
				args.platform = m.platform;
				args.readType = m.type;
				while (i + 1 < argc && argv[i + 1][0] != '-') args.reads.push_back(argv[++i]);
				if (args.reads.empty()) throw PipelineError(arg + ": expected at least one argument");
			}
			if (isMode) continue;

			if (arg == "-g" || arg == "--genome-size") genomeSize = needValue(i);
			else if (arg == "-o" || arg == "--out-dir") args.outDir = needValue(i);
			else if (arg == "-t" || arg == "--threads") args.threads = checkRange(arg, std::stoi(needValue(i)), 1, 128);
			else if (arg == "-i" || arg == "--iterations") args.numIters = checkRange(arg, std::stoi(needValue(i)), 0, 10);
			else if (arg == "-m" || arg == "--min-overlap") args.minOverlap = checkRange(arg, std::stoi(needValue(i)), 1000, 10000);
			else if (arg == "--asm-coverage") args.asmCoverage = std::stoi(needValue(i));
			else if (arg == "--hifi-error") hifiError = std::stod(needValue(i));
			else if (arg == "--read-error") args.readError = std::stod(needValue(i));
			else if (arg == "--extra-params") args.extraParams = needValue(i);
			else if (arg == "--meta") args.meta = true;
			else if (arg == "--keep-haplotypes") args.keepHaplotypes = true;
			else if (arg == "--no-alt-contigs") args.noAltContigs = true;
			else if (arg == "--scaffold") args.scaffold = true;
			else if (arg == "--debug") args.debug = true;
			else if (arg == "--deterministic") args.deterministic = true;
			else if (arg == "--stop-after") args.stopAfter = needValue(i);
			else if (arg == "--resume") args.resume = true;
			else if (arg == "--resume-from") args.resumeFrom = needValue(i);
			else if (arg == "--config-dir") args.configDir = needValue(i);
			else if (arg == "--plasmids" || arg == "--trestle") {}
			else if (arg == "-h" || arg == "--help") {printUsage(); return false;}
			else if (arg == "-v" || arg == "--version") {std::cout << FLYE_VERSION << std::endl; return false;}
			else throw PipelineError("unrecognized argument: " + arg);
		}

		if (args.readType.empty()) throw PipelineError("one of the read type options is required");
		if (args.outDir.empty()) throw PipelineError("the following arguments are required: -o/--out-dir");
		if (args.asmCoverage && genomeSize.empty())
		{
			throw PipelineError("--asm-coverage option requires genome size estimate (--genome-size)");
		}
		if (args.asmCoverage && args.meta) throw PipelineError("--asm-coverage is incompatible with --meta");
		if (hifiError != 0 && args.readError == 0) args.readError = hifiError;
		if (args.readError != 0 && args.readType == "raw")
		{
			throw PipelineError("--read-error can only be used with corr/hq/hifi modes");
		}
		if (args.readError > 1)
		{
			throw PipelineError("--read-error expressed as a decimal fraction, e.g. 0.01 or 0.03");
		}

		if (args.readError != 0)
		{
			std::ostringstream err;
			err << args.readError;	//Python formats floats with repr(); fine for typical values
			addExtraParam(args, "assemble_ovlp_divergence=" + err.str() +
						  ",repeat_graph_ovlp_divergence=" + err.str());
		}
		if (args.noAltContigs) addExtraParam(args, "remove_alt_edges=1");
		if (args.keepHaplotypes) addExtraParam(args, "aggressive_dup_filter=0");

		if (!genomeSize.empty())
		{
			bool digits = std::all_of(genomeSize.begin(), genomeSize.end(),
									  [](char c) {return c >= '0' && c <= '9';});
			std::string upper = genomeSize;
			for (char& c : upper) c = (char)toupper(c);
			args.genomeSize = digits ? std::stoll(genomeSize) : human2bytes(upper);
		}
		return true;
	}
}


namespace
{
	//reads a value from the flat params.json written by saveParams / Flye's Job.save
	bool jsonValue(const std::string& json, const std::string& key, std::string& value)
	{
		std::string pattern = "\"" + key + "\":";
		size_t pos = json.find(pattern);
		if (pos == std::string::npos) return false;
		pos += pattern.size();
		while (pos < json.size() && json[pos] == ' ') ++pos;
		if (pos < json.size() && json[pos] == '"')
		{
			size_t end = json.find('"', pos + 1);
			if (end == std::string::npos) return false;
			value = json.substr(pos + 1, end - pos - 1);
			return true;
		}
		size_t end = json.find_first_of(",}", pos);
		if (end == std::string::npos) return false;
		value = pycompat::pyStrip(json.substr(pos, end - pos));
		return true;
	}

	std::string readFile(const std::string& path)
	{
		std::ifstream in(path);
		std::stringstream ss;
		ss << in.rdbuf();
		return ss.str();
	}

	struct Stage
	{
		std::string name;
		std::vector<std::string> outFiles;
		std::function<void()> run;
	};
}

int pipeline_main(int argc, char** argv)
{
	Args args;
	try
	{
		if (!parseArgs(argc, argv, args)) return 0;
	}
	catch (const std::exception& e)
	{
		printUsage();
		std::cerr << "flye-modules pipeline: error: " << e.what() << std::endl;
		return 2;
	}

	try
	{
		makeDir(args.outDir);
		args.outDir = absPath(args.outDir);
		for (std::string& r : args.reads) r = absPath(r);
		args.logFile = joinPath(args.outDir, "flye.log");
		Logger::get().setOutputFile(args.logFile);
		Logger::get().setDebugging(args.debug);
		args.configDir = findConfigDir(args.configDir, argv[0]);
		const char* asmConfigs[][2] = {{"raw", "asm_raw_reads.cfg"},
			{"corrected", "asm_corrected_reads.cfg"}, {"hifi", "asm_hifi.cfg"},
			{"nano_hq", "asm_nano_hq.cfg"}, {"subasm", "asm_subasm.cfg"}};
		for (auto& c : asmConfigs)
		{
			if (args.readType == c[0]) args.asmConfig = joinPath(args.configDir, c[1]);
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}

	Logger::get().info() << "Starting Flye " << FLYE_VERSION << " (single-process pipeline)";
	std::string cmd;
	for (int i = 0; i < argc; ++i) cmd += (i ? " " : "") + std::string(argv[i]);
	Logger::get().debug() << "Cmd: " << cmd;

	try
	{
		for (const std::string& readFile : args.reads)
		{
			if (!fileExists(readFile)) throw PipelineError("Can't open " + readFile);
			if (readFile.find(' ') != std::string::npos)
			{
				throw PipelineError("Path to reads contain spaces: " + readFile);
			}
		}

		const std::string threads = std::to_string(args.threads);
		const std::string assemblyDir = joinPath(args.outDir, "00-assembly");
		const std::string draftAssembly = joinPath(assemblyDir, "draft_assembly.fasta");
		const std::string consensusDir = joinPath(args.outDir, "10-consensus");
		const std::string consensusFile = joinPath(consensusDir, "consensus.fasta");
		const std::string repeatDir = joinPath(args.outDir, "20-repeat");
		const std::string contiggerDir = joinPath(args.outDir, "30-contigger");
		const std::string polishingDir = joinPath(args.outDir, "40-polishing");
		const std::string rawContigs = joinPath(contiggerDir, "contigs.fasta");
		const std::string repeatStats = joinPath(contiggerDir, "contigs_stats.txt");
		const std::string saveFile = joinPath(args.outDir, "params.json");
		const bool withConsensus = args.readType != "subasm";
		const bool withPolishing = args.numIters > 0;
		const std::string disjointigs = withConsensus ? consensusFile : draftAssembly;
		const std::string contigsFile = withPolishing ?
			joinPath(polishingDir, "filtered_contigs.fasta") : rawContigs;
		const std::string polishedStats = withPolishing ?
			joinPath(polishingDir, "filtered_stats.txt") : "";
		const std::string polishedGfa = withPolishing ?
			joinPath(polishingDir, "polished_edges.gfa") : joinPath(contiggerDir, "graph_final.gfa");
		RunParams params;

		auto extraArgs = [&args](std::vector<std::string>& cmdline)
		{
			if (!args.extraParams.empty())
			{
				cmdline.push_back("--extra-params");
				cmdline.push_back(args.extraParams);
			}
		};

		//pipeline as a list of consecutive stages (main._create_job_list)
		std::vector<Stage> stages;
		stages.push_back({"configure", {}, [&]()
		{
			params = setupParams(args);
		}});

		stages.push_back({"assembly", {draftAssembly}, [&]()
		{
			makeDir(assemblyDir);
			Logger::get().info() << "Assembling disjointigs";
			std::vector<std::string> cmdline = {"assemble", "--reads", joinReads(args.reads),
				"--out-asm", draftAssembly, "--config", args.asmConfig, "--log", args.logFile,
				"--threads", args.deterministic ? "1" : threads};
			if (args.debug) cmdline.push_back("--debug");
			if (args.meta) cmdline.push_back("--meta");
			if (args.genomeSize)
			{
				cmdline.push_back("--genome-size");
				cmdline.push_back(std::to_string(args.genomeSize));
			}
			cmdline.push_back("--min-ovlp");
			cmdline.push_back(std::to_string(params.minOverlap));
			if (params.minReadLength > 0)
			{
				cmdline.push_back("--min-read");
				cmdline.push_back(std::to_string(params.minReadLength));
			}
			extraArgs(cmdline);
			if (runModule(assemble_main, cmdline) != 0) throw PipelineError("Assembly module failed");
			if (fileSize(draftAssembly) == 0)
			{
				throw PipelineError("No disjointigs were assembled - please check if the read type "
									"and genome size parameters are correct");
			}
			auto asmStats = finalize::shortStatistics(draftAssembly);
			Logger::get().debug() << "Disjointigs length: " << asmStats.first
				<< ", N50: " << asmStats.second;
		}});

		if (withConsensus)
		{
			stages.push_back({"consensus", {consensusFile}, [&]()
			{
				makeDir(consensusDir);
				Logger::get().info() << "Running Minimap2";
				aligner::AlignmentSet alignment =
					aligner::makeAlignment(draftAssembly, args.reads, args.threads,
										   args.platform, args.readType);
				samalign::SeqDict contigs = seqio::readSequenceDict(draftAssembly);
				Logger::get().info() << "Computing consensus";
				double meanAlnError = 0;
				samalign::SeqDict consensusFasta =
					consensus::getConsensus(alignment, contigs, args.threads, meanAlnError);
				char errBuf[64];
				snprintf(errBuf, sizeof(errBuf), "%f", meanAlnError);
				Logger::get().info() << "Alignment error rate: " << errBuf;
				seqio::writeFastaDict(consensusFasta, consensusFile);
			}});
		}

		stages.push_back({"repeat", {joinPath(repeatDir, "repeat_graph_dump"),
									 joinPath(repeatDir, "repeat_graph_edges.fasta"),
									 joinPath(repeatDir, "read_alignment_dump")}, [&]()
		{
			makeDir(repeatDir);
			Logger::get().info() << "Building and resolving repeat graph";
			std::vector<std::string> cmdline = {"repeat", "--disjointigs", disjointigs,
				"--reads", joinReads(args.reads), "--out-dir", repeatDir, "--config", args.asmConfig,
				"--log", args.logFile, "--threads", threads};
			if (args.debug) cmdline.push_back("--debug");
			if (args.meta) cmdline.push_back("--meta");
			if (args.keepHaplotypes) cmdline.push_back("--keep-haplotypes");
			cmdline.push_back("--min-ovlp");
			cmdline.push_back(std::to_string(params.minOverlap));
			extraArgs(cmdline);
			if (runModule(repeat_main, cmdline) != 0) throw PipelineError("Repeat module failed");
		}});

		stages.push_back({"contigger", {rawContigs, joinPath(contiggerDir, "graph_final.gv"),
										joinPath(contiggerDir, "graph_final.fasta"),
										joinPath(contiggerDir, "graph_final.gfa"), repeatStats,
										joinPath(contiggerDir, "scaffolds_links.txt")}, [&]()
		{
			makeDir(contiggerDir);
			Logger::get().info() << "Generating contigs";
			std::vector<std::string> cmdline = {"contigger",
				"--graph-edges", joinPath(repeatDir, "repeat_graph_edges.fasta"),
				"--reads", joinReads(args.reads), "--out-dir", contiggerDir,
				"--config", args.asmConfig,
				"--repeat-graph", joinPath(repeatDir, "repeat_graph_dump"),
				"--graph-aln", joinPath(repeatDir, "read_alignment_dump"),
				"--log", args.logFile, "--threads", threads};
			if (args.debug) cmdline.push_back("--debug");
			if (args.keepHaplotypes) cmdline.push_back("--no-scaffold");
			cmdline.push_back("--min-ovlp");
			cmdline.push_back(std::to_string(params.minOverlap));
			extraArgs(cmdline);
			if (runModule(contigger_main, cmdline) != 0) throw PipelineError("Contigger module failed");
			if (fileSize(rawContigs) == 0)
			{
				throw PipelineError("No contigs were assembled - pipeline stopped");
			}
		}});

		if (withPolishing)
		{
			stages.push_back({"polishing", {contigsFile, polishedStats, polishedGfa}, [&]()
			{
				makeDir(polishingDir);
				polish::PolishSettings settings;
				settings.platform = args.platform;
				settings.readType = args.readType;
				settings.configDir = args.configDir;
				settings.numThreads = args.threads;

				std::string statsFile;
				std::string polished = polish::polish(rawContigs, args.reads, polishingDir,
													  args.numIters, settings, statsFile);
				polish::filterByCoverage(args.meta, args.readType, statsFile, polished,
										 polishedStats, contigsFile);
				polish::generatePolishedEdges(joinPath(contiggerDir, "graph_final.fasta"),
											  joinPath(contiggerDir, "graph_final.gfa"),
											  contigsFile, polishingDir, settings, statsFile);
				std::remove(polished.c_str());
				if (fileSize(contigsFile) == 0)
				{
					throw PipelineError("No contigs were assembled - pipeline stopped");
				}
			}});
		}

		const std::string assembly = joinPath(args.outDir, "assembly.fasta");
		stages.push_back({"finalize", {assembly, joinPath(args.outDir, "assembly_info.txt"),
									   joinPath(args.outDir, "assembly_graph.gv"),
									   joinPath(args.outDir, "assembly_graph.gfa")}, [&]()
		{
			copyFile(joinPath(contiggerDir, "graph_final.gv"),
					 joinPath(args.outDir, "assembly_graph.gv"));
			copyFile(polishedGfa, joinPath(args.outDir, "assembly_graph.gfa"));
			const std::string links = args.scaffold ? joinPath(contiggerDir, "scaffolds_links.txt") : "";
			finalize::Scaffolds scaffolds = finalize::generateScaffolds(contigsFile, links, assembly);
			finalize::generateStats(repeatStats, polishedStats, scaffolds,
									joinPath(args.outDir, "assembly_info.txt"));
			Logger::get().info() << "Final assembly: " << assembly;
		}});

		bool knownStop = args.stopAfter.empty();
		for (const Stage& st : stages) knownStop |= st.name == args.stopAfter;
		if (!knownStop) throw PipelineError("Stop after: unknown stage '" + args.stopAfter + "'");

		//resuming (main._run)
		size_t currentJob = 0;
		if (args.resume || !args.resumeFrom.empty())
		{
			if (!fileExists(saveFile))
			{
				throw PipelineError("Can't resume, please check if the output directory contains "
									"intermediate results from previous Flye run. Otherwise, make "
									"a new run without --resume option");
			}
			Logger::get().info() << "Resuming previous run";
			const std::string saved = readFile(saveFile);
			std::string jobToResume = args.resumeFrom;
			if (jobToResume.empty() && !jsonValue(saved, "stage_name", jobToResume))
			{
				throw PipelineError("Can't resume: no stage name in " + saveFile);
			}

			bool canResume = false;
			for (size_t i = 0; i < stages.size(); ++i)
			{
				if (stages[i].name != jobToResume) continue;
				std::string value;
				if (!jsonValue(saved, "pipeline_version", value) ||
					value != std::to_string(PIPELINE_VERSION))
				{
					throw PipelineError("Inconsistent pipeline version");
				}
				if (jsonValue(saved, "min_overlap", value)) params.minOverlap = std::stoll(value);
				if (jsonValue(saved, "min_read_length", value)) params.minReadLength = std::stoll(value);
				params.configured = true;
				currentJob = i;
				//(Python indexes jobs[i - 1], i.e. the last stage for i = 0)
				const Stage& prev = stages[(i + stages.size() - 1) % stages.size()];
				for (const std::string& f : prev.outFiles)
				{
					if (!fileExists(f))
					{
						throw PipelineError("Can't resume: stage '" + prev.name + "' incomplete");
					}
				}
				canResume = true;
				break;
			}
			if (!canResume)
			{
				throw PipelineError("Can't resume: stage " + jobToResume + " does not exist");
			}
		}

		for (size_t i = currentJob; i < stages.size(); ++i)
		{
			saveParams(args, params, stages[i].name);
			Logger::get().info() << ">>>STAGE: " << stages[i].name;
			stages[i].run();
			if (args.stopAfter == stages[i].name)
			{
				if (i + 1 < stages.size()) saveParams(args, params, stages[i + 1].name);
				Logger::get().info() << "Pipeline stopped as requested by --stop-after";
				break;
			}
		}
	}
	catch (const std::exception& e)
	{
		Logger::get().error() << e.what();
		Logger::get().error() << "Pipeline aborted";
		return 1;
	}
	return 0;
}
