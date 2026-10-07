//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//Port of flye/polishing/polish.py

#pragma once

#include <string>
#include <vector>

namespace polish
{

struct PolishSettings
{
	std::string platform;		//"nano" / "pacbio"
	std::string readType;		//"raw", "corrected", "hifi", "nano_hq", "subasm"
	std::string configDir;		//flye/config/bin_cfg (substitution matrices)
	int numThreads = 1;
};

//polish.polish: returns the polished contigs file; writes contigs_stats.txt
//and base_coverage.bed.gz into workDir
std::string polish(const std::string& contigs, const std::vector<std::string>& reads,
				   const std::string& workDir, int numIters,
				   const PolishSettings& settings, std::string& statsFile);

//polish.filter_by_coverage
void filterByCoverage(bool meta, const std::string& readType,
					  const std::string& statsIn, const std::string& contigsIn,
					  const std::string& statsOut, const std::string& contigsOut);

//polish.generate_polished_edges
void generatePolishedEdges(const std::string& edgesFile, const std::string& gfaFile,
						   const std::string& polishedContigs, const std::string& workDir,
						   const PolishSettings& settings, const std::string& polishedStats);

}	//namespace polish
