//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//Port of flye/assembly/scaffolder.py

#pragma once

#include <string>
#include <vector>
#include <utility>

#include "py_compat.h"

namespace finalize
{

typedef pycompat::OrderedDict<std::vector<std::string>> Scaffolds;

//scaffolder.generate_scaffolds (linksFile empty = no scaffolding)
Scaffolds generateScaffolds(const std::string& contigsFile, const std::string& linksFile,
							const std::string& outScaffolds);

//scaffolder.generate_stats (polishedFile empty = no polishing stats)
void generateStats(const std::string& repeatFile, const std::string& polishedFile,
				   const Scaffolds& scaffolds, const std::string& outStats);

//scaffolder.short_statistics: total length and N50
std::pair<int64_t, int64_t> shortStatistics(const std::string& fastaFile);

}	//namespace finalize
