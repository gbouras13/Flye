//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//Port of flye/polishing/consensus.py

#pragma once

#include <string>

#include "aligner.h"
#include "samalign.h"

namespace consensus
{

//consensus.get_consensus: per-contig consensus of the aligned reads
samalign::SeqDict getConsensus(const aligner::AlignmentSet& alignment,
							   const samalign::SeqDict& contigs,
							   int numThreads, double& meanAlnError);

}	//namespace consensus
