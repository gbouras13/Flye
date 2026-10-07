//(c) 2026 by Authors
//This file is a part of the Flye package.
//Released under the BSD license (see LICENSE file)

//Port of flye/polishing/bubbles.py

#pragma once

#include <string>

#include "aligner.h"
#include "samalign.h"

namespace bubbles
{

//bubbles.make_bubbles: writes bubbles_out, returns per-contig mean coverage
pycompat::OrderedDict<int64_t> makeBubbles(const aligner::AlignmentSet& alignment,
										   const samalign::SeqDict& contigs,
										   const std::string& errMode, int numThreads,
										   const std::string& bubblesOut,
										   double& meanAlnError);

}	//namespace bubbles
