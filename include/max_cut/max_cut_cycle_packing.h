#pragma once

#include "max_cut_instance.hxx"
#include "cut_base/cut_base_packing.h"
#include <vector>

namespace LPMP {

// modification for max-cut of the multicut ICP algorithm from Lange et al's ICML18 algorithm.
void max_cut_cycle_packing(const max_cut_instance& input);
cycle_packing compute_max_cut_cycle_packing(const max_cut_instance& input);

triplet_max_cut_instance pack_max_cut_instance(const max_cut_instance& input, const cycle_packing& cp); 

} // end namespace LPMP
