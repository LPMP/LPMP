#pragma once

#include "multicut_instance.h"
#include "cut_base/cut_base_packing.h"
#include <vector>

namespace LPMP {

// implementation of the ICP algorithm from Lange et al's ICML18 algorithm.
void multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads);
cycle_packing compute_multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads);

triplet_multicut_instance pack_multicut_instance(const multicut_instance& input, const cycle_packing& cp); 

} // end namespace LPMP
