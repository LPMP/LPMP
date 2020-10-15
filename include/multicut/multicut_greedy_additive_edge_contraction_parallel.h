#pragma once

#include "multicut_instance.h"
#include "multicut/multicut_cycle_packing_parallel.h"

namespace LPMP {

    multicut_edge_labeling greedy_additive_edge_contraction_parallel(const multicut_instance& instance, const int nr_threads, const std::string option);

}
