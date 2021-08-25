#pragma once

#include "multicut_instance.h"

namespace LPMP {

    multicut_edge_labeling greedy_additive_edge_contraction(const multicut_instance& instance);
    std::tuple<multicut_edge_labeling, std::vector<int>> greedy_additive_edge_contraction_impl(const multicut_instance& instance);
}
