#pragma once

#include "multicut_instance.h"

namespace LPMP {

    multicut_edge_labeling multicut_greedy_edge_fixation(const multicut_instance& instance);
    std::tuple<multicut_edge_labeling, std::vector<int>> multicut_greedy_edge_fixation_impl(const multicut_instance& instance);

}
