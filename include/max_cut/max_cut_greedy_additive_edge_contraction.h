#pragma once

#include "max_cut_instance.hxx"

namespace LPMP {

    max_cut_edge_labeling greedy_additive_edge_contraction(const max_cut_instance& instance);

}

