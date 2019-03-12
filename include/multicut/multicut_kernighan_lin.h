#pragma once

#include "multicut_instance.h"
#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/kernighan-lin.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"
#include "andres/graph/multicut/greedy-fixation.hxx"

namespace LPMP {

   multicut_edge_labeling compute_gaec(const multicut_instance& instance);
   multicut_edge_labeling compute_multicut_kernighan_lin(const multicut_instance& instance, multicut_edge_labeling labeling = multicut_edge_labeling());
   multicut_edge_labeling compute_multicut_gaec_kernighan_lin(const multicut_instance& instance);
   multicut_edge_labeling compute_multicut_greedy_edge_fixation(const multicut_instance& instance);

}
