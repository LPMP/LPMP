#ifndef LPMP_MULTICUT_KERNIGHAN_LIN_H
#define LPMP_MULTICUT_KERNIGHAN_LIN_H

#include "multicut_instance.hxx"
#include "andres/graph/graph.hxx"
#include "andres/graph/multicut/kernighan-lin.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"

namespace LPMP {

   multicut_instance::edge_labeling compute_gaec(const multicut_instance& instance);
   multicut_instance::edge_labeling compute_multicut_kernighan_lin(const multicut_instance& instance, multicut_instance::edge_labeling labeling = multicut_instance::edge_labeling());
   multicut_instance::edge_labeling compute_multicut_gaec_kernighan_lin(const multicut_instance& instance);

}

#endif // LPMP_MULTICUT_KERNIGHAN_LIN_H
