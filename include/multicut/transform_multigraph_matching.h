#ifndef LPMP_TRANSFORM_MATCHING_PROBLEM_INPUT_H
#define LPMP_TRANSFORM_MATCHING_PROBLEM_INPUT_H

#include "graph_matching/matching_problem_input.h"
#include "multicut_instance.hxx"

namespace LPMP {

   namespace detail {
      struct multigraph_matching_nodes; // forward declaration
   }

   multicut_instance transform_multigraph_matching_to_correlation_clustering(const multigraph_matching_input& input, const detail::multigraph_matching_nodes& nodes);
   multicut_instance transform_multigraph_matching_to_correlation_clustering(const multigraph_matching_input& input);

   multigraph_matching_input::labeling transform_correlation_clustering_to_multigraph_matching(
      const multigraph_matching_input& mgm,
      const multicut_instance& cc,
      multicut_instance::edge_labeling l_input,
      const detail::multigraph_matching_nodes& nodes);
   multigraph_matching_input::labeling transform_correlation_clustering_to_multigraph_matching(
      const multigraph_matching_input& mgm,
      const multicut_instance& cc,
      multicut_instance::edge_labeling l_input);

   //multicut_instance::labeling transform_multigraph_matching_labeling_to_correlation_clustering(const multigraph_matching_input::labeling& input);

}

#endif // LPMP_TRANSFORM_MATCHING_PROBLEM_INPUT_H 
