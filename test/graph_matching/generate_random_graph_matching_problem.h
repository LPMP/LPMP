#ifndef LPMP_GENERATE_RANDOM_GRAPH_MATCHING_PROBLEM_H
#define LPMP_GENERATE_RANDOM_GRAPH_MATCHING_PROBLEM_H

#include "graph_matching/matching_problem_input.h"
#include <random>

namespace LPMP {

linear_assignment_problem_input generate_random_graph_matching_problem(const std::size_t no_nodes, std::random_device& rd)
{
   std::mt19937 gen{rd()};
   std::normal_distribution nd(-1.0,2.0);

   linear_assignment_problem_input output;

   for(std::size_t i=0; i<no_nodes; ++i)
      for(std::size_t j=0; j<no_nodes; ++j)
         output.add_assignment(i, j, nd(gen));

   return output;
}

multigraph_matching_input generate_random_multigraph_matching_problem(const std::size_t no_graphs, const std::size_t no_nodes, std::random_device& rd)
{
   multigraph_matching_input output;
   assert(no_graphs >= 3);

   for(std::size_t p=0; p<no_graphs; ++p) {
      for(std::size_t q=p+1; q<no_graphs; ++q) {
         output.push_back({p, q, generate_random_graph_matching_problem(no_nodes, rd)});
      }
   }

   return output;
}

}

#endif // LPMP_GENERATE_RANDOM_GRAPH_MATCHING_PROBLEM_H
