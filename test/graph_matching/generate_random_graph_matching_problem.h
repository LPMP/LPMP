#pragma once

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

graph_matching_input generate_random_graph_matching_problem(const std::size_t no_nodes, const double density, std::random_device& rd)
{
    std::mt19937 gen{rd()};
    std::normal_distribution nd(-1.0, 2.0);
    std::uniform_real_distribution ud(0.0, 1.0);

    graph_matching_input instance(generate_random_graph_matching_problem(no_nodes, rd));

    for(std::size_t a1=0; a1<instance.assignments.size(); ++a1) {
        // TODO: test a2=0; ...
        for(std::size_t a2=a1+1; a2<instance.assignments.size(); ++a2) {
            const std::size_t i1 = instance.assignments[a1].left_node;
            const std::size_t i2 = instance.assignments[a2].left_node;
            const std::size_t j1 = instance.assignments[a1].right_node;
            const std::size_t j2 = instance.assignments[a2].right_node;
            if(i1 == i2 || j1 == j2)
                continue;
            if(ud(gen) <= density) {
                instance.quadratic_terms.push_back({a1,a2,nd(gen)});
            } 
        }
    }

    return instance;
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

