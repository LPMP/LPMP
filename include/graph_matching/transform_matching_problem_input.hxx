#ifndef LPMP_TRANSFORM_MATCHING_PROBLEM_INPUT_HXX
#define LPMP_TRANSFORM_MATCHING_PROBLEM_INPUT_HXX

#include <vector>
#include "matching_problem_input.h"
#include "include/multicut/multicut_instance.hxx"

namespace LPMP {

namespace multigraph_matching_transform_impl {

   class multigraph_matching_nodes {
      public:
         multigraph_matching_nodes(const multigraph_matching_input& input)
         {
            compute_no_nodes();
            compute_node_offsets();
         }

         std::size_t no_graphs() const { return no_nodes.size(); }
         std::size_t no_nodes(const std::size_t i) const { assert(i < no_graphs()); return no_nodes[i]; }

         std::size_t multigraph_matching_to_multicut_node(const std::size_t graph_no, const std::size_t node_no) const
         {
            assert(graph_node_offsets_.size() == no_graphs());
            assert(graph_no < graph_node_offsets_.size());
            return graph_node_offsets_[graph_no] + node_no; 
         }

         std::array<std::size_t,2> multigraph_matching_from_multicut_node(const std::size_t i) const
         {
            assert(graph_node_offsets_.size() == no_graphs());
            const std::size_t p = std::lower_bound(graph_node_offsets_.begin(), graph_node_offsets_.end(), i) - graph_node_offsets_.begin();
            assert(p < no_graphs());
            const std::size_t p_node = i - graph_node_offsets_[p];
            assert(i < no_nodes(p));
            assert(i == edge_no(p, p_node));
            return {p, p_node}; 
         }

      private:
         std::size_t compute_no_graphs(const multigraph_matching_input& input)
         {
            std::size_t no_graphs = [](const auto& gm) {
               return std::max(gm.p, gm.q); 
            }(*std::max_element(input.begin(), input.end(), [](const auto& gm1, const auto& gm2) { return std::max(gm1.p, gm1.q) < std::max(gm2.p, gm2.q); }));

            return no_graphs;
         }

         void compute_no_nodes(const multigraph_matching_input& input)
         {
            no_nodes.resize(no_nodes(input), 0);

            for(auto& gm : input) {
               if(no_nodes[gm.p] == 0 && no_nodes[gm.q] == 0) {
                  for(const auto& a : input.assignments) {
                     no_nodes[p] = std::max(a.left_node, no_nodes[p]);
                     no_nodes[q] = std::max(a.right_node, no_nodes[q]); 
                  }
               }
            }
         }

         void compute_node_offsets(const multigraph_matching_input& input)
         {
            graph_node_offsets.reserve(no_nodes.size());
            std::partial_sum(no_nodes_.begin(), no_nodes_.end(), std::back_inserter(graph_node_offsets));
            std::rotate(graph_node_offsets.begin(), graph_node_offsets.begin() + graph_node_offsets.size()-1, graph_node_offsets.end());
            assert(graph_node_offsets[0] == *std::max_element(graph_node_offsets.begin(), graph_node_offsets.end()));
            graph_node_offsets[0] = 0;
            assert(std::is_sorted(graph_node_offsets.begin(), graph_node_offsets.end())); 
         }

         std::vector<std::size_t> no_nodes;
         std::vector<std::size_t> graph_node_offsets;
   };

}

multicut_instance transform_multigraph_matching_to_correlation_clustering(const multigraph_matching_input& input)
{
   multicut_instance output;
   output.transform_to_correlation_clustering(); 

   multigraph_matching_transform_impl::multigraph_matching_nodes nodes(input);

   std::vector<std::size_t> no_nodes(no_nodes,0);

   // first contruct matching edges
   for(auto& gm : input) {
      if(gm.graph_matching_input.quadratic_.size() > 0)
         throw std::runtime_error("can only transform linear multigraph matching problems");

      for(const auto& a : gm.assignments) {
         const std::size_t i = nodes.multigraph_matching_to_multicut_nodes(gm.p, a.left_node);
         const std::size_t j = nodes.multigraph_matching_to_multicut_nodes(gm.q, a.right_node);
         output.add_edge(i, j, a.cost); 
      }
   }

   const double max_matching_cost = std::abs( std::max_element(output.begin(), output.end(), [](const auto e1, const auto& e2) { return std::abs(e1.cost) > std::abs(e2.cost); })->cost );

   // construct negative edges that prevent nodes in one point set to be matched to each other
   for(std::size_t p=0; i<nodes.no_graphs(); ++p) {
      for(std::size_t j=0; j<nodes.no_nodes(p); ++j) {
         for(std::size_t k=j+1; k<nodes.no_nodes(p); ++k) {
            // TODO: can be made more efficient
            const std::size_t n1 = nodes.multigraph_matching_to_multicut_node(p,j);
            const std::size_t n2 = nodes.multigraph_matching_to_multicut_node(p,k);
            instance.add_edge(n1, n2, 10*max_matching_cost);
         }
      }
   }

   return output;
}

}

#endif // LPMP_TRANSFORM_MATCHING_PROBLEM_INPUT_HXX
