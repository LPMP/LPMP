#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "multicut/transform_multigraph_matching.h"
#include <iostream> // TODO: remove

namespace LPMP {

namespace detail {

   class multigraph_matching_nodes {
      public:
         multigraph_matching_nodes(const multigraph_matching_input& input)
         {
            compute_no_nodes(input);
            compute_node_offsets(input);
         }

         std::size_t no_graphs() const { return no_nodes_.size(); }
         std::size_t no_nodes(const std::size_t i) const { assert(i < no_graphs()); return no_nodes_[i]; }

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
            assert(i == multigraph_matching_to_multicut_node(p, p_node));
            return {p, p_node}; 
         }

      private:
         std::size_t compute_no_graphs(const multigraph_matching_input& input)
         {
            std::size_t no_graphs = 1 + [](const auto& gm) {
               return std::max(gm.left_graph_no, gm.right_graph_no); 
            }(*std::max_element(input.begin(), input.end(), [](const auto& gm1, const auto& gm2) { return std::max(gm1.left_graph_no, gm1.right_graph_no) < std::max(gm2.left_graph_no, gm2.right_graph_no); }));

            return no_graphs;
         }

         void compute_no_nodes(const multigraph_matching_input& input)
         {
            no_nodes_.resize(compute_no_graphs(input), 0);

            for(auto& gm : input) {
               for(const auto& a : gm.gm_input.assignments) {
                  no_nodes_[gm.left_graph_no] = std::max(a.left_node+1, no_nodes_[gm.left_graph_no]);
                  no_nodes_[gm.right_graph_no] = std::max(a.right_node+1, no_nodes_[gm.right_graph_no]); 
               }
            }
         }

         void compute_node_offsets(const multigraph_matching_input& input)
         {
            graph_node_offsets_.reserve(no_graphs());
            std::partial_sum(no_nodes_.begin(), no_nodes_.end(), std::back_inserter(graph_node_offsets_));
            std::rotate(graph_node_offsets_.begin(), graph_node_offsets_.begin() + graph_node_offsets_.size()-1, graph_node_offsets_.end());
            assert(graph_node_offsets_[0] == *std::max_element(graph_node_offsets_.begin(), graph_node_offsets_.end()));
            graph_node_offsets_[0] = 0;
            assert(std::is_sorted(graph_node_offsets_.begin(), graph_node_offsets_.end())); 
         }

         std::vector<std::size_t> no_nodes_;
         std::vector<std::size_t> graph_node_offsets_;
   };

} // namespace detail

multicut_instance transform_multigraph_matching_to_correlation_clustering(const multigraph_matching_input& input)
{
   detail::multigraph_matching_nodes nodes(input);
   return transform_multigraph_matching_to_correlation_clustering(input, nodes);
}

multicut_instance transform_multigraph_matching_to_correlation_clustering(const multigraph_matching_input& input, const detail::multigraph_matching_nodes& nodes)
{
   multicut_instance output;
   output.transform_to_correlation_clustering(); 

   std::vector<std::size_t> no_nodes(nodes.no_graphs(),0);

   std::vector<std::vector<char>> edge_present_;
   auto edge_present = [&](const std::size_t matching_no, const std::size_t i, const std::size_t j) -> bool {
      assert(matching_no < input.size());
      const auto& gm = input[matching_no];
      return edge_present_[matching_no][i * nodes.no_nodes(gm.right_graph_no) + j] == 1;
   };

   auto add_edge = [&](const std::size_t matching_no, const std::size_t i, const std::size_t j) {
      assert(matching_no < input.size());
      const auto& gm = input[matching_no];
      return edge_present_[matching_no][i * nodes.no_nodes(gm.right_graph_no) + j] = 1;
   };

   auto add_matching = [&](const std::size_t matching_no) {
      assert(matching_no < input.size());
      const auto& gm = input[matching_no];
      assert(matching_no == edge_present_.size());
      return edge_present_.push_back(std::vector<char>(nodes.no_nodes(gm.right_graph_no) * nodes.no_nodes(gm.left_graph_no),0));;
   };

   // first contruct matching edges
   for(std::size_t matching_no=0; matching_no<input.size(); ++matching_no) {
      const auto& gm = input[matching_no];
      if(gm.gm_input.quadratic_terms.size() > 0)
         throw std::runtime_error("can only transform linear multigraph matching problems");

      add_matching(matching_no);
      for(const auto& a : gm.gm_input.assignments) {
         const std::size_t i = nodes.multigraph_matching_to_multicut_node(gm.left_graph_no, a.left_node);
         const std::size_t j = nodes.multigraph_matching_to_multicut_node(gm.right_graph_no, a.right_node);
         output.add_edge(i, j, a.cost); 
         add_edge(matching_no, a.left_node, a.right_node);
      }
   }

   const double max_matching_cost = std::abs( std::max_element(output.begin(), output.end(), [](const auto e1, const auto& e2) { return std::abs(e1.cost) < std::abs(e2.cost); })->cost );
   const double non_matching_penalty = std::max(1e10, 1000*max_matching_cost);
   for(const auto& a : output.edges()) { assert(max_matching_cost >= std::abs(a.cost)); }

   // construct negative edges that prevent nodes in one point set to be matched to each other
   for(std::size_t p=0; p<nodes.no_graphs(); ++p) {
      for(std::size_t j=0; j<nodes.no_nodes(p); ++j) {
         for(std::size_t k=j+1; k<nodes.no_nodes(p); ++k) {
            // TODO: can be made more efficient
            const std::size_t n1 = nodes.multigraph_matching_to_multicut_node(p,j);
            const std::size_t n2 = nodes.multigraph_matching_to_multicut_node(p,k);
            output.add_edge(n1, n2, non_matching_penalty);
         }
      }
   }

   // if graph matching problem is not dense, add edges that disallow matchings for edges that are not present
   for(std::size_t matching_no=0; matching_no<input.size(); ++matching_no) {
      const auto& gm = input[matching_no];
      for(std::size_t i=0; i<nodes.no_nodes(gm.left_graph_no); ++i) {
         for(std::size_t j=0; j<nodes.no_nodes(gm.right_graph_no); ++j) {
            if(!edge_present(matching_no, i, j)) {
               const std::size_t n1 = nodes.multigraph_matching_to_multicut_node(gm.left_graph_no,i);
               const std::size_t n2 = nodes.multigraph_matching_to_multicut_node(gm.right_graph_no,j);
               output.add_edge(n1, n2, non_matching_penalty);
            }
         }
      }
   }

   return output;
}

multigraph_matching_input::labeling transform_correlation_clustering_to_multigraph_matching(
      const multigraph_matching_input& mgm,
      const multicut_instance& cc,
      multicut_instance::edge_labeling l_input)
{
   detail::multigraph_matching_nodes nodes(mgm);
   return transform_correlation_clustering_to_multigraph_matching(mgm, cc, l_input, nodes);
}

multigraph_matching_input::labeling transform_correlation_clustering_to_multigraph_matching(
      const multigraph_matching_input& mgm,
      const multicut_instance& cc,
      multicut_instance::edge_labeling l_input,
      const detail::multigraph_matching_nodes& nodes)
{
   l_input.transform_to_correlation_clustering();
   multigraph_matching_input::labeling l_output;
   l_output.reserve(mgm.size());

   // TODO: transform cc edge labeling to node labeling and then transform that to multigraph matching labeling
   std::size_t e = 0;
   for(auto& gm : mgm) {
      linear_assignment_problem_input::labeling l(nodes.no_nodes(gm.left_graph_no), std::numeric_limits<std::size_t>::max());
      l_output.push_back( {gm.left_graph_no, gm.right_graph_no, l} );
      for(const auto& a : gm.gm_input.assignments) {
         const std::size_t i = nodes.multigraph_matching_to_multicut_node(gm.left_graph_no, a.left_node);
         const std::size_t j = nodes.multigraph_matching_to_multicut_node(gm.right_graph_no, a.right_node);
         if(l_input[e] == 1) {
            l_output.back().labeling[a.left_node] = a.right_node;
            //std::cout << e << "\n";
         }
         ++e;
      }
   }

   // following edges should disallow infeasible pairwise matchings: they must all be zero.
   for(; e<l_input.size(); ++e)
      if(l_input[e] == 1)
         throw std::runtime_error("correlation clustering not corresponding to multigraph matching.");

   assert(e == l_input.size());
   assert(l_output.check_primal_consistency());

   return l_output;
}

/*
multicut_instance::labeling transform_multigraph_matching_labeling_to_correlation_clustering(
      const multigraph_matching_input::labeling& input
      const multigraph_matching_input& mgm,
      const multicut_instance& cc,
      const detail::multigraph_matching_nodes& nodes)
{
   multicut_instance::labeling output;
   // first add the matching edges
   for(const auto& e : cc.edges()) {

   }

   // now add intra-graph edges. Possibly input is not feasible. In this case not all edges are 0.
}
*/

} // namespace LPMP
