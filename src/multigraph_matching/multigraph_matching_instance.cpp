#include  "multigraph_matching/multigraph_matching_instance.h"
#include "union_find.hxx"
#include <unordered_set>

namespace LPMP {

         graph_size::graph_size(const multigraph_matching_instance& instance)
         {
            compute_no_nodes(instance);
            compute_node_offsets();
         }

         graph_size::graph_size(const multigraph_matching_labeling& labeling)
         {
            compute_no_nodes(labeling);
            compute_node_offsets();
         }

         std::size_t graph_size::no_graphs() const { return no_nodes_.size(); }
         std::size_t graph_size::total_no_nodes() const { return graph_node_offsets_.back() + no_nodes_.back(); }
         std::size_t graph_size::no_nodes(const std::size_t i) const { assert(i < no_graphs()); return no_nodes_[i]; }

         std::vector<std::size_t>& graph_size::no_nodes() { return no_nodes_; }
         const std::vector<std::size_t>& graph_size::no_nodes() const { return no_nodes_; }

         std::size_t graph_size::node_no(const std::size_t graph_no, const std::size_t node_no) const
         {
            assert(graph_node_offsets_.size() == no_graphs());
            assert(graph_no < graph_node_offsets_.size());
            return graph_node_offsets_[graph_no] + node_no; 
         }

         std::array<std::size_t,2> graph_size::graph_node_no(const std::size_t i) const
         {
            assert(graph_node_offsets_.size() == no_graphs());
            const std::size_t p = std::lower_bound(graph_node_offsets_.begin(), graph_node_offsets_.end(), i+1)-1 - graph_node_offsets_.begin();
            assert(p < no_graphs());
            const std::size_t p_node = i - graph_node_offsets_[p];
            assert(p_node < no_nodes(p));
            assert(i == node_no(p, p_node));
            return {p, p_node}; 
         }

         std::size_t graph_size::compute_no_graphs(const multigraph_matching_instance& instance)
         {
            std::size_t no_graphs = 1 + [](const auto& gm) {
               return std::max(gm.left_graph_no, gm.right_graph_no); 
            }(*std::max_element(instance.begin(), instance.end(), [](const auto& gm1, const auto& gm2) { return std::max(gm1.left_graph_no, gm1.right_graph_no) < std::max(gm2.left_graph_no, gm2.right_graph_no); }));

            return no_graphs;
         }

         void graph_size::compute_no_nodes(const multigraph_matching_labeling& labeling)
         {
             const std::size_t no_graphs = 1 + std::max_element(labeling.begin(), labeling.end(), [](const auto& gm1, const auto& gm2) { return gm1.right_graph_no < gm2.right_graph_no; })->right_graph_no;

             no_nodes_.resize(no_graphs, 0);
             for(const auto& gm : labeling) {
                 no_nodes_[gm.left_graph_no] = std::max(no_nodes_[gm.left_graph_no], gm.labeling.no_left_nodes());
                 no_nodes_[gm.right_graph_no] = std::max(no_nodes_[gm.right_graph_no], gm.labeling.highest_matched_node()+1);
             }

             compute_node_offsets();
         }

         void graph_size::compute_no_nodes(const multigraph_matching_instance& instance)
         {
            no_nodes_.resize(compute_no_graphs(instance), 0);

            for(auto& gm : instance) {
               for(const auto& a : gm.gm_instance.assignments()) {
                  no_nodes_[gm.left_graph_no] = std::max(a.left_node+1, no_nodes_[gm.left_graph_no]);
                  no_nodes_[gm.right_graph_no] = std::max(a.right_node+1, no_nodes_[gm.right_graph_no]); 
               }
            }
         }

         void graph_size::compute_node_offsets()
         {
            graph_node_offsets_.reserve(no_graphs());
            std::partial_sum(no_nodes_.begin(), no_nodes_.end(), std::back_inserter(graph_node_offsets_));
            std::rotate(graph_node_offsets_.begin(), graph_node_offsets_.begin() + graph_node_offsets_.size()-1, graph_node_offsets_.end());
            assert(graph_node_offsets_[0] == *std::max_element(graph_node_offsets_.begin(), graph_node_offsets_.end()));
            graph_node_offsets_[0] = 0;
            assert(std::is_sorted(graph_node_offsets_.begin(), graph_node_offsets_.end())); 
         }


   double multigraph_matching_instance::evaluate(const multigraph_matching_labeling& l) const
   {
      if(l.size() != this->size())
         throw std::runtime_error("multigraph matching labeling must have equal number of matchings as instance.");

      if(!l.check_primal_consistency())
         return std::numeric_limits<double>::infinity();

      double cost = 0.0;
      for(std::size_t i=0; i<this->size(); ++i) {
         if(l[i].left_graph_no != (*this)[i].left_graph_no || l[i].right_graph_no != (*this)[i].right_graph_no)
            throw std::runtime_error("labeling acts on different graph matching problems than instance.");

         cost += (*this)[i].gm_instance.evaluate(l[i].labeling);
      }

      return cost;
   }

   /*
   class node_mapping {
      public:
      node_mapping(const multigraph_matching_instance::labeling& mgm)
      {
         const std::size_t no_graphs = 1 + std::max_element(mgm.begin(), mgm.end(), [](const auto& gm1, const auto& gm2) { return gm1.right_graph_no < gm2.right_graph_no; })->right_graph_no;

         no_nodes_.resize(no_graphs, 0);
         for(const auto& gm : mgm) {
            no_nodes_[gm.left_graph_no] = std::max(no_nodes_[gm.left_graph_no], gm.labeling.no_left_nodes());
            no_nodes_[gm.right_graph_no] = std::max(no_nodes_[gm.right_graph_no], gm.labeling.highest_matched_node()+1);
         }

         sum_no_nodes_ = std::accumulate(no_nodes_.begin(), no_nodes_.end(), 0);

         graph_node_offsets_.resize(no_graphs);
         graph_node_offsets_[0] = 0;
         std::partial_sum(no_nodes_.begin(), no_nodes_.end()-1, graph_node_offsets_.begin()+1); 
      }

      std::size_t no_graphs() const { return no_nodes_.size(); }

      std::size_t no_nodes() const { return sum_no_nodes_; }

      std::size_t no_nodes(const std::size_t i) const
      {
         assert(i<no_graphs());
         return no_nodes_[i]; 
      }

      std::size_t node_no(const std::size_t graph_no, const std::size_t node) const
      {
         assert(graph_no < graph_node_offsets_.size());
         assert(node < no_nodes_[graph_no]);
         const std::size_t n = graph_node_offsets_[graph_no] + node;
         assert(n < std::accumulate(no_nodes_.begin(), no_nodes_.end(), 0));
         return n; 
      }

      std::array<std::size_t,2> node_no(const std::size_t i) const
      {
         assert(i < no_nodes());
         const std::size_t graph_no = std::upper_bound(graph_node_offsets_.begin(), graph_node_offsets_.end(), i)-1 - graph_node_offsets_.begin();
         assert(graph_no < graph_node_offsets_.size());
         assert(graph_node_offsets_[graph_no] <= i);
         const std::size_t node_no = i - graph_node_offsets_[graph_no];
         assert(node_no < no_nodes(graph_no));
         return {graph_no, node_no};
      }

      private:
      std::size_t sum_no_nodes_ = 0;
      std::vector<std::size_t> no_nodes_; 
      std::vector<std::size_t> graph_node_offsets_;
   };
   */

      bool multigraph_matching_labeling::check_primal_consistency() const
      {
         for(const auto& gm : *this)
            if(!gm.labeling.check_primal_consistency())
               return false;
               
         // check whether cycles are consistenct
         graph_size nn(*this);

         union_find uf(nn.no_nodes().size());
         for(const auto& gm : *this)
            for(std::size_t i=0; i<gm.labeling.size(); ++i)
               if(gm.labeling[i] != std::numeric_limits<std::size_t>::max())
                  uf.merge( nn.node_no(gm.left_graph_no, i), nn.node_no(gm.right_graph_no, gm.labeling[i]) );

         // check whether for all nodes of the same graph the identifiers of the union find datastructure are distinct
         for(std::size_t i=0; i<nn.no_graphs(); ++i) {
            std::unordered_set<std::size_t> identifiers;
            for(std::size_t l=0; l<nn.no_nodes(i); ++l)
               identifiers.insert( uf.find(nn.node_no(i,l)) );
            if(identifiers.size() < nn.no_nodes(i))
               return false;
         }

         // TODO: this can only be done when we have access to edges of the underlying pairwise graph matching problems
         // check whether edges that are not present in the pairwise graph matching are disconnected 

         return true;
      }


/*
        double multigraph_matching_instance::evaluate(const multigraph_matching_labeling& l) const
        {
            if(l.size() != this->size())
                throw std::runtime_error("multigraph matching labeling must have equal number of matchings as instance.");

            if(!l.check_primal_consistency())
                return std::numeric_limits<double>::infinity();

            double cost = 0.0;
            for(std::size_t i=0; i<this->size(); ++i) {
                if(l[i].left_graph_no != (*this)[i].left_graph_no || l[i].right_graph_no != (*this)[i].right_graph_no)
                    throw std::runtime_error("labeling acts on different graph matching problems than instance.");

                cost += (*this)[i].gm_instance.evaluate(l[i].labeling);
            }

            return cost;
        }
        */

} // namespace LPMP
