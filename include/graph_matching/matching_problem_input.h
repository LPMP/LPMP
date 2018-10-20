#ifndef LPMP_MATCHING_PROBLEM_INPUT_H
#define LPMP_MATCHING_PROBLEM_INPUT_H

#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <cassert>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include "union_find.hxx"

// TODO: split into header and cpp file

namespace LPMP {

struct linear_assignment_problem_input {

   void add_assignment(const std::size_t left_node, const std::size_t right_node, const double cost)
   {
     assignments.push_back({left_node, right_node, cost}); 
     no_left_nodes = std::max(no_left_nodes, left_node+1);
     no_right_nodes = std::max(no_right_nodes, right_node+1);
   }

   struct assignment {
      std::size_t left_node;
      std::size_t right_node;
      double cost;
   };
   std::vector<assignment> assignments;

   std::size_t no_left_nodes = 0;
   std::size_t no_right_nodes = 0; 

   std::size_t no_mcf_nodes() const { return no_left_nodes + no_right_nodes + 2; } 
   std::size_t no_mcf_edges() const { return assignments.size() + no_left_nodes + no_right_nodes + 1; }

   template<typename MCF>
   void initialize_mcf(MCF& mcf, const double scaling = 1.0) const
   {
      const std::size_t no_nodes = no_left_nodes + no_right_nodes;

      for(const auto a : assignments)
         mcf.add_edge(a.left_node, no_left_nodes + a.right_node, 0, 1, scaling*a.cost);

      for(std::size_t i=0; i<no_left_nodes; ++i) {
         mcf.add_edge(i, no_nodes+1, 0, 1, 0.0); // for non-assignment
         mcf.add_node_excess(i, 1);
      }

      for(std::size_t i=0; i<no_right_nodes; ++i) {
         mcf.add_edge(no_nodes, no_left_nodes + i, 0, 1, 0.0); // for non-assignment
         mcf.add_node_excess(no_left_nodes + i, -1);
      }

      mcf.add_edge(no_nodes, no_nodes + 1, 0, no_nodes, 0.0);
      mcf.add_node_excess(no_nodes, no_right_nodes);
      mcf.add_node_excess(no_nodes+1, -no_left_nodes);

      mcf.order();
      mcf.solve();
   }

   void normalize()
   {
      // merge parallel edges
      std::sort(assignments.begin(), assignments.end(), [](const auto& a1, const auto& a2) { 
         std::array<std::size_t,2> e1 {a1.left_node, a1.right_node};
         std::array<std::size_t,2> e2 {a2.left_node, a2.right_node};
         return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end());
         });

      std::vector<assignment> normalized_assignments; 

      // merge matching edge copies and add them to edges
      for(std::size_t k=0; k<assignments.size();) {
         normalized_assignments.push_back(assignments[k]);
         ++k; 
         while(k<assignments.size() && normalized_assignments.back().left_node == assignments[k].left_node && normalized_assignments.back().right_node == assignments[k].right_node) {
            normalized_assignments.back().cost += assignments[k].cost;
            ++k; 
         }
      }

      std::swap(normalized_assignments, assignments); 
   }

   class labeling : public std::vector<std::size_t> {
      public:
      using std::vector<std::size_t>::vector;

      template<typename STREAM>
      void WritePrimal(STREAM& s) const
      {
         for(std::size_t i=0; i<this->size(); ++i) {
            if((*this)[i] != std::numeric_limits<std::size_t>::max())
               s << i << " -> " << (*this)[i] << "\n";
            else
               s << i << " not matched\n";
         } 
      }

      bool check_primal_consistency() const
      {
         std::vector<char> labels_taken(this->highest_matched_node()+1, 0);
         for(const std::size_t l : *this) {
            if(l == std::numeric_limits<std::size_t>::max())
               continue;
            if(labels_taken[l] > 0)
               return false;
            labels_taken[l] = 1;
         }
         return true;
      }

      std::size_t no_left_nodes() const { return this->size(); }

      std::size_t highest_matched_node() const
      {
         std::size_t max_label = 0;
         for(std::size_t i=0; i<this->size(); ++i)
            if((*this)[i] != std::numeric_limits<std::size_t>::max())
               max_label = std::max(max_label, (*this)[i]);
         return max_label;
      }

   };


   double evaluate(const labeling& l) const
   {
      if(!l.no_left_nodes() == no_left_nodes)
         throw std::runtime_error("labeling must have equal number of left nodes as matching problem.");

      double cost = 0.0;

      std::vector<std::size_t> nodes_taken(no_right_nodes,0);
      for(const auto& a : assignments) {
         if(l[a.left_node] == a.right_node) {
            cost += a.cost;
            assert(a.left_node < nodes_taken.size());
            nodes_taken[a.left_node]++;
         }
      }

      if(*std::max_element(nodes_taken.begin(), nodes_taken.end()) > 1)
         return std::numeric_limits<double>::infinity();

      return cost;
   }

};

struct graph_matching_input : public linear_assignment_problem_input {
   using labeling = linear_assignment_problem_input::labeling;
    struct quadratic {
        std::size_t assignment_1, assignment_2;
        double cost;
    };
    std::vector<quadratic> quadratic_terms;

    double evaluate(const labeling& l) const
    {
       double cost = linear_assignment_problem_input::evaluate(l);
       for(const auto& q : quadratic_terms) {
          const auto a1 = this->assignments[q.assignment_1];
          const auto a2 = this->assignments[q.assignment_2];
          if(l[a1.left_node] == a1.right_node && l[a2.left_node] == a2.right_node)
             cost += q.cost;
       }
       return cost;
    }
};

struct multigraph_matching_input_entry {
    std::size_t left_graph_no, right_graph_no;
    graph_matching_input gm_input;

};


struct multigraph_matching_input : public std::vector<multigraph_matching_input_entry> {

   struct multigraph_matching_labeling_entry {
      std::size_t left_graph_no, right_graph_no;
      linear_assignment_problem_input::labeling labeling;
   };

   class labeling : public std::vector<multigraph_matching_labeling_entry> {
      public:
      using std::vector<multigraph_matching_labeling_entry>::vector;

      bool check_primal_consistency() const
      {
         for(const auto& gm : *this)
            if(!gm.labeling.check_primal_consistency())
               return false;
               
         // check whether cycles are consistenct
         node_mapping nn(*this);

         union_find uf(nn.no_nodes());
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

      template<typename STREAM>
      void write_primal_matching(STREAM& s) const
      {
         for(const auto& gm : *this) {
            s << "graph matching " << gm.left_graph_no << " -> " << gm.right_graph_no << "\n";
            gm.labeling.WritePrimal(s);
         }
      }

      // write our primal matching as table with three columns: graph number, node number, cluster identifier
      template<typename STREAM>
      void write_primal_clustering(STREAM& s) const
      {
         node_mapping nm(*this);
         union_find uf(nm.no_nodes());

         for(const auto& gm : *this)
            for(std::size_t i=0; i<gm.labeling.size(); ++i)
               if(gm.labeling[i] != std::numeric_limits<std::size_t>::max())
                  uf.merge( nm.node_no(gm.left_graph_no, i), nm.node_no(gm.right_graph_no, gm.labeling[i]) );

         for(std::size_t i=0; i<nm.no_nodes(); ++i) {
            const auto [graph_no, node_no] = nm.node_no(i);
            s << graph_no << ", " << node_no << ", " << uf.find(i) << "\n";
         }
      }
   };

   double evaluate(const labeling& l) const
   {
      if(l.size() != this->size())
         throw std::runtime_error("multigraph matching labeling must have equal number of matchings as input.");

      if(!l.check_primal_consistency())
         return std::numeric_limits<double>::infinity();

      double cost = 0.0;
      for(std::size_t i=0; i<this->size(); ++i) {
         if(l[i].left_graph_no != (*this)[i].left_graph_no || l[i].right_graph_no != (*this)[i].right_graph_no)
            throw std::runtime_error("labeling acts on different graph matching problems than instance.");

         cost += (*this)[i].gm_input.evaluate(l[i].labeling);
      }

      return cost;
   }

   private:
   class node_mapping {
      public:
      node_mapping(const multigraph_matching_input::labeling& mgm)
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
};

} // namespace LPMP

#endif // LPMP_MATCHING_PROBLEM_INPUT_H
