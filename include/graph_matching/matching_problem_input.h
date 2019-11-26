#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <cassert>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include "union_find.hxx"
#include <iostream>

// TODO: 
//      - split into header and cpp file and rename to graph_matching_instance.{h|cpp}
//      - rename input to instance
//      - make inner labeling class to standalone class

namespace LPMP {

struct linear_assignment_problem_input {

    constexpr static std::size_t no_assignment = std::numeric_limits<std::size_t>::max();

   void add_assignment(const std::size_t left_node, const std::size_t right_node, const double cost)
   {
       assert(left_node != no_assignment || right_node != no_assignment);
       assignments.push_back({left_node, right_node, cost}); 
       if(left_node != no_assignment)
           no_left_nodes = std::max(no_left_nodes, left_node+1);
       if(right_node != no_assignment)
           no_right_nodes = std::max(no_right_nodes, right_node+1);
   }

   struct assignment {
      std::size_t left_node;
      std::size_t right_node;
      double cost;
   };
   std::vector<assignment> assignments;
   double constant_ = 0.0;

   std::size_t no_left_nodes = 0;
   std::size_t no_right_nodes = 0; 

   std::size_t no_mcf_nodes() const { return no_left_nodes + no_right_nodes + 2; } 
   std::size_t no_mcf_edges() const { return assignments.size() + no_left_nodes + no_right_nodes + 1; }

   void add_to_constant(const double x) { constant_ += x; }
   double get_constant() const { return constant_; }

   template<typename MCF>
   void initialize_mcf(MCF& mcf, const double scaling = 1.0) const
   {
       const std::size_t no_nodes = no_left_nodes + no_right_nodes;
       const std::size_t no_mcf_nodes = no_left_nodes + no_right_nodes + 2;
       const std::size_t no_mcf_edges = assignments.size() + no_left_nodes + no_right_nodes + 1;

       mcf = MCF(no_mcf_nodes, no_mcf_edges);

       std::vector<std::size_t> non_assignment_left(no_left_nodes);
       for(std::size_t i=0; i<no_left_nodes; ++i) {
           non_assignment_left[i] = mcf.add_edge(i, no_nodes+1, 0, 1, 0.0); // for non-assignment
           mcf.add_node_excess(i, 1);
       }

       std::vector<std::size_t> non_assignment_right(no_right_nodes);
       for(std::size_t i=0; i<no_right_nodes; ++i) {
           non_assignment_right[i] = mcf.add_edge(no_nodes, no_left_nodes + i, 0, 1, 0.0); // for non-assignment
           mcf.add_node_excess(no_left_nodes + i, -1);
       }

       mcf.add_edge(no_nodes, no_nodes + 1, 0, no_nodes, 0.0);
       mcf.add_node_excess(no_nodes, no_right_nodes);
       mcf.add_node_excess(no_nodes+1, -no_left_nodes);

       for(const auto a : assignments) {
           if(a.left_node == no_assignment) {
               mcf.update_cost(non_assignment_right[a.right_node], scaling*a.cost);
           } else if(a.right_node == no_assignment) {
               mcf.update_cost(non_assignment_left[a.left_node], scaling*a.cost);
           } else {
               mcf.add_edge(a.left_node, no_left_nodes + a.right_node, 0, 1, scaling*a.cost);
           }
       } 

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

      template<typename MATRIX>
      labeling(const MATRIX& m)
      {
          assert(false); // shall last entry in each row/column denote non-assignment?
         // check that each row and columns has at most one 1 entry
         for(std::size_t i=0; i<m.cols(); ++i) {
            for(std::size_t j=0; j<m.rows(); ++j) {
               assert(m(i,j) == 0 || m(i,j) == 1);
            }
         }
         for(std::size_t i=0; i<m.cols(); ++i) {
            assert(m.col(i).sum() == 1);
         }
         for(std::size_t j=0; j<m.rows(); ++j) {
            assert(m.row(j).sum() == 1);
         }
         this->resize(m.cols(), std::numeric_limits<std::size_t>::max());
         for(std::size_t i=0; i<m.cols(); ++i) {
            for(std::size_t j=0; j<m.rows(); ++j) {
               if(m(i,j) == 1) {
                  assert((*this)[i] == std::numeric_limits<std::size_t>::max());
                  (*this)[i] = j;
               } 
            }
         }
         assert(this->size() == m.cols());
      }

      template<typename STREAM>
      void WritePrimal(STREAM& s) const
      {
         for(std::size_t i=0; i<this->size(); ++i) {
            if((*this)[i] != no_assignment)
               s << i << " -> " << (*this)[i] << "\n";
            else
               s << i << " not matched\n";
         } 
      }

      bool check_primal_consistency() const
      {
         std::vector<char> labels_taken(this->highest_matched_node()+1, 0);
         for(const std::size_t l : *this) {
            if(l == no_assignment)
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
            if((*this)[i] != no_assignment)
               max_label = std::max(max_label, (*this)[i]);
         return max_label;
      }

   };


   double evaluate(const labeling& l) const
   {
      if(l.no_left_nodes() != no_left_nodes)
         throw std::runtime_error("labeling must have equal number of left nodes as matching problem.");

      double cost = constant_;

      std::vector<char> left_node_assigned(no_left_nodes,0);
      std::vector<char> right_node_assigned(no_right_nodes,0);
      for(std::size_t i=0; i<l.size(); ++i) {
          if(l[i] != no_assignment) {
              assert(l[i] < no_right_nodes);
              if(left_node_assigned[i] == 1)
                  return std::numeric_limits<double>::infinity();
              left_node_assigned[i] = 1;

              if(right_node_assigned[l[i]] == 1)
                  return std::numeric_limits<double>::infinity();
              right_node_assigned[l[i]] = 1;
          }
      }

      for(const auto& a : assignments) {
          if(a.left_node != no_assignment && a.right_node != no_assignment) {
              if(l[a.left_node] == a.right_node)
                  cost += a.cost; 
          } else if(a.left_node == no_assignment) {
              assert(a.right_node != no_assignment);
              if(right_node_assigned[a.right_node] == 0)
                  cost += a.cost; 
          } else if(a.right_node == no_assignment) {
              assert(a.left_node != no_assignment);
              if(left_node_assigned[a.left_node] == 0)
                  cost += a.cost; 
          } else {
              assert(false);
          }
      }

      return cost;
   }

};

struct graph_matching_input : public linear_assignment_problem_input {

    graph_matching_input() {}

    graph_matching_input(const linear_assignment_problem_input& instance)
        : linear_assignment_problem_input(instance) {}

    using labeling = linear_assignment_problem_input::labeling;

    struct quadratic {
        std::size_t assignment_1, assignment_2;
        double cost;
        bool operator<(const quadratic& o) const 
        {
            if(assignment_1 != o.assignment_1) return assignment_1 < o.assignment_1;
            else return assignment_2 < o.assignment_2;
        } 
    };
    std::vector<quadratic> quadratic_terms;

    void add_quadratic_term(const std::size_t assignment_1, const std::size_t assignment_2, const double cost)
    {
        assert(assignment_1 != assignment_2);
        assert(std::max(assignment_1, assignment_2) < this->assignments.size());
        const std::size_t a1 = std::min(assignment_1, assignment_2);
        const std::size_t a2 = std::max(assignment_1, assignment_2);
        quadratic_terms.push_back({a1, a2, cost});
    }

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

    void normalize_quadratic_terms()
    {
        std::sort(quadratic_terms.begin(), quadratic_terms.end());

        std::vector<quadratic> normalized_quadratic_terms; 

        // merge matching copies and add them to normalized vector
        for(std::size_t k=0; k<quadratic_terms.size();) {
            normalized_quadratic_terms.push_back(quadratic_terms[k]);
            ++k; 
            while(k<quadratic_terms.size() && normalized_quadratic_terms.back().assignment_1 == quadratic_terms[k].assignment_1 && normalized_quadratic_terms.back().assignment_2 == quadratic_terms[k].assignment_2) {
                normalized_quadratic_terms.back().cost += quadratic_terms[k].cost;
                ++k; 
            }
        }

        std::swap(normalized_quadratic_terms, quadratic_terms); 
    }

   template<typename STREAM>
       void write(STREAM& s) const
       {
           // init line
           s << "p " << no_left_nodes << " " << no_right_nodes << " " << this->assignments.size() << " " << quadratic_terms.size() << "\n";
 
           // linear assignment costs
           for(std::size_t i=0; i<this->assignments.size(); ++i) {
               const auto& a = this->assignments[i];
               s << "a " << i << " " << a.left_node << " " << a.right_node << " " << a.cost << "\n";
           }

           // quadratic terms
           for(const auto& q : quadratic_terms) {
               s << "e " << q.assignment_1 << " " << q.assignment_2 << " " << q.cost << "\n";
           }
       }
};

struct multigraph_matching_input_entry {
    std::size_t left_graph_no, right_graph_no;
    graph_matching_input gm_input;
}; 

struct multigraph_matching_input : public std::vector<multigraph_matching_input_entry> {

   class graph_size {
      public:
         template<typename SIZE_ITERATOR>
         graph_size(SIZE_ITERATOR size_begin, SIZE_ITERATOR size_end)
         : no_nodes_(size_begin, size_end)
         {
            compute_node_offsets(); 
         }
         graph_size(const multigraph_matching_input& input)
         {
            compute_no_nodes(input);
            compute_node_offsets();
         }
         graph_size() {}

         std::size_t no_graphs() const { return no_nodes_.size(); }
         std::size_t total_no_nodes() const { return graph_node_offsets_.back() + no_nodes_.back(); }
         std::size_t no_nodes(const std::size_t i) const { assert(i < no_graphs()); return no_nodes_[i]; }

         auto& no_nodes() { return no_nodes_; }
         const auto& no_nodes() const { return no_nodes_; }

         std::size_t node_no(const std::size_t graph_no, const std::size_t node_no) const
         {
            assert(graph_node_offsets_.size() == no_graphs());
            assert(graph_no < graph_node_offsets_.size());
            return graph_node_offsets_[graph_no] + node_no; 
         }

         std::array<std::size_t,2> graph_node_no(const std::size_t i) const
         {
            assert(graph_node_offsets_.size() == no_graphs());
            const std::size_t p = std::lower_bound(graph_node_offsets_.begin(), graph_node_offsets_.end(), i+1)-1 - graph_node_offsets_.begin();
            assert(p < no_graphs());
            const std::size_t p_node = i - graph_node_offsets_[p];
            assert(p_node < no_nodes(p));
            assert(i == node_no(p, p_node));
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

         void compute_node_offsets()
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


   struct multigraph_matching_labeling_entry {
      std::size_t left_graph_no, right_graph_no;
      linear_assignment_problem_input::labeling labeling;
   };

   class labeling : public std::vector<multigraph_matching_labeling_entry> {
      public:
      using std::vector<multigraph_matching_labeling_entry>::vector;

      template<typename MATRIX>
      labeling(const MATRIX& m, const graph_size& mgm_size)
      {
          assert(false); // add treatment of non-assignment. Which matrix format should we take?
         assert(m.cols() == m.rows());
         //for(std::size_t i=0; i<m.cols(); ++i) {
         //   assert(m(i,i) == 1);
         //}
         for(std::size_t p=0; p<mgm_size.no_graphs(); ++p) {
            for(std::size_t q=p+1; q<mgm_size.no_graphs(); ++q) {
               // extract submatrix
               const std::size_t first_col = mgm_size.node_no(q,0);
               const std::size_t no_cols = mgm_size.no_nodes(q);
               const std::size_t first_row = mgm_size.node_no(p,0);
               const std::size_t no_rows = mgm_size.no_nodes(p);

               std::cout << first_row << "," << first_col << "," << no_rows << "," << no_cols << "\n";
               auto gm_matching = m.block(first_row, first_col, no_rows, no_cols);
               //auto gm_matching = m.block(first_col, first_row, no_cols, no_rows);
               linear_assignment_problem_input::labeling gm_labeling(gm_matching);
               this->push_back({p,q, gm_labeling});
            }
         }
      }

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

         s << "# ${graph_no}, ${node_no}, ${cluster_id}\n";
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
   // TODO: remove class, it is implemented already generally
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
