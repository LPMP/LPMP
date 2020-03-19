#pragma once

#include <vector>
#include <Eigen/Eigen>
#include <array>
#include <algorithm>
#include <limits>
#include <cassert>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include "union_find.hxx"
#include "hash_helper.hxx"
#include <iostream> // TODO: remove

// TODO: split into header and cpp file

namespace LPMP {

struct linear_assignment_problem_input {

    constexpr static std::size_t no_assignment = std::numeric_limits<std::size_t>::max();

   void add_assignment(const std::size_t left_node, const std::size_t right_node, const double cost)
   {
      // TODO: remove, only test
      //if(left_node == no_assignment || right_node == no_assignment)
      //   assignments.push_back({left_node, right_node, 0.0});
      //else
         assignments.push_back({left_node, right_node, cost});
      no_left_nodes = std::max(no_left_nodes, left_node + 1);
      no_right_nodes = std::max(no_right_nodes, right_node + 1);
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
      const std::size_t no_mcf_edges = assignments.size() + no_left_nodes + no_right_nodes + 1; // TODO: filter out the assignment that have non-assignment in one of their endpoints

       mcf = MCF(no_mcf_nodes, no_mcf_edges);

       std::vector<double> left_non_assignment_costs(no_left_nodes, 0.0);
       std::vector<double> right_non_assignment_costs(no_right_nodes, 0.0);
       for (const auto a : assignments)
       {
          if (a.left_node != no_assignment && a.right_node != no_assignment)
          {
             mcf.add_edge(a.left_node, no_left_nodes + a.right_node, 0, 1, scaling * a.cost);
          }
          else if (a.left_node == no_assignment)
          {
             assert(a.right_node < right_non_assignment_costs.size());
             right_non_assignment_costs[a.right_node] += scaling * a.cost;
          }
          else if (a.right_node == no_assignment)
          {
             assert(a.left_node < left_non_assignment_costs.size());
             left_non_assignment_costs[a.left_node] += scaling * a.cost;
          }
       }

       for (std::size_t i = 0; i < no_left_nodes; ++i)
       {
          mcf.add_edge(i, no_nodes + 1, 0, 1, left_non_assignment_costs[i]); // for non-assignment
          mcf.add_node_excess(i, 1);
       }

      for(std::size_t i=0; i<no_right_nodes; ++i) {
         mcf.add_edge(no_nodes, no_left_nodes + i, 0, 1, right_non_assignment_costs[i]); // for non-assignment
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

      template<typename Derived>
          struct is_matrix_expression
          : std::is_base_of<Eigen::MatrixBase<std::decay_t<Derived> >, std::decay_t<Derived> >
          {};

      template<typename M, typename = typename std::enable_if<is_matrix_expression<M>::value>::type >
      labeling(const M& m)
      {
         // check that each row and columns has at most one 1 entry
         std::cout << m << "\n";
         for(std::size_t i=0; i<m.cols(); ++i) {
            for(std::size_t j=0; j<m.rows(); ++j) {
               assert(m(i,j) == 0 || m(i,j) == 1);
            }
         }
         for(std::size_t i=0; i<m.cols(); ++i) {
            assert(m.col(i).sum() <= 1);
         }
         for(std::size_t j=0; j<m.rows(); ++j) {
            assert(m.row(j).sum() <= 1);
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

      std::vector<char> compute_right_node_taken() const
      {
         std::vector<char> right_node_taken(highest_matched_node()+1, 0);
         for(const std::size_t r : *this) 
         {
            if(r != no_assignment) 
            {
               assert(r < highest_matched_node()+1);
               right_node_taken[r] = 1;
            }
         }
         return right_node_taken;
      }

      bool assignment_active(const std::size_t left, const std::size_t right, const std::vector<char>& right_node_taken) const
      {
         assert(left != no_assignment || right != no_assignment);
         assert(left == no_assignment || left < no_left_nodes());
         //assert(right == no_assignment || right < highest_matched_node()+1);
         assert(right_node_taken.size() > highest_matched_node());
         if(left != no_assignment && right != no_assignment)
            return (*this)[left] == right;
         else if (right == no_assignment)
            return (*this)[left] == no_assignment;
         else
           return right >= right_node_taken.size() || right_node_taken[right] == 0; 
      }

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


   std::tuple<bool, std::vector<char>> feasible(const labeling& l) const
   {
       if(l.no_left_nodes() != no_left_nodes)
           throw std::runtime_error("labeling must have equal number of left nodes as matching problem.");

       std::vector<char> nodes_taken(no_right_nodes,0);
       for(std::size_t i=0; i<l.size(); ++i) {
           if(l[i] != no_assignment) {
               assert(l[i] < nodes_taken.size());
               if(nodes_taken[l[i]] == 1)
                   return {false, {}};
               nodes_taken[l[i]] = 1;
           }
       }
       return {true, nodes_taken};
       // TODO: check that all assignments are allowed w.r.t. costs
   }

   double evaluate_linear(const labeling &l) const
   {
      double cost = constant_;

      const std::vector<char> right_node_taken = l.compute_right_node_taken();

      for (const auto &a : assignments)
      {
         assert(a.left_node < no_left_nodes || a.left_node == no_assignment);
         assert(a.right_node < no_right_nodes || a.right_node == no_assignment);
         assert(a.left_node != no_assignment || a.right_node != no_assignment);

         if (l.assignment_active(a.left_node, a.right_node, right_node_taken))
            cost += a.cost;
       }

      return cost;
   }

   double evaluate(const labeling &l) const
   {
      return evaluate_linear(l);
       const auto [labeling_feasible, right_nodes_taken] = feasible(l);
       if(!labeling_feasible)
           return std::numeric_limits<double>::infinity();

       double cost = constant_;

       std::vector<char> right_node_taken(no_right_nodes, 0);
       for (const std::size_t r : l)
          if (r != no_assignment)
             right_node_taken[r] = 1;

       for (const auto &a : assignments)
       {
          assert(a.left_node < no_left_nodes || a.left_node == no_assignment);
          assert(a.right_node < no_right_nodes || a.right_node == no_assignment);
          assert(a.left_node != no_assignment || a.right_node != no_assignment);
          if (a.left_node != no_assignment && l[a.left_node] == a.right_node)
          {
             cost += a.cost;
          }
          else if (a.left_node == no_assignment && right_nodes_taken[a.right_node] == 0)
          {
             cost += a.cost;
          }
       }

      return cost;
   }

   std::string linear_identifier(const std::size_t l, const std::size_t r) const
   {
       assert(l < no_left_nodes);
       assert(r < no_right_nodes);
       return std::string("x_") + std::to_string(l) + "_" + std::to_string(r);
   }

   template<typename STREAM>
       void write_linear_objective(STREAM& s) const
       {
           for(const auto& a : assignments)
               s << (a.cost < 0.0 ? "-" : "+") << std::abs(a.cost) << " " << linear_identifier(a.left_node, a.right_node) << "\n";
       }
   template<typename STREAM>
       void write_linear_constraints(STREAM& s) const
       {
           std::vector<std::vector<std::size_t>> left(no_left_nodes);
           std::vector<std::vector<std::size_t>> right(no_right_nodes);
           for(const auto& a : assignments) {
               left[a.left_node].push_back(a.right_node);
               right[a.right_node].push_back(a.left_node);
           }
           for(std::size_t i=0; i<no_left_nodes; ++i) {
               if(left[i].size() > 0) {
                   for(const std::size_t j : left[i]) {
                       s << "+ " << linear_identifier(i,j) << " ";
                   }
                   s << " <= 1\n";//"x_" << i << "_no_assignment = 1\n"; 
               }
           }
           for(std::size_t j=0; j<no_right_nodes; ++j) {
               if(right[j].size() > 0) {
                   for(const std::size_t i : right[j]) {
                       s << "+ " << linear_identifier(i,j) << " ";
                   }
                   s << " <= 1\n";//"x_" << i << "_no_assignment = 1\n"; 
               }
           }
       }
   template<typename STREAM>
       void write_linear_variables(STREAM& s) const
       {
           for(const auto& a : assignments)
               s << linear_identifier(a.left_node, a.right_node) << "\n";
       }

   template<typename STREAM>
       void write_lp(STREAM& s) const
       {
           s << "Minimize\n";
           write_linear_objective(s);

           s << "Subject To\n";
           write_linear_constraints(s);

           s << "Bounds\n";
           s << "Binaries\n";
           write_linear_variables(s);

           s << "End\n";
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
    };
    std::vector<quadratic> quadratic_terms;

    void add_quadratic_term(const std::size_t assignment_1, const std::size_t assignment_2, const double cost)
    {
        assert(assignment_1 != assignment_2);
        assert(std::max(assignment_1, assignment_2) < this->assignments.size());
        const std::size_t a1 = std::min(assignment_1, assignment_2);
        const std::size_t a2 = std::max(assignment_1, assignment_2);
        //if(
        //   this->assignments[a1].left_node == no_assignment ||
        //   this->assignments[a1].right_node == no_assignment ||
        //   this->assignments[a2].left_node == no_assignment ||
        //   this->assignments[a2].right_node == no_assignment)
        //   quadratic_terms.push_back({a1, a2, 0.0});
        //else
           quadratic_terms.push_back({a1, a2, cost});
    }

    double evaluate_quadratic(const labeling &l) const
    {
       double cost = 0.0;
       const auto right_node_taken = l.compute_right_node_taken();
       for(const auto& q : quadratic_terms) {
          const auto a1 = this->assignments[q.assignment_1];
          const auto a2 = this->assignments[q.assignment_2];
          if (l.assignment_active(a1.left_node, a1.right_node, right_node_taken) && l.assignment_active(a2.left_node, a2.right_node, right_node_taken))
             cost += q.cost;
       }
       return cost; 
    }

    double evaluate(const labeling &l) const
    {
       return this->evaluate_linear(l) + evaluate_quadratic(l);
       double cost = linear_assignment_problem_input::evaluate(l);
       for(const auto& q : quadratic_terms) {
          const auto a1 = this->assignments[q.assignment_1];
          const auto a2 = this->assignments[q.assignment_2];
          if(l[a1.left_node] == a1.right_node && l[a2.left_node] == a2.right_node)
             cost += q.cost;
       }
       return cost;
    }

    std::string quadratic_identifier(const std::array<std::size_t,2> left_nodes, const std::array<std::size_t,2> right_nodes) const
    {
        return std::string("q_") + std::to_string(left_nodes[0]) + "_" + std::to_string(left_nodes[1]) + "_" + std::to_string(right_nodes[0]) + "_" + std::to_string(right_nodes[1]);
    }

    std::string quadratic_identifier(const std::size_t assignment_1, const std::size_t assignment_2) const
    {
        const std::size_t l1 = this->assignments[assignment_1].left_node;
        const std::size_t l2 = this->assignments[assignment_2].left_node;
        const std::size_t r1 = this->assignments[assignment_1].right_node;
        const std::size_t r2 = this->assignments[assignment_2].right_node;
        return quadratic_identifier({l1,l2}, {r1,r2});
    }

    template<typename STREAM>
        void write_quadratic_objective(STREAM& s) const
        {
            for(const auto& q : quadratic_terms)
                s << (q.cost < 0.0 ? "-" : "+") << std::abs(q.cost) << " " << quadratic_identifier(q.assignment_1, q.assignment_2) << "\n";
        }

    template<typename STREAM>
        void write_quadratic_constraints(STREAM& s) const
        {
           std::unordered_set<std::array<std::size_t,2>> left_pairwise_potentials;
           std::unordered_set<std::array<std::size_t,2>> right_pairwise_potentials;

           for(const auto& q : quadratic_terms) {
                const std::size_t a_1 = q.assignment_1;
                const std::size_t a_2 = q.assignment_2;
                const std::size_t left_node_1 = this->assignments[a_1].left_node;
                const std::size_t left_node_2 = this->assignments[a_2].left_node;
                const std::size_t right_node_1 = this->assignments[a_1].left_node;
                const std::size_t right_node_2 = this->assignments[a_2].left_node;
                left_pairwise_potentials.insert({left_node_1, left_node_2});
                right_pairwise_potentials.insert({left_node_1, left_node_2});
           }

           std::vector<std::vector<std::size_t>> left_labels(this->no_left_nodes);
           std::vector<std::vector<std::size_t>> right_labels(this->no_right_nodes);
           for(const auto& a : this->assignments) {
               left_labels[a.left_node].push_back(a.right_node);
               right_labels[a.right_node].push_back(a.left_node);
           }

           for(const auto [l_1,l_2] : left_pairwise_potentials) {
               for(const auto r_1 : left_labels[l_1]) {
                   for(const auto r_2 : left_labels[l_2]) {
                       s << " - " << quadratic_identifier({l_1,l_2}, {r_1,r_2});
                   }
                   s << " + " << this->linear_identifier(l_1,r_1) << " = 0\n";
               }

               for(const auto r_2 : left_labels[l_2]) {
                   for(const auto r_1 : left_labels[l_1]) {
                       s << " - " << quadratic_identifier({l_1,l_2}, {r_1,r_2});
                   }
                   s << " + " << this->linear_identifier(l_2,r_2) << " = 0\n";
               }
           }
        }
    template<typename STREAM>
        void write_quadratic_variables(STREAM& s) const
        {
            for(const auto& q : quadratic_terms)
                s << " " << quadratic_identifier(q.assignment_1, q.assignment_2) << "\n";
        }

    template<typename STREAM>
        void write_lp(STREAM& s) const
        {
            s << "Minimize\n";
            this->write_linear_objective(s);
            write_quadratic_objective(s);

            s << "Subject To\n";
            this->write_linear_constraints(s);
            write_quadratic_constraints(s);

            s << "Bounds\n";
            s << "Binaries\n";
            this->write_linear_variables(s);
            this->write_quadratic_variables(s);

            s << "End\n";
        }

    template<typename STREAM>
        void write_torresani_et_al(STREAM& s) const
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

   public:
   template<typename STREAM>
       void write_torresani_et_al(STREAM& s) const
       {
           for(std::size_t c=0; c<(*this).size(); ++c) {
               const auto gm_input = (*this)[c];
               s << "gm " << gm_input.left_graph_no << " " << gm_input.right_graph_no << "\n";
               gm_input.gm_input.write_torresani_et_al(s);
           } 
       } 
};

} // namespace LPMP
