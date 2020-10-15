#pragma once

#include "mrf/mrf_problem_construction.hxx"
#include "min_cost_flow_factor_ssp.hxx"
#include "tree_decomposition.hxx"
#include "graph_matching_input.h"
#include "graph_matching_frank_wolfe.h"

namespace LPMP {

template<typename MRF_CONSTRUCTOR>
class graph_matching_mrf_constructor {
public:
    using FMC = typename MRF_CONSTRUCTOR::FMC;

    template<typename SOLVER>
    graph_matching_mrf_constructor(SOLVER& solver)
        : lp_(&solver.GetLP()),
        left_mrf(solver),
        right_mrf(solver),
        construction_arg_("", "graphMatchingConstruction", "mode of constructing pairwise potentials for graph matching", false, "left", "{left|right|both_sides}", solver.get_cmd()),
        rounding_arg_("", "graphMatchingRounding", "method for graph matching rounding", false, "mcf", "{mcf/fw}", solver.get_cmd()),
        frank_wolfe_iterations_arg_("","graphMatchingFrankWolfeIterations", "how many iterations to run the Frank Wolfe method for graph matching rounding", false, 400, "integer > 0", solver.get_cmd()) 
    {}

    graph_matching_mrf_constructor(LP<FMC> *lp, const std::string &construction_method, const std::string &rounding_method, const std::size_t frank_wolfe_iterations)
        : lp_(lp),
          left_mrf(lp),
          right_mrf(lp),
          construction_arg_("", "", "", false, construction_method, ""),
          rounding_arg_("", "", "", false, rounding_method, ""),
          frank_wolfe_iterations_arg_("", "", "", false, frank_wolfe_iterations, "")
    {}

    void order_factors()
    {
       assert(no_left_nodes() > 0 && no_right_nodes() > 0);

       left_mrf.order_factors();
       right_mrf.order_factors();

       assert(left_mrf.get_number_of_variables() > 0);
       assert(right_mrf.get_number_of_variables() > 0);

       auto* f_l = left_mrf.get_unary_factor( left_mrf.get_number_of_variables()-1 );
       auto* f_r = right_mrf.get_unary_factor(0);

       lp_->add_factor_relation(f_l, f_r);
    }

    void begin()
    {
       // presolve w.r.t. linear costs with minimum cost flow solver and initialize unary mrf factors with solution
      //min_cost_flow_factor(assignments_);
    }

    void pre_iterate()
    {
        reparametrize_linear_assignment_problem();
    }

    std::vector<FactorTypeAdapter*> get_factors()
    {
       auto left_factors = left_mrf.get_factors();
       auto right_factors = right_mrf.get_factors();
       left_factors.insert(left_factors.end(), right_factors.begin(), right_factors.end());
       return left_factors; 
    }

    void send_messages_to_unaries()
    {
       left_mrf.send_messages_to_unaries();
       right_mrf.send_messages_to_unaries();
    }

    void ComputePrimal()
    {
       send_messages_to_unaries();
       const auto l = [&]() {
       if(rounding_arg_.getValue() == "mcf") {
           return compute_primal_mcf_solution();
       } else if(rounding_arg_.getValue() == "fw") {
          return compute_primal_fw_solution();
       } else {
           throw std::runtime_error("rounding method not recognized.");
       }
       }();
       read_in_labeling(l);
       const double labeling_cost = this->lp_->EvaluatePrimal(); // TODO: storing best solution should be done by the solver class.
       if(labeling_cost < best_labeling_cost_) {
          best_labeling_cost_ = labeling_cost;
          best_labeling_ = l; 
       }
    }

    bool CheckPrimalConsistency() const
    {
        if(debug()) 
            std::cout << "check graph matching assignment\n";

        std::vector<char> labels_taken(inverse_graph_.size(),false);
        for(std::size_t i=0; i<left_mrf.get_number_of_variables(); ++i) {
            const auto state = left_mrf.get_unary_factor(i)->get_factor()->primal();
            if(state < graph_[i].size()) {
                const std::size_t label = graph_[i][state];
                if(labels_taken[label]) { 
                    if(debug())
                       std::cout << "var " << i << ", state " << state << ", label " << label << " conflict\n";
                    return false; 
                }  
                labels_taken[ label ] = true;
            }  
        }  

        return true;
    }

    INDEX Tighten(const INDEX no_constraints_to_add)
    {
       return this->left_mrf.Tighten(no_constraints_to_add) + this->right_mrf.Tighten(no_constraints_to_add);
    }

    std::size_t no_left_nodes() const { return left_mrf.get_number_of_variables(); }
    std::size_t no_right_nodes() const { return right_mrf.get_number_of_variables(); }

    template<typename STREAM>
    void WritePrimal(STREAM& s) const
    {
       for(std::size_t i=0; i<graph_.size(); ++i) {
          const std::size_t primal_label = left_mrf.get_unary_factor(i)->get_factor()->primal();
          if(primal_label < graph_[i].size())
             s << i << " -> " << graph_[i][primal_label] << "\n";
          else
             s << i << " not matched\n";
       }
    }

    void construct(const graph_matching_input& gm_input)
    {
       construct_graphs(gm_input);
       construct_unary_factors(gm_input);
       construct_pairwise_factors(gm_input);
       mcf_ = std::make_unique<mcf_solver_type>(gm_input.no_mcf_nodes(), gm_input.no_mcf_edges());
       gm_input.initialize_mcf(*mcf_);
    }

    // solve underlying linear assignment problem with combinatorial problems. Use optimal dual node potentials to reparametrize unary factors in MRFs
    void reparametrize_linear_assignment_problem()
    {
       auto compute_lb = [&]() {
          double lb = 0.0;
          for(std::size_t i=0; i<left_mrf.get_number_of_variables(); ++i) {
             lb += left_mrf.get_unary_factor(i)->LowerBound();
          }
          for(std::size_t i=0; i<right_mrf.get_number_of_variables(); ++i) {
             lb += right_mrf.get_unary_factor(i)->LowerBound();
          }
          return lb;
       };

       const std::size_t no_left_nodes = left_mrf.get_number_of_variables();
       const std::size_t no_right_nodes = right_mrf.get_number_of_variables(); 

       auto compute_potential_bound = [&]() {
          double dual_cost = 0.0;
          for(std::size_t i=0; i<no_left_nodes; ++i) {
             dual_cost += mcf_->potential(i);
          }
          dual_cost += no_right_nodes * mcf_->potential(no_left_nodes + no_right_nodes);
          for(std::size_t i=0; i<no_right_nodes; ++i) {
             dual_cost -= mcf_->potential(no_left_nodes+i);
          }
          dual_cost -= no_left_nodes * mcf_->potential(no_left_nodes + no_right_nodes+1);
          return dual_cost;
       };

       read_in_mcf_costs();

       mcf_->solve(); 
       assert(std::abs(mcf_->potential(no_left_nodes+no_right_nodes) - mcf_->potential(no_left_nodes+no_right_nodes+1)) <= eps); 
       assert(compute_lb() <= mcf_->objective() + eps);

       write_back_mcf_costs(); 
       assert(std::abs(compute_lb() + compute_potential_bound() - mcf_->objective())/double(no_left_nodes + no_right_nodes) <= eps);

       left_mrf.get_lp()->add_to_constant(compute_potential_bound()); 
    }

    linear_assignment_problem_input::labeling compute_primal_mcf_solution()
    {
       if(debug())
          std::cout << "compute primal mcf solution\n";
       read_in_mcf_costs();
       mcf_->solve();

       linear_assignment_problem_input::labeling labeling(no_left_nodes(), std::numeric_limits<std::size_t>::max());

       for(std::size_t i=0; i<no_left_nodes(); ++i) {
          auto& u = *left_mrf.get_unary_factor(i)->get_factor();
          auto e = mcf_->first_outgoing_arc(i);
          for(std::size_t l=0; l<u.size()-1; ++l, ++e) {
             assert(mcf_->tail(e) == i);
             if(mcf_->flow(e) == 1) {
                assert(labeling[i] == std::numeric_limits<std::size_t>::max());
                labeling[i] = graph_[i][l];
             }
          }
          assert(e - mcf_->first_outgoing_arc(i) == graph_[i].size());
       }

       assert(labeling.check_primal_consistency());
       return labeling;
    }

    linear_assignment_problem_input::labeling compute_primal_fw_solution()
    {
       const graph_matching_input input = export_graph_matching_input();
       const graph_matching_input::labeling labeling = compute_primal_mcf_solution();

       graph_matching_frank_wolfe_options o;
       o.max_iter = frank_wolfe_iterations_arg_.getValue();
       graph_matching_frank_wolfe gm_fw(input, labeling, o);
       const auto l = gm_fw.solve();

       assert(input.evaluate(l) < std::numeric_limits<double>::infinity());
       assert(input.evaluate(l) >= this->lp_->LowerBound() - 1e-8);
       return l;
    }

    graph_matching_input export_linear_graph_matching_input() const
    {
       graph_matching_input instance;

       const std::size_t no_left_nodes = left_mrf.get_number_of_variables();
       const std::size_t no_right_nodes = right_mrf.get_number_of_variables();

       instance.add_to_constant(left_mrf.get_lp()->get_constant());

       // linear terms
       for(std::size_t i=0; i<no_left_nodes; ++i) {
          auto& u = *left_mrf.get_unary_factor(i)->get_factor();
          for(std::size_t l=0; l<u.size()-1; ++l) {
             instance.add_assignment(i, graph_[i][l], u[l]); // - u.back());
          }
          instance.add_assignment(i, graph_matching_input::no_assignment, u.back()); // non-assignment
       }

       for(std::size_t i=0; i<no_right_nodes; ++i) {
          auto& u = *right_mrf.get_unary_factor(i)->get_factor();
          for(std::size_t l=0; l<u.size()-1; ++l) {
             instance.add_assignment(inverse_graph_[i][l], i, u[l]); // - u.back());
          }
          instance.add_assignment(graph_matching_input::no_assignment, i, u.back()); // non-assignment
       }

       instance.normalize(); // merge parallel edges
       return instance;
    }

    graph_matching_input export_graph_matching_input() const
    {
       if(debug())
          std::cout << "export linear assignment problem solution based on current reparametrization to graph matching problem\n";
       graph_matching_input instance = export_linear_graph_matching_input();

       std::unordered_map<std::array<std::size_t,2>, std::size_t> assignment_nr;
       for(std::size_t idx=0; idx<instance.assignments.size(); ++idx) {
           const auto& a = instance.assignments[idx];
           assignment_nr.insert({{a.left_node, a.right_node}, idx});
       }

       // quadratic terms
       for(std::size_t pairwise_idx=0; pairwise_idx<left_mrf.get_number_of_pairwise_factors(); ++pairwise_idx) {
           const auto [i,j] = left_mrf.get_pairwise_variables(pairwise_idx);
           const auto& pot = *left_mrf.get_pairwise_factor(pairwise_idx)->get_factor();
           assert(pot.dim1() == left_mrf.get_number_of_labels(i) && pot.dim1() == graph_[i].size()+1);
           assert(pot.dim2() == left_mrf.get_number_of_labels(j) && pot.dim2() == graph_[j].size()+1);
           for(std::size_t l_i=0; l_i+1<pot.dim1(); ++l_i) {
               assert(assignment_nr.count({i, graph_[i][l_i]}) > 0);
               const std::size_t a1 = assignment_nr.find({i, graph_[i][l_i]})->second;
               for(std::size_t l_j=0; l_j+1<pot.dim2(); ++l_j) {
                   if(graph_[i][l_i] != graph_[j][l_j]) {
                       if(std::abs(pot(l_i,l_j)) >= 1e-8) {
                           assert(assignment_nr.count({j, graph_[j][l_j]}) > 0);
                           const std::size_t a2 = assignment_nr.find({j, graph_[j][l_j]})->second;
                           instance.add_quadratic_term(a1, a2, pot(l_i,l_j)); 
                       }
                   } else {
                       assert(pot(l_i,l_j) == std::numeric_limits<double>::infinity());
                   }
               }
           }

           const std::size_t i_non_assignment = assignment_nr.find({i, graph_matching_input::no_assignment})->second;
           const std::size_t j_non_assignment = assignment_nr.find({j, graph_matching_input::no_assignment})->second;

           for(std::size_t l_i=0; l_i+1<pot.dim1(); ++l_i) {
               const std::size_t a1 = assignment_nr.find({i, graph_[i][l_i]})->second;
               instance.add_quadratic_term(a1, j_non_assignment, pot(l_i, pot.dim2()-1)); 
           }
           for(std::size_t l_j=0; l_j+1<pot.dim2(); ++l_j) {
               const std::size_t a2 = assignment_nr.find({j, graph_[j][l_j]})->second;
               instance.add_quadratic_term(i_non_assignment, a2, pot(pot.dim1()-1, l_j)); 
           }

           instance.add_quadratic_term(i_non_assignment, j_non_assignment, pot(pot.dim1()-1, pot.dim2()-1));
       }

       for(std::size_t pairwise_idx=0; pairwise_idx<right_mrf.get_number_of_pairwise_factors(); ++pairwise_idx) {
           const auto [i,j] = right_mrf.get_pairwise_variables(pairwise_idx);
           const auto& pot = *right_mrf.get_pairwise_factor(pairwise_idx)->get_factor();
           assert(pot.dim1() == right_mrf.get_number_of_labels(i) && pot.dim1() == inverse_graph_[i].size()+1);
           assert(pot.dim2() == right_mrf.get_number_of_labels(j) && pot.dim2() == inverse_graph_[j].size()+1);
           for(std::size_t l_i=0; l_i+1<pot.dim1(); ++l_i) {
               assert(assignment_nr.count({inverse_graph_[i][l_i], i}) > 0);
               const std::size_t a1 = assignment_nr.find({inverse_graph_[i][l_i], i})->second;
               for(std::size_t l_j=0; l_j+1<pot.dim2(); ++l_j) {
                   if(inverse_graph_[i][l_i] != inverse_graph_[j][l_j]) {
                       if(std::abs(pot(l_i,l_j)) >= 1e-8) {
                           assert(assignment_nr.count({inverse_graph_[j][l_j], j}) > 0);
                           const std::size_t a2 = assignment_nr.find({inverse_graph_[j][l_j], j})->second;
                           instance.add_quadratic_term(a1, a2, pot(l_i,l_j)); 
                       }
                   } else {
                       assert(pot(l_i,l_j) == std::numeric_limits<double>::infinity());
                   }
               }
           }

           const std::size_t i_non_assignment = assignment_nr.find({graph_matching_input::no_assignment, i})->second;
           const std::size_t j_non_assignment = assignment_nr.find({graph_matching_input::no_assignment, j})->second;

           for(std::size_t l_i=0; l_i+1<pot.dim1(); ++l_i) {
               const std::size_t a1 = assignment_nr.find({inverse_graph_[i][l_i], i})->second;
               instance.add_quadratic_term(a1, j_non_assignment, pot(l_i, pot.dim2()-1)); 
           }
           for(std::size_t l_j=0; l_j+1<pot.dim2(); ++l_j) {
               const std::size_t a2 = assignment_nr.find({inverse_graph_[j][l_j], j})->second;
               instance.add_quadratic_term(i_non_assignment, a2, pot(pot.dim1()-1, l_j)); 
           }

           instance.add_quadratic_term(i_non_assignment, j_non_assignment, pot(pot.dim1()-1, pot.dim2()-1));
       }

       //instance.normalize_quadratic_terms();

       return instance;
    }

    void read_in_labeling(const linear_assignment_problem_input::labeling& l)
    {
       assert(l.check_primal_consistency());
       if(l.size() != no_left_nodes())
          throw std::runtime_error("labeling has different number of entries than graph matching problem.");
       if(l.highest_matched_node() > no_right_nodes())
          throw std::runtime_error("labeling has more right nodes than model.");

       for(std::size_t i=0; i<no_left_nodes(); ++i) {
          auto* left_factor = left_mrf.get_unary_factor(i)->get_factor();
          left_factor->primal() = get_left_index(i, l[i]);
          assert(left_factor->primal() < left_factor->size());
          if(l[i] < std::numeric_limits<std::size_t>::max()) {
             assert(graph_[i][ get_left_index(i,l[i]) ] == l[i]);
          }

          if(l[i] < no_right_nodes()) {
             auto* right_factor = right_mrf.get_unary_factor(l[i])->get_factor();
             right_factor->primal() = get_right_index(l[i], i);
             assert(inverse_graph_[l[i]][ get_right_index(l[i], i) ] == i);
             assert(right_factor->primal()+1 < right_factor->size());
          }
       }

       for(std::size_t i=0; i<no_right_nodes(); ++i) {
             auto* right_factor = right_mrf.get_unary_factor(i)->get_factor();
             if(right_factor->primal() >= right_factor->size())
                right_factor->primal() = right_factor->size()-1;
       }

       for(std::size_t i=0; i<no_left_nodes(); ++i)
          left_mrf.get_unary_factor(i)->propagate_primal_through_messages();
       for(std::size_t i=0; i<no_right_nodes(); ++i)
          right_mrf.get_unary_factor(i)->propagate_primal_through_messages(); 

       assert(l == write_out_labeling());
       assert(CheckPrimalConsistency());
    }

    linear_assignment_problem_input::labeling write_out_labeling() const
    {
       linear_assignment_problem_input::labeling output;
       output.reserve(no_left_nodes());

       for(std::size_t i=0; i<no_left_nodes(); ++i) {
          auto& left_factor = *left_mrf.get_unary_factor(i)->get_factor();
          if(left_factor.primal() >= left_factor.size()) 
             throw std::runtime_error("no valid solution for graph matching problem.");
          if(left_factor.primal() < left_factor.size()-1)
             output.push_back(graph_[i][left_factor.primal()]);
          else
             output.push_back(std::numeric_limits<std::size_t>::max());
       }

       assert(output.size() == no_left_nodes()); 
       return output;
    }

    graph_matching_input::labeling best_labeling() const
    {
       return best_labeling_;
    }

    bool has_edge(const std::size_t node_left, const std::size_t node_right) const
    {
       assert(node_left < no_left_nodes());
       assert(node_right < no_right_nodes());
       assert(std::binary_search(graph_[node_left].begin(), graph_[node_left].end(), node_right) == std::binary_search(inverse_graph_[node_right].begin(), inverse_graph_[node_right].end(), node_left));
       return std::binary_search(graph_[node_left].begin(), graph_[node_left].end(), node_right);
    } 

    std::size_t get_left_index(const std::size_t left_node, const std::size_t right_node) const
    {
       return get_index(left_node, right_node, graph_);
    }

    std::size_t get_right_index(const std::size_t right_node, const std::size_t left_node) const
    {
       return get_index(right_node, left_node, inverse_graph_);
    } 

protected:

    static std::size_t get_index(const std::size_t node_1, const std::size_t node_2, const std::vector<std::vector<std::size_t>>& graph)
    {
       assert(node_1 < graph.size());
       assert(std::is_sorted(graph[node_1].begin(), graph[node_1].end()));
       for (std::size_t i = 0; i + 1 < graph[node_1].size(); ++i)
       {
          assert(graph[node_1][i] < graph[node_1][i + 1]);
       }
       if (node_2 == graph_matching_input::no_assignment)
          return graph[node_1].size();

       assert(std::binary_search(graph[node_1].begin(), graph[node_1].end(), node_2));
       const std::size_t index = std::lower_bound(graph[node_1].begin(), graph[node_1].end(), node_2) - graph[node_1].begin();
       assert(index < graph[node_1].size() && graph[node_1][index] == node_2);
       return index;
    }


    void construct_graphs(const linear_assignment_problem_input& input)
    {
       graph_.resize(input.no_left_nodes);
       inverse_graph_.resize(input.no_right_nodes);
       for(const auto& a : input.assignments) {
          assert(a.left_node < graph_.size() || a.left_node == graph_matching_input::no_assignment);
          assert(a.right_node < inverse_graph_.size() || a.right_node == graph_matching_input::no_assignment);

          if (a.left_node != graph_matching_input::no_assignment && a.right_node != graph_matching_input::no_assignment) {
             graph_[a.left_node].push_back(a.right_node);
             inverse_graph_[a.right_node].push_back(a.left_node);
          }
       }

        for(auto& v : graph_) { std::sort(v.begin(), v.end()); }
        for(auto& v : inverse_graph_) { std::sort(v.begin(), v.end()); }
    }

   void construct_empty_unary_factors(MRF_CONSTRUCTOR& mrf, std::vector<std::vector<std::size_t>>& graph)
   {
      assert(mrf.get_number_of_variables() == 0);
      for(std::size_t i=0; i<graph.size(); ++i) {
         std::vector<REAL> costs(graph[i].size()+1, 0.0);
         mrf.add_unary_factor(costs);
      }
   }
   void construct_unary_factors(const linear_assignment_problem_input& input)
   {
       construct_empty_unary_factors(left_mrf, graph_);
       construct_empty_unary_factors(right_mrf, inverse_graph_);

       const auto [left_weight, right_weight] = [&]() -> std::array<double,2> {
           if(construction_arg_.getValue() == "left") {
               return {1.0, 0.0};
           } else if(construction_arg_.getValue() == "right") {
               return {0.0, 1.0};
           } else if(construction_arg_.getValue() == "both_sides") {
               return {0.5, 0.5};
           } else {
               throw std::runtime_error("graphMatchingConstruction argument not recognized");
           }
       }();

       for(const auto& a : input.assignments) {
          auto [left_index, left_factor] = [&]() -> std::tuple<std::size_t, typename MRF_CONSTRUCTOR::UnaryFactorContainer*> {
             if (a.left_node != graph_matching_input::no_assignment)
             {
                const std::size_t left_index = get_left_index(a.left_node, a.right_node);
                auto *left_factor = left_mrf.get_unary_factor(a.left_node);
                assert((*left_factor->get_factor())[left_index] == 0.0);
                return {left_index, left_factor};
             }
             else
             {
                return {graph_matching_input::no_assignment, nullptr};
             }
          }();

          auto [right_index, right_factor] = [&]() -> std::tuple<std::size_t, typename MRF_CONSTRUCTOR::UnaryFactorContainer*> {
             if (a.right_node != graph_matching_input::no_assignment)
             {
                const std::size_t right_index = get_right_index(a.right_node, a.left_node);
                auto *right_factor = right_mrf.get_unary_factor(a.right_node);
                assert((*right_factor->get_factor())[right_index] == 0.0);
                return {right_index, right_factor};
             }
             else
             {
                return {graph_matching_input::no_assignment, nullptr};
             }
          }();

          if(left_factor != nullptr && right_factor != nullptr) {
             (*left_factor->get_factor())[left_index] = left_weight * a.cost;
             (*right_factor->get_factor())[right_index] = right_weight * a.cost;
          }
          else if (left_factor != nullptr)
          {
             (*left_factor->get_factor())[left_index] = a.cost;
          }
          else if (right_factor != nullptr)
          {
             (*right_factor->get_factor())[right_index] = a.cost;
          }
       }
   }

   template<Chirality CHIRALITY>
   void construct_pairwise_factors(MRF_CONSTRUCTOR& mrf, const REAL scaling, std::vector<std::vector<std::size_t>>& graph, const graph_matching_input& input)
   {
      assert(scaling == 0.0 || scaling == 0.5 || scaling == 1.0);
      auto term_on_left_side_only = [&](const std::array<std::size_t,2> a1, const std::array<std::size_t,2> a2) {
         if(a1[1] == graph_matching_input::no_assignment || a2[1] == graph_matching_input::no_assignment)
            return true;
         return false;
      };
      auto term_on_right_side_only = [&](const std::array<std::size_t,2> a1, const std::array<std::size_t,2> a2) {
         if(a1[0] == graph_matching_input::no_assignment || a2[0] == graph_matching_input::no_assignment)
            return true;
         return false;
      };
      auto term_on_both_sides = [&](const std::array<std::size_t, 2> a1, const std::array<std::size_t, 2> a2) {
         if (a1[0] != graph_matching_input::no_assignment && a1[1] != graph_matching_input::no_assignment && a2[0] != graph_matching_input::no_assignment && a2[1] != graph_matching_input::no_assignment)
            return true;
         return false;
      };
      for (const auto &q : input.quadratic_terms)
      {
         const auto assignment_1 = input.assignments[q.assignment_1];
         const auto assignment_2 = input.assignments[q.assignment_2];
         const bool on_left_only = term_on_left_side_only({assignment_1.left_node, assignment_1.right_node}, {assignment_2.left_node, assignment_2.right_node});
         const bool on_right_only = term_on_right_side_only({assignment_1.left_node, assignment_1.right_node}, {assignment_2.left_node, assignment_2.right_node});
         const bool on_both_sides = term_on_both_sides({assignment_1.left_node, assignment_1.right_node}, {assignment_2.left_node, assignment_2.right_node});
         assert(int(on_left_only) + int(on_right_only) + int(on_both_sides) == 1);

         const auto [node_1, node_2, label_1, label_2] = [&]() {
            if (Chirality::left == CHIRALITY)
            {
               if (assignment_1.left_node < assignment_2.left_node)
               {
                  return std::make_tuple(assignment_1.left_node, assignment_2.left_node, assignment_1.right_node, assignment_2.right_node);
               }
               else
               {
                  return std::make_tuple(assignment_2.left_node, assignment_1.left_node, assignment_2.right_node, assignment_1.right_node);
               }
            }
            else
            {
               assert(CHIRALITY == Chirality::right);
               if (assignment_1.right_node < assignment_2.right_node)
               {
                  return std::make_tuple(assignment_1.right_node, assignment_2.right_node, assignment_1.left_node, assignment_2.left_node);
               }
               else
               {
                  return std::make_tuple(assignment_2.right_node, assignment_1.right_node, assignment_2.left_node, assignment_1.left_node);
               }
            }
         }();

         assert(node_1 < node_2 || node_1 == graph_matching_input::no_assignment);

         if (node_1 != graph_matching_input::no_assignment && node_2 != graph_matching_input::no_assignment)
         {
            if (!mrf.has_pairwise_factor(node_1, node_2))
            {
               auto *f = mrf.add_empty_pairwise_factor(node_1, node_2);
               // add infinities on diagonal
               for (std::size_t i1 = 0; i1 < graph[node_1].size(); ++i1)
               {
                  for (std::size_t i2 = 0; i2 < graph[node_2].size(); ++i2)
                  {
                     if (graph[node_1][i1] == graph[node_2][i2])
                     {
                        f->get_factor()->cost(i1, i2) = std::numeric_limits<REAL>::infinity();
                     }
                  }
               }
            }

            auto *f = mrf.get_pairwise_factor(node_1, node_2);

            const auto index_1 = get_index(node_1, label_1, graph);
            const auto index_2 = get_index(node_2, label_2, graph);

            if (on_left_only || on_right_only) {
               assert(index_1 + 1 == f->get_factor()->dim1() || index_2 + 1 == f->get_factor()->dim2());
               f->get_factor()->cost(index_1, index_2) += q.cost; // TODO: remove scaling for entry that is only present in left or right mrf. Correct?
            } else
               f->get_factor()->cost(index_1, index_2) += scaling * q.cost;
         }
       } 
   }
   void construct_pairwise_factors(const graph_matching_input& input)
   {
      // TODO: check if all quadratic entries are on correct side when using only left or right
       if(construction_arg_.getValue() == "left") {
           construct_pairwise_factors<Chirality::left>(left_mrf, 1.0, graph_, input);
       } else if(construction_arg_.getValue() == "right") {
           construct_pairwise_factors<Chirality::right>(right_mrf, 1.0, inverse_graph_, input);
       } else if(construction_arg_.getValue() == "both_sides") {
           construct_pairwise_factors<Chirality::left>(left_mrf, 0.5, graph_, input);
           construct_pairwise_factors<Chirality::right>(right_mrf, 0.5, inverse_graph_, input);
       } else {
           throw std::runtime_error("graphMatchingConstruction argument not recognized");
       }
   }

    void read_in_mcf_costs()
    {
       const std::size_t no_left_nodes = left_mrf.get_number_of_variables();
       const std::size_t no_right_nodes = right_mrf.get_number_of_variables();

       // read out costs from unary factors and write them to mcf
       mcf_->reset_costs();
       assert(std::abs(mcf_->solve()) <= eps);

       for(std::size_t i=0; i<left_mrf.get_number_of_variables(); ++i) {
          auto& u = *left_mrf.get_unary_factor(i)->get_factor();
          assert(mcf_->no_outgoing_arcs(i) == u.size());
          auto e = mcf_->first_outgoing_arc(i); 
          for(std::size_t l=0; l<u.size(); ++l, ++e) {
             assert(mcf_->cost(e) == 0.0);
             mcf_->update_cost(e, u[l]);
          }
       }

       for(std::size_t i=0; i<right_mrf.get_number_of_variables(); ++i) {
          auto& u = *right_mrf.get_unary_factor(i)->get_factor();
          assert(mcf_->no_outgoing_arcs(no_left_nodes+i) == u.size());
          auto e = mcf_->first_outgoing_arc(no_left_nodes+i);
          for(std::size_t l=0; l<u.size(); ++l, ++e) {
             mcf_->update_cost(e, -u[l]);
          }
       }
    }

    void write_back_mcf_costs()
    {
       const std::size_t no_left_nodes = left_mrf.get_number_of_variables();
       const std::size_t no_right_nodes = right_mrf.get_number_of_variables();

       for(std::size_t i=0; i<left_mrf.get_number_of_variables(); ++i) {
          auto& u = *left_mrf.get_unary_factor(i)->get_factor();
          auto e = mcf_->first_outgoing_arc(i);
          for(std::size_t l=0; l<u.size()-1; ++l, ++e) {
             u[l] = 0.5*mcf_->reduced_cost(e);
             assert(mcf_->tail(e) == i);
          }
          assert(e == mcf_->first_outgoing_arc(i)+u.size()-1);
          u[u.size()-1] = mcf_->reduced_cost( mcf_->first_outgoing_arc(i) + u.size()-1 );
          assert(std::count_if(u.begin(), u.end(), [&](const double& x) { return x < -eps; }) <= 1);
       }

       for(std::size_t i=0; i<right_mrf.get_number_of_variables(); ++i) {
          auto& u = *right_mrf.get_unary_factor(i)->get_factor();
          assert(mcf_->no_outgoing_arcs(no_left_nodes+i) == u.size());
          auto e = mcf_->first_outgoing_arc(no_left_nodes+i);
          for(std::size_t l=0; l<u.size()-1; ++l, ++e) {
             u[l] = -0.5*mcf_->reduced_cost(e);
             assert(mcf_->tail(e) == no_left_nodes+i);
          }
          assert(e == mcf_->first_outgoing_arc(no_left_nodes+i)+u.size()-1);
          u[u.size()-1] = -mcf_->reduced_cost( mcf_->first_outgoing_arc(no_left_nodes+i) + u.size()-1 );
          assert(std::count_if(u.begin(), u.end(), [&](const double& x) { return x < -eps; }) <= 1);
       }
    }

   LP<FMC>* lp_;
public:
   MRF_CONSTRUCTOR left_mrf, right_mrf;

   // TODO: make two_dim_variable_array<std::size_t> out of two structures below
   std::vector<std::vector<std::size_t>> graph_;
   std::vector<std::vector<std::size_t>> inverse_graph_;

   using mcf_solver_type = MCF::SSP<long,REAL>;
   std::unique_ptr<mcf_solver_type> mcf_;

private:
    TCLAP::ValueArg<std::string> construction_arg_;
    TCLAP::ValueArg<std::string> rounding_arg_;
    TCLAP::ValueArg<std::size_t> frank_wolfe_iterations_arg_;

protected: 
   // TODO: remove and change mcf construction in mcf solver
   struct assignment {const std::size_t left; const std::size_t right; const REAL cost;};
   std::vector<assignment> assignments_;
   struct quadratic {const std::size_t assignment_1; const std::size_t assignment_2; const REAL cost;};
   std::vector<quadratic> quadratic_; 
   graph_matching_input::labeling best_labeling_;
   double best_labeling_cost_ = std::numeric_limits<double>::infinity();
};

template<typename GRAPH_MATCHING_MRF_CONSTRUCTOR, typename ASSIGNMENT_MESSAGE>
class graph_matching_constructor : public GRAPH_MATCHING_MRF_CONSTRUCTOR {
public:
    using FMC = typename GRAPH_MATCHING_MRF_CONSTRUCTOR::FMC;
    using GRAPH_MATCHING_MRF_CONSTRUCTOR::GRAPH_MATCHING_MRF_CONSTRUCTOR;

    void construct(const graph_matching_input& gm_input)
    {
        GRAPH_MATCHING_MRF_CONSTRUCTOR::construct(gm_input);

        for(std::size_t left_node=0; left_node<this->graph_.size(); ++left_node) {
            auto* l = this->left_mrf.get_unary_factor(left_node);
            for(std::size_t left_index=0; left_index<this->graph_[left_node].size(); ++left_index) {
                const auto right_node = this->graph_[left_node][left_index];
                auto* r = this->right_mrf.get_unary_factor(right_node);
                this->lp_->template add_message<ASSIGNMENT_MESSAGE>(l, r, left_index, this->get_right_index(right_node, left_node));
            }
        }
    }
};

template<typename GRAPH_MATCHING_CONSTRUCTOR, typename MCF_FACTOR, typename MCF_MESSAGE>
class graph_matching_mcf_constructor : public GRAPH_MATCHING_CONSTRUCTOR {
public:
    using FMC = typename GRAPH_MATCHING_CONSTRUCTOR::FMC;

    using GRAPH_MATCHING_CONSTRUCTOR::GRAPH_MATCHING_CONSTRUCTOR;

    std::vector<FactorTypeAdapter*> get_factors()
    {
       auto factors = GRAPH_MATCHING_CONSTRUCTOR::get_factors();
       assert(mcf_factor != nullptr);
       factors.push_back(mcf_factor);
       return factors;
    }

    void construct(const graph_matching_input& input, factor_tree<FMC>* tree = nullptr)
    {
        GRAPH_MATCHING_CONSTRUCTOR::construct(input);

        // build assignment problem
        const auto no_left_nodes = this->graph_.size();
        const auto no_right_nodes = this->inverse_graph_.size();
        const auto no_nodes = no_left_nodes + no_right_nodes;
        const auto no_edges = this->assignments_.size() + no_left_nodes + no_right_nodes + 1;

        std::vector<typename min_cost_flow_factor::Edge> edges;
        edges.reserve(no_edges);
        for(const auto a : this->assignments_) {
            edges.push_back({a.left, no_left_nodes + a.right, 0, 1, 0.0});
        }
        std::vector<SIGNED_INDEX> demands;
        demands.reserve(no_left_nodes + no_right_nodes + 2);

        // slack edges
        for(std::size_t i=0; i<no_left_nodes; ++i) {
            edges.push_back({i,no_nodes + 1, 0, 1, 0.0}); // for non-assignment
            demands.push_back(1);
        }
        for(std::size_t i=0; i<no_right_nodes; ++i) {
            edges.push_back({no_nodes, no_left_nodes + i, 0, 1, 0.0}); // for non-assignment
            demands.push_back(-1);
        }
        edges.push_back({no_nodes, no_nodes + 1, 0, no_left_nodes + no_right_nodes, 0.0});
        demands.push_back(no_right_nodes);
        demands.push_back(-no_left_nodes);

        mcf_factor = this->lp_->template add_factor<MCF_FACTOR>(edges, demands);

        // connect assignment factor with unaries
        // left side
        {
            std::vector<std::vector<std::size_t>> edgeId(no_left_nodes);
            for(std::size_t i=0; i<no_left_nodes; ++i) {
                //assert(mrf_left.GetNumberOfLabels(i) == mcf->NoArcs(i));
                auto *u = this->left_mrf.get_unary_factor(i);
                const auto first_arc = mcf_factor->get_factor()->mcf_.first_outgoing_arc(i);
                const auto no_arcs = mcf_factor->get_factor()->mcf_.no_outgoing_arcs(i);
                auto* m = this->lp_->template add_message<MCF_MESSAGE>(u, mcf_factor, first_arc, no_arcs);

                if(tree) {
                    tree->add_message(m, Chirality::right);
                }
            }
        }
        // right side
        {
            for(std::size_t i=0; i<no_right_nodes; ++i) {
                //assert(mrf_right.GetNumberOfLabels(i) == mcf->NoArcs(no_left_nodes + i));
                auto *u = this->right_mrf.get_unary_factor(i);
                const auto first_arc = mcf_factor->get_factor()->mcf_.first_outgoing_arc(no_left_nodes+i);
                const auto no_arcs = mcf_factor->get_factor()->mcf_.no_outgoing_arcs(no_left_nodes+i);
                auto* m = this->lp_->template add_message<MCF_MESSAGE>(u, mcf_factor, first_arc, no_arcs);

                if(tree) {
                    tree->add_message(m,Chirality::right);
                }
            }
        }
    }

private:

    MCF_FACTOR* mcf_factor = nullptr;
};

template<typename GRAPH_MATCHING_CONSTRUCTOR, typename INTER_QUADRATIC_MESSAGE_CONTAINER>
class graph_matching_inter_quadratic_message_constructor : public GRAPH_MATCHING_CONSTRUCTOR {
public:
    using FMC = typename GRAPH_MATCHING_CONSTRUCTOR::FMC;
    using inter_quadratic_message_container = INTER_QUADRATIC_MESSAGE_CONTAINER;

    using GRAPH_MATCHING_CONSTRUCTOR::GRAPH_MATCHING_CONSTRUCTOR;
    /*
    template<typename SOLVER>
        graph_matching_inter_quadratic_message_constructor(SOLVER& s)
        //: GRAPH_MATCHING_CONSTRUCTOR(&s.GetLP(), "both_sides", "mcf")
        : GRAPH_MATCHING_CONSTRUCTOR(s)
        {}

        graph_matching_inter_quadratic_message_constructor(LP<FMC>* lp, const std::string& construction_method)
        : GRAPH_MATCHING_CONSTRUCTOR(lp, construction_method)
        {}
        */

    template<typename PAIRWISE_FACTOR, typename INDICES_ITERATOR>
    inter_quadratic_message_container* add_inter_quadratic_message(PAIRWISE_FACTOR* l, PAIRWISE_FACTOR* r, 
          INDICES_ITERATOR left_indices_begin, INDICES_ITERATOR left_indices_end,
          INDICES_ITERATOR right_indices_begin, INDICES_ITERATOR right_indices_end)
    {
       return this->lp_->template add_message<inter_quadratic_message_container>(l, r, left_indices_begin, left_indices_end, right_indices_begin, right_indices_end); 
    }

    template<typename PAIRWISE_FACTOR_CONTAINER>
    inter_quadratic_message_container* add_interquadratic_message(
          const std::size_t left_node_1, const std::size_t left_node_2, 
          const std::size_t right_node_1, const std::size_t right_node_2,
          PAIRWISE_FACTOR_CONTAINER* l, PAIRWISE_FACTOR_CONTAINER* r 
          )
    {
       assert(left_node_1 < left_node_2);
       assert(right_node_1 < right_node_2);

       const bool l1r1 = this->has_edge(left_node_1, right_node_1);
       const bool l1r2 = this->has_edge(left_node_1, right_node_2);
       const bool l2r1 = this->has_edge(left_node_2, right_node_1);
       const bool l2r2 = this->has_edge(left_node_2, right_node_2);

       const std::size_t l1r1_index = l1r1 ? this->get_left_index(left_node_1, right_node_1) : std::numeric_limits<std::size_t>::max();
       const std::size_t r1l1_index = l1r1 ? this->get_right_index(right_node_1, left_node_1) : std::numeric_limits<std::size_t>::max();
       const std::size_t l1r2_index = l1r2 ? this->get_left_index(left_node_1, right_node_2) : std::numeric_limits<std::size_t>::max();
       const std::size_t r2l1_index = l1r2 ? this->get_right_index(right_node_2, left_node_1) : std::numeric_limits<std::size_t>::max();
       const std::size_t l2r1_index = l2r1 ? this->get_left_index(left_node_2, right_node_1) : std::numeric_limits<std::size_t>::max();
       const std::size_t r1l2_index = l2r1 ? this->get_right_index(right_node_1, left_node_2) : std::numeric_limits<std::size_t>::max();
       const std::size_t l2r2_index = l2r2 ? this->get_left_index(left_node_2, right_node_2) : std::numeric_limits<std::size_t>::max();
       const std::size_t r2l2_index = l2r2 ? this->get_right_index(right_node_2, left_node_2) : std::numeric_limits<std::size_t>::max();

       if(l1r1) { assert(l1r1_index < std::numeric_limits<std::size_t>::max() && r1l1_index < std::numeric_limits<std::size_t>::max()); }
       if(l1r2) { assert(l1r2_index < std::numeric_limits<std::size_t>::max() && r2l1_index < std::numeric_limits<std::size_t>::max()); }
       if(l2r1) { assert(l2r1_index < std::numeric_limits<std::size_t>::max() && r1l2_index < std::numeric_limits<std::size_t>::max()); }
       if(l2r2) { assert(l2r2_index < std::numeric_limits<std::size_t>::max() && r2l2_index < std::numeric_limits<std::size_t>::max()); }

       std::array<std::size_t,2> matching_l_indices = {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()};
       std::array<std::size_t,2> matching_r_indices = {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()}; 
       std::size_t matching_counter = 0;
       if(l1r1 && l2r2) {
          matching_l_indices[matching_counter] = l1r1_index * (l->get_factor()->dim2()) + l2r2_index;
          matching_r_indices[matching_counter++] = r1l1_index * (r->get_factor()->dim2()) + r2l2_index;
       }
       if(l1r2_index < std::numeric_limits<std::size_t>::max() && l2r1_index < std::numeric_limits<std::size_t>::max()) {
          matching_l_indices[matching_counter] = l1r2_index * (l->get_factor()->dim2()) + l2r1_index;
          matching_r_indices[matching_counter++] = r1l2_index * (r->get_factor()->dim2()) + r2l1_index; 
       }
       if(matching_counter == 2 && matching_l_indices[0] > matching_l_indices[1]) {
          std::swap(matching_l_indices[0], matching_l_indices[1]);
          std::swap(matching_r_indices[0], matching_r_indices[1]); 
       }

       if(matching_counter > 0) 
          return add_inter_quadratic_message(l, r, matching_l_indices.begin(), matching_l_indices.begin() + matching_counter, matching_r_indices.begin(), matching_r_indices.begin() + matching_counter);
       else
          return nullptr;
    }

    void construct(const graph_matching_input& gm_input)
    {
       GRAPH_MATCHING_CONSTRUCTOR::construct(gm_input);

       for(std::size_t left_pairwise_factor_id=0; left_pairwise_factor_id<this->left_mrf.get_number_of_pairwise_factors(); ++left_pairwise_factor_id) {
          auto* l = this->left_mrf.get_pairwise_factor(left_pairwise_factor_id);
          const auto [l_node_1,l_node_2] = this->left_mrf.get_pairwise_variables(left_pairwise_factor_id);
          assert(l_node_1 < l_node_2);
          for(std::size_t l_index_1=0; l_index_1<l->get_factor()->dim1()-1; ++l_index_1) {
             const std::size_t r_node_1 = this->graph_[l_node_1][l_index_1];
             for(std::size_t l_index_2=0; l_index_2<l->get_factor()->dim2()-1; ++l_index_2) {
                const std::size_t r_node_2 = this->graph_[l_node_2][l_index_2];
                if(r_node_1 == r_node_2) continue;
                const std::size_t min_r_node = std::min(r_node_1, r_node_2);
                const std::size_t max_r_node = std::max(r_node_1, r_node_2);
                if(this->right_mrf.has_pairwise_factor(min_r_node, max_r_node)) {
                   auto* r = this->right_mrf.get_pairwise_factor(min_r_node, max_r_node);

                   add_interquadratic_message(l_node_1, l_node_2, min_r_node, max_r_node, l, r);


/*
                   const std::size_t l1r1_index = this->get_left_index(l_node_1, r_node_1);
                   const std::size_t l1r2_index = this->get_left_index(l_node_1, r_node_2);
                   const std::size_t l2r1_index = this->get_left_index(l_node_2, r_node_1);
                   const std::size_t l2r2_index = this->get_left_index(l_node_2, r_node_2);

                   assert(l1r1_index != l1r2_index);
                   assert(l2r1_index != l2r2_index);

                   const std::size_t r1l1_index = this->get_right_index(r_node_1, l_node_1);
                   const std::size_t r1l2_index = this->get_right_index(r_node_1, l_node_2);
                   const std::size_t r2l1_index = this->get_right_index(r_node_2, l_node_1);
                   const std::size_t r2l2_index = this->get_right_index(r_node_2, l_node_2);

                   assert(r1l1_index != r1l2_index);
                   assert(r2l1_index != r2l2_index);

                   std::array<std::size_t,2> matching_l_indices = {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()};
                   std::array<std::size_t,2> matching_r_indices = {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()}; 
                   std::size_t cur_index = 0;

                   if(l1r1_index < std::numeric_limits<std::size_t>::max() && l2r2_index < std::numeric_limits<std::size_t>::max()) {
                      matching_l_indices[cur_index] = l1r1_index * (l->get_factor()->dim2()) + l2r2_index;
                      matching_r_indices[cur_index] = r1l1_index * (r->get_factor()->dim2()) + r2l2_index;
                      assert(matching_l_indices[cur_index] < l->get_factor()->dim1()*l->get_factor()->dim2());
                      assert(matching_r_indices[cur_index] < r->get_factor()->dim1()*r->get_factor()->dim2());
                      ++cur_index;
                   }
                   if(l1r2_index < std::numeric_limits<std::size_t>::max() && l2r1_index < std::numeric_limits<std::size_t>::max()) {
                      matching_l_indices[cur_index] = l1r2_index * (l->get_factor()->dim2()) + l2r1_index;
                      matching_r_indices[cur_index] = r1l2_index * (r->get_factor()->dim2()) + r2l1_index; 
                      assert(matching_l_indices[cur_index] < l->get_factor()->dim1()*l->get_factor()->dim2());
                      assert(matching_r_indices[cur_index] < r->get_factor()->dim1()*r->get_factor()->dim2());
                      ++cur_index;
                   }
                   assert(cur_index > 0);

                   add_inter_quadratic_message(l, r, matching_l_indices.begin(), matching_l_indices.begin() + cur_index, matching_r_indices.begin(), matching_r_indices.begin() + cur_index);
                   */
                }
             }
          }
       }

       /*
       for(std::size_t left_pairwise_factor_id=0; left_pairwise_factor_id<this->left_mrf.get_number_of_pairwise_factors(); ++left_pairwise_factor_id) {
          for(std::size_t right_pairwise_factor_id=0; right_pairwise_factor_id<this->right_mrf.get_number_of_pairwise_factors(); ++right_pairwise_factor_id) {
             // check for overlap between assignment pairs
             auto* l = this->left_mrf.get_pairwise_factor(left_pairwise_factor_id);
             auto* r = this->right_mrf.get_pairwise_factor(right_pairwise_factor_id);

             const auto [l1,l2] = this->left_mrf.get_pairwise_variables(left_pairwise_factor_id);
             assert(l1<l2);
             const auto [r1,r2] = this->right_mrf.get_pairwise_variables(right_pairwise_factor_id);
             assert(r1<r2);

             const std::size_t l1r1_index = this->get_left_index(l1,r1);
             const std::size_t l1r2_index = this->get_left_index(l1,r2);
             const std::size_t l2r1_index = this->get_left_index(l2,r1);
             const std::size_t l2r2_index = this->get_left_index(l2,r2);

             assert(l1r1_index != l1r2_index);
             assert(l2r1_index != l2r2_index);

             const std::size_t r1l1_index = this->get_right_index(r1,l1);
             const std::size_t r1l2_index = this->get_right_index(r1,l2);
             const std::size_t r2l1_index = this->get_right_index(r2,l1);
             const std::size_t r2l2_index = this->get_right_index(r2,l2);

             assert(r1l1_index != r1l2_index);
             assert(r2l1_index != r2l2_index);

             std::vector<std::size_t> left_indices;
             std::vector<std::size_t> right_indices;

             if(l1r1_index < std::numeric_limits<std::size_t>::max() && l2r2_index < std::numeric_limits<std::size_t>::max()) {
                left_indices.push_back(l1r1_index * (l->get_factor()->dim2()) + l2r2_index);
                right_indices.push_back(r1l1_index * (r->get_factor()->dim2()) + r2l2_index); 
                assert(left_indices.back() < l->get_factor()->dim1()*l->get_factor()->dim2());
                assert(right_indices.back() < r->get_factor()->dim1()*r->get_factor()->dim2());
             }
             if(l1r2_index < std::numeric_limits<std::size_t>::max() && l2r1_index < std::numeric_limits<std::size_t>::max()) {
                left_indices.push_back(l1r2_index * (l->get_factor()->dim2()) + l2r1_index);
                right_indices.push_back(r1l2_index * (r->get_factor()->dim2()) + r2l1_index); 
                assert(left_indices.back() < l->get_factor()->dim1()*l->get_factor()->dim2());
                assert(right_indices.back() < r->get_factor()->dim1()*r->get_factor()->dim2());
             }

             if(left_indices.size() == 2) {
                assert(left_indices[0] < left_indices[1]);
             }
             if(right_indices.size() == 2) {
                assert(right_indices[0] < right_indices[1]);
             }

             if(left_indices.size() > 0)
                add_inter_quadratic_message(l,r, left_indices, right_indices);
          }
       }
       */
    }

};

} // namespace LPMP
