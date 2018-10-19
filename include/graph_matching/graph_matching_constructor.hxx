#ifndef LPMP_GRAPH_MATCHING_CONSTRUCTOR_HXX
#define LPMP_GRAPH_MATCHING_CONSTRUCTOR_HXX

#include "graph_matching.h"
#include "mrf/mrf_problem_construction.hxx"
#include "min_cost_flow_factor_ssp.hxx"
#include "tree_decomposition.hxx"
#include "graph_matching_input.h"

namespace LPMP {

template<typename MRF_CONSTRUCTOR>
class graph_matching_mrf_constructor {
public:
    using FMC = typename MRF_CONSTRUCTOR::FMC;

    template<typename SOLVER>
    graph_matching_mrf_constructor(SOLVER& solver)
        : lp_(&solver.GetLP()),
        left_mrf(solver),
        right_mrf(solver) 
    {}

    graph_matching_mrf_constructor(LP<FMC>* lp)
        : lp_(lp),
        left_mrf(lp),
        right_mrf(lp)
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

    void ComputePrimal()
    {
       auto l = compute_primal_mcf_solution();
       read_in_labeling(l);
    }

    bool CheckPrimalConsistency() const
    {
        if(debug()) std::cout << "check graph matching assignment\n";

        std::vector<bool> labels_taken(inverse_graph_.size(),false);
        for(std::size_t i=0; i<left_mrf.get_number_of_variables(); ++i) {
            const auto state = left_mrf.get_unary_factor(i)->get_factor()->primal();
            if(state < graph_[i].size()) {
                const auto label = graph_[i][state];
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

    graph_matching_input export_graph_matching_input(const bool include_quadratic = true) const
    {
       if(debug())
          std::cout << "export linear assignment problem solution based on current reparametrization to graph matching problem\n";
       assert(include_quadratic == false);
       graph_matching_input output;

       const std::size_t no_left_nodes = left_mrf.get_number_of_variables();
       const std::size_t no_right_nodes = right_mrf.get_number_of_variables();

       for(std::size_t i=0; i<no_left_nodes; ++i) {
          auto& u = *left_mrf.get_unary_factor(i)->get_factor();
          for(std::size_t l=0; l<u.size()-1; ++l) {
             output.add_assignment(i, graph_[i][l], u[l] - u.back());
          }
       }

       for(std::size_t i=0; i<no_right_nodes; ++i) {
          auto& u = *right_mrf.get_unary_factor(i)->get_factor();
          for(std::size_t l=0; l<u.size()-1; ++l) {
             output.add_assignment(inverse_graph_[i][l], i, u[l] - u.back());
          }
       }

       output.normalize(); // merge parallel edges
       return output;
    }

    void read_in_labeling(const linear_assignment_problem_input::labeling& l)
    {
       assert(l.check_primal_consistency());
       if(l.size() != no_left_nodes())
          throw std::runtime_error("labeling has different number of entries than graph matching problem.");

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
       for(std::size_t i=0; i+1<graph[node_1].size(); ++i) {
          assert(graph[node_1][i] < graph[node_1][i+1]);
       }
       if(node_2 == std::numeric_limits<std::size_t>::max())
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
          assert(a.left_node < graph_.size());
          assert(a.right_node < inverse_graph_.size());

          graph_[a.left_node].push_back(a.right_node);
          inverse_graph_[a.right_node].push_back(a.left_node);
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

       for(const auto& a : input.assignments) {
           const std::size_t left_index = get_left_index(a.left_node, a.right_node);
           const std::size_t right_index = get_right_index(a.right_node, a.left_node);

           auto* left_factor = left_mrf.get_unary_factor(a.left_node);
           assert((*left_factor->get_factor())[left_index] == 0.0);
           (*left_factor->get_factor())[left_index] = 0.5*a.cost;

           auto* right_factor = right_mrf.get_unary_factor(a.right_node);
           assert((*right_factor->get_factor())[right_index] == 0.0);
           (*right_factor->get_factor())[right_index] = 0.5*a.cost;
       }
   }

   template<Chirality CHIRALITY>
   void construct_pairwise_factors(MRF_CONSTRUCTOR& mrf, const REAL scaling, std::vector<std::vector<std::size_t>>& graph, const graph_matching_input& input)
   {
       for(const auto& q : input.quadratic_terms) {
           const auto assignment_1 = input.assignments[q.assignment_1];
           const auto assignment_2 = input.assignments[q.assignment_2];

           const auto [node_1, node_2, label_1, label_2] = [&]()  {
               if(Chirality::left == CHIRALITY) {
                   if(assignment_1.left_node < assignment_2.left_node) { 
                       return std::make_tuple(assignment_1.left_node, assignment_2.left_node, assignment_1.right_node, assignment_2.right_node);
                   } else {
                       return std::make_tuple(assignment_2.left_node, assignment_1.left_node, assignment_2.right_node, assignment_1.right_node);
                   }
               } else {
                   assert(CHIRALITY == Chirality::right);
                   if(assignment_1.right_node < assignment_2.right_node) { 
                       return std::make_tuple(assignment_1.right_node, assignment_2.right_node, assignment_1.left_node, assignment_2.left_node);
                   } else {
                       return std::make_tuple(assignment_2.right_node, assignment_1.right_node, assignment_2.left_node, assignment_1.left_node);
                   }
               }
           }();

           assert(node_1 < node_2);

           if(!mrf.has_pairwise_factor(node_1, node_2)) {
               auto* f = mrf.add_empty_pairwise_factor(node_1, node_2);
               // add infinities on diagonal
               for(std::size_t i1=0; i1<graph[node_1].size(); ++i1) {
                   for(std::size_t i2=0; i2<graph[node_2].size(); ++i2) {
                       if(graph[node_1][i1] == graph[node_2][i2]) {
                           f->get_factor()->cost(i1, i2) = std::numeric_limits<REAL>::infinity();
                       }
                   }
               }


           } 

           auto* f = mrf.get_pairwise_factor(node_1, node_2);

           const auto index_1 = get_index(node_1, label_1, graph);
           const auto index_2 = get_index(node_2, label_2, graph);

           f->get_factor()->cost(index_1, index_2) = scaling*q.cost; 
       } 
   }
   void construct_pairwise_factors(const graph_matching_input& input)
   {
       construct_pairwise_factors<Chirality::left>(left_mrf, 0.5, graph_, input);
       construct_pairwise_factors<Chirality::right>(right_mrf, 0.5, inverse_graph_, input);
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

protected: 
   // TODO: remove and change mcf construction in mcf solver
   struct assignment {const std::size_t left; const std::size_t right; const REAL cost;};
   std::vector<assignment> assignments_;
   struct quadratic {const std::size_t assignment_1; const std::size_t assignment_2; const REAL cost;};
   std::vector<quadratic> quadratic_; 
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

} // namespace LPMP

#endif // LPMP_GRAPH_MATCHING_CONSTRUCTOR_HXX

