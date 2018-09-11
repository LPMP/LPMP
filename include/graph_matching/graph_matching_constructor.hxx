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

    void begin()
    {
       std::cout << "order factors\n";
       left_mrf.order_factors();
       right_mrf.order_factors();

       assert(left_mrf.get_number_of_variables() > 0);
       assert(right_mrf.get_number_of_variables() > 0);

       auto* f_l = left_mrf.get_unary_factor( left_mrf.get_number_of_variables()-1 );
       auto* f_r = right_mrf.get_unary_factor(0);

       lp_->add_factor_relation(f_l, f_r);
    }

    bool CheckPrimalConsistency() const
    {
        if(debug()) { std::cout << "check assignment\n"; }

        std::vector<bool> labels_taken(inverse_graph_.size(),false);
        for(std::size_t i=0; i<left_mrf.get_number_of_variables(); ++i) {
            const auto state = left_mrf.get_unary_factor(i)->get_factor()->primal();
            if(state < graph_[i].size()) {
                const auto label = graph_[i][state];
                if(labels_taken[label]) { 
                    if(debug()) { std::cout << "var " << i << ", state " << state << ", label " << label << " conflict\n"; }
                    return false; 
                }  
                labels_taken[ label ] = true;
            }  
        }  

        return true;
    }  


    void set_no_left_nodes(const std::size_t i)
    {
        graph_.resize(i);
    }
    void set_no_right_nodes(const std::size_t i)
    {
        inverse_graph_.resize(i); 
    }

    std::size_t add_assignment(const std::size_t left, const std::size_t right, const REAL cost)
    {
        assert(left < graph_.size());
        assert(right < inverse_graph_.size());

        graph_[left].push_back(right);
        inverse_graph_[right].push_back(left);

        assignments_.push_back({left,right,cost});
        return assignments_.size()-1;
    }

    void add_quadratic_cost(const std::size_t assignment_1, const std::size_t assignment_2, const REAL cost)
    {
       assert(assignment_1 != assignment_2);
       assert(assignment_1 < assignments_.size()); 
       assert(assignment_2 < assignments_.size()); 

       quadratic_.push_back({assignment_1, assignment_2, cost});
    }

    void construct()
    {
        for(const auto& v : graph_) { assert(std::is_sorted(v.begin(), v.end())); }
        for(const auto& v : inverse_graph_) { assert(std::is_sorted(v.begin(), v.end())); }
        construct_unary_factors();
        construct_pairwise_factors();
    }

    void read_input(const graph_matching_input& gm_input)
    {
        assert(std::is_sorted( gm_input.assignment_.begin(), gm_input.assignment_.end()));
        set_no_left_nodes(gm_input.no_left_nodes_);
        set_no_right_nodes(gm_input.no_right_nodes_);
        for(const auto& a : gm_input.assignment_) { add_assignment(a.left_node_, a.right_node_, a.cost_); }
        for(const auto& q : gm_input.quadratic_) { add_quadratic_cost(q.assignment_1, q.assignment_2, q.cost); }

        //constexpr PairwiseConstruction pc = FmcConstruction(typename SOLVER::FMC{});
    }

protected:
   void construct_empty_unary_factors(MRF_CONSTRUCTOR& mrf, const std::vector<std::vector<std::size_t>>& graph)
   {
       assert(mrf.get_number_of_variables() == 0);
       for(std::size_t i=0; i<graph.size(); ++i) {
           std::vector<REAL> costs(graph[i].size()+1, 0.0);
           mrf.add_unary_factor(costs);
       }
   }
   void construct_unary_factors()
   {
       construct_empty_unary_factors(left_mrf, graph_);
       construct_empty_unary_factors(right_mrf, inverse_graph_);

       std::vector<std::size_t> left_counter(graph_.size());
       std::vector<std::size_t> right_counter(inverse_graph_.size());
       for(const auto& assignment : assignments_) {
           const auto left_node = assignment.left;
           const auto right_node = assignment.right;
           const auto cost = assignment.cost;

           auto* left_factor = left_mrf.get_unary_factor(left_node);
           (*left_factor->get_factor())[left_counter[left_node]++] = 0.5*cost;

           auto* right_factor = right_mrf.get_unary_factor(right_node);
           (*right_factor->get_factor())[right_counter[right_node]++] = 0.5*cost;
       }
   }

   template<Chirality CHIRALITY>
   void construct_pairwise_factors(MRF_CONSTRUCTOR& mrf, const REAL scaling, std::vector<std::vector<std::size_t>>& graph)
   {
       for(const auto& q : quadratic_) {
           const auto assignment_1 = assignments_[q.assignment_1];
           const auto assignment_2 = assignments_[q.assignment_2];

           const auto [node_1, node_2, label_1, label_2] = [&]()  {
               if(Chirality::left == CHIRALITY) {
                   if(assignment_1.left < assignment_2.left) { 
                       return std::make_tuple(assignment_1.left, assignment_2.left, assignment_1.right, assignment_2.right);
                   } else {
                       return std::make_tuple(assignment_2.left, assignment_1.left, assignment_2.right, assignment_1.right);
                   }
               } else {
                   assert(CHIRALITY == Chirality::right);
                   if(assignment_1.right < assignment_2.right) { 
                       return std::make_tuple(assignment_1.right, assignment_2.right, assignment_1.left, assignment_2.left);
                   } else {
                       return std::make_tuple(assignment_2.right, assignment_1.right, assignment_2.left, assignment_1.left);
                   }
               }
           }();

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

           const auto index_1 = std::find(graph[node_1].begin(), graph[node_1].end(), label_1) - graph[node_1].begin();
           assert(label_1 < graph_[node_1].size());
           const auto index_2 = std::find(graph[node_2].begin(), graph[node_2].end(), label_2) - graph[node_2].begin();
           assert(label_2 < graph_[node_2].size());

           f->get_factor()->cost(index_1, index_2) = scaling*q.cost; 
       } 
   }
   void construct_pairwise_factors()
   {
       construct_pairwise_factors<Chirality::left>(left_mrf, 0.5, graph_);
       construct_pairwise_factors<Chirality::right>(right_mrf, 0.5, inverse_graph_);
   }

   LP<FMC>* lp_;
public:
   MRF_CONSTRUCTOR left_mrf, right_mrf;

   std::vector<std::vector<std::size_t>> graph_;
   std::vector<std::vector<std::size_t>> inverse_graph_;

protected: 
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

    void construct()
    {
        GRAPH_MATCHING_MRF_CONSTRUCTOR::construct();

        std::vector<std::size_t> right_label_counter(this->right_mrf.get_number_of_variables(), 0);
        for(std::size_t i=0; i<this->graph_.size(); ++i) {
            auto* l = this->left_mrf.get_unary_factor(i);
            for(std::size_t xi=0; xi<this->graph_[i].size(); ++xi) {
                const auto state = this->graph_[i][xi];
                auto* r = this->right_mrf.get_unary_factor(state);
                this->lp_->template add_message<ASSIGNMENT_MESSAGE>(l, r, xi, right_label_counter[state]);
                right_label_counter[state]++;
            }
        }
    }
};

template<typename GRAPH_MATCHING_CONSTRUCTOR, typename MCF_FACTOR, typename MCF_MESSAGE>
class graph_matching_mcf_constructor : public GRAPH_MATCHING_CONSTRUCTOR {
public:
    using FMC = typename GRAPH_MATCHING_CONSTRUCTOR::FMC;

    using GRAPH_MATCHING_CONSTRUCTOR::GRAPH_MATCHING_CONSTRUCTOR;

    void construct(factor_tree<FMC>* tree = nullptr)
    {
        GRAPH_MATCHING_CONSTRUCTOR::construct();

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

        auto* f = this->lp_->template add_factor<MCF_FACTOR>(edges, demands);

        // connect assignment factor with unaries
        // left side
        {
            std::vector<std::vector<std::size_t>> edgeId(no_left_nodes);
            for(std::size_t i=0; i<no_left_nodes; ++i) {
                //assert(mrf_left.GetNumberOfLabels(i) == mcf->NoArcs(i));
                auto *u = this->left_mrf.get_unary_factor(i);
                const auto first_arc = f->get_factor()->mcf_.first_outgoing_arc(i);
                const auto no_arcs = f->get_factor()->mcf_.no_outgoing_arcs(i);
                auto* m = this->lp_->template add_message<MCF_MESSAGE>(u, f, first_arc, no_arcs);

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
                const auto first_arc = f->get_factor()->mcf_.first_outgoing_arc(no_left_nodes+i);
                const auto no_arcs = f->get_factor()->mcf_.no_outgoing_arcs(no_left_nodes+i);
                auto* m = this->lp_->template add_message<MCF_MESSAGE>(u, f, first_arc, no_arcs);

                if(tree) {
                    tree->add_message(m,Chirality::right);
                }
            }
        }
    }
};

template<typename GRAPH_MATCHING_CONSTRUCTOR>
class graph_matching_constructor_tightening : public GRAPH_MATCHING_CONSTRUCTOR {
public:
    using FMC = typename GRAPH_MATCHING_CONSTRUCTOR::FMC;

    using GRAPH_MATCHING_CONSTRUCTOR::GRAPH_MATCHING_CONSTRUCTOR;

    INDEX Tighten(const INDEX no_constraints_to_add)
    {
        return this->left_mrf.Tighten(no_constraints_to_add+1) + this->right_mrf.Tighten(no_constraints_to_add+1);
    }
};

} // namespace LPMP

#endif // LPMP_GRAPH_MATCHING_CONSTRUCTOR_HXX

