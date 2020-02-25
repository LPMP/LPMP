#pragma once

#include "graph_matching/graph_matching_constructor.hxx"
#include "multigraph_matching_input.h"
#include <unordered_map>
#include <memory>
#include <variant>
#include "LP.h"
#include "config.hxx"
#include <omp.h>
#include <atomic>
#include <future>
#include "multicut/multicut_instance.h"
#include "multicut/multicut_kernighan_lin.h"
#include "multicut/transform_multigraph_matching.h"
#include "multigraph_matching_synchronization.h"

namespace LPMP {

// works for three matching problems forming a triangle
template<typename GRAPH_MATCHING_CONSTRUCTOR,
   typename TRIPLET_CONSISTENCY_FACTOR, typename TRIPLET_CONSISTENCY_FACTOR_ZERO, 
   typename PQ_ROW_TRIPLET_CONSISTENCY_MESSAGE, typename QR_COLUMN_TRIPLET_CONSISTENCY_MESSAGE, typename PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE,
   typename PQ_ROW_TRIPLET_CONSISTENCY_ZERO_MESSAGE, typename QR_COLUMN_TRIPLET_CONSISTENCY_ZERO_MESSAGE
   >
class multigraph_matching_constructor {

   using ptr_to_triplet_consistency_factor = std::variant<TRIPLET_CONSISTENCY_FACTOR*, TRIPLET_CONSISTENCY_FACTOR_ZERO*>;

    struct triplet_consistency_factor {
        std::size_t p,q,r; // graphs
        std::size_t p_node, r_node;

        bool operator==(const triplet_consistency_factor& o) const 
        {
            return p == o.p && q == o.q && r == o.r && p_node == o.p_node && r_node == o.r_node;
        }
    };
    struct triplet_consistency_factor_hash {
        std::size_t operator()(const triplet_consistency_factor& t) const
        {
            assert(t.p < t.r);
            size_t hash = std::hash<std::size_t>()(t.p);
            hash = hash::hash_combine(hash, std::hash<std::size_t>()(t.q));
            hash = hash::hash_combine(hash, std::hash<std::size_t>()(t.r));
            hash = hash::hash_combine(hash, std::hash<std::size_t>()(t.p_node));
            hash = hash::hash_combine(hash, std::hash<std::size_t>()(t.r_node));
            return hash; 
        }
    };

    // this struct caches constructors of given graph matching triplet, so they need not be looked up repeatedly in triplet enumeration
    // given a triplet_consistency_matching_factor, it retrieves data needed for constructing consistency factors from the graph matching prolems.
    struct graph_matching_triplet {
       std::size_t p,q,r; // graphs
       GRAPH_MATCHING_CONSTRUCTOR* pr_constructor = nullptr;
       GRAPH_MATCHING_CONSTRUCTOR* pq_constructor = nullptr;
       GRAPH_MATCHING_CONSTRUCTOR* qr_constructor = nullptr;

        bool triplet_consistency_factor_matches(const triplet_consistency_factor& t) const
        {
           std::array<std::size_t,3> graphs {p,q,r};
           std::sort(graphs.begin(), graphs.end());

           std::array<std::size_t,3> t_graphs {t.p,t.q,t.r};
           std::sort(t_graphs.begin(), t_graphs.end());

           return graphs == t_graphs;
        }

        auto get_pr_factor_left(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           auto* pr_constructor = get_pr_constructor(t);
           auto* f_pr_left = pr_constructor->left_mrf.get_unary_factor(t.p_node);
           const auto& p_matching_labels = pr_constructor->graph_[t.p_node];
           const auto p_index = std::find(p_matching_labels.begin(), p_matching_labels.end(), t.r_node) - p_matching_labels.begin(); 

           return std::make_tuple(f_pr_left, p_index);
        }
        auto get_pr_factor_right(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           auto* pr_constructor = get_pr_constructor(t);
           auto* f_pr_right = pr_constructor->right_mrf.get_unary_factor(t.r_node);
           const auto& r_matching_labels = pr_constructor->inverse_graph_[t.r_node];
           const auto r_index = std::find(r_matching_labels.begin(), r_matching_labels.end(), t.p_node) - r_matching_labels.begin();
           return std::make_tuple(f_pr_right, r_index); 
        }

        auto* get_pq_factor(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           auto* pq_constructor = get_pq_constructor(t);
           if(t.p < t.q)
              return pq_constructor->left_mrf.get_unary_factor(t.p_node);
           else
              return pq_constructor->right_mrf.get_unary_factor(t.p_node);
        }
        const auto& get_pq_factor_labels(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           auto* pq_constructor = get_pq_constructor(t);
           auto* f_pq = get_pq_factor(t);
           if(t.p < t.q)
              return pq_constructor->graph_[t.p_node]; 
           else
              return pq_constructor->inverse_graph_[t.p_node]; 
        }

        auto* get_qr_factor(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           auto* qr_constructor = get_qr_constructor(t);
           if(t.q < t.r)
              return qr_constructor->right_mrf.get_unary_factor(t.r_node);
           else
              return qr_constructor->left_mrf.get_unary_factor(t.r_node);
        }
        const auto& get_qr_factor_labels(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           auto* qr_constructor = get_qr_constructor(t);
           auto* f_qr = get_qr_factor(t);
           // TODO: correct way around?
           if(t.q < t.r)
              return qr_constructor->inverse_graph_[t.r_node]; 
           else
              return qr_constructor->graph_[t.r_node]; 
        }

        auto get_pq_qr_indices(const triplet_consistency_factor& t) const
        {
           // TODO: remove
           assert(triplet_consistency_factor_matches(t));

           auto* f_pq = get_pq_factor(t);
           const auto& pq_labels = get_pq_factor_labels(t);
           auto* f_qr = get_qr_factor(t);
           const auto& qr_labels = get_qr_factor_labels(t);

           // record which entries in f_pq and f_qr point towards the same nodes in q
           std::vector<std::size_t> common_nodes;
           std::set_intersection(pq_labels.begin(), pq_labels.end(), qr_labels.begin(), qr_labels.end(), std::back_inserter(common_nodes));

           std::vector<std::size_t> pq_factor_indices;
           pq_factor_indices.reserve(common_nodes.size());
           auto pq_iterator = pq_labels.begin();

           std::vector<std::size_t> qr_factor_indices;
           qr_factor_indices.reserve(common_nodes.size());
           auto qr_iterator = qr_labels.begin();

           for(const auto common_node : common_nodes) {
              pq_iterator = std::find(pq_iterator, pq_labels.end(), common_node);
              assert(pq_iterator != pq_labels.end());
              const auto pq_idx = pq_iterator - pq_labels.begin();
              pq_factor_indices.push_back(pq_idx);

              qr_iterator = std::find(qr_iterator, qr_labels.end(), common_node);
              assert(qr_iterator != qr_labels.end());
              const auto qr_idx = qr_iterator - qr_labels.begin();
              qr_factor_indices.push_back(qr_idx); 
           }

           return std::make_tuple(std::move(pq_factor_indices), std::move(qr_factor_indices));
        }

        GRAPH_MATCHING_CONSTRUCTOR* get_pr_constructor(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));
           if(t.p == p && t.r == r) return pr_constructor;
           if(t.p == p && t.r == q) return pq_constructor;
           if(t.p == q && t.r == r) return qr_constructor;
           assert(false);
           return nullptr;
        }

        GRAPH_MATCHING_CONSTRUCTOR* get_pq_constructor(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));

           if(t.p == p && t.q == q) return pq_constructor;
           if(t.p == q && t.q == p) return pq_constructor;

           if(t.p == p && t.q == r) return pr_constructor;
           if(t.p == r && t.q == p) return pr_constructor;

           if(t.p == q && t.q == r) return qr_constructor;
           if(t.p == r && t.q == q) return qr_constructor;

           assert(false);
           return nullptr;
        }

        GRAPH_MATCHING_CONSTRUCTOR* get_qr_constructor(const triplet_consistency_factor& t) const
        {
           assert(triplet_consistency_factor_matches(t));

           if(t.q == p && t.r == q) return pq_constructor;
           if(t.q == q && t.r == p) return pq_constructor;

           if(t.q == p && t.r == r) return pr_constructor;
           if(t.q == r && t.r == p) return pr_constructor;

           if(t.q == q && t.r == r) return qr_constructor;
           if(t.q == r && t.r == q) return qr_constructor;

           assert(false);
           return nullptr;
        }
    };

    graph_matching_triplet get_graph_matching_triplet(const std::size_t p, const std::size_t q, const std::size_t r) const
    {
       assert(p<q && q<r && r<no_graphs());
       graph_matching_triplet t;
       t.p = p;
       t.q = q;
       t.r = r;
       t.pr_constructor = get_graph_matching_constructor(p,r);
       t.pq_constructor = get_graph_matching_constructor(p,q);
       t.qr_constructor = get_graph_matching_constructor(q,r);
       return t; 
    }

    struct graph_matching {
        std::size_t p; // left graph
        std::size_t q; // right graph

        bool operator==(const graph_matching& o) const 
        {
            return p == o.p && q == o.q;
        } 

        bool operator<(const graph_matching& o) const
        {
           std::array<std::size_t,2> idx1 {p,q};
           std::array<std::size_t,2> idx2 {o.p,o.q};
           return std::lexicographical_compare(idx1.begin(), idx1.end(), idx2.begin(), idx2.end());
        }
    };
    struct graph_matching_hash  {
        std::size_t operator()(const graph_matching& g) const
        {
            assert(g.p<g.q);
            return hash::hash_combine(std::hash<std::size_t>()(g.p), std::hash<std::size_t>()(g.q));
        }
    };

public:
    using FMC = typename GRAPH_MATCHING_CONSTRUCTOR::FMC;

    template<typename SOLVER>
        multigraph_matching_constructor(SOLVER& s)
        : lp_(&s.GetLP()),
        presolve_iterations_arg_("", "presolveIterations", "number of iterations for presolving pairwise graph matching instances", false, 0, &positiveIntegerConstraint_, s.get_cmd()), 
        presolve_reparametrization_arg_("", "presolveReparametrization", "mode of reparametrization for presolving pairwise graph matching instances", false, "anisotropic", "{anisotropic|uniform}:${leave_percentage}", s.get_cmd()),
        primal_checking_triplets_arg_("", "primalCheckingTriplets", "number of triplet consistency factors to include during primal feasibility check", false, 0, &positiveIntegerConstraint_, s.get_cmd()),
        mcf_reparametrization_arg_("", "mcfReparametrization", "enable reparametrization by solving a linear assignment problem with a minimimum cost flow solver", s.get_cmd(), true),
        mcf_primal_rounding_arg_("", "mcfRounding", "enable runding by solving a linear assignment problem with a minimimum cost flow solver", s.get_cmd(), true),
        graph_matching_construction_arg_("", "graphMatchingConstruction", "mode of constructing pairwise potentials for graph matching", false, "both_sides", "{left|right|both_sides}", s.get_cmd()),
        output_format_arg_("", "multigraphMatchingOutputFormat", "output format for multigraph matching", false, "matching", "{matching|clustering}", s.get_cmd()),
        primal_rounding_algorithms_arg_("", "multigraphMatchingRoundingMethod", "correlation clustering algorithms that are used for rounding a solution", false, "gaec_KL", "{gaec_KL|KL|MCF_KL|MCF_PS}", s.get_cmd())
        
        {
           omp_set_nested(0); // there is parallelism in the graph matching constructors, which we suppress hereby. TODO: should be better set somewhere else
        }

    void order_factors()
    {
       std::vector<graph_matching> graph_matching_indices;

       assert(graph_matching_constructors.size() > 0);
       for(auto& c : graph_matching_constructors) {
          c.second->order_factors();
          graph_matching_indices.push_back(c.first); 
       }

       // order graph matching problems lexicographically
       std::sort(graph_matching_indices.begin(), graph_matching_indices.end());
       auto* graph_matching_constructor_prev = get_graph_matching_constructor(graph_matching_indices[0].p, graph_matching_indices[0].q);
       for(std::size_t i=1; i<graph_matching_indices.size(); ++i) {
          auto* graph_matching_constructor_next = get_graph_matching_constructor(graph_matching_indices[i].p, graph_matching_indices[i].q);
          auto* f_prev = graph_matching_constructor_prev->right_mrf.get_unary_factor( graph_matching_constructor_prev->right_mrf.get_number_of_variables()-1 ); 
          auto* f_next = graph_matching_constructor_next->left_mrf.get_unary_factor(0);
          lp_->add_factor_relation(f_prev, f_next);

          graph_matching_constructor_prev = graph_matching_constructor_next;
       }
    }

    // preoptimize in parallel pairwise graph matching problems
    void begin()
    {
       const std::size_t presolve_iterations = presolve_iterations_arg_.getValue(); 
       if(presolve_iterations == 0) return; 

       if(debug()) { std::cout << "presolve pairwise graph matching problems\n"; }

#pragma omp parallel for schedule(dynamic)
       for(std::size_t i=0; i<graph_matching_constructors.size(); ++i) {
          auto* c = graph_matching_constructors[i].second.get();
          auto gm_factors = c->get_factors();
          std::vector<FactorTypeAdapter*> gm_update_factors;
          std::copy_if(gm_factors.begin(), gm_factors.end(), std::back_inserter(gm_update_factors), [](FactorTypeAdapter* f) { return f->FactorUpdated(); } );

          auto [omega_forward, receive_mask_forward] = compute_anisotropic_weights(gm_factors.begin(), gm_factors.end(), 0.0);
          auto [omega_backward, receive_mask_backward] = compute_anisotropic_weights(gm_factors.rbegin(), gm_factors.rend(), 0.0);
          for(std::size_t iter=0; iter<presolve_iterations; ++iter) {
             c->reparametrize_linear_assignment_problem();
             compute_pass(gm_update_factors.begin(), gm_update_factors.end(), omega_forward.begin(), receive_mask_forward.begin());
             compute_pass(gm_update_factors.rbegin(), gm_update_factors.rend(), omega_backward.begin(), receive_mask_backward.begin()); 
          } 
       } 
    }

    void End()
    {
       // read out last mgm solution
       if(primal_result_handle_.valid()) {
          primal_result_handle_.wait();
          if(debug()) 
             std::cout << "read in primal mgm solution\n";
          auto mgm_sol = primal_result_handle_.get();
          read_in_labeling(mgm_sol); 
       }
    }

    void pre_iterate()
    {
#pragma omp parallel for schedule(guided)
       for(std::size_t i=0; i<graph_matching_constructors.size(); ++i) {
          graph_matching_constructors[i].second->pre_iterate();
       }
    }

    GRAPH_MATCHING_CONSTRUCTOR* add_graph_matching_problem(const std::size_t p, const std::size_t q, LP<FMC>* lp)
    {
        assert(!has_graph_matching_problem(p,q));
        auto ptr = std::make_unique<GRAPH_MATCHING_CONSTRUCTOR>(lp, graph_matching_construction_arg_.getValue(), "mcf");
        auto* c = ptr.get();
        graph_matching_constructors.push_back(std::make_pair(graph_matching({p,q}), std::move(ptr))); 
        graph_matching_constructors_map.insert(std::make_pair(graph_matching({p,q}), c)); 
        no_graphs_ = std::max(no_graphs_, q+1);
        return c;
    }

    bool has_graph_matching_problem(const std::size_t p, const std::size_t q) const
    {
        assert(p<q);
        return graph_matching_constructors_map.count({p,q}) > 0;
    }

    GRAPH_MATCHING_CONSTRUCTOR* get_graph_matching_constructor(const std::size_t p, const std::size_t q) const
    {
        assert(has_graph_matching_problem(p,q));
        return graph_matching_constructors_map.find(graph_matching({p,q}))->second;
    }

    bool has_triplet_consistency_factor(const triplet_consistency_factor& t) const
    {
        assert(t.p < t.r);
        return triplet_consistency_factors.count(t) > 0; 
    }

    ptr_to_triplet_consistency_factor add_triplet_consistency_factor(const triplet_consistency_factor& t)
    {
       assert(!has_triplet_consistency_factor(t));
       auto* pr_c = get_graph_matching_constructor(t.p, t.r);

       auto* f_pq = get_matching_factor(t.p, t.q, t.p_node);
       auto* f_qr = get_matching_factor(t.r, t.q, t.r_node);

       const auto pq_labels = get_triplet_consistency_labels(t.p, t.q, t.p_node);
       assert(pq_labels.size() == f_pq->get_factor()->size()-1);
       const auto qr_labels = get_triplet_consistency_labels(t.r, t.q, t.r_node);
       assert(qr_labels.size() == f_qr->get_factor()->size()-1);

       if(pr_c->has_edge(t.p_node, t.r_node)) { 
          auto* f_pr_left = get_matching_factor(t.p, t.r, t.p_node);
          auto* f_pr_right = get_matching_factor(t.r, t.p, t.r_node);

          auto* f_triplet = lp_->template add_factor<TRIPLET_CONSISTENCY_FACTOR>(pq_labels, qr_labels);
          triplet_consistency_factors.insert(std::make_pair(t, f_triplet));

          lp_->template add_message<PQ_ROW_TRIPLET_CONSISTENCY_MESSAGE>(f_pq, f_triplet);
          lp_->template add_message<QR_COLUMN_TRIPLET_CONSISTENCY_MESSAGE>(f_qr, f_triplet);

          const std::size_t p_index = get_matching_index(t.p, t.r, t.p_node, t.r_node);
          assert(p_index < f_pr_left->get_factor()->size()-1);
          const std::size_t r_index = get_matching_index(t.r, t.p, t.r_node, t.p_node);
          assert(r_index < f_pr_right->get_factor()->size()-1);

          lp_->template add_message<PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE>(f_pr_left, f_triplet, p_index);
          lp_->template add_message<PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE>(f_pr_right, f_triplet, r_index); 

          return f_triplet;
       } else {
          auto* f_triplet = lp_->template add_factor<TRIPLET_CONSISTENCY_FACTOR_ZERO>(pq_labels, qr_labels);
          triplet_consistency_factors.insert(std::make_pair(t, f_triplet));

          lp_->template add_message<PQ_ROW_TRIPLET_CONSISTENCY_ZERO_MESSAGE>(f_pq, f_triplet);
          lp_->template add_message<QR_COLUMN_TRIPLET_CONSISTENCY_ZERO_MESSAGE>(f_qr, f_triplet);

          return f_triplet;
       }
    }

    // (i) compute lower bound before reparametrizing
    // (ii) add triplet consistency factor
    // (iii) reparametrize, i.e. send messages to triplet consistency factor
    // (iv) compute new lower bound
    // (v) restore previous factor state
    double triplet_consistency_dual_increase(const triplet_consistency_factor& t, const graph_matching_triplet& gm_t) const
    {
       auto* f_pq = gm_t.get_pq_factor(t);
       auto& f_pq_factor = *f_pq->get_factor();
       const auto& pq_labels = gm_t.get_pq_factor_labels(t);

       auto* f_qr = gm_t.get_qr_factor(t);
       auto& f_qr_factor = *f_qr->get_factor();
       const auto& qr_labels = gm_t.get_qr_factor_labels(t);

       simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left> m_pq;
       simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right> m_qr;

       auto* pr_c = gm_t.get_pr_constructor(t);

       if(pr_c->has_edge(t.p_node, t.r_node)) { // check increase for triplet consistency factor
          auto [f_pr_left, p_index] = gm_t.get_pr_factor_left(t);
          auto& f_pr_left_factor = *f_pr_left->get_factor();

          auto [f_pr_right, r_index] = gm_t.get_pr_factor_right(t);
          auto& f_pr_right_factor = *f_pr_right->get_factor();

          const double prev_lb = f_pq_factor.LowerBound() + f_qr_factor.LowerBound() + f_pr_left_factor.LowerBound() + f_pr_right_factor.LowerBound();

          multigraph_matching_triplet_consistency_factor f(pq_labels, qr_labels);

          simplex_multigraph_matching_triplet_scalar_consistency_message m_pr_left(p_index);
          simplex_multigraph_matching_triplet_scalar_consistency_message m_pr_right(r_index);

          vector<double> msg_val_pq(f_pq_factor.size()-1, 0.0);
          m_pq.send_message_to_right(f_pq_factor, msg_val_pq, 1.0);
          m_pq.RepamRight(f, -1.0*msg_val_pq);
          m_pq.RepamLeft(f_pq_factor, msg_val_pq);

          vector<double> msg_val_qr(f_qr_factor.size()-1, 0.0);
          m_qr.send_message_to_right(f_qr_factor, msg_val_qr, 1.0);
          m_qr.RepamRight(f, -1.0*msg_val_qr);
          m_qr.RepamLeft(f_qr_factor, msg_val_qr);

          array<double,1> msg_val_pr_left = {0.0};
          m_pr_left.send_message_to_right(f_pr_left_factor, msg_val_pr_left, 1.0);
          m_pr_left.RepamRight(f, -msg_val_pr_left[0], 0);
          m_pr_left.RepamLeft(f_pr_left_factor, msg_val_pr_left[0], 0);

          array<double,1> msg_val_pr_right = {0.0};
          m_pr_right.send_message_to_right(f_pr_right_factor, msg_val_pr_right, 1.0);
          m_pr_right.RepamRight(f, -msg_val_pr_right[0], 0);
          m_pr_right.RepamLeft(f_pr_right_factor, msg_val_pr_right[0], 0);

          const double after_lb = f.LowerBound() + f_pq_factor.LowerBound() + f_qr_factor.LowerBound() + f_pr_left_factor.LowerBound() + f_pr_right_factor.LowerBound();

          // revert changes
          m_pq.RepamLeft(f_pq_factor, -1.0*msg_val_pq);
          m_qr.RepamLeft(f_qr_factor, -1.0*msg_val_qr);
          m_pr_left.RepamLeft(f_pr_left_factor, -msg_val_pr_left[0], 0);
          m_pr_right.RepamLeft(f_pr_right_factor, -msg_val_pr_right[0], 0);

          assert(std::abs(prev_lb - (f_pq->LowerBound() + f_qr->LowerBound() + f_pr_left->LowerBound() + f_pr_right->LowerBound())) <= eps);
          assert(after_lb >= prev_lb - eps);

          return after_lb - prev_lb; 

       } else { // check increase for zero factor

          const double prev_lb = f_pq_factor.LowerBound() + f_qr_factor.LowerBound();

          multigraph_matching_triplet_consistency_factor_zero f(pq_labels, qr_labels);

          simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left> m_pq;
          simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right> m_qr;

          vector<double> msg_val_pq(f_pq_factor.size()-1, 0.0);
          m_pq.send_message_to_right(f_pq_factor, msg_val_pq, 1.0);
          m_pq.RepamRight(f, -1.0*msg_val_pq);
          m_pq.RepamLeft(f_pq_factor, msg_val_pq);

          vector<double> msg_val_qr(f_qr_factor.size()-1, 0.0);
          m_qr.send_message_to_right(f_qr_factor, msg_val_qr, 1.0);
          m_qr.RepamRight(f, -1.0*msg_val_qr);
          m_qr.RepamLeft(f_qr_factor, msg_val_qr);

          const double after_lb = f.LowerBound() + f_pq_factor.LowerBound() + f_qr_factor.LowerBound();

          m_pq.RepamLeft(f_pq_factor, -1.0*msg_val_pq);
          m_qr.RepamLeft(f_qr_factor, -1.0*msg_val_qr);

          assert(std::abs(prev_lb - (f_pq->LowerBound() + f_qr->LowerBound())) <= eps);
          assert(after_lb >= prev_lb - eps);

          return after_lb - prev_lb; 
       } 
    }

    // enumerate all graph matching triplets
    // process graph matching triplets such that pairwise graph matching problems do not overlap
    template<typename FUNC>
    void for_each_triplet_consistency_factor(FUNC&& func) const
    {
       // enumerate all graphs
       std::deque<std::array<std::size_t,3>> gms;
       const auto n = no_graphs();
       for(std::size_t r=0; r<n; ++r) {
          for(std::size_t q=0; q<r; ++q) {
             for(std::size_t p=0; p<q; ++p) {
                gms.push_back({p,q,r});
             }
          }
       }

       std::mutex matchings_processed_mutex;
       // TODO: using stack allocator is not admissible here: use std::vector instead
       matrix<std::size_t> matchings_currently_processed(no_graphs(), no_graphs(), 0);
       auto get_graph_matching_problems = [&]() -> std::array<std::size_t,3> {
          std::lock_guard<std::mutex> lock(matchings_processed_mutex); 
          std::size_t no_tries = 0;
          while(!gms.empty()) {
             const auto [p,q,r] = gms.front();
             gms.pop_front();
             std::array<std::size_t,3> indices {p,q,r};
             std::sort(indices.begin(), indices.end());
             if(! (matchings_currently_processed(indices[0], indices[1]) == 0 && matchings_currently_processed(indices[0], indices[2]) == 0 && matchings_currently_processed(indices[1], indices[2]) == 0) ) {
                gms.push_back({p,q,r});
                ++no_tries;
                if(no_tries >= gms.size())
                   break;
                else
                   continue;
             }
             matchings_currently_processed(indices[0], indices[1]) = 1;
             matchings_currently_processed(indices[0], indices[2]) = 1;
             matchings_currently_processed(indices[1], indices[2]) = 1;
             return {p,q,r};
          }
          return {0,0,0};
       };

       auto release_graph_matching_problems = [&](const std::size_t p, const std::size_t q, const std::size_t r) {
          std::lock_guard<std::mutex> lock(matchings_processed_mutex);
          std::array<std::size_t,3> indices {p,q,r};
          std::sort(indices.begin(), indices.end());
          assert(matchings_currently_processed(indices[0], indices[1]) == 1);
          assert(matchings_currently_processed(indices[0], indices[2]) == 1);
          assert(matchings_currently_processed(indices[1], indices[2]) == 1);
          matchings_currently_processed(indices[0], indices[1]) = 0;
          matchings_currently_processed(indices[0], indices[2]) = 0;
          matchings_currently_processed(indices[1], indices[2]) = 0;
       };

#pragma omp parallel shared(gms, matchings_currently_processed, matchings_processed_mutex)
       {
          while(!gms.empty()) {
             const auto [p,q,r] = get_graph_matching_problems();
             if(p == 0 && q == 0 && r == 0) break;

             graph_matching_triplet gm_t = get_graph_matching_triplet(p,q,r);
             const auto p_no_nodes = gm_t.pq_constructor->left_mrf.get_number_of_variables();
             const auto q_no_nodes = gm_t.pq_constructor->right_mrf.get_number_of_variables();
             const auto r_no_nodes = gm_t.pr_constructor->right_mrf.get_number_of_variables();

             // enumerate all triplet_consistency_factor for given triplet of graphs
             triplet_consistency_factor t;

             t.p = p; t.q = q; t.r = r;
             for(std::size_t p_node=0; p_node<p_no_nodes; ++p_node) {
                for(std::size_t r_node=0; r_node<r_no_nodes; ++r_node) {
                   t.p_node = p_node; t.r_node = r_node;
                   func(t, gm_t);
                }
             }

             t.p = p; t.q = r; t.r = q;
             for(std::size_t p_node=0; p_node<p_no_nodes; ++p_node) {
                for(std::size_t q_node=0; q_node<q_no_nodes; ++q_node) {
                   t.p_node = p_node; t.r_node = q_node;
                   func(t, gm_t);
                }
             }

             t.p = q; t.q = p; t.r = r;
             for(std::size_t q_node=0; q_node<q_no_nodes; ++q_node) {
                for(std::size_t r_node=0; r_node<r_no_nodes; ++r_node) {
                   t.p_node = q_node; t.r_node = r_node;
                   func(t, gm_t);
                }
             }

             release_graph_matching_problems(p,q,r);
          }
       }
    }

    INDEX Tighten(const INDEX no_constraints_to_add)
    {
        // iterate over all triplets of graphs and enumerate all possible triplet consistency factors that can be added. 
        // Record guaranteed dual increase of adding the triplet consistency factor
        std::vector< std::vector<std::pair<triplet_consistency_factor, double>> > triplet_consistency_candidates_local(omp_get_max_threads());

        auto compute_dual_increase = [&](const triplet_consistency_factor& t, const graph_matching_triplet& gm_t) {
           const double guaranteed_dual_increase = this->triplet_consistency_dual_increase(t, gm_t); 
           if(guaranteed_dual_increase >= eps) {
              triplet_consistency_candidates_local[omp_get_thread_num()].push_back( std::make_pair(t, guaranteed_dual_increase) );
           } 
        };
        for_each_triplet_consistency_factor(compute_dual_increase);

        std::vector<std::pair<triplet_consistency_factor, double>> triplet_consistency_candidates = std::move(triplet_consistency_candidates_local[0]);
        for(std::size_t i=1; i<triplet_consistency_candidates_local.size(); ++i)
           triplet_consistency_candidates.insert(triplet_consistency_candidates.end(), triplet_consistency_candidates_local[i].begin(), triplet_consistency_candidates_local[i].end());

        std::sort(triplet_consistency_candidates.begin(), triplet_consistency_candidates.end(), [](const auto& t1, const auto& t2) { return t1.second > t2.second; });

        std::size_t no_constraints_added = 0;
        for(auto [t,cost] : triplet_consistency_candidates) {
            if(!has_triplet_consistency_factor(t)) {
                add_triplet_consistency_factor(t);
                no_constraints_added++;
                if(no_constraints_added >= no_constraints_to_add)
                    break;
            } 
        }

        if(diagnostics())
            std::cout << "Added " << no_constraints_added << " triplet consistency factor for multigraph matching\n";

        // tighten for graph matching constructors
        // TODO: do in parallel. But: constructors will add factors in parallel, leading to data corruption. On the other hand cycle search for mrfs is already parallelized
        const std::size_t no_constraints_to_add_per_gm = no_constraints_to_add / (((graph_matching_constructors.size() * ( graph_matching_constructors.size() -1))/2)) + 1;
        for(auto& c : graph_matching_constructors) {
           no_constraints_added += c.second->Tighten(no_constraints_to_add_per_gm);
        }

        return no_constraints_added;
    }

    void construct(const multigraph_matching_input& mgm_instance)
    {
        for(const auto& gm : mgm_instance) {
            auto* gm_constructor = add_graph_matching_problem(gm.left_graph_no, gm.right_graph_no, lp_);
            gm_constructor->construct(gm.gm_input);
        }
    }

    bool CheckPrimalConsistency() const
    {
       if(debug())
          std::cout << "check primal consistency in multigraph matching constructor\n";

       const auto labeling = write_out_labeling();
       return labeling.check_primal_consistency();

       /*
       std::atomic<bool> gm_feasible = true;
#pragma omp parallel for
       for(std::size_t i=0; i<graph_matching_constructors.size(); ++i) {
          if(gm_feasible && !graph_matching_constructors[i].second->CheckPrimalConsistency()) {
             gm_feasible = false;
          } 
       }

       if(debug())
          std::cout << "graph matching problem " << (gm_feasible ? "feasible" : "infeasible") << "\n";

       if(!gm_feasible) 
          return false;

       // check cycle consistency
       bool feasible = true;
       std::vector< std::vector<std::pair<triplet_consistency_factor, double>> > triplet_consistency_candidates_local(omp_get_max_threads());

       auto check_primal_consistency_func = [&](const triplet_consistency_factor& t, const graph_matching_triplet& gm_t) {
          auto [f_pr_left, p_index] = gm_t.get_pr_factor_left(t);
          auto [f_pr_right, r_index] = gm_t.get_pr_factor_right(t);

          auto* f_pq = gm_t.get_pq_factor(t);
          auto* f_qr = gm_t.get_qr_factor(t);

          const bool pr_match = (f_pr_left->get_factor()->primal() == p_index);
          assert(pr_match == (f_pr_right->get_factor()->primal() == r_index));
          for(std::size_t i=0; i<pq_indices.size(); ++i) {
             const bool pq_match = (f_pq->get_factor()->primal() == pq_indices[i]);
             const bool qr_match = (f_qr->get_factor()->primal() == qr_indices[i]);
             const std::size_t no_matches = std::size_t(pr_match) + std::size_t(pq_match) + std::size_t(qr_match);
             if(no_matches == 2) {
                feasible = false;
                if(!has_triplet_consistency_factor(t)) {
                   const double guaranteed_dual_increase = this->triplet_consistency_dual_increase(t, gm_t); 
                   triplet_consistency_candidates_local[omp_get_thread_num()].push_back({t, guaranteed_dual_increase});
                } 
             } 
          }
       };
       for_each_triplet_consistency_factor(check_primal_consistency_func);

       std::vector<std::pair<triplet_consistency_factor, double>> triplet_consistency_candidates = std::move(triplet_consistency_candidates_local[0]);
       for(std::size_t i=1; i<triplet_consistency_candidates_local.size(); ++i) {
          triplet_consistency_candidates.insert(triplet_consistency_candidates.end(), triplet_consistency_candidates_local[i].begin(), triplet_consistency_candidates_local[i].end());
       } 
       
       std::sort(triplet_consistency_candidates.begin(), triplet_consistency_candidates.end(), [](const auto& t1, const auto& t2) { return t1.second > t2.second; });
       std::size_t no_factors_added = 0;
       for(auto t : triplet_consistency_candidates) {
          if(no_factors_added >= primal_checking_triplets_arg_.getValue()) 
             break;
          if(!has_triplet_consistency_factor(t.first)) {
             add_triplet_consistency_factor(t.first);
             ++no_factors_added;
          }
       }

       if(!feasible && debug()) {
          std::cout << "multigraph matching constraints violated\n";
          if(primal_checking_triplets_arg_.getValue() > 0)
             std::cout << "Added " << no_factors_added << " triplet consistency factors with violated primal assignments\n";
       }

       return feasible;
       */
    }

    template<typename STREAM>
    void WritePrimal(STREAM& s) const
    {
       const multigraph_matching_input::labeling l = write_out_labeling();
       const std::string mode(output_format_arg_.getValue());

       if(std::string("matching") == mode)
          l.write_primal_matching(s);
       else if(std::string("clustering") == mode)
          l.write_primal_clustering(s);
       else
          throw std::runtime_error("output format must be {matching|clustering}");
    }

    void send_messages_to_unaries()
    {
       for(auto& c : graph_matching_constructors)
          c.second->send_messages_to_unaries();

       for(auto& t : triplet_consistency_factors)
          std::visit([&](auto&& t) { this->send_messages_to_unaries(t); }, t.second);
    }

    matrix<int> compute_allowed_matching_matrix() const
    {
       auto mgm = export_linear_multigraph_matching_input();
       multigraph_matching_input::graph_size mgm_size(mgm);

       matrix<int> m(mgm_size.total_no_nodes(), mgm_size.total_no_nodes(), 0);

       for(std::size_t i=0; i<m.dim1(); ++i) 
          for(std::size_t j=0; j<m.dim2(); ++j) 
             m(i,j) = 0;

       for(std::size_t p=0; p<mgm_size.no_graphs(); ++p) {
          for(std::size_t q=p+1; q<mgm_size.no_graphs(); ++q) {
             auto* c = get_graph_matching_constructor(p, q);
             const auto& graph = c->graph_;
             for(std::size_t i=0; i<graph.size(); ++i) {
                for(std::size_t j : graph[i]) {
                   m(mgm_size.node_no(p,i), mgm_size.node_no(q,j)) = 1;
                   m(mgm_size.node_no(q,j), mgm_size.node_no(p,i))  = 1;
                }
             }
          }
       } 

       return m;
    }

    // start with possibly inconsistent primal labeling obtained by individual graph matching roundings.
    // remote cycles that are inconsistent through a multicut solver
    void ComputePrimal()
    {
       if(!mcf_primal_rounding_arg_.getValue())
          return;

       if(primal_result_handle_.valid() && primal_result_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
          if(debug()) 
             std::cout << "read in primal mgm solution\n";
          auto mgm_sol = primal_result_handle_.get();
          read_in_labeling(mgm_sol); 
       }

       if(!primal_result_handle_.valid() || primal_result_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::deferred) {
          if(debug()) 
             std::cout << "construct mgm rounding problem\n";

          if(debug()) std::cout << "send messages to unaries\n";
          send_messages_to_unaries(); 

          enum class rounding_method {gaec_KL, KL, mcf_KL, mcf_ps};
          rounding_method rm = [&]() {
          if(primal_rounding_algorithms_arg_.getValue() == "gaec_KL") 
             return rounding_method::gaec_KL;
          else if(primal_rounding_algorithms_arg_.getValue() == "KL")
             return rounding_method::KL;
          else if(primal_rounding_algorithms_arg_.getValue() == "MCF_KL") // K&L on infeasible partial minimum cost flow matchings
             return rounding_method::mcf_KL;
          else if(primal_rounding_algorithms_arg_.getValue() == "MCF_PS") // permutation synchronization
             return rounding_method::mcf_ps;
          else
             throw std::runtime_error("could not recognize correlation clustering rounding method");
          }();

          multigraph_matching_input::labeling labeling_to_improve;
          if(rm == rounding_method::mcf_KL || rm == rounding_method::mcf_ps)
             for(auto& c : graph_matching_constructors)
                labeling_to_improve.push_back( {c.first.p, c.first.q, c.second->compute_primal_mcf_solution()} );

          // TODO: split round_primal_async based on rounding method into multiple lambdas
          auto mgm = std::make_shared<multigraph_matching_input>(export_linear_multigraph_matching_input());
          const auto allowed_matchings = compute_allowed_matching_matrix();
          auto round_primal_async = ([=](std::shared_ptr<multigraph_matching_input> mgm, multigraph_matching_input::labeling labeling_to_improve) {
                  multigraph_matching_correlation_clustering_transform mgm_cc_trafo(mgm);
                  auto cc = mgm_cc_trafo.get_correlatino_clustering_instance();
                  auto mc = cc.transform_to_multicut();
                  if(rm == rounding_method::gaec_KL) {
                     auto mc_sol = compute_multicut_gaec_kernighan_lin(mc); // seems better than just computing with Kernighan&Lin
                     auto cc_sol = mc_sol.transform_to_correlation_clustering();
                     auto mgm_sol = mgm_cc_trafo.transform(cc_sol); 
                     return mgm_sol;
                  } else if(rm == rounding_method::KL) {
                     auto mc_sol = compute_multicut_kernighan_lin(mc);
                     auto cc_sol = mc_sol.transform_to_correlation_clustering();
                     auto mgm_sol = mgm_cc_trafo.transform(cc_sol); 
                     return mgm_sol;
                  } else if(rm == rounding_method::mcf_KL) {
                     // try fixing infeasible matchings computed by individual graph matching solvers with Kernighan&Lin
                     auto cc_sol_to_improve = mgm_cc_trafo.transform(labeling_to_improve);
                     auto mc_sol_to_improve = cc_sol_to_improve.transform_to_multicut();
                     auto mc_sol = compute_multicut_kernighan_lin(mc, mc_sol_to_improve);
                     auto cc_sol = mc_sol.transform_to_correlation_clustering();
                     auto mgm_sol = mgm_cc_trafo.transform(cc_sol); 
                     return mgm_sol;
                  } else if(rm == rounding_method::mcf_ps) {
                     multigraph_matching_input::graph_size gs(*mgm);
                     synchronize_multigraph_matching(gs, labeling_to_improve, allowed_matchings); // for sparse assignment problems
                     //synchronize_multigraph_matching(gs, labeling_to_improve); // when all edges are present in pairwise matching subproblems
                     return labeling_to_improve; 
                  } else {
                     throw std::runtime_error("rounding method not supported");
                  }

                });
          primal_result_handle_ = std::async(std::launch::async, round_primal_async, mgm, labeling_to_improve);
       }
    }

    void read_in_labeling(const multigraph_matching_input::labeling& l)
    {
       if(l.size() != graph_matching_constructors.size())
          throw std::runtime_error("number of graph matching labelings and graph mathcing problems disagrees.");

       for(const auto& l_gm : l) {
          auto* c = get_graph_matching_constructor(l_gm.left_graph_no, l_gm.right_graph_no);
          c->read_in_labeling(l_gm.labeling);
       }
    }

    multigraph_matching_input::labeling get_primal() const
    {
       multigraph_matching_input::labeling output;
       output.reserve(graph_matching_constructors.size());

       for(auto& c : graph_matching_constructors)
          output.push_back( {c.first.p, c.first.q, c.second->write_out_labeling()} );

       return output; 
    }

    // TODO: deprecated
    multigraph_matching_input::labeling write_out_labeling() const
    {
        return get_primal();
    }

    multigraph_matching_input export_multigraph_matching_input() const
    {
        multigraph_matching_input mgm;
        for(const auto& gm : graph_matching_constructors) {
            mgm.push_back({gm.first.p, gm.first.q, gm.second->export_graph_matching_input()});
        }
        return mgm;
    }


private:
    std::size_t no_graphs_ = 0;

    std::size_t no_graphs() const
    {
       assert(no_graphs_ == std::max_element(graph_matching_constructors.begin(), graph_matching_constructors.end(), [](const auto& g1, const auto& g2) { return g1.first.q < g2.first.q; })->first.q + 1);
       return no_graphs_;
    }

    multigraph_matching_input export_linear_multigraph_matching_input() const
    {
       multigraph_matching_input mgm;
       for(const auto& gm : graph_matching_constructors)
          mgm.push_back({gm.first.p, gm.first.q, gm.second->export_linear_graph_matching_input()});
       return mgm;
    }

    std::vector<std::size_t> get_triplet_consistency_labels(const std::size_t p, const std::size_t q, const std::size_t p_node) const
    {
       if(p < q) {
          auto* c = get_graph_matching_constructor(p,q);
          assert(p_node < c->left_mrf.get_number_of_variables());
          return c->graph_[p_node];
       } else {
          auto* c = get_graph_matching_constructor(q,p);
          assert(p_node < c->right_mrf.get_number_of_variables());
          return c->inverse_graph_[p_node];
       }
    }

    auto get_matching_factor(const std::size_t p, const std::size_t q, const std::size_t p_node) const
    {
       if(p < q) {
          auto* c = get_graph_matching_constructor(p,q);
          return c->left_mrf.get_unary_factor(p_node);
       } else {
          auto* c = get_graph_matching_constructor(q,p);
          return c->right_mrf.get_unary_factor(p_node); 
       } 
    }

    std::size_t get_matching_index(const std::size_t p, const std::size_t q, const std::size_t p_node, const std::size_t q_node) const
    {
       if(p < q)
          return get_graph_matching_constructor(p,q)->get_left_index(p_node,q_node);
       else
          return get_graph_matching_constructor(q,p)->get_right_index(p_node,q_node);
    }

   std::tuple< PQ_ROW_TRIPLET_CONSISTENCY_MESSAGE*, QR_COLUMN_TRIPLET_CONSISTENCY_MESSAGE*, PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE*, PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE* >
   get_unary_messages(TRIPLET_CONSISTENCY_FACTOR* t) const
   {
      auto pq_row_msgs = t->template get_messages<PQ_ROW_TRIPLET_CONSISTENCY_MESSAGE>();
      auto qr_column_msgs = t->template get_messages<QR_COLUMN_TRIPLET_CONSISTENCY_MESSAGE>();
      auto pr_scalar_msgs = t->template get_messages<PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE>();
      assert(pq_row_msgs.size() == 1 && qr_column_msgs.size() == 1 && pr_scalar_msgs.size() == 2);
      return {pq_row_msgs[0], qr_column_msgs[0], pr_scalar_msgs[0], pr_scalar_msgs[1]};
   }
    void send_messages_to_unaries(TRIPLET_CONSISTENCY_FACTOR* t)
    {
       auto msgs = get_unary_messages(t);
       std::get<0>(msgs)->send_message_to_right(1.0);
       std::get<1>(msgs)->send_message_to_right(1.0);
       std::get<2>(msgs)->send_message_to_right(1.0);
       std::get<3>(msgs)->send_message_to_right(1.0);

       std::get<0>(msgs)->send_message_to_left(0.25);
       std::get<1>(msgs)->send_message_to_left(0.25);
       std::get<2>(msgs)->send_message_to_left(0.25);
       std::get<3>(msgs)->send_message_to_left(1.0);
       std::get<2>(msgs)->send_message_to_left(1.0/3.0);
       std::get<1>(msgs)->send_message_to_left(1.0/3.0);
       std::get<0>(msgs)->send_message_to_left(1.0);
       std::get<1>(msgs)->send_message_to_left(0.5);
       std::get<2>(msgs)->send_message_to_left(1.0);
       std::get<1>(msgs)->send_message_to_left(1.0); 
    }

   std::tuple<PQ_ROW_TRIPLET_CONSISTENCY_ZERO_MESSAGE*, QR_COLUMN_TRIPLET_CONSISTENCY_ZERO_MESSAGE*>
   get_unary_messages(TRIPLET_CONSISTENCY_FACTOR_ZERO* t) const
   {
      auto pq_row_msgs = t->template get_messages<PQ_ROW_TRIPLET_CONSISTENCY_ZERO_MESSAGE>();
      auto qr_column_msgs = t->template get_messages<QR_COLUMN_TRIPLET_CONSISTENCY_ZERO_MESSAGE>();
      assert(pq_row_msgs.size() == 1 && qr_column_msgs.size() == 1);
      return {pq_row_msgs[0], qr_column_msgs[0]}; 
   }

   void send_messages_to_unaries(TRIPLET_CONSISTENCY_FACTOR_ZERO* t)
   {
       auto msgs = get_unary_messages(t);
       std::get<0>(msgs)->send_message_to_right(1.0);
       std::get<1>(msgs)->send_message_to_right(1.0);

       std::get<0>(msgs)->send_message_to_left(0.5);
       std::get<1>(msgs)->send_message_to_left(1.0);
       std::get<0>(msgs)->send_message_to_left(1.0);
   } 

    LP<FMC>* lp_;
    //std::unordered_map<triplet_consistency_factor, TRIPLET_CONSISTENCY_FACTOR*, triplet_consistency_factor_hash> triplet_consistency_factors;
    std::unordered_map<triplet_consistency_factor, ptr_to_triplet_consistency_factor, triplet_consistency_factor_hash> triplet_consistency_factors;

    std::vector<std::pair<graph_matching, std::unique_ptr<GRAPH_MATCHING_CONSTRUCTOR>>> graph_matching_constructors;
    std::unordered_map<graph_matching, GRAPH_MATCHING_CONSTRUCTOR*, graph_matching_hash> graph_matching_constructors_map;

    PositiveIntegerConstraint positiveIntegerConstraint_;
    TCLAP::ValueArg<std::size_t> presolve_iterations_arg_;
    TCLAP::ValueArg<std::string> presolve_reparametrization_arg_;
    TCLAP::ValueArg<std::size_t> primal_checking_triplets_arg_; // no triplets to include during primal checking 
    TCLAP::SwitchArg mcf_reparametrization_arg_; // TODO: this should be part of graph matching constructor
    TCLAP::SwitchArg mcf_primal_rounding_arg_; // TODO: this should be part of graph matching constructor as well
    TCLAP::ValueArg<std::string> graph_matching_construction_arg_; // same as construction_arg_ in the graph matching constructor
    mutable TCLAP::ValueArg<std::string> output_format_arg_; // mutable should not be necessary, but TCLAP's getValue is not const.
    TCLAP::ValueArg<std::string> primal_rounding_algorithms_arg_; // which multicut algorithms to run on

    std::future<multigraph_matching_input::labeling> primal_result_handle_;
}; 

} // namespace LPMP
