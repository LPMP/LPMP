#ifndef LPMP_MULTIGRAPH_MATCHING_CONSTRUCTOR_HXX
#define LPMP_MULTIGRAPH_MATCHING_CONSTRUCTOR_HXX

#include "graph_matching_constructor.hxx"
#include "multigraph_matching_input.h"
#include <unordered_map>
#include <memory>
#include "LP.h"

namespace LPMP {

// works for three matching problems forming a triangle
template<typename GRAPH_MATCHING_CONSTRUCTOR, typename TRIPLET_CONSISTENCY_FACTOR, typename PQ_ROW_TRIPLET_CONSISTENCY_MESSAGE, typename QR_COLUMN_TRIPLET_CONSISTENCY_MESSAGE, typename PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE>
class multigraph_matching_constructor {

    struct triplet_consistency_factor {
        std::size_t p,q,r; // graphs
        std::size_t p_node, r_node;
        bool operator==(const triplet_consistency_factor& o) const 
        {
            return p == o.p && q == o.q && r == o.r && p_node == o.p_node && r_node == o.r_node;
        }
    }; 
    struct triplet_consistency_factor_hash  {
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

    struct graph_matching{
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
        : lp_(&s.GetLP())
        {}

    void begin()
    {
       std::vector<graph_matching> graph_matching_indices;

       assert(graph_matching_constructors.size() > 0);
       for(auto& c : graph_matching_constructors) {
          c.second->begin();
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

    template<typename SOLVER>
    GRAPH_MATCHING_CONSTRUCTOR* add_graph_matching_problem(const std::size_t p, const std::size_t q, SOLVER& s)
    {
        assert(!has_graph_matching_problem(p,q));
        auto ptr = std::make_unique<GRAPH_MATCHING_CONSTRUCTOR>(s);
        auto* c = ptr.get();
        graph_matching_constructors.insert(std::make_pair(graph_matching({p,q}), std::move(ptr))); 
        return c;
    }

    bool has_graph_matching_problem(const std::size_t p, const std::size_t q) const
    {
        assert(p<q);
        return graph_matching_constructors.count({p,q}) > 0;
    }

    GRAPH_MATCHING_CONSTRUCTOR* get_graph_matching_constructor(const std::size_t p, const std::size_t q) const
    {
        assert(has_graph_matching_problem(p,q));
        return graph_matching_constructors.find(graph_matching({p,q}))->second.get();
    }

    bool has_triplet_consistency_factor(const triplet_consistency_factor& t) const
    {
        assert(t.p < t.r);
        return triplet_consistency_factors.count(t) > 0; 
    }

    // indices for p->q-unary, q->r-unary, p->r unary, r->q unary partaking in triplet consistency factor.
    // associated types are std::vector<std::size_t>, std::vector<std::size_t>, std::size_t, std::size_t
    // this is followed by the unaries that are used to that end: matching factor p->q-unary, q->r-unary, p->r unary, r->q unary
    auto
    get_triplet_consistency_factor_data(const triplet_consistency_factor& t) const
    {
        const auto [p,q,r,p_node,r_node] = std::make_tuple(t.p, t.q, t.r, t.p_node, t.r_node);
        assert(p<r);
        triplet_consistency_factor triplet({p,q,r, p_node, r_node});

        // compute p->q and p->r matching data for triplet consistency factor
        auto get_matching_data = [&](const std::size_t p, const std::size_t q, const std::size_t p_node) {
            if(p<q) {
                auto* pq_constructor = get_graph_matching_constructor(p,q);
                auto* f_pq = pq_constructor->left_mrf.get_unary_factor(p_node);
                const auto& pq_labels = pq_constructor->graph_[p_node]; 
                return std::make_tuple(pq_constructor, f_pq, pq_labels);
            } else {
                auto* pq_constructor = get_graph_matching_constructor(q,p);
                auto* f_pq = pq_constructor->right_mrf.get_unary_factor(p_node);
                const auto& pq_labels = pq_constructor->inverse_graph_[p_node]; 
                return std::make_tuple(pq_constructor, f_pq, pq_labels);
            }
        };

        auto [pq_constructor, f_pq, pq_labels] = get_matching_data(p,q, p_node);
        auto [qr_constructor, f_qr, qr_labels] = get_matching_data(r,q, r_node);

        // record which entries in f_pq and f_qr point towards the same nodes in q
        std::vector<std::size_t> common_nodes;
        std::set_intersection(pq_labels.begin(), pq_labels.end(), qr_labels.begin(), qr_labels.end(), std::back_inserter(common_nodes));
        assert(common_nodes.size() > 0);

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

        // find the two unaries in p->r graph matching
        auto pr_constructor = get_graph_matching_constructor(p,r);

        auto* f_pr_left = pr_constructor->left_mrf.get_unary_factor(p_node);
        const auto& p_matching_labels = pr_constructor->graph_[p_node];
        const auto p_index = std::find(p_matching_labels.begin(), p_matching_labels.end(), r_node) - p_matching_labels.begin();

        auto* f_pr_right = pr_constructor->right_mrf.get_unary_factor(r_node);
        const auto& r_matching_labels = pr_constructor->inverse_graph_[r_node];
        const auto r_index = std::find(r_matching_labels.begin(), r_matching_labels.end(), p_node) - r_matching_labels.begin();

        return std::make_tuple(pq_factor_indices, qr_factor_indices, p_index, r_index, f_pq, f_qr, f_pr_left, f_pr_right);
    }

    TRIPLET_CONSISTENCY_FACTOR* add_triplet_consistency_factor(const triplet_consistency_factor& t)
    {
        const auto [p,q,r,p_node,r_node] = std::make_tuple(t.p, t.q, t.r, t.p_node, t.r_node);
        assert(p<r);
        triplet_consistency_factor triplet({p,q,r, p_node, r_node});
        assert(triplet_consistency_factors.count(triplet) == 0);

        const auto [pq_indices, qr_indices, p_index, r_index, f_pq, f_qr, f_pr_left, f_pr_right] = get_triplet_consistency_factor_data(t);
        assert(pq_indices.size() == qr_indices.size());

        auto* f_triplet = lp_->template add_factor<TRIPLET_CONSISTENCY_FACTOR>(pq_indices.size());
        triplet_consistency_factors.insert(std::make_pair(t, f_triplet));

        lp_->template add_message<PQ_ROW_TRIPLET_CONSISTENCY_MESSAGE>(f_pq, f_triplet, pq_indices.begin(), pq_indices.end());
        lp_->template add_message<QR_COLUMN_TRIPLET_CONSISTENCY_MESSAGE>(f_qr, f_triplet, qr_indices.begin(), qr_indices.end());

        lp_->template add_message<PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE>(f_pr_left, f_triplet, p_index);
        lp_->template add_message<PR_SCALAR_TRIPLET_CONSISTENCY_MESSAGE>(f_pr_right, f_triplet, r_index); 

        triplet_consistency_factors.insert(std::make_pair(t,f_triplet));

        return f_triplet;
    }

    REAL triplet_consistency_dual_increase(triplet_consistency_factor t,
            GRAPH_MATCHING_CONSTRUCTOR* pq_constructor, GRAPH_MATCHING_CONSTRUCTOR* qr_constructor, GRAPH_MATCHING_CONSTRUCTOR* pr_constructor) const
    {
        auto [pq_indices, qr_indices, p_index, r_index, f_pq, f_qr, f_pr_left, f_pr_right] = get_triplet_consistency_factor_data(t);
        assert(pq_indices.size() == qr_indices.size());

        // compute lower bound before reparametrizing
        const REAL prev_lb = f_pq->LowerBound() + f_qr->LowerBound() + f_pr_left->LowerBound() + f_pr_right->LowerBound();

        // send messages into triplet consistency factor
        simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left> m_pq(pq_indices.begin(), pq_indices.end());
        simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right> m_qr(qr_indices.begin(), qr_indices.end());
        simplex_multigraph_matching_triplet_scalar_consistency_message m_pr_left(p_index);
        simplex_multigraph_matching_triplet_scalar_consistency_message m_pr_right(r_index);
        multigraph_matching_triplet_consistency_factor f(pq_indices.size());

        vector<REAL> msg_val_pq(pq_indices.size(), 0.0);
        m_pq.send_message_to_right(*f_pq->get_factor(), msg_val_pq, 1.0);
        m_pq.RepamRight(f, -1.0*msg_val_pq);
        m_pq.RepamLeft(*f_pq->get_factor(), msg_val_pq);

        vector<REAL> msg_val_qr(qr_indices.size(), 0.0);
        m_qr.send_message_to_right(*f_qr->get_factor(), msg_val_qr, 1.0);
        m_qr.RepamRight(f, -1.0*msg_val_qr);
        m_qr.RepamLeft(*f_qr->get_factor(), msg_val_qr);

        array<REAL,1> msg_val_pr_left = {0.0};
        m_pr_left.send_message_to_right(*f_pr_left->get_factor(), msg_val_pr_left, 1.0);
        m_pr_left.RepamRight(f, -msg_val_pr_left[0], 0);
        m_pr_left.RepamLeft(*f_pr_left->get_factor(), msg_val_pr_left[0], 0);

        array<REAL,1> msg_val_pr_right = {0.0};
        m_pr_right.send_message_to_right(*f_pr_right->get_factor(), msg_val_pr_right, 1.0);
        m_pr_right.RepamRight(f, -msg_val_pr_right[0], 0);
        m_pr_right.RepamLeft(*f_pr_right->get_factor(), msg_val_pr_right[0], 0);

        const REAL after_lb = f.LowerBound() + f_pq->LowerBound() + f_qr->LowerBound() + f_pr_left->LowerBound() + f_pr_right->LowerBound();

        // revert changes
        m_pq.RepamLeft(*f_pq->get_factor(), -1.0*msg_val_pq);
        m_qr.RepamLeft(*f_qr->get_factor(), -1.0*msg_val_qr);
        m_pr_left.RepamLeft(*f_pr_left->get_factor(), -msg_val_pr_left[0], 0);
        m_pr_right.RepamLeft(*f_pr_right->get_factor(), -msg_val_pr_right[0], 0);

        assert(prev_lb == f_pq->LowerBound() + f_qr->LowerBound() + f_pr_left->LowerBound() + f_pr_right->LowerBound());
        assert(after_lb >= prev_lb - eps);

        return after_lb - prev_lb; 
    }

    INDEX Tighten(const INDEX no_constraints_to_add)
    {
        // iterate over all triplets of graphs and enumerate all possible triplet consistency factors that can be added. 
        // Record guaranteed dual increase of adding the triplet consistency factor
        std::vector<std::pair<triplet_consistency_factor, REAL>> triplet_consistency_candidates;

        // enumerate all graph matching triplets
        for(auto& g : graph_matching_constructors) {
            std::cout << g.first.p << "," << g.first.q << "\n";
        }
        const auto no_graphs = std::max_element(graph_matching_constructors.begin(), graph_matching_constructors.end(), [](const auto& g1, const auto& g2) { return g1.first.q < g2.first.q; })->first.q + 1;
        for(std::size_t r=0; r<no_graphs; ++r) {
            for(std::size_t q=0; q<r; ++q) {
                for(std::size_t p=0; p<q; ++p) {
                    auto* pq_constructor = get_graph_matching_constructor(p,q);
                    auto* qr_constructor = get_graph_matching_constructor(q,r);
                    auto* pr_constructor = get_graph_matching_constructor(p,r);
                    const auto p_no_nodes = pq_constructor->left_mrf.get_number_of_variables();
                    const auto q_no_nodes = pq_constructor->right_mrf.get_number_of_variables();
                    const auto r_no_nodes = pr_constructor->right_mrf.get_number_of_variables();

                    // enumerate all triplet_consistency_factor
                    triplet_consistency_factor t;

                    t.p = p; t.q = q; t.r = r;
                    for(std::size_t p_node=0; p_node<p_no_nodes; ++p_node) {
                        for(std::size_t r_node=0; r_node<r_no_nodes; ++r_node) {
                            t.p_node = p_node; t.r_node = r_node;
                            const auto guaranteed_dual_increase = triplet_consistency_dual_increase(t, pq_constructor, qr_constructor, pr_constructor); 
                            if(guaranteed_dual_increase >= eps) {
                                triplet_consistency_candidates.push_back( std::make_pair(t, guaranteed_dual_increase) );
                            }
                        }
                    }

                    t.p = p; t.q = r; t.r = q;
                    for(std::size_t p_node=0; p_node<p_no_nodes; ++p_node) {
                        for(std::size_t q_node=0; q_node<q_no_nodes; ++q_node) {
                            t.p_node = p_node; t.r_node = q_node;
                            const auto guaranteed_dual_increase = triplet_consistency_dual_increase(t, pq_constructor, qr_constructor, pr_constructor); 
                            if(guaranteed_dual_increase >= eps) {
                                triplet_consistency_candidates.push_back( std::make_pair(t, guaranteed_dual_increase) );
                            } 
                        }
                    }

                    t.p = q; t.q = p; t.r = r;
                    for(std::size_t q_node=0; q_node<q_no_nodes; ++q_node) {
                        for(std::size_t r_node=0; r_node<r_no_nodes; ++r_node) {
                            t.p_node = q_node; t.r_node = r_node;
                            const auto guaranteed_dual_increase = triplet_consistency_dual_increase(t, pq_constructor, qr_constructor, pr_constructor); 
                            if(guaranteed_dual_increase >= eps) {
                                triplet_consistency_candidates.push_back( std::make_pair(t, guaranteed_dual_increase) );
                            }

                        }
                    }

                }
            }
        }

        std::sort(triplet_consistency_candidates.begin(), triplet_consistency_candidates.end(), [](const auto& t1, const auto& t2) { return t1.second > t2.second; });

        std::size_t no_constraints_added = 0;
        for(auto [t,cost] : triplet_consistency_candidates) {
            if(!has_triplet_consistency_factor(t)) {
                add_triplet_consistency_factor(t);
                no_constraints_added++;
                if(no_constraints_added >= no_constraints_to_add) {
                    break;
                }
            } 
        }

        if(debug()) {
            std::cout << "Added " << no_constraints_added << " triplet consistency factor for multigraph matching\n";
        }

        return no_constraints_added;
    }

    void construct(const mgm_input& input)
    {
        for(const auto& gm : input) {
            auto* gm_constructor = add_graph_matching_problem(gm.left_graph_no, gm.right_graph_no, lp_);
            gm_constructor->read_input(gm.gm_input);
            gm_constructor->construct();
        }
    }

private:
    LP<FMC>* lp_;
    std::unordered_map<triplet_consistency_factor, TRIPLET_CONSISTENCY_FACTOR*, triplet_consistency_factor_hash> triplet_consistency_factors;
    std::unordered_map<graph_matching, std::unique_ptr<GRAPH_MATCHING_CONSTRUCTOR>, graph_matching_hash> graph_matching_constructors; 
};


} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_CONSTRUCTOR_HXX
