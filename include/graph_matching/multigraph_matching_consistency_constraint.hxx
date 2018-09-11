#ifndef LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX
#define LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX

#include "LP.h"
#include "vector"

namespace LPMP {

// Z_ij = \sum_p X_ip Y_pj
class multigraph_matching_triplet_consistency_factor {
public:
    // cost vectors
    REAL cost_z;
    vector<REAL> cost_x;
    vector<REAL> cost_y;

    // primal
    std::size_t x, y;

    multigraph_matching_triplet_consistency_factor(const std::size_t dim)
    : cost_z(0.0),
    cost_x(dim),
    cost_y(dim)
    {
        std::fill(cost_x.begin(), cost_x.end(), 0.0);
        std::fill(cost_y.begin(), cost_y.end(), 0.0);
        assert(dim > 1);
    }

    // put in general place somewhere
    template<typename ITERATOR>
    std::array<std::size_t, 2> min_indices(ITERATOR begin, ITERATOR end) const
    {
        REAL smallest_taken = std::numeric_limits<REAL>::infinity();
        REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
        std::size_t min_idx, second_min_idx;
        for(auto it=begin; it!=end; ++it) {
            assert(false); 
        } 
        return {min_idx, second_min_idx};
    }

    REAL EvaluatePrimal() const
    {
        if(x > cost_x.size() || y > cost_y.size()) { return std::numeric_limits<REAL>::infinity(); }
        else if(x == y) { return cost_z + cost_x[x] + cost_y[y]; }
        else { return cost_x[x] + cost_y[y]; }
    }

    std::array<REAL,2> z_marginals() const
    {
        // z = 0
        REAL min_different_component = std::numeric_limits<REAL>::infinity();
        REAL min_x = cost_x[0];
        REAL min_y = cost_y[0];

        // z = 1
        REAL min_same_component = cost_x[0] + cost_y[0];

        for(std::size_t i=1; i<cost_x.size(); ++i) {
            min_different_component = std::min({min_different_component, cost_x[i] + min_y, min_x + cost_y[i]}); 
            min_x = std::min(min_x, cost_x[i]);
            min_y = std::min(min_y, cost_y[i]);

            min_same_component = std::min(min_same_component, cost_x[i] + cost_y[i]);
        }

        const REAL cost_0 = std::min({min_different_component, min_x, min_y, 0.0});
        const REAL cost_1 = std::min(min_same_component + cost_z, cost_z);

        return {cost_0, cost_1};
    }

    vector<REAL> marginals(const vector<REAL>& v1, const vector<REAL>& v2) const
    {
        assert(v1.size() == v2.size());
        vector<REAL> marginals(v1.size());
        for(std::size_t i=0; i<marginals.size(); ++i) { marginals[i] = v1[i]; }

        for(std::size_t i=0; i<v1.size(); ++i) {
            for(std::size_t j=0; j<v2.size(); ++j) {
                if(i == j) {
                    marginals[i] = std::min(marginals[i], v1[i] + v2[i] + cost_z);
                } else {
                    marginals[i] = std::min(marginals[i], v1[i] + v2[j]); 
                }
            }
        }

        const REAL cost_not_taken = std::min({0.0, cost_z, v2.min()});
        for(std::size_t i=0; i<marginals.size(); ++i) { marginals[i] -= cost_not_taken; }

        return marginals;
    }

    vector<REAL> x_marginals() const { return marginals(cost_x, cost_y); }
    vector<REAL> y_marginals() const { return marginals(cost_y, cost_x); }

    REAL LowerBound() const
    {
        const auto z_marg = z_marginals();
        return std::min(z_marg[0], z_marg[1]);
        // z = 0
        const auto x_min = cost_x.min();
        const auto y_min = cost_y.min();
        REAL cost_0 = std::min({0.0, cost_x.min(), cost_y.min()});
        for(std::size_t i=0; i<cost_x.size(); ++i) {
            for(std::size_t j=0; j<cost_y.size(); ++j) {
                if(i != j) {
                    cost_0 = std::min(cost_0, cost_x[i] + cost_y[j]);
                } 
            }
        }

        // z = 1
        REAL cost_1 = std::numeric_limits<REAL>::infinity();
        for(std::size_t i=0; i<cost_x.size(); ++i) {
            cost_1 = std::min(cost_1, cost_x[i] + cost_y[i]); 
        }
        cost_1 += cost_z;

        return std::min(cost_0, cost_1);
    }

    void MaximizePotentialAndComputePrimal() 
    {
    }
    

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(x,y); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cost_z, cost_x, cost_y ); }

   auto export_variables() { return std::tie(cost_z, cost_y, cost_y); }

   void init_primal() { x = std::numeric_limits<std::size_t>::max(); y = std::numeric_limits<std::size_t>::max(); }

   template<typename EXTERNAL_SOLVER>
   void construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::scalar z, typename EXTERNAL_SOLVER::vector x, typename EXTERNAL_SOLVER::vector y) const
   {
   }
   template<typename EXTERNAL_SOLVER>
   void convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::scalar z, typename EXTERNAL_SOLVER::vector x, typename EXTERNAL_SOLVER::vector y) const
   {
   }
};

class simplex_multigraph_matching_triplet_scalar_consistency_message {
public:
    simplex_multigraph_matching_triplet_scalar_consistency_message(const std::size_t idx) : idx_(idx) {} 

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX msg_dim)
   {
       assert(msg_dim == 0);
       l[idx_] += msg;
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX msg_dim)
   {
       assert(msg_dim == 0);
       r.cost_z += msg;
   }
   
    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
    {
        const auto cost = l[idx_];
        const auto cost_other = l.min_except(idx_);
        msg[0] -= omega*(cost - cost_other);
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
    {
        const auto z_marginals = r.z_marginals();
        msg[0] -= omega*(z_marginals[1] - z_marginals[0]);
    }

private:
    const std::size_t idx_;
};

// left factor: simplex factor
// right factor multigraph_matching_triplet_consistency_factor
template<Chirality C> // does the message connect to x (left) or y (right) variable
class simplex_multigraph_matching_triplet_vector_consistency_message {
public:
    template<typename ITERATOR>
    simplex_multigraph_matching_triplet_vector_consistency_message(ITERATOR idx_begin, ITERATOR idx_end)
    : idx_(idx_begin, idx_end)
    {}

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX dim)
   {
       assert(dim < idx_.size());
       l[idx_[dim]] += msg;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
   {
       assert(msg.size() == idx_.size());
       for(std::size_t i=0; i<idx_.size(); ++i) {
           l[idx_[i]] += msg[i];
       }
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX dim)
   {
       assert(dim < idx_.size());
       assert(idx_.size() == r.cost_x.size() && idx_.size() == r.cost_y.size());

       if(C == Chirality::left) {
           r.cost_x[dim] += msg;
       } else {
           assert(C == Chirality::right);
           r.cost_y[dim] += msg;
       }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
   {
       assert(msg.size() == idx_.size());
       if(C == Chirality::left) {
           assert(msg.size() == r.cost_x.size());
           for(std::size_t i=0; i<r.cost_x.size(); ++i) {
               r.cost_x[i] += msg[i]; 
           }
       } else {
           assert(C == Chirality::right);
           assert(msg.size() == r.cost_y.size());
           for(std::size_t i=0; i<r.cost_y.size(); ++i) {
               r.cost_y[i] += msg[i]; 
           } 
       }
   }
   
    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
    {
        std::vector<bool> edge_taken(l.size(), false);
        for(auto i : idx_) edge_taken[i] = true;

        REAL smallest_taken = std::numeric_limits<REAL>::infinity();
        REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
        REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
        for(std::size_t i=0; i<l.size(); ++i) {
            const auto val = l[i];
            if(edge_taken[i]) {
                const REAL min = std::min(smallest_taken, val);
                const REAL max = std::max(smallest_taken, val);
                smallest_taken = min;
                second_smallest_taken = std::min(max, second_smallest_taken);
            } else {
                smallest_not_taken = std::min(val, smallest_not_taken);
            }
        }

        const auto set_to_cost = std::min(second_smallest_taken, smallest_not_taken);

        for(std::size_t i=0; i<idx_.size(); ++i) {
            msg[i] -= omega*(l[idx_[i]] - set_to_cost);
        }
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
    {
       if(C == Chirality::left) {
           auto msg_val =  r.x_marginals();
           msg -= omega*msg_val; 
       } else {
           assert(C == Chirality::right);
           auto msg_val =  r.y_marginals();
           msg -= omega*msg_val; 
       }
    }

private:
    vector<std::size_t> idx_;
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX
