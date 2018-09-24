#ifndef LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX
#define LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX

#include "LP.h"
#include "vector"
#include <array>
#include <limits>

namespace LPMP {

// Z_ij = \sum_p X_ip Y_pj
class multigraph_matching_triplet_consistency_factor {
public:
   constexpr static std::size_t primal_not_set = std::numeric_limits<std::size_t>::max(); // primal is not set to any value
   constexpr static std::size_t primal_inactive = std::numeric_limits<std::size_t>::max()-1;

    // cost vectors
    REAL cost_z;
    vector<REAL> cost_x;
    vector<REAL> cost_y;

    // primal
    std::size_t x, y, z;

    multigraph_matching_triplet_consistency_factor(const std::size_t dim)
    : cost_z(0.0),
    cost_x(dim, 0.0),
    cost_y(dim, 0.0)
    {
        assert(dim >= 1);
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

    bool primal_feasible(const std::size_t _x, const std::size_t _y, const std::size_t _z) const
    {
       if(_x == primal_not_set || _y == primal_not_set || _z == primal_not_set) return false;
       const std::size_t no_inactive = std::size_t(_x == primal_inactive) + std::size_t(_y == primal_inactive) + std::size_t(_z == primal_inactive);
       if(no_inactive >= 2) return true;
       if(_y != _x && _z == primal_inactive) return true;
       if(_y == _x && _z == 1) return true;
       if(debug())
          std::cout << "triplet consistency constraint infeasible\n";
       return false;
    }

    bool primal_feasible() const { return primal_feasible(x,y,z); }

    REAL EvaluatePrimal() const
    {
       if(primal_feasible()) {
          return evaluate(x,y,z);
       } else {
          return std::numeric_limits<REAL>::infinity();
       }
    }

    REAL evaluate(const std::size_t _x, const std::size_t _y, const std::size_t _z) const
    {
       double cost = 0.0;
       if(_x < cost_x.size()) cost += cost_x[_x];
       if(_y < cost_y.size()) cost += cost_y[_y];
       if(_z == 1) cost += cost_z;
       return cost;
    } 

    template<typename FUNC>
    void for_each_labeling(FUNC&& f) const
    {
       // z == primal_inactive
       f(primal_inactive,primal_inactive,primal_inactive);

       for(std::size_t i=0; i<cost_x.size(); ++i) {
          f(i,primal_inactive,primal_inactive);
       }

       for(std::size_t i=0; i<cost_x.size(); ++i) {
          f(primal_inactive,i,primal_inactive);
       }

       for(std::size_t i=0; i<cost_x.size(); ++i) {
          for(std::size_t j=0; j<cost_y.size(); ++j) {
             if(i != j) {
                f(i,j,primal_inactive);
             } 
          }
       }

       // z == 1
       f(primal_inactive, primal_inactive, 1);

       for(std::size_t i=0; i<cost_x.size(); ++i) {
          f(i,i,1);
       }
    }

    void MaximizePotentialAndComputePrimal()
    {
       double best_cost = std::numeric_limits<double>::infinity();
       std::size_t best_x = primal_not_set;
       std::size_t best_y = primal_not_set;
       std::size_t best_z = primal_not_set;

       auto update_primal = [&](const std::size_t _x, const std::size_t _y, const std::size_t _z) {
          assert(primal_feasible(_x,_y,_z));
          if((x == primal_not_set || x == _x) && (y == primal_not_set || y == _y) && (z == primal_not_set || z == _z)) {
             const double cost = evaluate(_x,_y,_z);
             if(cost < best_cost) {
                best_cost = cost;
                best_x = _x;
                best_y = _y;
                best_z = _z;
             }
          }
       };

       for_each_labeling(update_primal);

       assert(primal_feasible(best_x,best_y,best_z)); 

       x = best_x;
       y = best_y;
       z = best_z; 
    }

    std::array<REAL,2> z_marginals() const
    {
       // z = multigraph_matching_triplet_consistency_factor::primal_inactive
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
        // z = multigraph_matching_triplet_consistency_factor::primal_inactive
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

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(x,y); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cost_z, cost_x, cost_y ); }

   auto export_variables() { return std::tie(cost_z, cost_y, cost_y); }

   void init_primal() 
   { 
      x = primal_not_set;
      y = primal_not_set;
      z = primal_not_set;
   }

   // if two primal variables are set, we must set the third one accordingly
   void PropagatePrimal()
   {
      const std::size_t no_variables_set = std::size_t(x != primal_not_set) + std::size_t(x != primal_not_set) + std::size_t(z != primal_not_set);
      if(no_variables_set == 2) {
         if(x == std::numeric_limits<std::size_t>::max()) {
            if(z == 1)
               x = y;
         } else if(y == std::numeric_limits<std::size_t>::max()) {
            if(z == 1)
               y = x; 
         } else {
            assert(z == std::numeric_limits<std::size_t>::max());
            if(x == y)
               z = 1;
            else
               z = multigraph_matching_triplet_consistency_factor::primal_inactive;
         }
      }
   } 

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
   void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX msg_dim) const
   {
       assert(msg_dim == 0);
       l[idx_] += msg;
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX msg_dim) const
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

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR> 
    bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
    {
       if(l.primal() < l.size()) {
          if(l.primal() == idx_) {
             const bool changed = (r.z != 1);
             r.z = 1;
             return changed;
          } else {
             const bool changed = (r.z != multigraph_matching_triplet_consistency_factor::primal_inactive);
             r.z = multigraph_matching_triplet_consistency_factor::primal_inactive;
             return changed; 
          }
       }
       return false;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       if(r.z != std::numeric_limits<std::size_t>::max()) {
          if(r.z == multigraph_matching_triplet_consistency_factor::primal_inactive) {
             const bool changed = (l.primal() == idx_);
             if(changed)
                l.init_primal(); 
             return changed;
          } else {
             assert(r.z == 1);
             const bool changed = (l.primal() != idx_);
             l.primal() = idx_;
             return changed;
          }
       }
       return false;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       const bool left_active = l.primal() == idx_;
       const bool right_active = (r.z == 1);
       if(left_active != right_active) {
             std::cout << "inconsistency in scalar consistency message\n";
       }
       return left_active == right_active;
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
   void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX dim) const
   {
       assert(dim < idx_.size());
       l[idx_[dim]] += msg;
   } 

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg) const
   {
       assert(msg.size() == idx_.size());
       for(std::size_t i=0; i<idx_.size(); ++i) {
           l[idx_[i]] += msg[i];
       }
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX dim) const
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
   void RepamRight(RIGHT_FACTOR& r, const MSG& msg) const
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

    // given simplex labeling, return corresponding label in triplet consistency factor.
    std::size_t get_triplet_consistency_factor_index(const std::size_t l) const
    {
       for(std::size_t i=0; i<idx_.size(); ++i)
          if(l == idx_[i]) 
             return i;
       return multigraph_matching_triplet_consistency_factor::primal_inactive; 
    } 

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR> 
    void receive_restricted_message_from_right(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       if(C == Chirality::left) {
          if(r.x != multigraph_matching_triplet_consistency_factor::primal_inactive) {
             for(std::size_t i=0; i<idx_.size(); ++i) {
                l[idx_[i]] = std::numeric_limits<REAL>::infinity();
             } 
          } else {
             for(std::size_t i=0; i<l.size(); ++i) {
                if(r.x != idx_[i])
                l[i] = std::numeric_limits<REAL>::infinity();
             } 
          }
       } else {

       }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR> 
    bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
    {
       if(l.primal() < l.size()) {
          const std::size_t i = get_triplet_consistency_factor_index(l.primal());

          if(C == Chirality::left) {
             const bool changed = (r.x != i);
             r.x = i;
             return changed;
          } else {
             const bool changed = (r.y != i);
             r.y = i;
             return changed;
          }
       }
       return false;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       const std::size_t label = C == Chirality::left ? r.x : r.y;

       if(label != multigraph_matching_triplet_consistency_factor::primal_inactive && label != multigraph_matching_triplet_consistency_factor::primal_not_set) {
          const bool changed = (l.primal() != idx_[label]);
          l.primal() = idx_[label];
          return changed;
       } else {
          return false;
       }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       for(std::size_t i=0; i<idx_.size(); ++i) {
          const std::size_t label = idx_[i];
          const bool left_active = (l.primal() == label);
          const bool right_active = [&]() {
             if(C == Chirality::left)
                return r.x == label;
             else 
                return r.y == label;
          }();
          if(left_active != right_active) {
             std::cout << "inconsistency in vector consistency message\n";
             return false;
          }
       }
       return true;
    }

private:
    vector<std::size_t> idx_;
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX
