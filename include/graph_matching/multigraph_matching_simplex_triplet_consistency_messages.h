#ifndef LPMP_MULTIGRAPH_MATCHING_SIMPLEX_TIRPLET_CONSISTENCY_MESSAGES_H
#define LPMP_MULTIGRAPH_MATCHING_SIMPLEX_TIRPLET_CONSISTENCY_MESSAGES_H

#include "config.hxx"
#include "multigraph_matching_triplet_consistency_factor.h"

namespace LPMP {

class simplex_multigraph_matching_triplet_scalar_consistency_message {
public:
    simplex_multigraph_matching_triplet_scalar_consistency_message(const std::size_t idx) : idx_(idx) {} 

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& l, const double msg, const std::size_t msg_dim) const
   {
       assert(msg_dim == 0);
       assert(idx_ < l.size()-1);
       l[idx_] += msg;
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const double msg, const std::size_t msg_dim) const
   {
       assert(msg_dim == 0);
       r.cost_z += msg;
   }
   
    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const double omega = 1.0)
    {
        const auto cost = l[idx_];
        const auto cost_other = l.min_except(idx_);
        msg[0] -= omega*(cost - cost_other);
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const double omega = 1.0)
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
             const bool changed = (r.z != multigraph_matching_primal_inactive);
             r.z = multigraph_matching_primal_inactive;
             return changed; 
          }
       }
       return false;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal_deactivated(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       if(r.z != std::numeric_limits<std::size_t>::max()) {
          if(r.z == multigraph_matching_primal_inactive) {
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
       assert(idx_ < l.size()-1);
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

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& l, const double msg, const std::size_t dim) const
   {
       assert(dim < l.size()-1);
       l[dim] += msg;
   } 

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg) const
   {
       assert(msg.size() == l.size()-1);
       for(std::size_t i=0; i<l.size()-1; ++i) {
           l[i] += msg[i];
       }
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const double msg, const std::size_t dim) const
   {
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
       if(C == Chirality::left) {
           assert(msg.size() == r.cost_x.size());
           for(std::size_t i=0; i<r.cost_x.size(); ++i)
               r.cost_x[i] += msg[i]; 
       } else {
           assert(C == Chirality::right);
           assert(msg.size() == r.cost_y.size());
           for(std::size_t j=0; j<r.cost_y.size(); ++j)
               r.cost_y[j] += msg[j]; 
       }
   }
   
    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const double omega = 1.0)
    {
        for(std::size_t i=0; i<l.size()-1; ++i) {
            msg[i] -= omega*(l[i] - l.back());
        }
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const double omega = 1.0)
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

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR> 
    void receive_restricted_message_from_right(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       if(C == Chirality::left) {
          if(r.x != multigraph_matching_primal_inactive) {
             for(std::size_t i=0; i<l.size()-1; ++i) {
                l[i] = std::numeric_limits<double>::infinity();
             } 
          } else {
             for(std::size_t i=0; i<l.size(); ++i) {
                if(r.x != multigraph_matching_primal_inactive)
                   l[i] = std::numeric_limits<double>::infinity();
             } 
          }
       } else {
          assert(false); 
       }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR> 
    bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
    {
       const std::size_t right_label = l.primal() == l.size()-1 ? multigraph_matching_primal_inactive : l.primal();
       std::size_t& right_primal = C == Chirality::left ? r.x : r.y;
       const bool changed = right_primal != right_label;
       right_primal = right_label;
       return changed;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal_deactivated(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       const std::size_t right_label = C == Chirality::left ? r.x : r.y;
       const std::size_t left_label = right_label == multigraph_matching_primal_inactive ? l.size()-1 : right_label;

       const bool changed = l.primal() != left_label;
       l.primal() = left_label;
       return changed;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       const std::size_t right_label = C == Chirality::left ? r.x : r.y;
       const std::size_t left_label = right_label == multigraph_matching_primal_inactive ? l.size()-1 : right_label;
       assert(l.size() - 1 == (C == Chirality::left ? r.cost_x.size() : r.cost_y.size()));

       return l.primal() == left_label;
    }
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_SIMPLEX_TIRPLET_CONSISTENCY_MESSAGES_H
