#ifndef LPMP_GRAPH_MATCHING_INTER_QUADRATIC_MESSAGE_HXX
#define LPMP_GRAPH_MATCHING_INTER_QUADRATIC_MESSAGE_HXX

#include <array>

namespace LPMP {

class graph_matching_inter_quadratic_message {
public:
   template<typename INDICES_ITERATOR>
   graph_matching_inter_quadratic_message(INDICES_ITERATOR left_indices_begin, INDICES_ITERATOR left_indices_end, INDICES_ITERATOR right_indices_begin, INDICES_ITERATOR right_indices_end)
   {
      assert(std::distance(left_indices_begin, left_indices_end) == std::distance(right_indices_begin, right_indices_end));
      assert(std::distance(left_indices_begin, left_indices_end) >= 1);
      assert(std::distance(left_indices_begin, left_indices_end) <= 2);

      for(std::size_t i=0; i<std::distance(left_indices_begin, left_indices_end); ++i)
         left_indices_[i] = left_indices_begin[i];
      for(std::size_t i=0; i<std::distance(right_indices_begin, right_indices_end); ++i)
         right_indices_[i] = right_indices_begin[i];

      assert(left_indices_[0] < left_indices_[1]);
      assert(right_indices_[0] < right_indices_[1]);
   }


   template<typename FACTOR, typename MSG>
   void send_message_to_left(const FACTOR& rightPot, MSG& msg, const double omega)
   {
      send_message(rightPot, msg, right_indices_, omega);
   }
   template<typename FACTOR, typename MSG>
   void send_message_to_right(const FACTOR& leftPot, MSG& msg, const double omega)
   {
      send_message(leftPot, msg, left_indices_, omega);
   }

   template<typename FACTOR, typename MSG, typename INDICES_ARRAY>
   void send_message(const FACTOR& pot, MSG& msg, const INDICES_ARRAY& indices, const double omega = 1.0)
   {
      std::array<double,2> msg_vars;
      msg_vars[0] = pot[indices[0]];
      std::size_t no_idxs = 1;
      if(indices[1] < std::numeric_limits<std::size_t>::max()) {
         msg_vars[1] = pot[indices[1]];
         ++no_idxs;
      }

      const double min_not_taken = [&]() {
         if(no_idxs == 1) { 
            return pot.lower_bound_except( std::array<std::size_t,1>({indices[0]}) );
         } else {
            assert(no_idxs == 2);
            return pot.lower_bound_except(indices);
         }
      }();

      msg_vars[0] -= min_not_taken;
      msg_vars[1] -= min_not_taken;

      msg[0] -= omega*msg_vars[0];
      if(indices[1] < std::numeric_limits<std::size_t>::max()) {
         msg[1] -= omega*msg_vars[1];
      }
   }

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& l, const double msg_val, const std::size_t msg_idx) const
   {
      assert(left_indices_[msg_idx] < std::numeric_limits<std::size_t>::max());
      l.cost(left_indices_[msg_idx]) += msg_val;
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const double msg_val, const std::size_t msg_idx) const
   {
      assert(right_indices_[msg_idx] < std::numeric_limits<std::size_t>::max());
      r.cost(right_indices_[msg_idx]) += msg_val;
   }


private:
   std::array<std::size_t,2> left_indices_ = {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()};
   std::array<std::size_t,2> right_indices_ = {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()};
};

} // namespace LPMP

#endif // LPMP_GRAPH_MATCHING_INTER_QUADRATIC_MESSAGE_HXX

