#ifndef LPMP_MESSAGE_PASSING_SCHEDULE_HXX
#define LPMP_MESSAGE_PASSING_SCHEDULE_HXX

#include <cstddef>

namespace LPMP {

namespace message_passing_schedule {

   using type = std::size_t;

   constexpr static std::size_t send_message_to_left = 1<<0;
   constexpr static std::size_t send_message_to_right = 1<<1;
   constexpr static std::size_t receive_message_from_left = 1<<2;
   constexpr static std::size_t receive_message_from_right = 1<<3;

   constexpr static std::size_t always_send_message_to_left = 1<<4 | send_message_to_left;
   constexpr static std::size_t always_send_message_to_right = 1<<5 | send_message_to_right;
   constexpr static std::size_t always_receive_message_from_left = 1<<6 | receive_message_from_left;
   constexpr static std::size_t always_receive_message_from_right = 1<<7 | receive_message_from_right;

   constexpr static std::size_t left = receive_message_from_right | send_message_to_right;
   constexpr static std::size_t right = receive_message_from_left | send_message_to_left; 
   constexpr static std::size_t full = left | right;
   constexpr static std::size_t only_send = send_message_to_left | send_message_to_right;
   constexpr static std::size_t none = 0;

   constexpr static std::size_t always_left = always_receive_message_from_right | always_send_message_to_right;
   constexpr static std::size_t always_right = always_receive_message_from_left | always_send_message_to_left; 

   constexpr bool sends_message_to_left(const std::size_t mps) { return (mps & send_message_to_left) == send_message_to_left; }
   constexpr bool sends_message_to_right(const std::size_t mps) { return (mps & send_message_to_right) == send_message_to_right; }
   constexpr bool receives_message_from_left(const std::size_t mps) { return (mps & receive_message_from_left) == receive_message_from_left; }
   constexpr bool receives_message_from_right(const std::size_t mps) { return (mps & receive_message_from_right) == receive_message_from_right; }

   constexpr bool always_sends_message_to_left(const std::size_t mps) { return (mps & always_send_message_to_left) == always_send_message_to_left; }
   constexpr bool always_sends_message_to_right(const std::size_t mps) { return (mps & always_send_message_to_right) == always_send_message_to_right; }
   constexpr bool always_receives_message_from_left(const std::size_t mps) { return (mps & always_receive_message_from_left) == always_receive_message_from_left; }
   constexpr bool always_receives_message_from_right(const std::size_t mps) { return (mps & always_receive_message_from_right) == always_receive_message_from_right; }

} // namespace message_passing_schedule

namespace message_passing_schedule_factor_view {

   using type = std::size_t;

   constexpr static std::size_t send_message_to_adjacent_factor = 1<<0;
   constexpr static std::size_t send_message_from_adjacent_factor = 1<<1;
   constexpr static std::size_t receive_message_from_adjacent_factor = 1<<2;
   constexpr static std::size_t receive_message_to_adjacent_factor = 1<<3;

   constexpr static std::size_t always_send_message_to_adjacent_factor = 1<<4 | send_message_to_adjacent_factor;
   constexpr static std::size_t always_send_message_from_adjacent_factor = 1<<5 | send_message_from_adjacent_factor;
   constexpr static std::size_t always_receive_message_from_adjacent_factor = 1<<6 | receive_message_from_adjacent_factor;
   constexpr static std::size_t always_receive_message_to_adjacent_factor = 1<<7 | receive_message_to_adjacent_factor;

   constexpr bool sends_message_to_adjacent_factor(const std::size_t mps) { return (mps & send_message_to_adjacent_factor) == send_message_to_adjacent_factor; }
   constexpr bool sends_message_from_adjacent_factor(const std::size_t mps) { return (mps & send_message_from_adjacent_factor) == send_message_from_adjacent_factor; }
   constexpr bool receives_message_from_adjacent_factor(const std::size_t mps) { return (mps & receive_message_from_adjacent_factor) == receive_message_from_adjacent_factor; }
   constexpr bool receives_message_to_adjacent_factor(const std::size_t mps) { return (mps & receive_message_to_adjacent_factor) == receive_message_to_adjacent_factor; }

   constexpr bool always_sends_message_to_adjacent_factor(const std::size_t mps) { return (mps & always_send_message_to_adjacent_factor) == always_send_message_to_adjacent_factor; }
   constexpr bool always_sends_message_from_adjacent_factor(const std::size_t mps) { return (mps & always_send_message_from_adjacent_factor) == always_send_message_from_adjacent_factor; }
   constexpr bool always_receives_message_from_adjacent_factor(const std::size_t mps) { return (mps & always_receive_message_from_adjacent_factor) == always_receive_message_from_adjacent_factor; }
   constexpr bool always_receives_message_to_adjacent_factor(const std::size_t mps) { return (mps & always_receive_message_to_adjacent_factor) == always_receive_message_to_adjacent_factor; }

   constexpr std::size_t right_factor_view(const std::size_t mps) { return mps; }

   constexpr std::size_t left_factor_view(const std::size_t mps)
   {
      // flip every two bits as follows
      // mask 0101....01 and move one left
      const std::size_t mask01 = 1 | (1 << 2) | (1 << 4) | (1 << 6) | (1 << 8);
      const std::size_t left = (mps & mask01) << 1;

      // mask 1010....10 and move one right
      const std::size_t mask10 = (1 << 1) | (1 << 3) | (1 << 5) | (1 << 7) | (1 << 9);
      const std::size_t right = (mps & mask10) >> 1;

      return left | right; 
   }

} // namespace message_passing_schedule_factor_view

} // namespace LPMP

#endif // LPMP_MESSAGE_PASSING_SCHEDULE_HXX
