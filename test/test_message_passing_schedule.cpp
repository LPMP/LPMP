#include "test.h"
#include "message_passing_schedule.hxx"

using namespace LPMP;

int main(int argc, char** argv)
{
   // left
   {
      message_passing_schedule::type mps = message_passing_schedule::left;
      test(message_passing_schedule::sends_message_to_right(mps) && message_passing_schedule::receives_message_from_right(mps) && !message_passing_schedule::sends_message_to_left(mps) && !message_passing_schedule::receives_message_from_left(mps));
      test(!message_passing_schedule::always_sends_message_to_right(mps) && !message_passing_schedule::always_receives_message_from_right(mps) && !message_passing_schedule::always_sends_message_to_left(mps) && !message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(!message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      }
   }

   // right
   {
      message_passing_schedule::type mps = message_passing_schedule::right;
      test(!message_passing_schedule::sends_message_to_right(mps) && !message_passing_schedule::receives_message_from_right(mps) && message_passing_schedule::sends_message_to_left(mps) && message_passing_schedule::receives_message_from_left(mps));
      test(!message_passing_schedule::always_sends_message_to_right(mps) && !message_passing_schedule::always_receives_message_from_right(mps) && !message_passing_schedule::always_sends_message_to_left(mps) && !message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(!message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      }
   }

   // full
   {
      message_passing_schedule::type mps = message_passing_schedule::full;
      test(message_passing_schedule::sends_message_to_right(mps) && message_passing_schedule::receives_message_from_right(mps) && message_passing_schedule::sends_message_to_left(mps) && message_passing_schedule::receives_message_from_left(mps));
      test(!message_passing_schedule::always_sends_message_to_right(mps) && !message_passing_schedule::always_receives_message_from_right(mps) && !message_passing_schedule::always_sends_message_to_left(mps) && !message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      } 
   }

   // none
   {
      message_passing_schedule::type mps = message_passing_schedule::none;
      test(!message_passing_schedule::sends_message_to_right(mps) && !message_passing_schedule::receives_message_from_right(mps) && !message_passing_schedule::sends_message_to_left(mps) && !message_passing_schedule::receives_message_from_left(mps));
      test(!message_passing_schedule::always_sends_message_to_right(mps) && !message_passing_schedule::always_receives_message_from_right(mps) && !message_passing_schedule::always_sends_message_to_left(mps) && !message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(!message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(!message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      }
   }

   // only_send
   {
      message_passing_schedule::type mps = message_passing_schedule::only_send;
      test(message_passing_schedule::sends_message_to_right(mps) && !message_passing_schedule::receives_message_from_right(mps) && message_passing_schedule::sends_message_to_left(mps) && !message_passing_schedule::receives_message_from_left(mps));
      test(!message_passing_schedule::always_sends_message_to_right(mps) && !message_passing_schedule::always_receives_message_from_right(mps) && !message_passing_schedule::always_sends_message_to_left(mps) && !message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      }
   }

   // always_left
   {
      message_passing_schedule::type mps = message_passing_schedule::always_left;
      test(message_passing_schedule::sends_message_to_right(mps) && message_passing_schedule::receives_message_from_right(mps) && !message_passing_schedule::sends_message_to_left(mps) && !message_passing_schedule::receives_message_from_left(mps));
      test(message_passing_schedule::always_sends_message_to_right(mps) && message_passing_schedule::always_receives_message_from_right(mps) && !message_passing_schedule::always_sends_message_to_left(mps) && !message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(!message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      }
   }

   // always_right
   {
      message_passing_schedule::type mps = message_passing_schedule::always_right;
      test(!message_passing_schedule::sends_message_to_right(mps) && !message_passing_schedule::receives_message_from_right(mps) && message_passing_schedule::sends_message_to_left(mps) && message_passing_schedule::receives_message_from_left(mps));
      test(!message_passing_schedule::always_sends_message_to_right(mps) && !message_passing_schedule::always_receives_message_from_right(mps) && message_passing_schedule::always_sends_message_to_left(mps) && message_passing_schedule::always_receives_message_from_left(mps));

      {
         message_passing_schedule_factor_view::type mps_left = message_passing_schedule_factor_view::left_factor_view(mps);
         test(!message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_left) && message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_left));
         test(!message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_left) && message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_left) && !message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_left) && message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_left));
      }

      {
         message_passing_schedule_factor_view::type mps_right = message_passing_schedule_factor_view::right_factor_view(mps);
         test(message_passing_schedule_factor_view::sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::sends_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::receives_message_to_adjacent_factor(mps_right));
         test(message_passing_schedule_factor_view::always_sends_message_to_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_sends_message_from_adjacent_factor(mps_right) && message_passing_schedule_factor_view::always_receives_message_from_adjacent_factor(mps_right) && !message_passing_schedule_factor_view::always_receives_message_to_adjacent_factor(mps_right));
      }
   }
}
