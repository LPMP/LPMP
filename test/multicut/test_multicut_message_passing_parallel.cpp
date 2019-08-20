#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "multicut/multicut_text_input.h"
#include "multicut/multicut_message_passing_parallel.h"
#include "multicut/multicut_cycle_packing_parallel.h"
#include "test.h"

using namespace LPMP;

int main()
{
   {
      multicut_instance test_instance;
      test_instance.add_edge(0,1,1);
      test_instance.add_edge(0,2,1);
      test_instance.add_edge(1,2,-1.5);
      test_instance.add_edge(1,3,1);
      test_instance.add_edge(2,3,1);

      {
         auto [final_instance, lower_bound] = multicut_message_passing_parallel(test_instance, false, 1);
         const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(final_instance, 1, "non-blocking");
         test(final_instance.evaluate(sol) == 0.0);
         test(lower_bound == 0.0);
      }
      {
         auto [final_instance, lower_bound] = multicut_message_passing_parallel(test_instance, false, 2);
         const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(final_instance, 2, "non-blocking");
         test(final_instance.evaluate(sol) == 0.0);
         test(lower_bound == 0.0);
      }
   }
   {
      multicut_instance test_instance;
      test_instance.add_edge(0,1,1);
      test_instance.add_edge(0,2,1);
      test_instance.add_edge(1,2,-3);
      test_instance.add_edge(1,3,1);
      test_instance.add_edge(2,3,1);

      {
         auto [final_instance, lower_bound] = multicut_message_passing_parallel(test_instance, false, 1);
         const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(final_instance, 1, "non-blocking");
         test(final_instance.evaluate(sol) - (-1) < 1e-10);
         test(lower_bound - (-1) < 1e-10);
      }
      {
         auto [final_instance, lower_bound] = multicut_message_passing_parallel(test_instance, false, 2);
         const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(final_instance, 2, "non-blocking");
         test(final_instance.evaluate(sol) - (-1) < 1e-10);
         test(lower_bound - (-1) < 1e-10);
      }
   }
} 
