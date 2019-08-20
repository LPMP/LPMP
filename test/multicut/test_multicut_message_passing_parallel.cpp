#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "../src/multicut/multicut_message_passing_parallel.cpp"
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
   {
      tf::Taskflow taskflow;
      tf::Executor executor(1);
      atomic_edge_container edges = {
         {{0,1},1},
         {{1,2},-1.5},
         {{0,2},1},
         {{1,3},1},
         {{2,3},1}
      };
      edge_to_triangle_map M = {
         {{0,1},{2}},
         {{0,2},{1}},
         {{1,2},{0,3}},
         {{1,3},{2}},
         {{2,3},{1}},
      };
      multicut_triangle_factor T = {
         {{0,1,2},{0,0,0}},
         {{1,2,3},{0,0,0}}
      };

      LPMP::send_weights_to_triplets_parallel(taskflow, edges, T, M, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      auto lb1 = compute_lower_bound(edges, T);
      std::cout << "Lower bound after MP step 1: " << lb1 << std::endl;
      taskflow.clear();


      LPMP::send_triplets_to_edge_parallel(taskflow, edges, T, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      auto lb2 = compute_lower_bound(edges, T);
      std::cout << "Lower bound after MP step 2:" << lb2 << std::endl;
   }
   {
      tf::Taskflow taskflow;
      tf::Executor executor(1);
      atomic_edge_container edges = {
         {{0,1},1},
         {{1,2},-3},
         {{0,2},1},
         {{1,3},1},
         {{2,3},1}
      };
      edge_to_triangle_map M = {
         {{0,1},{2}},
         {{0,2},{1}},
         {{1,2},{0,3}},
         {{1,3},{2}},
         {{2,3},{1}},
      };
      multicut_triangle_factor T = {
         {{0,1,2},{0,0,0}},
         {{1,2,3},{0,0,0}}
      };

      LPMP::send_weights_to_triplets_parallel(taskflow, edges, T, M, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      auto lb1 = compute_lower_bound(edges, T);
      std::cout << "Lower bound after MP step 1: " << lb1 << std::endl;
      taskflow.clear();


      LPMP::send_triplets_to_edge_parallel(taskflow, edges, T, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      auto lb2 = compute_lower_bound(edges, T);
      std::cout << "Lower bound after MP step 2:" << lb2 << std::endl;
   }
} 
