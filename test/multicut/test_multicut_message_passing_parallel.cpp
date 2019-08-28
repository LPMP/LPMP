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
      std::vector<edge_item> edge_to_triangle;
      edge_to_triangle.push_back(edge_item{{0,1},1,{{2,0}}});
      edge_to_triangle.push_back(edge_item{{1,2},-1.5,{{0,0},{3,1}}});
      edge_to_triangle.push_back(edge_item{{0,2},1,{{1,0}}});
      edge_to_triangle.push_back(edge_item{{1,3},1,{{2,1}}});
      edge_to_triangle.push_back(edge_item{{2,3},1,{{1,1}}});
      std::vector<triangle_item> triangle_to_edge;
      triangle_to_edge.push_back(triangle_item{{0,1,2},{0,0,0},{0,1,2}});
      triangle_to_edge.push_back(triangle_item{{0,1,2},{0,0,0},{1,4,3}});

      LPMP::send_weights_to_triplets_parallel(taskflow, edge_to_triangle, triangle_to_edge, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      std::vector<edge_t> other_edges = {};
      auto lb1 = compute_lower_bound(other_edges, edge_to_triangle, triangle_to_edge);
      std::cout << "Lower bound after MP step 1: " << lb1 << std::endl;
      taskflow.clear();


      LPMP::send_triplets_to_edge_parallel(taskflow, edge_to_triangle, triangle_to_edge, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      auto lb2 = compute_lower_bound(other_edges, edge_to_triangle, triangle_to_edge);
      std::cout << "Lower bound after MP step 2:" << lb2 << std::endl;
   }
   {
      tf::Taskflow taskflow;
      tf::Executor executor(1);
      std::vector<edge_item> edge_to_triangle;
      edge_to_triangle.push_back(edge_item{{0,1},1,{{2,0}}});
      edge_to_triangle.push_back(edge_item{{1,2},-3,{{0,0},{3,1}}});
      edge_to_triangle.push_back(edge_item{{0,2},1,{{1,0}}});
      edge_to_triangle.push_back(edge_item{{1,3},1,{{2,1}}});
      edge_to_triangle.push_back(edge_item{{2,3},1,{{1,1}}});
      std::vector<triangle_item> triangle_to_edge;
      triangle_to_edge.push_back(triangle_item{{0,1,2},{0,0,0},{0,1,2}});
      triangle_to_edge.push_back(triangle_item{{0,1,2},{0,0,0},{1,4,3}});

      LPMP::send_weights_to_triplets_parallel(taskflow, edge_to_triangle, triangle_to_edge, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      std::vector<edge_t> other_edges = {};
      auto lb1 = compute_lower_bound(other_edges, edge_to_triangle, triangle_to_edge);
      std::cout << "Lower bound after MP step 1: " << lb1 << std::endl;
      taskflow.clear();


      LPMP::send_triplets_to_edge_parallel(taskflow, edge_to_triangle, triangle_to_edge, 1);
      executor.run(taskflow);
      executor.wait_for_all();
      auto lb2 = compute_lower_bound(other_edges, edge_to_triangle, triangle_to_edge);
      std::cout << "Lower bound after MP step 2:" << lb2 << std::endl;
   }
} 
