#include "graph_matching/multigraph_matching_triplet_consistency_factor.h"
#include "../src/multigraph_matching/multigraph_matching_triplet_consistency_factor.cpp"
#include "generate_random_label_set.hxx"
#include "test/test.h"
#include <array>
#include <vector>
#include <random>
#include <numeric>

using namespace LPMP;

int main(int argc, char** argv)
{
   std::random_device rd{};
   std::mt19937 gen{rd()};
   std::normal_distribution nd(-1.0,2.0);

   // first test with identical labels
   for(std::size_t no_labels=1; no_labels<100; ++no_labels) {
      std::vector<double> cost_x(no_labels);
      std::vector<double> cost_y(no_labels);
      for(auto& x : cost_x) x = nd(gen);
      for(auto& y : cost_y) y = nd(gen);

      std::vector<std::size_t> labels_x(no_labels);
      std::vector<std::size_t> labels_y(no_labels);
      std::iota(labels_x.begin(), labels_x.end(), 0);
      std::iota(labels_y.begin(), labels_y.end(), 0);

      test(detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y) == detail::min_diagonal_off_diagonal_sum_naive(cost_x, cost_y, labels_x, labels_y));
   } 

   // now with varying label sets
   for(std::size_t no_labels_x=1; no_labels_x<100; ++no_labels_x) {
      for(std::size_t no_labels_y=1; no_labels_y<100; ++no_labels_y) {
         std::vector<double> cost_x(no_labels_x);
         std::vector<double> cost_y(no_labels_y);
         for(auto& x : cost_x) x = nd(gen);
         for(auto& y : cost_y) y = nd(gen);

         std::vector<std::size_t> labels_x = generate_random_label_set(no_labels_x,100);
         std::vector<std::size_t> labels_y = generate_random_label_set(no_labels_y,100);

         auto method = detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y);
         auto method_slow = detail::min_diagonal_off_diagonal_sum_naive(cost_x, cost_y, labels_x, labels_y);
         test(method == method_slow);
      }
   }

   // test two min indices computation
   for(std::size_t i=2; i<100; ++i) {
      std::vector<double> vec(i);
      for(double& x : vec) x = nd(gen);
      const auto [min_idx, second_min_idx] = detail::two_min_indices(vec);
      test(vec[min_idx] == *std::min_element(vec.begin(), vec.end()));
      vec[min_idx] = std::numeric_limits<double>::infinity();
      test(vec[second_min_idx] == *std::min_element(vec.begin(), vec.end()));
   }
}
