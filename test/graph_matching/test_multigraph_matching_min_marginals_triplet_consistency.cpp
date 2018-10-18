#include "graph_matching/multigraph_matching_consistency_constraint.hxx"
#include "generate_random_label_set.hxx"
#include "test/test.h"
#include <array>
#include <vector>
#include <random>

using namespace LPMP;

template<typename COST_X, typename COST_Y, typename LABELS_X, typename LABELS_Y>
std::array<double,2> min_diagonal_off_diagonal_sum_slow(const COST_X& cost_x, const COST_Y& cost_y, const LABELS_X& labels_x, const LABELS_Y& labels_y)
{
   double min_same_labels = std::numeric_limits<double>::infinity();
   double min_different_labels = std::numeric_limits<double>::infinity();
   for(std::size_t i=0; i<cost_x.size(); ++i) {
      for(std::size_t j=0; j<cost_y.size(); ++j) {
         if(labels_x[i] == labels_y[j])
            min_same_labels = std::min(min_same_labels, cost_x[i] + cost_y[j]);
         else
            min_different_labels = std::min(min_different_labels, cost_x[i] + cost_y[j]);
      }
   }

   return {min_same_labels, min_different_labels}; 
}

int main(int argc, char** argv)
{
   std::random_device rd{};
   std::mt19937 gen{rd()};
   std::normal_distribution nd(-1.0,2.0);

   // first test with identical labels
   for(std::size_t no_labels=2; no_labels<100; ++no_labels) {
      std::vector<double> cost_x(no_labels);
      std::vector<double> cost_y(no_labels);
      for(auto& x : cost_x) x = nd(gen);
      for(auto& y : cost_y) y = nd(gen);

      std::vector<std::size_t> labels_x(no_labels);
      std::vector<std::size_t> labels_y(no_labels);
      std::iota(labels_x.begin(), labels_x.end(), 0);
      std::iota(labels_y.begin(), labels_y.end(), 0);

      test(detail::min_diagonal_off_diagonal_sum(cost_x, cost_y) == min_diagonal_off_diagonal_sum_slow(cost_x, cost_y, labels_x, labels_y));
      test(detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y) == min_diagonal_off_diagonal_sum_slow(cost_x, cost_y, labels_x, labels_y));
   } 

   // now with varying label sets
   for(std::size_t no_labels_x=2; no_labels_x<100; ++no_labels_x) {
      for(std::size_t no_labels_y=2; no_labels_y<100; ++no_labels_y) {
         std::vector<double> cost_x(no_labels_x);
         std::vector<double> cost_y(no_labels_y);
         for(auto& x : cost_x) x = nd(gen);
         for(auto& y : cost_y) y = nd(gen);

         std::vector<std::size_t> labels_x = generate_random_label_set(no_labels_x,100);
         std::vector<std::size_t> labels_y = generate_random_label_set(no_labels_y,100);

         auto method = detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y);
         auto method_slow = min_diagonal_off_diagonal_sum_slow(cost_x, cost_y, labels_x, labels_y);
         test(method == method_slow);
      }
   }

   for(std::size_t i=2; i<100; ++i) {
      std::vector<double> vec(i);
      for(double& x : vec) x = nd(gen);
      const auto [min_idx, second_min_idx] = detail::two_min_indices(vec);
      test(vec[min_idx] == *std::min_element(vec.begin(), vec.end()));
      vec[min_idx] = std::numeric_limits<double>::infinity();
      test(vec[second_min_idx] == *std::min_element(vec.begin(), vec.end()));
   }
}
