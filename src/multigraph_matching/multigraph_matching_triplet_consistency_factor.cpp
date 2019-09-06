#include "multigraph_matching/multigraph_matching_triplet_consistency_factor.h"
#include "vector"
#include <array>
#include <limits>
#include <cassert>

namespace LPMP {

// TODO: put detail into separate file
namespace detail {

template<typename COST_X, typename COST_Y, typename LABELS_X, typename LABELS_Y>
std::array<double,2> min_diagonal_off_diagonal_sum_naive(const COST_X& cost_x, const COST_Y& cost_y, const LABELS_X& labels_x, const LABELS_Y& labels_y)
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

// compute minimum over configuration with same labels and minimum over sum with different labels
template<typename COST_X, typename COST_Y, typename LABELS_X, typename LABELS_Y>
std::array<double,2> min_diagonal_off_diagonal_sum(const COST_X& cost_x, const COST_Y& cost_y, const LABELS_X& labels_x, const LABELS_Y& labels_y)
{
   // TODO: possibly add fast version for small label space
   assert(std::is_sorted(labels_x.begin(), labels_x.end()));
   assert(std::is_sorted(labels_y.begin(), labels_y.end()));
   assert(cost_x.size() == labels_x.size());
   assert(cost_y.size() == labels_y.size());

   if(cost_x.size() == 1) {
      return min_diagonal_off_diagonal_sum_naive(cost_x, cost_y, labels_x, labels_y); 
   } else if(cost_y.size() == 1) {
      return min_diagonal_off_diagonal_sum_naive(cost_y, cost_x, labels_y, labels_x); 
   }

   assert(cost_x.size() >= 2);
   assert(cost_y.size() >= 2);

   std::size_t min_index_x = std::numeric_limits<std::size_t>::max();
   std::size_t second_min_index_x = std::numeric_limits<std::size_t>::max();
   double min_x = std::numeric_limits<double>::infinity();
   double second_min_x = std::numeric_limits<double>::infinity();

   std::size_t min_index_y = std::numeric_limits<std::size_t>::max();
   std::size_t second_min_index_y = std::numeric_limits<std::size_t>::max();
   double min_y = std::numeric_limits<double>::infinity();
   double second_min_y = std::numeric_limits<double>::infinity();

   double min_same_labels = std::numeric_limits<double>::infinity();

   auto update_minima = [](std::size_t& min_index, std::size_t& second_min_index, double& min_value, double& second_min_value, const double value, const std::size_t index)
   {
      if(min_value >= value) {
         second_min_value = min_value;
         min_value = value;
         second_min_index = min_index;
         min_index = index;
      } else if(second_min_value >= value) {
         second_min_value = value;
         second_min_index = index;
      }
   };

   std::size_t i=0;
   std::size_t j=0;
   while(i<cost_x.size() && j<cost_y.size()) {
      if(labels_x[i] == labels_y[j]) {
         min_same_labels = std::min(min_same_labels, cost_x[i] + cost_y[j]);
         update_minima(min_index_x, second_min_index_x, min_x, second_min_x, cost_x[i], i);
         update_minima(min_index_y, second_min_index_y, min_y, second_min_y, cost_y[j], j);
         ++i;
         ++j;
      } else if(labels_x[i] < labels_y[j]) {
         update_minima(min_index_x, second_min_index_x, min_x, second_min_x, cost_x[i], i);
         ++i; 
      }  else {
         assert(labels_x[i] > labels_y[j]);
         update_minima(min_index_y, second_min_index_y, min_y, second_min_y, cost_y[j], j);
         ++j;
      } 
   }
   assert(i == cost_x.size() || j == cost_y.size());

   for(;i<cost_x.size(); ++i) {
      update_minima(min_index_x, second_min_index_x, min_x, second_min_x, cost_x[i], i);
   }
   for(;j<cost_y.size(); ++j) {
      update_minima(min_index_y, second_min_index_y, min_y, second_min_y, cost_y[j], j);
   }

   assert(cost_x[min_index_x] == *std::min_element(cost_x.begin(), cost_x.end()));
   assert(cost_x[min_index_x] <= cost_x[second_min_index_x]);
   assert(cost_y[min_index_y] == *std::min_element(cost_y.begin(), cost_y.end()));
   assert(cost_y[min_index_y] <= cost_y[second_min_index_y]);

   double min_different_labels = std::numeric_limits<double>::infinity();
   if(labels_x[min_index_x] != labels_y[min_index_y])
      min_different_labels = cost_x[min_index_x] + cost_y[min_index_y];
   else
      min_different_labels = std::min({cost_x[min_index_x] + cost_y[second_min_index_y], cost_x[second_min_index_x] + cost_y[min_index_y]});

   return {min_same_labels, min_different_labels};
}

// efficient version of above when labels of x and y are identical
template<typename COST_X, typename COST_Y>
std::array<double,2> min_diagonal_off_diagonal_sum(const COST_X& cost_x, const COST_Y& cost_y)
{
   assert(cost_x.size() == cost_y.size());
   assert(cost_x.size() > 0);

   double min_different_labels = std::numeric_limits<double>::infinity();
   double min_x = cost_x[0];
   double min_y = cost_y[0];

   double min_same_labels = cost_x[0] + cost_y[0];

   for(std::size_t i=1; i<cost_x.size(); ++i) {
      min_different_labels = std::min({min_different_labels, cost_x[i] + min_y, min_x + cost_y[i]}); 
      min_x = std::min(min_x, cost_x[i]);
      min_y = std::min(min_y, cost_y[i]);

      min_same_labels = std::min(min_same_labels, cost_x[i] + cost_y[i]);
   }

   return {min_same_labels, min_different_labels};
} 

template<typename VECTOR>
std::array<std::size_t,2> two_min_indices(const VECTOR& vec)
{
   assert(vec.size() >= 2);
   std::size_t min_index = std::numeric_limits<std::size_t>::max();
   double min_value = std::numeric_limits<double>::infinity();
   std::size_t second_min_index = std::numeric_limits<std::size_t>::max();
   double second_min_value = std::numeric_limits<double>::infinity();

   for(std::size_t i=0; i<vec.size(); ++i) {
      if(min_value >= vec[i]) {
         second_min_value = min_value;
         min_value = vec[i];
         second_min_index = min_index;
         min_index = i;
      } else if(second_min_value >= vec[i]) {
         second_min_value = vec[i];
         second_min_index = i;
      }
   }

   assert(min_index < vec.size());
   assert(second_min_index < vec.size());
   return {min_index, second_min_index};
}

} // namespace detail

// multigraph_matching_triplet_consistency_factor

bool multigraph_matching_triplet_consistency_factor::primal_feasible(const std::size_t _x, const std::size_t _y, const std::size_t _z) const
{
   if(_x == multigraph_matching_primal_not_set || _y == multigraph_matching_primal_not_set || _z == multigraph_matching_primal_not_set) return false;
   if(_y < labels_y.size() && _x < labels_x.size() && _z == 1 && labels_x[_x] != labels_y[_y]) return false;
   const std::size_t no_inactive = std::size_t(_x == multigraph_matching_primal_inactive) + std::size_t(_y == multigraph_matching_primal_inactive) + std::size_t(_z == multigraph_matching_primal_inactive);
   if(no_inactive >= 2) return true;
   if(labels_y[_y] != labels_x[_x] && _z == multigraph_matching_primal_inactive) return true;
   if(labels_y[_y] == labels_x[_x] && _z == 1) return true;
   return false;
}

double multigraph_matching_triplet_consistency_factor::EvaluatePrimal() const
{
   if(primal_feasible())
      return evaluate(x,y,z);
   else
      return std::numeric_limits<double>::infinity();
}

double multigraph_matching_triplet_consistency_factor::evaluate(const std::size_t _x, const std::size_t _y, const std::size_t _z) const
{
   assert(primal_feasible(_x,_y,_z));
   double cost = 0.0;
   if(_x < cost_x.size()) cost += cost_x[_x];
   if(_y < cost_y.size()) cost += cost_y[_y];
   if(_z == 1) cost += cost_z;
   return cost;
} 

void multigraph_matching_triplet_consistency_factor::MaximizePotentialAndComputePrimal()
{
   double best_cost = std::numeric_limits<double>::infinity();
   std::size_t best_x = multigraph_matching_primal_not_set;
   std::size_t best_y = multigraph_matching_primal_not_set;
   std::size_t best_z = multigraph_matching_primal_not_set;

   auto update_primal = [&](const std::size_t _x, const std::size_t _y, const std::size_t _z) {
      assert(primal_feasible(_x,_y,_z));
      if((x == multigraph_matching_primal_not_set || x == _x) && (y == multigraph_matching_primal_not_set || y == _y) && (z == multigraph_matching_primal_not_set || z == _z)) {
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

std::array<double,2> multigraph_matching_triplet_consistency_factor::z_marginals() const
{
   const auto [min_same_labels, min_different_labels] = detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y);
   const double cost_0 = std::min({min_different_labels, cost_x.min(), cost_y.min(), 0.0});
   const double cost_1 = std::min({cost_z, min_same_labels + cost_z});
   return {cost_0, cost_1};
}

vector<double> multigraph_matching_triplet_consistency_factor::marginals(const vector<double>& v1, const vector<double>& v2) const
{
   assert(false);
   // TODO: remove
   assert(v1.size() == v2.size());
   vector<double> marginals(v1.size());
   for(std::size_t i=0; i<marginals.size(); ++i) { marginals[i] = v1[i]; }

   for(std::size_t i=0; i<v1.size(); ++i) {
      for(std::size_t j=0; j<v2.size(); ++j) {
         if(i == j)
            marginals[i] = std::min(marginals[i], v1[i] + v2[i] + cost_z);
         else
            marginals[i] = std::min(marginals[i], v1[i] + v2[j]); 
      }
   }

   const double cost_not_taken = std::min({0.0, cost_z, v2.min()});
   for(std::size_t i=0; i<marginals.size(); ++i) { marginals[i] -= cost_not_taken; }

   return marginals;
   //
}

double multigraph_matching_triplet_consistency_factor::LowerBound() const
{
   const auto z_marg = z_marginals();
#ifndef NDEBUG
   const auto x_min = cost_x.min();
   const auto y_min = cost_y.min();
   // cost for z ==0/1
   double cost_0 = std::min({0.0, cost_x.min(), cost_y.min()});
   double cost_1 = 0.0;
   for(std::size_t i=0; i<cost_x.size(); ++i) {
      for(std::size_t j=0; j<cost_y.size(); ++j) {
         if(labels_x[i] != labels_y[j]) {
            cost_0 = std::min(cost_0, cost_x[i] + cost_y[j]);
         } else {
            cost_1 = std::min(cost_1, cost_x[i] + cost_y[j]); 
         }
      }
   }
   cost_1 += cost_z;

   assert(std::min(cost_0, cost_1) == std::min(z_marg[0], z_marg[1]));
#endif
   return std::min(z_marg[0], z_marg[1]);
}

void multigraph_matching_triplet_consistency_factor::init_primal() 
{
   x = multigraph_matching_primal_not_set;
   y = multigraph_matching_primal_not_set;
   z = multigraph_matching_primal_not_set;
}

// if two primal variables are set, we must set the third one accordingly
void multigraph_matching_triplet_consistency_factor::PropagatePrimal()
{
   return; //TODO: activate again
   const std::size_t no_variables_set = std::size_t(x != multigraph_matching_primal_not_set) + std::size_t(x != multigraph_matching_primal_not_set) + std::size_t(z != multigraph_matching_primal_not_set);
   if(no_variables_set == 2) {
      if(x == multigraph_matching_primal_not_set) {
         if(z == 1)
            x = y;
      } else if(y == multigraph_matching_primal_not_set) {
         if(z == 1)
            y = x; 
      } else {
         assert(z == multigraph_matching_primal_not_set);
         if(labels_x[x] == labels_y[y])
            z = 1;
         else
            z = multigraph_matching_primal_inactive;
      }
   }
}

vector<double> multigraph_matching_triplet_consistency_factor::marginals(const vector<double>& cost_1, const vector<double>& cost_2, const vector<std::size_t>& labels_1, const vector<std::size_t>& labels_2) const
{
   vector<double> marginals(cost_1.size());
   if(cost_2.size() == 1) {

      const double cost_0 = std::min({0.0, cost_z, cost_2[0]});
      for(std::size_t i=0; i<cost_1.size(); ++i) {
         if(labels_1[i] == labels_2[0])
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[0] + cost_z);
         else
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[0]);
         marginals[i] -= cost_0; 
      }

   } else {

      const auto [min_idx_2, second_min_idx_2] = detail::two_min_indices(cost_2);

      double cost_0 = std::min({0.0, cost_z, cost_2[min_idx_2]});

      std::size_t j=0;
      for(std::size_t i=0; i<cost_1.size(); ++i) {
         if(labels_1[i] == labels_2[min_idx_2]) 
            marginals[i] = std::min({cost_1[i], cost_1[i] + cost_2[second_min_idx_2]});
         else
            marginals[i] = std::min({cost_1[i], cost_1[i] + cost_2[min_idx_2]});

         while(j<labels_2.size() && labels_2[j] < labels_1[i])
            ++j; 

         if(j<labels_2.size() && labels_1[i] == labels_2[j]) {
            marginals[i] = std::min(marginals[i], cost_1[i] + cost_2[j] + cost_z); 
            ++j;
         }
         marginals[i] -= cost_0; 
      }

   }
   return marginals;
}

// multigraph_matching_triplet_consistency_factor_zero

double multigraph_matching_triplet_consistency_factor_zero::LowerBound() const
{
   const double off_diagonal = detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y)[1];
   return std::min({0.0, cost_x.min(), cost_y.min(), off_diagonal});
}

double multigraph_matching_triplet_consistency_factor_zero::EvaluatePrimal() const
{
   return evaluate(x,y);
}

void multigraph_matching_triplet_consistency_factor_zero::init_primal() 
{ 
   x = multigraph_matching_primal_not_set;
   y = multigraph_matching_primal_not_set;
}

bool multigraph_matching_triplet_consistency_factor_zero::primal_feasible(const std::size_t _x, const std::size_t _y) const
{
   if(_x == multigraph_matching_primal_not_set || _y == multigraph_matching_primal_not_set) return false;
   if(_x == multigraph_matching_primal_inactive || _y == multigraph_matching_primal_inactive) return true;
   assert(_x < cost_x.size() && _y < cost_y.size());
   return labels_x[_x] != labels_y[_y];
}

double multigraph_matching_triplet_consistency_factor_zero::evaluate(const std::size_t _x, const std::size_t _y) const
{
   if(!primal_feasible(_x, _y)) 
      return std::numeric_limits<double>::infinity();

   assert(_x < labels_x.size() || _x == multigraph_matching_primal_inactive);
   assert(_y < labels_y.size() || _y == multigraph_matching_primal_inactive);

   if(_x < labels_x.size() && _y < labels_y.size())
      return cost_x[_x] + cost_y[_y];
   else if(_x == multigraph_matching_primal_inactive && _y < labels_y.size())
      return cost_y[_y];
   else if(_y == multigraph_matching_primal_inactive && _x < labels_x.size())
      return cost_x[_x];
   else
      return 0.0;
}

void multigraph_matching_triplet_consistency_factor_zero::MaximizePotentialAndComputePrimal()
{
   double best_cost = std::numeric_limits<double>::infinity();
   std::size_t best_x = multigraph_matching_primal_not_set;
   std::size_t best_y = multigraph_matching_primal_not_set;

   auto update_primal = [&](const std::size_t _x, const std::size_t _y) {
      assert(primal_feasible(_x,_y));
      if((x == multigraph_matching_primal_not_set || x == _x) && (y == multigraph_matching_primal_not_set || y == _y)) {
         const double cost = evaluate(_x,_y);
         if(cost < best_cost) {
            best_cost = cost;
            best_x = _x;
            best_y = _y;
         }
      }
   };

   for_each_labeling(update_primal);

   assert(primal_feasible(best_x,best_y)); 

   x = best_x;
   y = best_y;
}

vector<double> multigraph_matching_triplet_consistency_factor_zero::marginals(const vector<double>& cost_1, const vector<double>& cost_2, const vector<std::size_t>& labels_1, const vector<std::size_t>& labels_2)
{
   vector<double> marginals(cost_1.size());
   if(cost_2.size() == 1) {

      const double cost_zero = std::min(0.0, cost_2[0]);
      for(std::size_t i=0; i<cost_1.size(); ++i) {
         if(labels_1[i] == labels_2[0])
            marginals[i] = cost_1[i];
         else
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[0]);
         marginals[i] -= cost_zero;
      }

   } else {

      const auto [min_index_2, second_min_index_2] = detail::two_min_indices(cost_2);

      const double cost_zero = std::min(0.0, cost_2[min_index_2]);
      for(std::size_t i=0; i<cost_1.size(); ++i) {
         if(labels_1[i] == labels_2[min_index_2])
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[second_min_index_2]);
         else
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[min_index_2]);
         marginals[i] -= cost_zero; 
      }

   }
   return marginals;
}

// TODO: use this function to accelerate for the dense assignment case
// in this case we can use faster version for marginal computation
bool multigraph_matching_triplet_consistency_factor_zero::labels_equal() const
{
   return labels_x.size() == 0 && labels_y.size() == 0;
}

} // namespace LPMP
