#ifndef LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX
#define LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX

#include "LP.h"
#include "vector"
#include <array>
#include <limits>

namespace LPMP {

namespace {
   constexpr static std::size_t primal_not_set = std::numeric_limits<std::size_t>::max(); // primal is not set to any value, i.e. invalid
   constexpr static std::size_t primal_inactive = std::numeric_limits<std::size_t>::max()-1; 
}

// TODO: put detail into separate file
namespace detail {

// compute minimum over configuration with same labels and minimum over sum with different labels
template<typename COST_X, typename COST_Y, typename LABELS_X, typename LABELS_Y>
std::array<double,2> min_diagonal_off_diagonal_sum(const COST_X& cost_x, const COST_Y& cost_y, const LABELS_X& labels_x, const LABELS_Y& labels_y)
{
   // TODO: possibly add fast version for small label space

   assert(std::is_sorted(labels_x.begin(), labels_x.end()));
   assert(std::is_sorted(labels_y.begin(), labels_y.end()));
   assert(cost_x.size() == labels_x.size());
   assert(cost_y.size() == labels_y.size());
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

// Z_ij = \sum_p X_ip Y_pj
class multigraph_matching_triplet_consistency_factor {
public:

   // cost vectors
   REAL cost_z = 0.0;
   vector<REAL> cost_x;
   vector<REAL> cost_y;

   // primal
   std::size_t x, y, z;

   template<typename LABEL_X, typename LABEL_Y>
   multigraph_matching_triplet_consistency_factor(LABEL_X labels_x, LABEL_Y labels_y)
   : multigraph_matching_triplet_consistency_factor(labels_x.begin(), labels_x.end(), labels_y.begin(), labels_y.end())
   {}

   template<typename LABEL_ITERATOR_X, typename LABEL_ITERATOR_Y>
   multigraph_matching_triplet_consistency_factor(LABEL_ITERATOR_X labels_x_begin, LABEL_ITERATOR_X labels_x_end, LABEL_ITERATOR_Y labels_y_begin, LABEL_ITERATOR_Y labels_y_end)
   : 
      cost_z(0.0),
      cost_x(std::distance(labels_x_begin, labels_x_end), 0.0),
      cost_y(std::distance(labels_y_begin, labels_y_end), 0.0),
      labels_x(labels_x_begin, labels_x_end),
      labels_y(labels_y_begin, labels_y_end)
   {
      assert(std::is_sorted(labels_x.begin(), labels_x.end()));
      assert(std::is_sorted(labels_y.begin(), labels_y.end()));
   }

    bool primal_feasible(const std::size_t _x, const std::size_t _y, const std::size_t _z) const
    {
       if(_x == primal_not_set || _y == primal_not_set || _z == primal_not_set) return false;
       if(_y < labels_y.size() && _x < labels_x.size() && _z == 1 && labels_x[_x] != labels_y[_y]) return false;
       const std::size_t no_inactive = std::size_t(_x == primal_inactive) + std::size_t(_y == primal_inactive) + std::size_t(_z == primal_inactive);
       if(no_inactive >= 2) return true;
       if(labels_y[_y] != labels_x[_x] && _z == primal_inactive) return true;
       if(labels_y[_y] == labels_x[_x] && _z == 1) return true;
       return false;
    }

    bool primal_feasible() const { return primal_feasible(x,y,z); }

    REAL EvaluatePrimal() const
    {
       if(primal_feasible())
          return evaluate(x,y,z);
       else
          return std::numeric_limits<REAL>::infinity();
    }

    REAL evaluate(const std::size_t _x, const std::size_t _y, const std::size_t _z) const
    {
       double cost = 0.0;
       if(_x < cost_x.size()) cost += cost_x[_x];
       if(_y < cost_y.size()) cost += cost_y[_y];
       if(_z == 1) cost += cost_z;
       return cost;
    } 

    template<typename FUNC>
    void for_each_labeling(FUNC&& f) const
    {
       f(primal_inactive,primal_inactive,primal_inactive);
       f(primal_inactive, primal_inactive, 1);

       for(std::size_t i=0; i<cost_x.size(); ++i) {
          f(i,primal_inactive,primal_inactive);
       }

       for(std::size_t j=0; j<cost_y.size(); ++j) {
          f(primal_inactive,j,primal_inactive);
       }

       for(std::size_t i=0; i<cost_x.size(); ++i) {
          for(std::size_t j=0; j<cost_y.size(); ++j) {
             if(labels_x[i] != labels_y[j]) {
                f(i,j,primal_inactive);
             } else {
                f(i,j,1);
             } 
          }
       }
    }

    void MaximizePotentialAndComputePrimal()
    {
       double best_cost = std::numeric_limits<double>::infinity();
       std::size_t best_x = primal_not_set;
       std::size_t best_y = primal_not_set;
       std::size_t best_z = primal_not_set;

       auto update_primal = [&](const std::size_t _x, const std::size_t _y, const std::size_t _z) {
          assert(primal_feasible(_x,_y,_z));
          if((x == primal_not_set || x == _x) && (y == primal_not_set || y == _y) && (z == primal_not_set || z == _z)) {
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

    std::array<REAL,2> z_marginals() const
    {
       const auto [min_same_labels, min_different_labels] = detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y);
       const double cost_0 = std::min({min_different_labels, cost_x.min(), cost_y.min(), 0.0});
       const double cost_1 = std::min({cost_z, min_same_labels + cost_z});
       return {cost_0, cost_1};
    }

    vector<REAL> marginals(const vector<REAL>& v1, const vector<REAL>& v2) const
    {
       assert(false);
       // TODO: remove
       assert(v1.size() == v2.size());
       vector<REAL> marginals(v1.size());
       for(std::size_t i=0; i<marginals.size(); ++i) { marginals[i] = v1[i]; }

       for(std::size_t i=0; i<v1.size(); ++i) {
          for(std::size_t j=0; j<v2.size(); ++j) {
             if(i == j)
                marginals[i] = std::min(marginals[i], v1[i] + v2[i] + cost_z);
             else
                marginals[i] = std::min(marginals[i], v1[i] + v2[j]); 
          }
       }

        const REAL cost_not_taken = std::min({0.0, cost_z, v2.min()});
        for(std::size_t i=0; i<marginals.size(); ++i) { marginals[i] -= cost_not_taken; }

        return marginals;
        //
    }

    vector<REAL> x_marginals() const { return marginals(cost_x, cost_y, labels_x, labels_y); }
    vector<REAL> y_marginals() const { return marginals(cost_y, cost_x, labels_y, labels_x); }

    REAL LowerBound() const
    {
        const auto z_marg = z_marginals();
        return std::min(z_marg[0], z_marg[1]);
        // z = multigraph_matching_triplet_consistency_factor::primal_inactive
        const auto x_min = cost_x.min();
        const auto y_min = cost_y.min();
        REAL cost_0 = std::min({0.0, cost_x.min(), cost_y.min()});
        for(std::size_t i=0; i<cost_x.size(); ++i) {
            for(std::size_t j=0; j<cost_y.size(); ++j) {
                if(i != j) {
                    cost_0 = std::min(cost_0, cost_x[i] + cost_y[j]);
                } 
            }
        }

        // z = 1
        REAL cost_1 = std::numeric_limits<REAL>::infinity();
        for(std::size_t i=0; i<cost_x.size(); ++i) {
            cost_1 = std::min(cost_1, cost_x[i] + cost_y[i]); 
        }
        cost_1 += cost_z;

        return std::min(cost_0, cost_1);
    }

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(x,y); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cost_z, cost_x, cost_y ); }

   auto export_variables() { return std::tie(cost_x, cost_y, cost_z); }

   void init_primal() 
   { 
      x = primal_not_set;
      y = primal_not_set;
      z = primal_not_set;
   }

   // if two primal variables are set, we must set the third one accordingly
   void PropagatePrimal()
   {
      const std::size_t no_variables_set = std::size_t(x != primal_not_set) + std::size_t(x != primal_not_set) + std::size_t(z != primal_not_set);
      if(no_variables_set == 2) {
         if(x == std::numeric_limits<std::size_t>::max()) {
            if(z == 1)
               x = y;
         } else if(y == std::numeric_limits<std::size_t>::max()) {
            if(z == 1)
               y = x; 
         } else {
            assert(z == std::numeric_limits<std::size_t>::max());
            if(x == y)
               z = 1;
            else
               z = primal_inactive;
         }
      }
   }

private:

    template<typename COST_1, typename COST_2, typename LABELS_1, typename LABELS_2>
    vector<double> marginals(const COST_1& cost_1, const COST_2& cost_2, const LABELS_1& labels_1, const LABELS_2& labels_2) const
    {
        const auto [min_idx_2, second_min_idx_2] = detail::two_min_indices(cost_2);
        vector<double> marginals(cost_1.size());

        double cost_0 = std::min({0.0, cost_z, cost_2[min_idx_2]});

        std::size_t j=0;
        for(std::size_t i=0; i<cost_1.size(); ++i) {
           if(labels_1[i] == labels_2[min_idx_2]) {
              marginals[i] = std::min({cost_1[i], cost_1[i] + cost_2[second_min_idx_2]});
           } else {
              marginals[i] = std::min({cost_1[i], cost_1[i] + cost_2[min_idx_2]});
           }

           while(j<labels_2.size() && labels_2[j] < labels_1[i]) {
              ++j; 
           }

           if(j<labels_2.size() && labels_1[i] == labels_2[j]) {
              marginals[i] = std::min(marginals[i], cost_1[i] + cost_2[j] + cost_z); 
              ++j;
           }
        }

        for(std::size_t i=0; i<marginals.size(); ++i)
           marginals[i] -= cost_0; 

        return marginals;
    }

   // TODO: add const
   vector<std::size_t> labels_x;
   vector<std::size_t> labels_y; 
};

// when edge is not present, we require X_ip Y_pj = 0
class multigraph_matching_triplet_consistency_factor_zero {
public:

   template<typename LABEL_X, typename LABEL_Y>
   multigraph_matching_triplet_consistency_factor_zero(LABEL_X labels_x, LABEL_Y labels_y)
   : multigraph_matching_triplet_consistency_factor_zero(labels_x.begin(), labels_x.end(), labels_y.begin(), labels_y.end())
   {}

   template<typename LABEL_ITERATOR_X, typename LABEL_ITERATOR_Y>
   multigraph_matching_triplet_consistency_factor_zero(LABEL_ITERATOR_X labels_x_begin, LABEL_ITERATOR_X labels_x_end, LABEL_ITERATOR_Y labels_y_begin, LABEL_ITERATOR_Y labels_y_end)
   : 
      cost_x(std::distance(labels_x_begin, labels_x_end), 0.0),
      cost_y(std::distance(labels_y_begin, labels_y_end), 0.0),
      labels_x(labels_x_begin, labels_x_end),
      labels_y(labels_y_begin, labels_y_end)
   {
      assert(std::is_sorted(labels_x.begin(), labels_x.end()));
      assert(std::is_sorted(labels_y.begin(), labels_y.end()));
   }

   double LowerBound() const
   {
      const double off_diagonal = detail::min_diagonal_off_diagonal_sum(cost_x, cost_y, labels_x, labels_y)[1];
      return std::min({0.0, cost_x.min(), cost_y.min(), off_diagonal});
   }

   double EvaluatePrimal() const
   {
      return evaluate(x,y);
   }

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(x,y); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(cost_x, cost_y); }
   auto export_variables() { return std::tie(cost_x, cost_y); }

   void init_primal() 
   { 
      x = primal_not_set;
      y = primal_not_set;
   }

   std::size_t x, y;
   vector<double> cost_x;
   vector<double> cost_y;

   bool primal_feasible(const std::size_t _x, const std::size_t _y) const
   {
      if(_x == primal_not_set || _y == primal_not_set) return false;
      if(_x == primal_inactive || _y == primal_inactive) return true;
      assert(_x < cost_x.size() && _y < cost_y.size());
      return labels_x[_x] != labels_y[_y];
   }

   bool primal_feasible() const { return primal_feasible(x,y); }

   double evaluate(const std::size_t _x, const std::size_t _y) const
   {
      if(!primal_feasible(_x, _y)) 
         return std::numeric_limits<double>::infinity();

      assert(_x < labels_x.size() || _x == primal_inactive);
      assert(_y < labels_y.size() || _y == primal_inactive);

      if(_x < labels_x.size() && _y < labels_y.size())
         return cost_x[_x] + cost_y[_y];
      else if(_x == primal_inactive && _y < labels_y.size())
         return cost_y[_y];
      else if(_y == primal_inactive && _x < labels_x.size())
         return cost_x[_x];
      else
         return 0.0;
   }

   template<typename FUNC>
   void for_each_labeling(FUNC&& f) const
   { 
      f(primal_inactive,primal_inactive);

      for(std::size_t i=0; i<cost_x.size(); ++i)
         f(i,primal_inactive);

      for(std::size_t j=0; j<cost_y.size(); ++j)
         f(primal_inactive,j);

      for(std::size_t i=0; i<cost_x.size(); ++i)
         for(std::size_t j=0; j<cost_y.size(); ++j)
            if(labels_x[i] != labels_y[j])
               f(i,j);
    }

    void MaximizePotentialAndComputePrimal()
    {
       double best_cost = std::numeric_limits<double>::infinity();
       std::size_t best_x = primal_not_set;
       std::size_t best_y = primal_not_set;

       auto update_primal = [&](const std::size_t _x, const std::size_t _y) {
          assert(primal_feasible(_x,_y));
          if((x == primal_not_set || x == _x) && (y == primal_not_set || y == _y)) {
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

   vector<double> x_marginals() const { return marginals(cost_x, cost_y, labels_x, labels_y); }
   vector<double> y_marginals() const { return marginals(cost_y, cost_x, labels_y, labels_x); }

private:
   template<typename COST_1, typename COST_2, typename LABELS_1, typename LABELS_2>
   static vector<REAL> marginals(const COST_1& cost_1, const COST_2& cost_2, const LABELS_1& labels_1, const LABELS_2& labels_2)
   {
      vector<double> marginals(cost_1.size());
      const auto [min_index_2, second_min_index_2] = detail::two_min_indices(cost_2);

      const double cost_zero = std::min(0.0, cost_2[min_index_2]);
      for(std::size_t i=0; i<cost_1.size(); ++i) {
         if(labels_1[i] == labels_2[min_index_2])
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[second_min_index_2]);
         else
            marginals[i] = std::min(cost_1[i], cost_1[i] + cost_2[min_index_2]);
         marginals[i] -= cost_zero; 
      }
      return marginals;
   }

   // TODO: use this function to accelerate for the dense assignment case
   // in this case we can use faster version for marginal computation
   bool labels_equal() const
   {
      return labels_x.size() == 0 && labels_y.size() == 0;
   }

   const vector<std::size_t> labels_x;
   const vector<std::size_t> labels_y;
};

class simplex_multigraph_matching_triplet_scalar_consistency_message {
public:
    simplex_multigraph_matching_triplet_scalar_consistency_message(const std::size_t idx) : idx_(idx) {} 

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX msg_dim) const
   {
       assert(msg_dim == 0);
       assert(idx_ < l.size()-1);
       l[idx_] += msg;
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX msg_dim) const
   {
       assert(msg_dim == 0);
       r.cost_z += msg;
   }
   
    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
    {
        const auto cost = l[idx_];
        const auto cost_other = l.min_except(idx_);
        msg[0] -= omega*(cost - cost_other);
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
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
             const bool changed = (r.z != primal_inactive);
             r.z = primal_inactive;
             return changed; 
          }
       }
       return false;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal_deactivated(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       if(r.z != std::numeric_limits<std::size_t>::max()) {
          if(r.z == primal_inactive) {
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
   void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX dim) const
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
   void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX dim) const
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
           for(std::size_t i=0; i<r.cost_y.size(); ++i)
               r.cost_y[i] += msg[i]; 
       }
   }
   
    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
    {
        for(std::size_t i=0; i<l.size()-1; ++i) {
            msg[i] -= omega*(l[i] - l.back());
        }
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
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
          if(r.x != primal_inactive) {
             for(std::size_t i=0; i<l.size()-1; ++i) {
                l[i] = std::numeric_limits<REAL>::infinity();
             } 
          } else {
             for(std::size_t i=0; i<l.size(); ++i) {
                if(r.x != primal_inactive)
                   l[i] = std::numeric_limits<REAL>::infinity();
             } 
          }
       } else {
          assert(false); 
       }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR> 
    bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
    {
       const std::size_t right_label = l.primal() == l.size()-1 ? primal_inactive : l.primal();
       std::size_t& right_primal = C == Chirality::left ? r.x : r.y;
       const bool changed = right_primal != right_label;
       right_primal = right_label;
       return changed;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool ComputeLeftFromRightPrimal_deactivated(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       const std::size_t right_label = C == Chirality::left ? r.x : r.y;
       const std::size_t left_label = right_label == primal_inactive ? l.size()-1 : right_label;

       const bool changed = l.primal() != left_label;
       l.primal() = left_label;
       return changed;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
    {
       const std::size_t right_label = C == Chirality::left ? r.x : r.y;
       const std::size_t left_label = right_label == primal_inactive ? l.size()-1 : right_label;
       assert(l.size() - 1 == (C == Chirality::left ? r.cost_x.size() : r.cost_y.size()));

       return l.primal() == left_label;
    }
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_CONSISTENCY_CONSTRAINT_HXX
