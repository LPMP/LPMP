#ifndef LPMP_MULTIGRAPH_MATCHING_TRIPLET_CONSISTENCY_FACTOR_H 
#define LPMP_MULTIGRAPH_MATCHING_TRIPLET_CONSISTENCY_FACTOR_H 

#include "config.hxx"
#include "vector.hxx"

namespace LPMP {

constexpr static std::size_t multigraph_matching_primal_not_set = std::numeric_limits<std::size_t>::max(); // primal is not set to any value, i.e. invalid
constexpr static std::size_t multigraph_matching_primal_inactive = std::numeric_limits<std::size_t>::max()-1; 

// Z_ij = \sum_p X_ip Y_pj
class multigraph_matching_triplet_consistency_factor {
public:

   // cost vectors
   double cost_z = 0.0;
   vector<double> cost_x;
   vector<double> cost_y;

   // primal
   std::size_t x, y, z;

   template<typename LABEL_X, typename LABEL_Y>
   multigraph_matching_triplet_consistency_factor(LABEL_X labels_x, LABEL_Y labels_y)
   : multigraph_matching_triplet_consistency_factor(labels_x.begin(), labels_x.end(), labels_y.begin(), labels_y.end())
   {}

   template<typename LABEL_ITERATOR_X, typename LABEL_ITERATOR_Y>
   multigraph_matching_triplet_consistency_factor(LABEL_ITERATOR_X labels_x_begin, LABEL_ITERATOR_X labels_x_end, LABEL_ITERATOR_Y labels_y_begin, LABEL_ITERATOR_Y labels_y_end);

    bool primal_feasible(const std::size_t _x, const std::size_t _y, const std::size_t _z) const;
    bool primal_feasible() const { return primal_feasible(x,y,z); }

    double EvaluatePrimal() const;
    double evaluate(const std::size_t _x, const std::size_t _y, const std::size_t _z) const;

    template<typename FUNC>
    void for_each_labeling(FUNC&& f) const;

    void MaximizePotentialAndComputePrimal();

    double LowerBound() const;

    std::array<double,2> z_marginals() const;
    vector<double> marginals(const vector<double>& v1, const vector<double>& v2) const;
    vector<double> x_marginals() const { return marginals(cost_x, cost_y, labels_x, labels_y); }
    vector<double> y_marginals() const { return marginals(cost_y, cost_x, labels_y, labels_x); }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(x,y); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cost_z, cost_x, cost_y ); }

    auto export_variables() { return std::tie(cost_x, cost_y, cost_z); }

    void init_primal();
    // if two primal variables are set, we must set the third one accordingly
   void PropagatePrimal();

private:

    vector<double> marginals(const vector<double>& cost_1, const vector<double>& cost_2, const vector<std::size_t>& labels_1, const vector<std::size_t>& labels_2) const;

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
   multigraph_matching_triplet_consistency_factor_zero(LABEL_ITERATOR_X labels_x_begin, LABEL_ITERATOR_X labels_x_end, LABEL_ITERATOR_Y labels_y_begin, LABEL_ITERATOR_Y labels_y_end);

   double LowerBound() const;
   double EvaluatePrimal() const;

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(x,y); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(cost_x, cost_y); }
   auto export_variables() { return std::tie(cost_x, cost_y); }

   void init_primal();

   std::size_t x, y;
   vector<double> cost_x;
   vector<double> cost_y;

   bool primal_feasible(const std::size_t _x, const std::size_t _y) const;
   bool primal_feasible() const { return primal_feasible(x,y); }

   double evaluate(const std::size_t _x, const std::size_t _y) const;

   template<typename FUNC>
   void for_each_labeling(FUNC&& f) const;

    void MaximizePotentialAndComputePrimal();

   vector<double> x_marginals() const { return marginals(cost_x, cost_y, labels_x, labels_y); }
   vector<double> y_marginals() const { return marginals(cost_y, cost_x, labels_y, labels_x); }

private:
   static vector<double> marginals(const vector<double>& cost_1, const vector<double>& cost_2, const vector<std::size_t>& labels_1, const vector<std::size_t>& labels_2);

   // TODO: use this function to accelerate for the dense assignment case
   // in this case we can use faster version for marginal computation
   bool labels_equal() const;

   const vector<std::size_t> labels_x;
   const vector<std::size_t> labels_y;
};

template<typename LABEL_ITERATOR_X, typename LABEL_ITERATOR_Y>
multigraph_matching_triplet_consistency_factor::multigraph_matching_triplet_consistency_factor(LABEL_ITERATOR_X labels_x_begin, LABEL_ITERATOR_X labels_x_end, LABEL_ITERATOR_Y labels_y_begin, LABEL_ITERATOR_Y labels_y_end)
   : 
      cost_z(0.0),
      cost_x(std::distance(labels_x_begin, labels_x_end), 0.0),
      cost_y(std::distance(labels_y_begin, labels_y_end), 0.0),
      labels_x(labels_x_begin, labels_x_end),
      labels_y(labels_y_begin, labels_y_end)
{
   assert(*std::max_element(cost_x.begin(), cost_x.end()) == 0.0);
   assert(*std::min_element(cost_x.begin(), cost_x.end()) == 0.0);
   assert(*std::max_element(cost_y.begin(), cost_y.end()) == 0.0);
   assert(*std::min_element(cost_y.begin(), cost_y.end()) == 0.0);
   assert(cost_x.size() == labels_x.size());
   assert(cost_y.size() == labels_y.size());
   assert(std::is_sorted(labels_x.begin(), labels_x.end()));
   assert(std::is_sorted(labels_y.begin(), labels_y.end()));
}

template<typename FUNC>
void multigraph_matching_triplet_consistency_factor::for_each_labeling(FUNC&& f) const
{
   f(multigraph_matching_primal_inactive,multigraph_matching_primal_inactive,multigraph_matching_primal_inactive);
   f(multigraph_matching_primal_inactive, multigraph_matching_primal_inactive, 1);

   for(std::size_t i=0; i<cost_x.size(); ++i) {
      f(i,multigraph_matching_primal_inactive,multigraph_matching_primal_inactive);
   }

   for(std::size_t j=0; j<cost_y.size(); ++j) {
      f(multigraph_matching_primal_inactive,j,multigraph_matching_primal_inactive);
   }

   for(std::size_t i=0; i<cost_x.size(); ++i) {
      for(std::size_t j=0; j<cost_y.size(); ++j) {
         if(labels_x[i] != labels_y[j]) {
            f(i,j,multigraph_matching_primal_inactive);
         } else {
            f(i,j,1);
         } 
      }
   }
} 


template<typename LABEL_ITERATOR_X, typename LABEL_ITERATOR_Y>
multigraph_matching_triplet_consistency_factor_zero::multigraph_matching_triplet_consistency_factor_zero(LABEL_ITERATOR_X labels_x_begin, LABEL_ITERATOR_X labels_x_end, LABEL_ITERATOR_Y labels_y_begin, LABEL_ITERATOR_Y labels_y_end)
   : 
      cost_x(std::distance(labels_x_begin, labels_x_end), 0.0),
      cost_y(std::distance(labels_y_begin, labels_y_end), 0.0),
      labels_x(labels_x_begin, labels_x_end),
      labels_y(labels_y_begin, labels_y_end)
{
   assert(cost_x.size() == labels_x.size());
   assert(cost_y.size() == labels_y.size());

   assert(std::is_sorted(labels_x.begin(), labels_x.end()));
   assert(std::is_sorted(labels_y.begin(), labels_y.end()));
}

template<typename FUNC>
void multigraph_matching_triplet_consistency_factor_zero::for_each_labeling(FUNC&& f) const
{ 
   f(multigraph_matching_primal_inactive,multigraph_matching_primal_inactive);

   for(std::size_t i=0; i<cost_x.size(); ++i)
      f(i,multigraph_matching_primal_inactive);

   for(std::size_t j=0; j<cost_y.size(); ++j)
      f(multigraph_matching_primal_inactive,j);

   for(std::size_t i=0; i<cost_x.size(); ++i)
      for(std::size_t j=0; j<cost_y.size(); ++j)
         if(labels_x[i] != labels_y[j])
            f(i,j);
}

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_TRIPLET_CONSISTENCY_FACTOR_H 

