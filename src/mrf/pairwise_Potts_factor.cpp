#include "mrf/pairwise_Potts_factor.h"

namespace LPMP {

pairwise_potts_factor::pairwise_potts_factor(const std::size_t dim1, const std::size_t dim2)
: pairwise_potts_factor(dim1, double(0.0))
{
   assert(dim1 == dim2);
}

pairwise_potts_factor::pairwise_potts_factor(const std::size_t dim1, const std::size_t dim2, const double diff_cost)
   : pairwise_potts_factor(dim1, diff_cost)
{
   assert(dim1 == dim2); 
} 

pairwise_potts_factor::pairwise_potts_factor(const std::size_t dim, const double diff_cost)
   : 
      vector<double>(2*dim,0.0), // holds reparametrisation
      diff_cost_(diff_cost) 
{}

double pairwise_potts_factor::operator()(const std::size_t x1, const std::size_t x2) const
{
   const double msg_val = msg1(x1) + msg2(x2);
   const double potts_val = x1 != x2 ? diff_cost() : 0.0;
   return msg_val + potts_val;
}

// min cost for same label and different label
std::array<double,2> pairwise_potts_factor::min_values() const
{
   const auto smallest2 = two_smallest_elements<double>(msg2_begin(), msg2_end());

   double min_same_label = std::numeric_limits<double>::infinity();
   double min_diff_label = std::numeric_limits<double>::infinity();
   for(std::size_t i=0; i<dim(); ++i) {
      const double same_label = (*this)[i] + (*this)[i+dim()];
      min_same_label = std::min(min_same_label, same_label);
      const double diff_label = (*this)[i] + diff_cost() + ((*this)[i+dim()] == smallest2[0] ? smallest2[1] : smallest2[0]);
      min_diff_label = std::min(min_diff_label, diff_label);
   }
   return {min_same_label, min_diff_label}; 
}

double pairwise_potts_factor::LowerBound() const {
   const auto v = min_values();
   return std::min(v[0], v[1]);
}

double pairwise_potts_factor::EvaluatePrimal() const
{
   if(primal_[0] < dim() && primal_[1] < dim()) {
      const double pairwise_cost = primal_[0] == primal_[1] ? 0.0 : diff_cost_;
      const double msg_cost = (*this)[primal_[0]] + (*this)[primal_[1] + dim()];
      return pairwise_cost + msg_cost; 
   } else {
      return std::numeric_limits<double>::infinity();
   }
}

void pairwise_potts_factor::MaximizePotentialAndComputePrimal()
{
   assert(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] == std::numeric_limits<std::size_t>::max());
   // TODO: make efficient
   if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] == std::numeric_limits<std::size_t>::max()) {
      double min_cost = std::numeric_limits<double>::infinity();
      for(std::size_t i=0; i<dim(); ++i) {
         for(std::size_t j=0; j<dim(); ++j) {
            if((*this)(i,j) < min_cost) {
               min_cost = (*this)(i,j);
               primal_ = {i,j};
            }
         }
      }
   } else if(primal_[0] < dim1() && primal_[1] == std::numeric_limits<std::size_t>::max()) {
      double min_cost = std::numeric_limits<double>::infinity();
      for(std::size_t j=0; j<dim2(); ++j) {
         if((*this)(primal_[0],j) < min_cost) {
            min_cost = (*this)(primal_[0],j);
            primal_[1] = j;
         }
      } 
   } else if(primal_[1] < dim2() && primal_[0] == std::numeric_limits<std::size_t>::max()) {
      double min_cost = std::numeric_limits<double>::infinity();
      for(std::size_t i=0; i<dim1(); ++i) {
         if((*this)(i,primal_[1]) < min_cost) {
            min_cost = (*this)(i,primal_[1]);
            primal_[0] = i;
         }
      }
   }
}

vector<double> pairwise_potts_factor::min_marginal_1() const
{
   vector<double> m(dim());

   const auto smallest2 = two_smallest_elements<double>(msg2_begin(), msg2_end());

   for(std::size_t i=0; i<dim(); ++i) {
      const double same_label = (*this)[i+dim()];
      const double diff_label = diff_cost() + ((*this)[i+dim()] == smallest2[0] ? smallest2[1] : smallest2[0]);
      m[i] = (*this)[i] + std::min(same_label, diff_label); 
   } 

   return m;
}

vector<double> pairwise_potts_factor::min_marginal_2() const
{
   vector<double> m(dim());

   const auto smallest2 = two_smallest_elements<double>(msg1_begin(), msg1_end());

   for(std::size_t i=0; i<dim(); ++i) {
      const double same_label = (*this)[i];
      const double diff_label = diff_cost() + ((*this)[i] == smallest2[0] ? smallest2[1] : smallest2[0]);
      m[i] = (*this)[i+dim()] + std::min(same_label, diff_label); 
   } 
   return m;
}

double pairwise_potts_factor::min_marginal_cut() const
{
   const auto v = min_values();
   return v[1] - v[0];
   //return v[0] - v[1];
}

} // namespace LPMP
