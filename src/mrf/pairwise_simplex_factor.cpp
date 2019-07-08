#include "mrf/pairwise_simplex_factor.h"
#include <array>
#include <cassert>

namespace LPMP {

PairwiseSimplexFactor::PairwiseSimplexFactor(const std::size_t _dim1, const std::size_t _dim2) 
   : pairwise_(_dim1, _dim2),
   left_msg_(_dim1),
   right_msg_(_dim2)
{
   std::fill(pairwise_.begin(), pairwise_.end(), 0.0);
   std::fill(left_msg_.begin(), left_msg_.end(), 0.0);
   std::fill(right_msg_.begin(), right_msg_.end(), 0.0);
}

PairwiseSimplexFactor::PairwiseSimplexFactor(const PairwiseSimplexFactor& o) 
   : pairwise_(o.dim1(), o.dim2()),
   left_msg_(o.dim1()),
   right_msg_(o.dim2())
{
   pairwise_ = o.pairwise_;
   left_msg_ = o.left_msg_;
   right_msg_ = o.right_msg_;
}

void PairwiseSimplexFactor::operator=(const PairwiseSimplexFactor& o) {
   assert(dim1() == o.dim1() && dim2() == o.dim2());
   pairwise_ = o.pairwise_;
   left_msg_ = o.left_msg_;
   right_msg_ = o.right_msg_;
}

double PairwiseSimplexFactor::operator[](const std::size_t idx) const 
{
   const auto [x,y] = get_indices(idx);
   return pairwise_(x,y) + left_msg_[x] + right_msg_[y];
}

// below is not nice: two different values, only differ by const!
double PairwiseSimplexFactor::operator()(const std::size_t x1, const std::size_t x2) const 
{
   assert(x1 < dim1() && x2 < dim2());
   return pairwise_(x1,x2) + left_msg_[x1] + right_msg_[x2];
}

double& PairwiseSimplexFactor::cost(const std::size_t x1, const std::size_t x2) 
{
   assert(x1 < dim1() && x2 < dim2());
   return pairwise_(x1,x2);
}

double& PairwiseSimplexFactor::cost(const std::size_t idx) 
{
   const auto [x,y] = get_indices(idx);
   return pairwise_(x,y);
}

// TODO: make more efficient based on SIMD instructions
double PairwiseSimplexFactor::LowerBound() const
{
   double lb = std::numeric_limits<double>::infinity();
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         lb = std::min(lb, (*this)(x1,x2));
      }
   }
   assert(std::isfinite(lb));
   return lb;
}

double PairwiseSimplexFactor::sensitivity() const
{
   double smallest = std::numeric_limits<double>::infinity();
   double second_smallest = std::numeric_limits<double>::infinity();
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         const double val = (*this)(x1,x2);
         const double min = std::min(smallest, val);
         const double max = std::max(smallest, val);
         smallest = min;
         second_smallest = std::min(max, second_smallest);
      }
   }
   assert(std::isfinite(second_smallest));
   return second_smallest - smallest;
} 

std::size_t PairwiseSimplexFactor::dim(const std::size_t d) const 
{ 
   assert(d<2);
   if(d == 0) { return dim1(); }
   else { return dim2(); }
}

void PairwiseSimplexFactor::init_primal() 
{
   primal_[0] = std::numeric_limits<std::size_t>::max();
   primal_[1] = std::numeric_limits<std::size_t>::max();
}

double 
PairwiseSimplexFactor::EvaluatePrimal() const { 
   if(primal_[0] >= dim1() || primal_[1] >= dim2()) {
      return std::numeric_limits<double>::infinity();
   }
   assert(primal_[0] < dim1());
   assert(primal_[1] < dim2());
   return (*this)(primal_[0], primal_[1]); 
}

void 
PairwiseSimplexFactor::MaximizePotentialAndComputePrimal() 
{
   if(primal_[0] >= dim1() && primal_[1] >= dim2()) {
      double min_val = std::numeric_limits<double>::infinity();
      for(std::size_t x1=0; x1<dim1(); ++x1) {
         for(std::size_t x2=0; x2<dim2(); ++x2) {
            if(min_val >= (*this)(x1,x2)) {
               min_val = (*this)(x1,x2);
               primal_[0] = x1;
               primal_[1] = x2;
            }
         }
      }
   } else if(primal_[0] >= dim1() && primal_[1] < dim2()) {
      double min_val = std::numeric_limits<double>::infinity();
      for(std::size_t x1=0; x1<dim1(); ++x1) {
         if(min_val >= (*this)(x1,primal_[1])) {
            min_val = (*this)(x1,primal_[1]);
            primal_[0] = x1;
         }
      } 
   } else if(primal_[1] >= dim2() && primal_[0] < dim1()) {
      double min_val = std::numeric_limits<double>::infinity();
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         if(min_val >= (*this)(primal_[0],x2)) {
            min_val = (*this)(primal_[0],x2);
            primal_[1] = x2;
         }
      }
   } else {
      assert(primal_[0] < dim1() && primal_[1] < dim2());
   }
}

vector<double> 
PairwiseSimplexFactor::min_marginal_1() const
{
   auto min = pairwise_.min1(right_msg_);
   min += left_msg_;
#ifndef NDEBUG
   for(std::size_t x1=0; x1<dim1(); ++x1) {
       double msg_test = std::numeric_limits<double>::infinity();
       for(std::size_t x2=0; x2<dim2(); ++x2)
           msg_test = std::min(msg_test, (*this)(x1,x2));
       assert(std::abs(msg_test - min[x1]) <= eps);
   } 
#endif
   return min;
}

vector<double> 
PairwiseSimplexFactor::min_marginal_2() const
{
   auto min = pairwise_.min2(left_msg_);
   min += right_msg_;
#ifndef NDEBUG
   for(std::size_t x2=0; x2<dim2(); ++x2) {
       double msg_test = std::numeric_limits<double>::infinity();
       for(std::size_t x1=0; x1<dim1(); ++x1)
           msg_test = std::min(msg_test, (*this)(x1,x2));
       assert(std::abs(msg_test - min[x2]) <= eps);
   } 
#endif 
   return min; 
}

std::array<std::size_t,2> PairwiseSimplexFactor::get_indices(const std::size_t idx) const
{
   assert(idx < size());
   const std::size_t x = idx / dim2();
   const std::size_t y = idx % dim2();
   return {x,y};
}
std::size_t PairwiseSimplexFactor::get_index(const std::size_t x, const std::size_t y) const
{
   assert(x < dim1());
   assert(y < dim2());
   return x*dim2() + y;
}

} // namespace LPMP
