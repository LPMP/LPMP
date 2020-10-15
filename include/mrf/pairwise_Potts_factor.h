#ifndef LPMP_PAIRWISE_POTTS_FACTOR_H
#define LPMP_PAIRWISE_POTTS_FACTOR_H

#include <array>
#include "vector.hxx"

namespace LPMP {

class pairwise_potts_factor : public vector<double> {
   template<typename MATRIX>
   static bool is_potts(const MATRIX& c);

public:
   template<typename VEC>
   pairwise_potts_factor(const std::size_t dim1, const std::size_t dim2, const VEC& cost);
   pairwise_potts_factor(const std::size_t dim1, const std::size_t dim2);
   pairwise_potts_factor(const std::size_t dim1, const std::size_t dim2, const double diff_cost);
   pairwise_potts_factor(const std::size_t dim, const double diff_cost);
   double operator()(const std::size_t x1, const std::size_t x2) const;

   // min cost for same label and different label
   std::array<double,2> min_values() const;

   double LowerBound() const;
   double EvaluatePrimal() const;
   void MaximizePotentialAndComputePrimal();

   vector<double> min_marginal_1() const;
   vector<double> min_marginal_2() const;
   double min_marginal_cut() const;

   std::size_t dim() const { return this->size()/2; }
   std::size_t dim1() const { return dim(); }
   std::size_t dim2() const { return dim(); }

   double* msg1_begin() const { return this->begin(); }
   double* msg1_end() const { return this->begin() + dim(); }
   double msg1(const std::size_t x1) const { assert(x1 < dim()); return (*this)[x1]; }
   double& msg1(const std::size_t x1) { assert(x1 < dim()); return (*this)[x1]; }

   double* msg2_begin() const { return this->begin() + dim(); }
   double* msg2_end() const { return this->end(); }
   double msg2(const std::size_t x2) const { assert(x2 < dim()); return (*this)[dim() + x2]; }
   double& msg2(const std::size_t x2) { assert(x2 < dim()); return (*this)[dim() + x2]; }

   double diff_cost() const { return diff_cost_; }
   double& diff_cost() { return diff_cost_; }

   void init_primal() { primal_[0] = std::numeric_limits<std::size_t>::max(); primal_[1] = std::numeric_limits<std::size_t>::max(); }
   auto& primal() { return primal_; }
   const auto& primal() const { return primal_; }

   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<double>*>(this), diff_cost_ ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

   auto export_variables() { return std::tie( *static_cast<vector<double>*>(this), diff_cost_ ); }

protected:
   double diff_cost_;
   std::array<std::size_t,2> primal_;
};


template<typename MATRIX>
bool pairwise_potts_factor::is_potts(const MATRIX& c)
{
   if(c.dim1() != c.dim2()) { return false; }
   if(c.dim1() <= 1) { return false; }
   for(std::size_t x1=0; x1<c.dim1(); ++x1) {
      for(std::size_t x2=0; x2<c.dim2(); ++x2) {
         if(x1 == x2) {
            if(c(x1,x2) != 0.0) return false;
         } else {
            if(c(x1,x2) != c(0,1)) return false;
         }
      }
   }
   return true;
}

template<typename VEC>
pairwise_potts_factor::pairwise_potts_factor(const std::size_t dim1, const std::size_t dim2, const VEC& cost)
   : pairwise_potts_factor(dim1, cost(0,1))
{
   assert(dim1 == dim2);
   assert(is_potts(cost));
   assert(cost(0,0) == 0.0);
}

} // namespace LPMP

#endif // LPMP_PAIRWISE_POTTS_FACTOR_H
