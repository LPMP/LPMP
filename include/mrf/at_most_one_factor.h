#pragma once

#include <cassert>
#include "vector.hxx"

namespace LPMP {

class at_most_one_factor : public vector<double> {
public:
   at_most_one_factor(const std::vector<double>& cost) : vector<double>(cost.begin(), cost.end()) {}
   at_most_one_factor(const std::size_t n) : vector<double>(n, 0.0) {}

   double LowerBound() const ;
   double EvaluatePrimal() const;
   void MaximizePotentialAndComputePrimal();

   // load/store function for the primal value
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<double>*>(this) ); }

   auto export_variables() { return std::tie(*static_cast<vector<double>*>(this)); }

   void init_primal() { primal_ = std::numeric_limits<std::size_t>::max(); }
   std::size_t primal() const { return primal_; }
   std::size_t& primal() { return primal_; }
   void primal(const std::size_t p) { primal_ = p; }

   double sensitivity() const;

   template<typename EXTERNAL_SOLVER>
   void construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v) const;
   template<typename EXTERNAL_SOLVER>
   void convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v);

private:
   std::size_t primal_;
};

inline double at_most_one_factor::LowerBound() const 
{ 
   const double lb = this->min();
   assert(std::isfinite(lb));
   return std::min(lb,double(0.0));
}

inline double at_most_one_factor::EvaluatePrimal() const 
{ 
   if(primal_ > size()) {
      return std::numeric_limits<double>::infinity();
   } else if(primal_ == size()) {
      return 0.0;
   } else {
      return (*this)[primal_]; 
   }
}

inline void at_most_one_factor::MaximizePotentialAndComputePrimal() 
{
   if(primal_ > size()) {
      auto min = std::min_element(this->begin(), this->end());
      if(*min < 0.0) {
         primal_ = min - this->begin();
      } else {
         primal_ = this->size();
      }
      assert(primal_ <= size());
   }
}

inline double at_most_one_factor::sensitivity() const
{
   auto minima = two_smallest_elements<double>(this->begin(), this->end());
   return minima[1] - minima[0];
}

template<typename EXTERNAL_SOLVER>
void at_most_one_factor::construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v) const
{
   s.add_at_most_one_constraint(v.begin(), v.end());
}

template<typename EXTERNAL_SOLVER>
void at_most_one_factor::convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v)
{
   primal_ = s.first_active(v.begin(), v.end());
}

} // namespace LPMP
