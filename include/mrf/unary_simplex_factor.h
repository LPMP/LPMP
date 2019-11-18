#pragma once

#include "vector.hxx"

namespace LPMP {

class UnarySimplexFactor : public vector<double> {
public:
    using vector<double>::vector;

    UnarySimplexFactor(std::size_t dim) : vector(dim, 0.0) {}

    template<typename VECTOR>
    UnarySimplexFactor(const VECTOR& vec);

   void print_potential();
   double LowerBound() const;
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

   // on which dimensions of cost does the current solution act? Out of this action the subgradient and the dot product between cost and subgradient can be computed
   // possibly this can be used for primal evaluation, making EvaluatePrimal superfluous.
   template<typename ARRAY>
   void apply(ARRAY& a) const;
   template<typename EXTERNAL_SOLVER>
   void construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v) const;
   template<typename EXTERNAL_SOLVER>
   void convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v);

private:
   std::size_t primal_;
};

template<typename VECTOR>
UnarySimplexFactor::UnarySimplexFactor(const VECTOR& vec)
   : vector(vec.size())
{
   for(std::size_t i=0; i<vec.size(); ++i) { (*this)[i] = vec[i]; }
}

inline void UnarySimplexFactor::print_potential()
{
   for(std::size_t i=0; i<size(); ++i) { std::cout << (*this)[i] << ", "; } 
   std::cout << "\n";
}

inline double UnarySimplexFactor::LowerBound() const { 
   const double lb = this->min();
   assert(std::isfinite(lb));
   return lb;
}

inline double UnarySimplexFactor::EvaluatePrimal() const 
{ 
   if(primal_ >= size()) {
      return std::numeric_limits<double>::infinity();
   }
   return (*this)[primal_]; 
}

inline void UnarySimplexFactor::MaximizePotentialAndComputePrimal() 
{
   if(primal_ >= size()) {
      primal_ = std::min_element(this->begin(), this->end()) - this->begin();
      assert(primal_ < size());
   }
}

template<typename ARRAY>
void UnarySimplexFactor::apply(ARRAY& a) const 
{ 
   assert(primal_ < this->size());
   a[primal_];
}

// example for constructing constraints
// note: this function can also be used to reduce sat constraints: Implement a second solver that implements this in terms of external_solver_interface functions!
// the same holds for LP-interface
template<typename EXTERNAL_SOLVER>
void UnarySimplexFactor::construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v) const
{
   s.add_simplex_constraint(v.begin(), v.end());
}

template<typename EXTERNAL_SOLVER>
void UnarySimplexFactor::convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v)
{
   primal_ = s.first_active(v.begin(), v.end());
}

} // namespace LPMP
