#pragma once

#include <array>
#include "vector.hxx"

namespace LPMP {

// factor assumes that triplet potentials is empty and holds only messages to pairwise factors, i.e. is latently factorizable
class SimpleTighteningTernarySimplexFactor : public tensor3_expression<double, SimpleTighteningTernarySimplexFactor> {
public:
   SimpleTighteningTernarySimplexFactor(const std::size_t _dim1, const std::size_t _dim2, const std::size_t _dim3);

   double LowerBound() const;
   double EvaluatePrimal() const;
   void MaximizePotentialAndComputePrimal();
   double operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) const;

   matrix<double> min_marginal12() const;
   matrix<double> min_marginal13() const;
   matrix<double> min_marginal23() const;

   std::size_t size() const { return dim1()*dim2()*dim3(); }
   std::size_t dim(const std::size_t d) const;
   std::size_t dim1() const { return msg12.dim1(); }
   std::size_t dim2() const { return msg23.dim1(); }
   std::size_t dim3() const { return msg13.dim2(); }

   void init_primal();
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( msg12, msg13, msg23 ); }
   auto export_variables() { return std::tie( msg12, msg13, msg23 ); }

   template<typename SOLVER>
   void construct_constraints(SOLVER& s, typename SOLVER::matrix msg12_vars, typename SOLVER::matrix msg13_vars, typename SOLVER::matrix msg23_vars) const;
   template<typename SOLVER>
   void convert_primal(SOLVER& s, typename SOLVER::matrix msg12_vars, typename SOLVER::matrix msg13_vars, typename SOLVER::matrix msg23_vars);

   const std::array<std::size_t,3>& primal() const { return primal_; }
   std::array<std::size_t,3>& primal() { return primal_; }

   matrix<double> msg12;
   matrix<double> msg13;
   matrix<double> msg23;

protected:
   std::array<std::size_t,3> primal_;
};

template<typename SOLVER>
void SimpleTighteningTernarySimplexFactor::construct_constraints(SOLVER& s, typename SOLVER::matrix msg12_vars, typename SOLVER::matrix msg13_vars, typename SOLVER::matrix msg23_vars) const
{
   auto joint_vars = s.add_tensor(dim1(),dim2(),dim3());

   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         auto slice = joint_vars.slice12(x1,x2);
         auto slice_sum = s.add_at_most_one_constraint(slice.begin(), slice.end());
         s.make_equal(msg12_vars(x1,x2), slice_sum);
      }
   }

   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x3=0; x3<dim3(); ++x3) {
         auto slice = joint_vars.slice13(x1,x3);
         auto slice_sum = s.add_at_most_one_constraint(slice.begin(), slice.end());
         s.make_equal(msg13_vars(x1,x3), slice_sum);
      }
   }

   for(std::size_t x2=0; x2<dim2(); ++x2) {
      for(std::size_t x3=0; x3<dim3(); ++x3) {
         auto slice = joint_vars.slice23(x2,x3);
         auto slice_sum = s.add_at_most_one_constraint(slice.begin(), slice.end());
         s.make_equal(msg23_vars(x2,x3), slice_sum);
      }
   }
}

template<typename SOLVER>
void SimpleTighteningTernarySimplexFactor::convert_primal(SOLVER& s, typename SOLVER::matrix msg12_vars, typename SOLVER::matrix msg13_vars, typename SOLVER::matrix msg23_vars)
{
   auto msg12_sol = s.first_active(msg12_vars);
   auto msg13_sol = s.first_active(msg13_vars);
   auto msg23_sol = s.first_active(msg23_vars);
   assert(msg12_sol[0] == msg12_sol[0]);
   assert(msg13_sol[1] == msg23_sol[1]);
   assert(msg12_sol[1] == msg23_sol[0]);

   primal_[0] = msg12_sol[0];
   primal_[1] = msg12_sol[1];
   primal_[2] = msg13_sol[1];
} 

} // namespace LPMP
