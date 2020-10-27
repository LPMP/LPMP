#ifndef LPMP_PAIRWISE_SIMPLEX_FACTOR_H
#define LPMP_PAIRWISE_SIMPLEX_FACTOR_H

#include "vector.hxx"

namespace LPMP {

// do zrobienia: if pairwise was supplied to us (e.g. external factor, then reflect this in constructor and only allocate space for messages.
// When tightening, we can simply replace pairwise pointer to external factor with an explicit copy. Reallocate left_msg_ and right_msg_ to make memory contiguous? Not sure, depends whether we use block_allocator, which will not acually release the memory
// when factor is copied, then pairwise_ must only be copied if it is actually modified. This depends on whether we execute SMRP or MPLP style message passing. Templatize for this possibility
class PairwiseSimplexFactor : public matrix_expression<double, PairwiseSimplexFactor> {
public:
   PairwiseSimplexFactor(const std::size_t _dim1, const std::size_t _dim2);
   template<typename MATRIX>
   PairwiseSimplexFactor(const std::size_t dim1, const std::size_t dim2, const MATRIX& m);
   PairwiseSimplexFactor(const PairwiseSimplexFactor& o);

   void operator=(const PairwiseSimplexFactor& o);
   double operator[](const std::size_t x) const;
   double operator()(const std::size_t x1, const std::size_t x2) const;
   double& cost(const std::size_t x1, const std::size_t x2);
   double& cost(const std::size_t idx);
   double LowerBound() const;
   template<std::size_t N>
   double lower_bound_except(const std::array<std::size_t,N> indices) const;
   double sensitivity() const;
   std::size_t size() const { return dim1()*dim2(); }
   std::size_t dim1() const { return left_msg_.size(); }
   std::size_t dim2() const { return right_msg_.size(); }
   std::size_t dim(const std::size_t d) const;
   double& pairwise(const std::size_t x1, const std::size_t x2) { assert(x1<dim1() && x2<dim2()); return pairwise_(x1,x2); }
   double& msg1(const std::size_t x1) { assert(x1<dim1()); return left_msg_[x1]; }
   double& msg2(const std::size_t x2) { assert(x2<dim2()); return right_msg_[x2]; }
   void init_primal();
   double EvaluatePrimal() const;
   void MaximizePotentialAndComputePrimal();

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_[0], primal_[1] ); }
   //template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cereal::binary_data( pairwise_, sizeof(double)*(size()+dim1()+dim2()) ) ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( left_msg_, right_msg_, pairwise_ ); }

   auto export_variables() { return std::tie(left_msg_, right_msg_, pairwise_); }

   vector<double> min_marginal_1() const;
   vector<double> min_marginal_2() const;

   template<typename ARRAY>
   void apply(ARRAY& a) const;
   template<typename SOLVER>
   void construct_constraints(SOLVER& s, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables) const;
   template<typename SOLVER>
   void convert_primal(SOLVER& s, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables);

   const std::array<std::size_t,2>& primal() const { return primal_; }
   std::array<std::size_t,2>& primal() { return primal_; }

private:
   mutable matrix<double> pairwise_;
   std::array<std::size_t,2> get_indices(const std::size_t idx) const;
   std::size_t get_index(const std::size_t x, const std::size_t y) const;
   vector<double> left_msg_;
   vector<double> right_msg_;
   std::array<std::size_t,2> primal_;
};

template<typename MATRIX>
PairwiseSimplexFactor::PairwiseSimplexFactor(const std::size_t dim1, const std::size_t dim2, const MATRIX& m)
: PairwiseSimplexFactor(dim1,dim2)
{
   for(std::size_t x1=0; x1<this->dim1(); ++x1) {
      for(std::size_t x2=0; x2<this->dim2(); ++x2) {
         this->cost(x1,x2) = m(x1,x2);
      }
   }
}

template<typename ARRAY>
void PairwiseSimplexFactor::apply(ARRAY& a) const
{
   assert(primal_[0] < dim1() && primal_[1] < dim2());
   a[primal_[0]];
   a[dim1() + primal_[1]]; 
   a[dim1() + dim2() + primal_[0]*dim2() + primal_[1]];
}

template<typename SOLVER>
void PairwiseSimplexFactor::construct_constraints(SOLVER& s, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables) const
{
   s.add_simplex_constraint(pairwise_variables.begin(), pairwise_variables.end());
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      auto slice = pairwise_variables.slice_left(x1);
      auto c = s.add_at_most_one_constraint(slice.begin(), slice.end());
      s.make_equal(c, left_unary_variables[x1]);
   }

   for(std::size_t x2=0; x2<dim2(); ++x2) {
      auto slice = pairwise_variables.slice_right(x2);
      auto c = s.add_at_most_one_constraint(slice.begin(), slice.end());
      s.make_equal(c, right_unary_variables[x2]);
   }
}

template<typename SOLVER>
void PairwiseSimplexFactor::convert_primal(SOLVER& s, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables)
{
   primal_[0] = s.first_active(left_unary_variables);
   primal_[1] = s.first_active(right_unary_variables);
   auto pairwise_sol = s.first_active(pairwise_variables);
   assert(pairwise_sol[0] == primal_[0]);
   assert(pairwise_sol[1] == primal_[1]);
}

template<std::size_t N>
double PairwiseSimplexFactor::lower_bound_except(const std::array<std::size_t,N> indices) const
{
   std::array<double,N> vars;
   for(std::size_t i=0; i<N;++i) {
      const auto [x,y] = get_indices(indices[i]);
      vars[i] = pairwise_(x,y);
      pairwise_(x,y) = std::numeric_limits<double>::infinity();
   }
   const double lb = LowerBound();
   for(std::size_t i=0; i<N;++i) {
      const auto [x,y] = get_indices(indices[i]);
      pairwise_(x,y) = vars[i];
   }
   return lb;
}

} // namespace LPMP

#endif // LPMP_PAIRWISE_SIMPLEX_FACTOR_H
