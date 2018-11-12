#ifndef LPMP_SIMPLEX_FACTOR_HXX
#define LPMP_SIMPLEX_FACTOR_HXX

#include "config.hxx"
#include "vector.hxx"
#include "serialization.hxx"

// Investigate contingency tables (Knuth) for more general tabular structures.

namespace LPMP {

// the polytope {x >= 0 : x_1 + ... + x_n = 1 }
// dual problem:
//      max_{z} z 
//      s.t.    z <= repamPot[i]

class UnarySimplexFactor : public vector<REAL> {
public:
    using vector<REAL>::vector;

    UnarySimplexFactor(std::size_t dim)
    : vector(dim, 0.0)
    {}

    template<typename VECTOR>
    UnarySimplexFactor(const VECTOR& vec)
    : vector(vec.size())
    {
        for(std::size_t i=0; i<vec.size(); ++i) { (*this)[i] = vec[i]; }
    }

   void print_potential()
   {
       for(std::size_t i=0; i<size(); ++i) { std::cout << (*this)[i] << ", "; } 
       std::cout << "\n";
   }

   REAL LowerBound() const { 
       //std::cout << "remove me " << this << ": " << (*this)[0] << "," << (*this)[1] << "\n";
      const REAL lb = this->min();
      assert(std::isfinite(lb));
      return lb;
   }

   REAL EvaluatePrimal() const 
   { 
      if(primal_ >= size()) {
         return std::numeric_limits<REAL>::infinity();
      }
      return (*this)[primal_]; 
   }

   void MaximizePotentialAndComputePrimal() 
   {
      if(primal_ >= size()) {
         primal_ = std::min_element(this->begin(), this->end()) - this->begin();
         assert(primal_ < size());
      }
   }

   // load/store function for the primal value
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }

   auto export_variables() { return std::tie(*static_cast<vector<REAL>*>(this)); }

   void init_primal() { primal_ = std::numeric_limits<INDEX>::max(); }
   INDEX primal() const { return primal_; }
   INDEX& primal() { return primal_; }
   void primal(const INDEX p) { primal_ = p; }

   // on which dimensions of cost does the current solution act? Out of this action the subgradient and the dot product between cost and subgradient can be computed
   // possibly this can be used for primal evaluation, making EvaluatePrimal superfluous.
   template<typename ARRAY>
   void apply(ARRAY& a) const 
   { 
      assert(primal_ < this->size());
      a[primal_];
   }

   INDEX subgradient(double* w) const
   {
      assert(primal_ < size());
      std::fill(w, w+this->size(), 0.0);
      w[primal_] = 1.0;
      return this->size();
   }
   REAL dot_product(double* w) const
   {
      return w[primal_];
   }

   // example for constructing constraints
   // note: this function can also be used to reduce sat constraints: Implement a second solver that implements this in terms of external_solver_interface functions!
   // the same holds for LP-interface
   template<typename EXTERNAL_SOLVER>
   void construct_constraints(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v) const
   {
      s.add_simplex_constraint(v.begin(), v.end());
   }
   template<typename EXTERNAL_SOLVER>
   void convert_primal(EXTERNAL_SOLVER& s, typename EXTERNAL_SOLVER::vector v)
   {
      primal_ = s.first_active(v.begin(), v.end());
   }

private:
   INDEX primal_;
};

// the polytope {x >= 0 : x_1 + ... + x_n <= 1 }

class at_most_one_factor : public vector<REAL> {
public:
   at_most_one_factor(const std::vector<REAL>& cost) : vector<REAL>(cost.begin(), cost.end()) {}
   at_most_one_factor(const INDEX n) : vector<REAL>(n, 0.0) {}

   REAL LowerBound() const { 
      const REAL lb = this->min();
      assert(std::isfinite(lb));
      return std::min(lb,REAL(0.0));
   }

   REAL EvaluatePrimal() const 
   { 
      if(primal_ > size()) {
         return std::numeric_limits<REAL>::infinity();
      } else if(primal_ == size()) {
        return 0.0;
      } else {
        return (*this)[primal_]; 
      }
   }
   void MaximizePotentialAndComputePrimal() 
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

   // load/store function for the primal value
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }

   auto export_variables() { return std::tie(*static_cast<vector<REAL>*>(this)); }

   void init_primal() { primal_ = std::numeric_limits<INDEX>::max(); }
   INDEX primal() const { return primal_; }
   INDEX& primal() { return primal_; }
   void primal(const INDEX p) { primal_ = p; }

   REAL sensitivity() const
   {
      auto minima = two_smallest_elements<REAL>(this->begin(), this->end());
      return minima[1] - minima[0];
   }

#ifdef WITH_SAT
   template<typename SAT_SOLVER>
   void construct_sat_clauses(SAT_SOLVER& s) const
   {
      auto literals = s.add_literal_vector(size() + 1);
      add_simplex_constraint_sat(s, literals.begin(), literals.end());
   }

   template<typename VEC>
   void reduce_sat(VEC& assumptions, const REAL th, sat_literal begin) const
   {
      const REAL lb = LowerBound();
      for(INDEX i=0; i<this->size()+1; ++i) {
         if((*this)[i] > lb + th) { 
            assumptions.push_back(-(begin+i));
         }
      } 
   }

   template<typename SAT_SOLVER>
   void convert_primal(SAT_SOLVER& s, sat_literal first)
   {
      for(sat_literal i=first; i<first+this->size()+1; ++i) {
         if(s.solution(i)) {
            primal_ = i-first;
         }
      }
   }
#endif

private:
   INDEX primal_;
};

// do zrobienia: if pairwise was supplied to us (e.g. external factor, then reflect this in constructor and only allocate space for messages.
// When tightening, we can simply replace pairwise pointer to external factor with an explicit copy. Reallocate left_msg_ and right_msg_ to make memory contiguous? Not sure, depends whether we use block_allocator, which will not acually release the memory
// when factor is copied, then pairwise_ must only be copied if it is actually modified. This depends on whether we execute SMRP or MPLP style message passing. Templatize for this possibility
class PairwiseSimplexFactor : public matrix_expression<REAL, PairwiseSimplexFactor> {
public:
   PairwiseSimplexFactor(const INDEX _dim1, const INDEX _dim2) 
     : pairwise_(_dim1, _dim2),
       left_msg_(_dim1),
       right_msg_(_dim2)
   {
      std::fill(pairwise_.begin(), pairwise_.end(), 0.0);
      std::fill(left_msg_.begin(), left_msg_.end(), 0.0);
      std::fill(right_msg_.begin(), right_msg_.end(), 0.0);
   }

   template<typename MATRIX>
   PairwiseSimplexFactor(const INDEX dim1, const INDEX dim2, const MATRIX& m) 
   : PairwiseSimplexFactor(dim1,dim2)
   {
      for(INDEX x1=0; x1<this->dim1(); ++x1) {
         for(INDEX x2=0; x2<this->dim2(); ++x2) {
            this->cost(x1,x2) = m(x1,x2);
         }
      }
   }

   PairwiseSimplexFactor(const PairwiseSimplexFactor& o) 
     : pairwise_(o.dim1(), o.dim2()),
       left_msg_(o.dim1()),
       right_msg_(o.dim2())
  {
      pairwise_ = o.pairwise_;
      left_msg_ = o.left_msg_;
      right_msg_ = o.right_msg_;
   }
   void operator=(const PairwiseSimplexFactor& o) {
      assert(dim1() == o.dim1() && dim2() == o.dim2());
      pairwise_ = o.pairwise_;
      left_msg_ = o.left_msg_;
      right_msg_ = o.right_msg_;
   }

   REAL operator[](const INDEX x) const {
      const INDEX x1 = x/dim2();
      const INDEX x2 = x%dim2();
      assert(x1 < dim1() && x2 < dim2());
      return pairwise_(x1,x2) + left_msg_[x1] + right_msg_[x2];
   }
   // below is not nice: two different values, only differ by const!
   REAL operator()(const INDEX x1, const INDEX x2) const {
      assert(x1 < dim1() && x2 < dim2());
      return pairwise_(x1,x2) + left_msg_[x1] + right_msg_[x2];
   }
   REAL& cost(const INDEX x1, const INDEX x2) {
      assert(x1 < dim1() && x2 < dim2());
      return pairwise_(x1,x2);
   }
   REAL LowerBound() const {
      REAL lb = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            lb = std::min(lb, (*this)(x1,x2));
         }
      }
      assert(std::isfinite(lb));
      return lb;
   }
   REAL sensitivity() const
   {
      REAL smallest = std::numeric_limits<REAL>::infinity();
      REAL second_smallest = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            const REAL val = (*this)(x1,x2);
            const REAL min = std::min(smallest, val);
            const REAL max = std::max(smallest, val);
            smallest = min;
            second_smallest = std::min(max, second_smallest);
         }
      }
      assert(std::isfinite(second_smallest));
      return second_smallest - smallest;
   } 
   
   INDEX size() const { return dim1()*dim2(); }
   INDEX dim1() const { return left_msg_.size(); }
   INDEX dim2() const { return right_msg_.size(); }
   INDEX dim(const INDEX d) const 
   { 
     assert(d<2);
     if(d == 0) { return dim1(); }
     else { return dim2(); }
   }
     
   REAL& pairwise(const INDEX x1, const INDEX x2) { assert(x1<dim1() && x2<dim2()); return pairwise_(x1,x2); }
   REAL& msg1(const INDEX x1) { assert(x1<dim1()); return left_msg_[x1]; }
   REAL& msg2(const INDEX x2) { assert(x2<dim2()); return right_msg_[x2]; }

   void init_primal() 
   {
      primal_[0] = std::numeric_limits<INDEX>::max();
      primal_[1] = std::numeric_limits<INDEX>::max();
   }

   REAL EvaluatePrimal() const { 
      if(primal_[0] >= dim1() || primal_[1] >= dim2()) {
         return std::numeric_limits<REAL>::infinity();
      }
      assert(primal_[0] < dim1());
      assert(primal_[1] < dim2());
      return (*this)(primal_[0], primal_[1]); 
   }
   void MaximizePotentialAndComputePrimal() 
   {
      if(primal_[0] >= dim1() && primal_[1] >= dim2()) {
         REAL min_val = std::numeric_limits<REAL>::infinity();
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
         REAL min_val = std::numeric_limits<REAL>::infinity();
         for(std::size_t x1=0; x1<dim1(); ++x1) {
            if(min_val >= (*this)(x1,primal_[1])) {
               min_val = (*this)(x1,primal_[1]);
               primal_[0] = x1;
            }
         } 
      } else if(primal_[1] >= dim2() && primal_[0] < dim1()) {
         REAL min_val = std::numeric_limits<REAL>::infinity();
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

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_[0], primal_[1] ); }
   //template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cereal::binary_data( pairwise_, sizeof(REAL)*(size()+dim1()+dim2()) ) ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( left_msg_, right_msg_, pairwise_ ); }

   auto export_variables() { return std::tie(left_msg_, right_msg_, pairwise_); }

   vector<REAL> min_marginal_1() const
   {
      auto min = pairwise_.min1(right_msg_);
      min += left_msg_;
#ifndef NDEBUG
      vector<REAL> msg_test(dim1(), std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<dim1(); ++x1) {
          for(INDEX x2=0; x2<dim2(); ++x2) {
             msg_test[x1] = std::min(msg_test[x1],(*this)(x1,x2));
          }
      } 
      for(INDEX i=0; i<dim1(); ++i) {
        assert(std::abs(msg_test[i] - min[i]) <= eps || msg_test[i] == min[i]);
      }
#endif
      return min;
   }

   vector<REAL> min_marginal_2() const
   {
      auto min = pairwise_.min2(left_msg_);
      min += right_msg_;
#ifndef NDEBUG
      vector<REAL> msg_test(dim2(), std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<dim1(); ++x1) {
        for(INDEX x2=0; x2<dim2(); ++x2) {
          msg_test[x2] = std::min(msg_test[x2],(*this)(x1,x2));
        }
      } 
      for(INDEX i=0; i<dim2(); ++i) {
        assert(std::abs(msg_test[i] - min[i]) <= eps || msg_test[i] == min[i]);
      }
#endif 
      return min; 
   }

   template<typename ARRAY>
   void apply(ARRAY& a) const
   {
      assert(primal_[0] < dim1() && primal_[1] < dim2());
      a[primal_[0]];
      a[dim1() + primal_[1]]; 
      a[dim1() + dim2() + primal_[0]*dim2() + primal_[1]];
   }

   template<typename SOLVER>
   void construct_constraints(SOLVER& s, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables) const
   {
     s.add_simplex_constraint(pairwise_variables.begin(), pairwise_variables.end());
      for(INDEX x1=0; x1<dim1(); ++x1) {
         auto slice = pairwise_variables.slice_left(x1);
         auto c = s.add_at_most_one_constraint(slice.begin(), slice.end());
         s.make_equal(c, left_unary_variables[x1]);
      }

      for(INDEX x2=0; x2<dim2(); ++x2) {
         auto slice = pairwise_variables.slice_right(x2);
         auto c = s.add_at_most_one_constraint(slice.begin(), slice.end());
         s.make_equal(c, right_unary_variables[x2]);
      }
   }

   //template<typename VEC>
   //void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
   //{
   //   sat_literal_vector left_unaries(dim1());
   //   sat_literal_vector right_unaries(dim2());
   //   sat_literal_matrix pairwise(dim1(),dim2());
   //   load_sat_literals(begin, left_unaries, right_unaries, pairwise);
//
//      const REAL lb = LowerBound();
//      for(INDEX x1=0; x1<this->dim1(); ++x1) {
//         for(INDEX x2=0; x2<this->dim2(); ++x2) {
//            if((*this)(x1,x2) > lb + th) { 
//               assumptions.push_back(-pairwise(x1,x2));
//            }
//         }
//      } 
//   }

   template<typename SOLVER>
   void convert_primal(SOLVER& s, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables)
   {
     primal_[0] = s.first_active(left_unary_variables);
     primal_[1] = s.first_active(right_unary_variables);
     auto pairwise_sol = s.first_active(pairwise_variables);
     assert(pairwise_sol[0] == primal_[0]);
     assert(pairwise_sol[1] == primal_[1]);
   }

   const std::array<INDEX,2>& primal() const { return primal_; }
   std::array<INDEX,2>& primal() { return primal_; }

   //INDEX primal_[0], primal_[1]; // not so nice: make getters and setters!
private:
   matrix<REAL> pairwise_;
   vector<REAL> left_msg_;
   vector<REAL> right_msg_;
   std::array<INDEX,2> primal_;
};

// factor assumes that triplet potentials is empty and holds only messages to pairwise factors, i.e. is latently factorizable
class SimpleTighteningTernarySimplexFactor : public tensor3_expression<REAL, SimpleTighteningTernarySimplexFactor> {
public:
   SimpleTighteningTernarySimplexFactor(const INDEX _dim1, const INDEX _dim2, const INDEX _dim3) : 
     msg12(_dim1, _dim2, 0.0),
     msg13(_dim1, _dim3, 0.0),
     msg23(_dim2, _dim3, 0.0)
   {
      std::fill(msg12.begin(), msg12.end(), 0.0);
      std::fill(msg13.begin(), msg13.end(), 0.0);
      std::fill(msg23.begin(), msg23.end(), 0.0);
   }
   /*
   SimpleTighteningTernarySimplexFactor(const SimpleTighteningTernarySimplexFactor& o) : 
     msg12_(o.msg12_),
     msg13_(o.msg13_),
     msg23_(o.msg23_)
  {}
   void operator=(const SimpleTighteningTernarySimplexFactor& o) {
      assert(dim1() == o.dim1() && dim2() == o.dim2() && dim3() == o.dim3());
      for(INDEX i=0; i<dim1()*dim2(); ++i) { msg12_[i] = o.msg12_[i]; }
      for(INDEX i=0; i<dim1()*dim3(); ++i) { msg13_[i] = o.msg13_[i]; }
      for(INDEX i=0; i<dim2()*dim3(); ++i) { msg23_[i] = o.msg23_[i]; }
   }
   */

   REAL LowerBound() const {
      REAL lb = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            for(INDEX x3=0; x3<dim3(); ++x3) {
               lb = std::min(lb, (*this)(x1,x2,x3));
            }
         }
      }
      assert(std::isfinite(lb));
      return lb;
   }

   REAL EvaluatePrimal() const
   {
      if(primal_[0] >= dim1() || primal_[1] >= dim2() || primal_[2] >= dim3()) {
         return std::numeric_limits<REAL>::infinity();
      }
      //assert((*this)(primal_[0], primal_[1], primal_[2]) < std::numeric_limits<REAL>::infinity());
      return (*this)(primal_[0], primal_[1], primal_[2]);
   }

   void MaximizePotentialAndComputePrimal()
   {
      if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x1=0; x1<dim1(); ++x1) {
            for(INDEX x2=0; x2<dim2(); ++x2) {
               for(INDEX x3=0; x3<dim3(); ++x3) {
                  if((*this)(x1,x2,x3) < cost) {
                     cost = (*this)(x1,x2,x3);
                     primal_ = {x1,x2,x3};
                  }
               }
            }
         } 
      } else if(primal_[0] < dim1() && primal_[1] < dim2() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x3=0; x3<dim3(); ++x3) {
            if((*this)(primal_[0],primal_[1],x3) < cost) {
               cost = (*this)(primal_[0],primal_[1],x3);
               primal_[2] = x3;
            }
         }
      } else if(primal_[0] < dim1() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] < dim3()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x2=0; x2<dim2(); ++x2) {
            if((*this)(primal_[0],x2,primal_[2]) < cost) {
               cost = (*this)(primal_[0],x2,primal_[2]);
               primal_[1] = x2;
            }
         }
      } else if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] < dim2() && primal_[2] < dim3()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x1=0; x1<dim1(); ++x1) {
            if((*this)(x1,primal_[1],primal_[2]) < cost) {
               cost = (*this)(x1,primal_[1],primal_[2]);
               primal_[0] = x1;
            }
         }
      } else if(primal_[0] < dim1() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x2=0; x2<dim2(); ++x2) {
            for(INDEX x3=0; x3<dim3(); ++x3) {
               if((*this)(primal_[0],x2,x3) < cost) {
                  cost = (*this)(primal_[0],x2,x3);
                  primal_ = {primal_[0],x2,x3};
               }
            }
         } 
      } else if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] < dim2() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x1=0; x1<dim1(); ++x1) {
            for(INDEX x3=0; x3<dim3(); ++x3) {
               if((*this)(x1,primal_[1],x3) < cost) {
                  cost = (*this)(x1,primal_[1],x3);
                  primal_ = {x1,primal_[1],x3};
               }
            }
         } 
      } else if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] < dim3()) {
         REAL cost = std::numeric_limits<REAL>::infinity();
         for(INDEX x1=0; x1<dim1(); ++x1) {
            for(INDEX x2=0; x2<dim2(); ++x2) {
               if((*this)(x1,x2,primal_[2]) < cost) {
                  cost = (*this)(x1,x2,primal_[2]);
                  primal_ = {x1,x2,primal_[2]};
               }
            }
         } 
      } else {
         assert(false);
      }
   }

   REAL operator()(const INDEX x1, const INDEX x2, const INDEX x3) const {
      return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }

   /*
   REAL operator[](const INDEX x) const {
      const INDEX x1 = x / (dim2()*dim3());
      const INDEX x2 = ( x % (dim2()*dim3()) ) / dim3();
      const INDEX x3 = x % dim3();
      return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }
   */

   template<typename MSG>
   void min_marginal12(MSG& msg) const {
      assert(msg.dim1() == dim1());
      assert(msg.dim2() == dim2());
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            msg(x1,x2) = std::numeric_limits<REAL>::infinity();
         }
      }
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
             for(INDEX x3=0; x3<dim3(); ++x3) {
                 msg(x1,x2) = std::min( msg(x1,x2), (*this)(x1,x2,x3));
             }
         }
      }
   }

   template<typename MSG>
   void min_marginal13(MSG& msg) const {
      assert(msg.dim1() == dim1());
      assert(msg.dim2() == dim3());

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            msg(x1,x3) = std::numeric_limits<REAL>::infinity();
         }
      }
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            for(INDEX x3=0; x3<dim3(); ++x3) {
               msg(x1,x3) = std::min(msg(x1,x3), (*this)(x1,x2,x3));
            }
         }
      }
   }

   template<typename MSG>
   void min_marginal23(MSG& msg) const {
      assert(msg.dim1() == dim2());
      assert(msg.dim2() == dim3());
      for(INDEX x2=0; x2<dim2(); ++x2) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            msg(x2,x3) = std::numeric_limits<REAL>::infinity();
         }
      }

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            for(INDEX x3=0; x3<dim3(); ++x3) {
               msg(x2,x3) = std::min(msg(x2,x3),  (*this)(x1,x2,x3));
            }
         }
      }
   }

   INDEX size() const { return dim1()*dim2()*dim3(); }

   INDEX dim(const INDEX d) const 
   { 
     assert(d<3);   
     switch (d) {
       case 0 : return dim1(); break;
       case 1 : return dim2(); break;
       case 2 : return dim3(); break;
       default: return 0;
     }
   }
   INDEX dim1() const { return msg12.dim1(); }
   INDEX dim2() const { return msg23.dim1(); }
   INDEX dim3() const { return msg13.dim2(); }

   void init_primal() {
      primal_[0] = std::numeric_limits<INDEX>::max();
      primal_[1] = std::numeric_limits<INDEX>::max();
      primal_[2] = std::numeric_limits<INDEX>::max();
   }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) 
   { 
      ar( msg12, msg13, msg23 );
   }

   auto export_variables() { return std::tie( msg12, msg13, msg23 ); }

   template<typename SOLVER>
   void construct_constraints(SOLVER& s, typename SOLVER::matrix msg12_vars, typename SOLVER::matrix msg13_vars, typename SOLVER::matrix msg23_vars) const
   {
      auto joint_vars = s.add_tensor(dim1(),dim2(),dim3());

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            auto slice = joint_vars.slice12(x1,x2);
            auto slice_sum = s.add_at_most_one_constraint(slice.begin(), slice.end());
            s.make_equal(msg12_vars(x1,x2), slice_sum);
         }
      }

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            auto slice = joint_vars.slice13(x1,x3);
            auto slice_sum = s.add_at_most_one_constraint(slice.begin(), slice.end());
            s.make_equal(msg13_vars(x1,x3), slice_sum);
         }
      }

      for(INDEX x2=0; x2<dim2(); ++x2) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            auto slice = joint_vars.slice23(x2,x3);
            auto slice_sum = s.add_at_most_one_constraint(slice.begin(), slice.end());
            s.make_equal(msg23_vars(x2,x3), slice_sum);
         }
      }
   }

//   template<typename VEC>
//   void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
//   {
//      sat_literal_matrix literals_12(dim1(),dim2());
//      sat_literal_matrix literals_13(dim1(),dim3());
//      sat_literal_matrix literals_23(dim2(),dim3());
//      sat_literal_tensor literals_123(dim1(),dim2(),dim3());
//      load_sat_literals(begin, literals_12, literals_13, literals_23, literals_123);
//
//      const REAL lb = LowerBound();
//      for(INDEX x1=0; x1<this->dim1(); ++x1) {
//         for(INDEX x2=0; x2<this->dim2(); ++x2) {
//            for(INDEX x3=0; x3<this->dim3(); ++x3) {
//               if((*this)(x1,x2,x3) > lb + th) { 
//                  assumptions.push_back(-literals_123(x1,x2,x3));
//               }
//            }
//         }
//      } 
//   }

   template<typename SOLVER>
   void convert_primal(SOLVER& s, typename SOLVER::matrix msg12_vars, typename SOLVER::matrix msg13_vars, typename SOLVER::matrix msg23_vars)
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

   const std::array<INDEX,3>& primal() const { return primal_; }
   std::array<INDEX,3>& primal() { return primal_; }

   matrix<REAL> msg12;
   matrix<REAL> msg13;
   matrix<REAL> msg23;

protected:
   std::array<INDEX,3> primal_;
};


/*
class TighteningTernarySimplexFactor {
public:
   TighteningTernarySimplexFactor(const INDEX dim1, const INDEX dim2, const INDEX dim3) : dim1_(dim1), dim2_(dim2), dim3_(dim3) 
   {
      const INDEX size = dim1_*dim2_ + dim1_*dim3_ + dim2_*dim3_;
      
      msg12_ = global_real_block_allocator.allocate(size + std::ceil((1.0*sizeof(INDEX))/(1.0*sizeof(REAL))*size));
      assert(msg12_ != nullptr);
      msg13_ = msg12_ + dim1_*dim2_;
      msg23_ = msg13_ + dim1_*dim3_;
      std::fill(msg12_, msg12_ + size, 0.0);

      msg12_sorted_ = reinterpret_cast<INDEX*>(msg12_ + size);
      msg13_sorted_ = msg12_sorted_ + dim1_*dim2_;
      msg23_sorted_ = msg13_sorted_ + dim1_*dim3_;
      std::iota(msg12_sorted_, msg12_sorted_ + dim1_*dim2_, 0);
      std::iota(msg13_sorted_, msg13_sorted_ + dim1_*dim3_, 0);
      std::iota(msg23_sorted_, msg23_sorted_ + dim1_*dim2_, 0);
   }

   // slightly extended algorithms as used by Caetano et al for faster message passing
   void sort_indices() const { // do zrobienia: use insertion sort as it is faster on nearly sorted data
      std::sort(msg12_sorted_, msg12_sorted_ + dim1_*dim2_, [&](const INDEX a, const INDEX b) { return msg12_[msg12_sorted_[a]] < msg12_[msg12_sorted_[b]]; });
      std::sort(msg13_sorted_, msg13_sorted_ + dim1_*dim3_, [&](const INDEX a, const INDEX b) { return msg13_[msg13_sorted_[a]] < msg13_[msg13_sorted_[b]]; });
      std::sort(msg23_sorted_, msg23_sorted_ + dim1_*dim2_, [&](const INDEX a, const INDEX b) { return msg23_[msg23_sorted_[a]] < msg23_[msg23_sorted_[b]]; });
   }

   REAL FindBest3(const INDEX x1, const INDEX x2) const {
      REAL min = std::numeric_limits<REAL>::infinity();
      INDEX end13 = msg13_sorted_inverse(x1,msg23_sorted(x2,0));
      INDEX end23 = msgmsg23_sorted_inverse(msg13_sorted(x1,0));
      for(INDEX l3=0; l3<dim3_; ++l3) {
         if(l3 >= std::min(end13, end23)) break;
         {
            const INDEX x3 = msg13_sorted(x1,l3);
            min_val = std::min(min_val, msg_13(x1,x3) + msg_23(x2,x3));
            end12 = std::min(end12,msg23_sorted_inverse(x2, x3) );
         }
         {
            const INDEX x3 = msg23_sorted(x2,l3);
            min_val = std::min(min_val, msg_13(x1,x3) + msg_23(x2,x3));
            end23 = std::min(end23,msg13_sorted_inverse(x2, x3) );
         }
      }
      return min_val;
   }

   REAL LowerBound() const {
      REAL val = std::numeric_limits<REAL>::infinity();
      sort_indices();
      for(INDEX l1=0; l1<dim1_*dim2_; ++l1) { // add early stopping by checking whether val < msg12_sorted(l1)
         const INDEX x1 = msg12_sorted_1(l1);
         const INDEX x2 = msg12_sorted_2(l1);
         // find best index x3 such that
         val = std::min(val, msg12(x1,x2), FindBest3(x1,x2));
      }

      return val;
   }

   REAL& msg12(const INDEX x1, const INDEX x2) { return msg12_[x1*dim2_ + x2]; }
   REAL& msg13(const INDEX x1, const INDEX x3) { return msg13_[x1*dim3_ + x3]; }
   REAL& msg23(const INDEX x2, const INDEX x3) { return msg12_[x2*dim3_ + x3]; }
   INDEX msg12_sorted_1(const INDEX x) const { return msg12_sorted_[x]/dim1_; }
   INDEX msg12_sorted_2(const INDEX x) const { return msg12_sorted_[x]%dim1_; }
private:
   const INDEX dim1_, dim2_, dim3_;
   REAL *msg12_, *msg13_, *msg23_;
   mutable INDEX *msg12_sorted_, *msg13_sorted_, *msg23_sorted_;

public:

};
*/

class pairwise_potts_factor : public vector<REAL> {
   template<typename MATRIX>
   static bool is_potts(const MATRIX& c)
   {
      if(c.dim1() != c.dim2()) { return false; }
      if(c.dim1() <= 1) { return false; }
      for(INDEX x1=0; x1<c.dim1(); ++x1) {
         for(INDEX x2=0; x2<c.dim2(); ++x2) {
            if(x1 == x2) {
              if(c(x1,x2) != 0.0) return false;
            } else {
               if(c(x1,x2) != c(0,1)) return false;
            }
         }
      }
      return true;
   }
public:
   template<typename VEC>
   pairwise_potts_factor(const INDEX dim1, const INDEX dim2, const VEC& cost)
   : pairwise_potts_factor(dim1, cost(0,1))
   {
      assert(dim1 == dim2);
      assert(is_potts(cost));
      assert(cost(0,0) == 0.0);
   }

   pairwise_potts_factor(const INDEX dim1, const INDEX dim2)
      : pairwise_potts_factor(dim1, REAL(0.0))
   {
      assert(dim1 == dim2);
   }

   pairwise_potts_factor(const INDEX dim1, const INDEX dim2, const REAL diff_cost)
      : pairwise_potts_factor(dim1, diff_cost)
   {
      assert(dim1 == dim2); 
   } 

   pairwise_potts_factor(const INDEX dim, const REAL diff_cost)
      : 
         vector<REAL>(2*dim,0.0), // holds reparametrisation
         diff_cost_(diff_cost) 
   {}

   REAL operator()(const INDEX x1, const INDEX x2) const
   {
      const REAL msg_val = msg1(x1) + msg2(x2);
      const REAL potts_val = x1 != x2 ? diff_cost() : 0.0;
      return msg_val + potts_val;
   }

   // min cost for same label and different label
   std::array<REAL,2> min_values() const
   {
      const auto smallest2 = two_smallest_elements<REAL>(msg2_begin(), msg2_end());

      REAL min_same_label = std::numeric_limits<REAL>::infinity();
      REAL min_diff_label = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<dim(); ++i) {
         const REAL same_label = (*this)[i] + (*this)[i+dim()];
         min_same_label = std::min(min_same_label, same_label);
         const REAL diff_label = (*this)[i] + diff_cost() + ((*this)[i+dim()] == smallest2[0] ? smallest2[1] : smallest2[0]);
         min_diff_label = std::min(min_diff_label, diff_label);
      }
      return {min_same_label, min_diff_label}; 
   }

   REAL LowerBound() const {
      const auto v = min_values();
      return std::min(v[0], v[1]);
   }

   REAL EvaluatePrimal() const
   {
      if(primal_[0] < dim() && primal_[1] < dim()) {
         const REAL pairwise_cost = primal_[0] == primal_[1] ? 0.0 : diff_cost_;
         const REAL msg_cost = (*this)[primal_[0]] + (*this)[primal_[1] + dim()];
         return pairwise_cost + msg_cost; 
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

   void MaximizePotentialAndComputePrimal()
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

   vector<REAL> min_marginal_1() const
   {
      vector<REAL> m(dim());

      const auto smallest2 = two_smallest_elements<REAL>(msg2_begin(), msg2_end());

      for(INDEX i=0; i<dim(); ++i) {
         const REAL same_label = (*this)[i+dim()];
         const REAL diff_label = diff_cost() + ((*this)[i+dim()] == smallest2[0] ? smallest2[1] : smallest2[0]);
         m[i] = (*this)[i] + std::min(same_label, diff_label); 
      } 

      return m;
   }

   vector<REAL> min_marginal_2() const
   {
      vector<REAL> m(dim());

      const auto smallest2 = two_smallest_elements<REAL>(msg1_begin(), msg1_end());

      for(INDEX i=0; i<dim(); ++i) {
         const REAL same_label = (*this)[i];
         const REAL diff_label = diff_cost() + ((*this)[i] == smallest2[0] ? smallest2[1] : smallest2[0]);
         m[i] = (*this)[i+dim()] + std::min(same_label, diff_label); 
      } 
      return m;
   }

   REAL min_marginal_cut() const
   {
      const auto v = min_values();
      return v[1] - v[0];
      //return v[0] - v[1];
   }

   INDEX dim() const { return this->size()/2; }
   INDEX dim1() const { return dim(); }
   INDEX dim2() const { return dim(); }

   REAL* msg1_begin() const { return this->begin(); }
   REAL* msg1_end() const { return this->begin() + dim(); }
   REAL msg1(const INDEX x1) const { assert(x1 < dim()); return (*this)[x1]; }
   REAL& msg1(const INDEX x1) { assert(x1 < dim()); return (*this)[x1]; }

   REAL* msg2_begin() const { return this->begin() + dim(); }
   REAL* msg2_end() const { return this->end(); }
   REAL msg2(const INDEX x2) const { assert(x2 < dim()); return (*this)[dim() + x2]; }
   REAL& msg2(const INDEX x2) { assert(x2 < dim()); return (*this)[dim() + x2]; }

   REAL diff_cost() const { return diff_cost_; }
   REAL& diff_cost() { return diff_cost_; }

   void init_primal() { primal_[0] = std::numeric_limits<INDEX>::max(); primal_[1] = std::numeric_limits<INDEX>::max(); }
   auto& primal() { return primal_; }
   const auto& primal() const { return primal_; }

   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this), diff_cost_ ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

   auto export_variables() { return std::tie( *static_cast<vector<REAL>*>(this), diff_cost_ ); }

protected:
   REAL diff_cost_;
   std::array<INDEX,2> primal_;

};

} // end namespace LPMP

#endif // LPMP_SIMPLEX_FACTOR_HXX

