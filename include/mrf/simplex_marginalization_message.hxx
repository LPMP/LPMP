#ifndef LPMP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
#define LPMP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX

#include "LP.h"
#include "factors_messages.hxx"
#include "vector.hxx"
#include "memory_allocator.hxx"
#include <cmath>
#include "config.hxx"

namespace LPMP {

// specialized messages between UnarySimplexFactor and PairwiseSimplexFactor
template<Chirality CHIRALITY, bool SUPPORT_INFINITY = true>
class UnaryPairwiseMessage {
public:
   UnaryPairwiseMessage(const INDEX i1, const INDEX i2) {} // obsolete
   UnaryPairwiseMessage() {} 

   static constexpr INDEX pairwise_index_ = CHIRALITY == Chirality::left ? 0 : 1;

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
   {
     // we assume that only r.right_primal was assigned, r.left_primal not
     //assert(r.primal_[0] == i1_);
     vector<REAL> msgs(r.dim(pairwise_index_));
     if(CHIRALITY == Chirality::left) {
       if(r.primal()[1] < r.dim2() && r.primal()[0] >= r.dim1()) {
         for(INDEX x=0; x<r.dim1(); ++x) {
           msgs[x] = r(x,r.primal()[1]);
         }
         msg -= msgs;
       } 
     } else {
       assert(CHIRALITY == Chirality::right);
       if(r.primal()[0] < r.dim1() && r.primal()[1] >= r.dim2()) {
         for(INDEX x=0; x<r.dim2(); ++x) {
           msgs[x] = r(r.primal()[0],x);
         }
         msg -= msgs;
       } 
     }
   }

   template<typename G>
   void RepamLeft(G& r, const REAL msg, const INDEX msg_dim)
   {
      assert(!std::isnan(msg));
      if(SUPPORT_INFINITY) {
         r[msg_dim] += normalize( msg );
      } else {
         r[msg_dim] += msg;
      }
      assert(!std::isnan(r[msg_dim]));
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, MSG& msg)
   {
      for(std::size_t i=0; i<l.size(); ++i) {
          RepamLeft(l, msg[i], i);
      }
   }

   template<typename A1, typename A2>
   void RepamRight(A1& r, const A2& msgs)
   {
      for(INDEX x=0; x<r.dim(pairwise_index_); ++x) {
         assert(!std::isnan(msgs[x]));
         const REAL val = SUPPORT_INFINITY ? normalize(msgs[x]) : msgs[x];
         if(CHIRALITY == Chirality::left) {
            r.msg1(x) += val;
            assert(!std::isnan(r.msg1(x)));
         } else {
            r.msg2(x) += val;
            assert(!std::isnan(r.msg2(x)));
         }
      }
   }
   template<typename G>
   void RepamRight(G& r, const REAL msg, const INDEX dim)
   {
      assert(!std::isnan(msg));
      const REAL val = SUPPORT_INFINITY ? normalize(msg) : msg;
      if(CHIRALITY == Chirality::left) {
         r.msg1(dim) += val;
         assert(!std::isnan(r.msg1(dim)));
      } else {
         r.msg2(dim) += val;
         assert(!std::isnan(r.msg2(dim)));
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      //assert(l.primal() < l.size());
      if(l.primal() < l.size()) {
         const bool changed = (l.primal() != r.primal()[pairwise_index_]);
         r.primal()[pairwise_index_] = l.primal();
         return changed;
      } else {
         return false;
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      //assert(r.primal()[pairwise_index_] < l.size());
      if(r.primal()[pairwise_index_] < l.size()) {
         const bool changed = (l.primal() != r.primal()[pairwise_index_]);
         l.primal() = r.primal()[pairwise_index_];
         return changed; 
      } else {
         return false;
      }
   }

    template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void construct_constraints(SOLVER& s, const LEFT_FACTOR& l, typename SOLVER::vector left_variables, const RIGHT_FACTOR& r, typename SOLVER::vector left_unary_variables, typename SOLVER::vector right_unary_variables, typename SOLVER::matrix pairwise_variables) const
    {
      if(CHIRALITY == Chirality::left) {
        s.make_equal(left_variables.begin(), left_variables.end(), left_unary_variables.begin(), left_unary_variables.end());
      } else {
        assert(CHIRALITY == Chirality::right);
        s.make_equal(left_variables.begin(), left_variables.end(), right_unary_variables.begin(), right_unary_variables.end());
      }
    }

    template<typename LEFT_FACTOR, typename G2>
    void send_message_to_right(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
      const REAL min = l.LowerBound();
      for(INDEX x=0; x<l.size(); ++x) {
        if(!SUPPORT_INFINITY) {
          assert(!std::isnan(l[x]) && l[x] != std::numeric_limits<REAL>::infinity());
        }
        msg[x] -= omega*(l[x] - min);
      }
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    REAL send_message_to_right_improvement(LEFT_FACTOR& l, RIGHT_FACTOR& r)
    {
       const auto before_right_lb = r.LowerBound(); 

       const REAL min = l.LowerBound();
       for(INDEX x=0; x<l.size(); ++x) {
           if(!SUPPORT_INFINITY) { assert(!std::isnan(l[x]) && l[x] != std::numeric_limits<REAL>::infinity()); }
           RepamRight(r, l[x] - min, x);
       }

       const auto after_right_lb = r.LowerBound();

       for(INDEX x=0; x<l.size(); ++x) {
           if(!SUPPORT_INFINITY) { assert(!std::isnan(l[x]) && l[x] != std::numeric_limits<REAL>::infinity()); }
           RepamRight(r, (l[x] - min), x);
       }

       return after_right_lb - before_right_lb; 
    }

    template<typename RIGHT_FACTOR, typename G2>
    void send_message_to_left(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
#ifndef NDEBUG
       const REAL before_lb = r.LowerBound();
#endif
       vector<REAL> msgs;
       if(CHIRALITY == Chirality::left) {
          msgs = r.min_marginal_1();
       } else {
          msgs = r.min_marginal_2(); 
       }
       if(!SUPPORT_INFINITY) {
         for(INDEX x=0; x<msgs.size(); ++x) {
           assert(!std::isnan(msgs[x]) && msgs[x] != std::numeric_limits<REAL>::infinity());
         }
       }
       const REAL min = msgs.min();
       for(INDEX i=0; i<msgs.size(); ++i) { msgs[i] -= min; }
       msg -= omega*msgs;
#ifndef NDEBUG
       const REAL after_lb = r.LowerBound();
       //assert(before_lb <= after_lb); 
#endif
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    REAL send_message_to_left_improvement(LEFT_FACTOR& l, RIGHT_FACTOR& r)
    {
       vector<REAL> msg;
       if(CHIRALITY == Chirality::left) {
          msg = r.min_marginal_1();
       } else {
          msg = r.min_marginal_2(); 
       }
        
       const REAL min = msg.min();
       for(INDEX i=0; i<msg.size(); ++i) { msg[i] -= min; }

       const auto before_left_lb = l.LowerBound();

       RepamLeft(l, msg);

       const auto after_left_lb = l.LowerBound();

       for(auto& x : msg) { x *= -1.0; }
       RepamLeft(l, msg);

       return after_left_lb - before_left_lb;
    }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      return l.primal() == r.primal()[pairwise_index_];
   } 
};

template<INDEX I1, INDEX I2>
class PairwiseTripletMessage {
public:
   PairwiseTripletMessage() {} 
   PairwiseTripletMessage(const INDEX, const INDEX, const INDEX) {} 
   ~PairwiseTripletMessage() {
      static_assert(I1 < I2 && I2 < 3,""); 
   } 

   // for primal computation as in TRW-S, we need to compute restricted messages as well
   template<typename RIGHT_FACTOR, typename G2>
      void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
      {
         throw std::runtime_error("rounding on pairwise factors is not currently supported");
         assert(false);
      }

   template<typename A1, typename A2>
      void RepamLeft(A1& l, const A2& msgs)
      {
         // do zrobienia: possibly use counter
         for(INDEX x1=0; x1<l.dim1(); ++x1) {
            for(INDEX x2=0; x2<l.dim2(); ++x2) {
               l.cost(x1,x2) += normalize( msgs(x1,x2) );
               assert(!std::isnan(l(x1,x2)));
            }
         }
      }
   template<typename A1, typename A2>
      void RepamRight(A1& r, const A2& msgs)
      {
         // do zrobienia: possibly use counter
         if(I1 == 0 && I2 == 1) {
            for(INDEX x1=0; x1<r.dim1(); ++x1) {
               for(INDEX x2=0; x2<r.dim2(); ++x2) {
                  assert(!std::isnan(msgs(x1,x2)));
                  r.msg12(x1,x2) += normalize( msgs(x1,x2) );
                  assert(!std::isnan(r.msg12(x1,x2)));
               }
            }
         } else
            if(I1 == 0 && I2 == 2) {
               for(INDEX x1=0; x1<r.dim1(); ++x1) {
                  for(INDEX x2=0; x2<r.dim3(); ++x2) {
                     assert(!std::isnan(msgs(x1,x2)));
                     r.msg13(x1,x2) += normalize( msgs(x1,x2) );
                     assert(!std::isnan(r.msg13(x1,x2)));
                  }
               }
            } else
               if(I1 == 1 && I2 == 2) {
                  for(INDEX x1=0; x1<r.dim2(); ++x1) {
                     for(INDEX x2=0; x2<r.dim3(); ++x2) {
                        assert(!std::isnan(msgs(x1,x2)));
                        r.msg23(x1,x2) += normalize( msgs(x1,x2) );
                        assert(!std::isnan(r.msg23(x1,x2)));
                     }
                  }
               } else {
                  assert(false);
               }
      }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
      void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
      {
        if(I1 == 0 && I2 == 1) {

          if(l.primal()[0] < l.dim1()) {
            r.primal()[0] = l.primal()[0];
          }
          if(l.primal()[1] < l.dim2()) {
            r.primal()[1] = l.primal()[1];
          }

        } else if(I1 == 0 && I2 == 2) {

          if(l.primal()[0] < l.dim1()) {
            r.primal()[0] = l.primal()[0];
          }
          if(l.primal()[1] < l.dim2()) {
            r.primal()[2] = l.primal()[1];
          }

        } else if(I1 == 1 && I2 == 2) {

          if(l.primal()[0] < l.dim1()) {
            r.primal()[1] = l.primal()[0];
          }
          if(l.primal()[1] < l.dim2()) {
            r.primal()[2] = l.primal()[1];
          }


        } else {
          assert(false);
        }
      }

   template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
     void construct_constraints(SOLVER& s, 
         const LEFT_FACTOR& l, typename SOLVER::vector l_left_msg_variables, typename SOLVER::vector l_right_msg_variables, typename SOLVER::matrix l_pairwise_variables, 
         const RIGHT_FACTOR& r, typename SOLVER::matrix r_msg_12_variables, typename SOLVER::matrix r_msg_13_variables, typename SOLVER::matrix r_msg_23_variables) const
     {
       if(I1 == 0 && I2 == 1) {
         s.make_equal(l_pairwise_variables.begin(), l_pairwise_variables.end(), r_msg_12_variables.begin(), r_msg_12_variables.end());
       } else if(I1 == 0 && I2 == 2) {
         s.make_equal(l_pairwise_variables.begin(), l_pairwise_variables.end(), r_msg_13_variables.begin(), r_msg_13_variables.end());
       } else if(I1 == 1 && I2 == 2) {
         s.make_equal(l_pairwise_variables.begin(), l_pairwise_variables.end(), r_msg_23_variables.begin(), r_msg_23_variables.end());
       } else {
         assert(false); // not possible
       }
     }

   template<typename LEFT_FACTOR, typename G2>
     void send_message_to_right(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
     {
         const auto dim1 = l.dim1();
         const auto dim2 = l.dim2();

         for(std::size_t x1=0; x1<dim1; ++x1) {
             for(std::size_t x2=0; x2<dim2; ++x2) {
                 assert(!std::isnan(l(x1,x2))); 
             } 
         }

         matrix<REAL> m(dim1, dim2);
         auto min = std::numeric_limits<REAL>::infinity();
         for(std::size_t x1=0; x1<dim1; ++x1) {
             for(std::size_t x2=0; x2<dim2; ++x2) {
                 const auto val = l(x1,x2);
                 m(x1,x2) = val;
                 min = std::min(min, val);
             } 
         }

         for(std::size_t x1=0; x1<dim1; ++x1) {
             for(std::size_t x2=0; x2<dim2; ++x2) {
                 m(x1,x2) -= min;
             }
         }

         msg -= omega*m;
     }

   template<typename RIGHT_FACTOR, typename G2>
     void send_message_to_left(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
     {
         matrix<REAL> msgs(r.dim(I1), r.dim(I2));
         if(I1 == 0 && I2 == 1) {
             r.min_marginal12(msgs);
         } else if(I1 == 0 && I2 == 2) {
             r.min_marginal13(msgs);
         } else if(I1 == 1 && I2 == 2) {
             r.min_marginal23(msgs);
         } else {
             assert(false);
         }

         for(auto it=msgs.begin(); it!=msgs.end(); ++it) {
             assert(!std::isnan(*it));
         }

         const auto min = msgs.min();
         for(std::size_t x1=0; x1<msgs.dim1(); ++x1) {
             for(std::size_t x2=0; x2<msgs.dim2(); ++x2) {
                 msgs(x1,x2) -= min;
             }
         }
         msg -= omega*msgs;
     }
};

} // end namespace LPMP

#endif // LPMP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
