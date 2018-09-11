#include "test.h"
#include "mrf/simplex_factor.hxx"

using namespace LPMP;
using namespace std;

void set_potts(PairwiseSimplexFactor& p, const REAL diff)
{
   for(INDEX x1=0; x1<p.dim1(); ++x1) { 
      for(INDEX x2=0; x2<p.dim2(); ++x2) { 
         if(x1 == x2) {
            p.cost(x1,x2) = 0.0;
         } else {
            p.cost(x1,x2) = diff;
         }
      }
   }
}

void test_factor_equal(pairwise_potts_factor& p1, PairwiseSimplexFactor p2)
{
   test( p1.LowerBound() == p2.LowerBound() ); 

   auto msg = p1.min_marginal_1();
   auto msg2 =  p2.min_marginal_1();
   test(msg == msg2);

   msg = p1.min_marginal_2(); 
   msg2 = p2.min_marginal_2();
   test(msg == msg2);
}

int main()
{
   pairwise_potts_factor potts(3,1.0);
   PairwiseSimplexFactor potts2(3,3);
   set_potts(potts2,1.0);


   // "lower bound, positive coupling"
   test_factor_equal(potts, potts2);

   // "lower bound, negative coupling" ) i
   {
     auto potts_neg = potts;
     potts_neg.diff_cost() = -1.0;
     auto potts2_neg = potts2;
     set_potts(potts2_neg,-1.0);
     test_factor_equal(potts_neg, potts2_neg);
   }

   potts.msg1(0) = -0.1; potts2.msg1(0) = -0.1;
   potts.msg1(1) = 0.5;  potts2.msg1(1) = 0.5;
   potts.msg1(2) = 0.8;  potts2.msg1(2) = 0.8;

   potts.msg2(0) = 1.5;  potts2.msg2(0) = 1.5;
   potts.msg2(1) = 1.0;  potts2.msg2(1) = 1.0;

   // "lower bound, positive coupling, with messages"
   test_factor_equal(potts, potts2);

   // "lower bound, negative coupling, with messages"
   potts.diff_cost() = -1.0;
   set_potts(potts2,-1.0);
   test_factor_equal(potts, potts2);

}
