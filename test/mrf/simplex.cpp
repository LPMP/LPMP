#include "test.h"
#include "../test_message.hxx"
#include "mrf/simplex_factor.hxx"
#include <vector>
#include <random>

using namespace LPMP;

int main() {

   std::random_device rd;

  { // unary simplex
    std::vector<double> cost {0.1, 0.2, 0.05, 1};
    UnarySimplexFactor simplex(cost);

    test( simplex.LowerBound() == 0.05 );

    simplex.init_primal();
    simplex.MaximizePotentialAndComputePrimal();
    test(simplex.primal() == 2);
    test(simplex.EvaluatePrimal() ==simplex.LowerBound());
  }

  {
     for(std::size_t i=1; i<100; ++i) {
        UnarySimplexFactor s(i);
        test_factor(s, rd);
     }
  }


  { // pairwise simplex 
    PairwiseSimplexFactor simplex(3,3);
    for(INDEX x1=0; x1<simplex.dim1(); ++x1) {
      for(INDEX x2=0; x2<simplex.dim2(); ++x2) {
        if(x1 != x2) {
          simplex.cost(x1,x2) = 0;
        } else {
          simplex.cost(x1,x2) = -REAL(x1)-1.0;
        }
      }
    }

    test( simplex.LowerBound() == -3.0 );

    simplex.init_primal();
    simplex.MaximizePotentialAndComputePrimal();
    test(simplex.primal()[0] == 2);
    test(simplex.primal()[1] == 2);
    test(simplex.EvaluatePrimal() ==simplex.LowerBound());
  }

  {
     for(std::size_t i=1; i<100; ++i) {
        for(std::size_t j=1; j<100; ++j) {
           PairwiseSimplexFactor p(i,j);
           test_factor(p, rd);
        }
     }
  }

  {
     for(std::size_t i=1; i<20; ++i) {
        for(std::size_t j=1; j<20; ++j) {
           for(std::size_t k=1; k<20; ++k) {
              SimpleTighteningTernarySimplexFactor t(i,j,k);
              test_factor(t, rd);
           }
        }
     }
  }
}
