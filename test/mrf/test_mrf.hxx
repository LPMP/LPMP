#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "test.h"

namespace LPMP {

// UAI test input. Note: not all unaries are present, hence zero unaries must be added.
std::string uai_test_input = 
R"(MARKOV
3
2 2 3
3
1 0
2 0 1
2 1 2

2
 0.436 0.564

4
 0.128 0.872
 0.920 0.080

6
 0.210 0.333 0.457
 0.811 0.000 0.189 
)";

std::vector<std::string> solver_options = {
   {"graphical model test"},
   {"--maxIter"}, {"500"},
   {"--timeout"}, {"60"}, // one minute
   {"--lowerBoundComputationInterval"}, {"1"},
   {"--primalComputationInterval"}, {"5"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"--roundingReparametrization"}, {"anisotropic"},
   {"-v"}, {"2"},
   {"--inputFile"}, uai_test_input
};

std::string uai_test_input_2 = 
R"(MARKOV
3
2 2 2 
5
1 0 
1 1 
1 2 
2 0 1 
2 1 2 

2
0.000000 0.000000 

2
0.000000 0.000000 

2
0.000000 0.000000 

4
30.000000 9.000000
10.000000 19.000000 

4
9.000000 17.000000 
38.000000 8.000000 
)";

std::vector<std::string> solver_options_2 = {
   {"graphical model test"},
   {"--maxIter"}, {"500"},
   {"--timeout"}, {"60"}, // one minute
   {"--lowerBoundComputationInterval"}, {"1"},
   {"--primalComputationInterval"}, {"5"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"--roundingReparametrization"}, {"anisotropic"},
   {"-v"}, {"0"},
   {"--inputFile"}, uai_test_input_2
};

matrix<REAL> construct_potts(const std::size_t n1, const std::size_t n2, const REAL diag_val, const REAL off_diag_val)
{
    matrix<REAL> p(n1,n2);
    for(std::size_t i=0; i<n1; ++i) {
        for(std::size_t j=0; j<n2; ++j) {
            if(i == j) {
                p(i,j) = diag_val;
            } else {
                p(i,j) = off_diag_val;
            }
        }
    }

    for(std::size_t i=std::min(n1,n2); i<n1; ++i) {
        for(std::size_t j=0; j<n2; ++j) {
            p(i,j) = 100.0;
        }
    }

    for(std::size_t i=0; i<n1; ++i) {
        for(std::size_t j=std::min(n1,n2); j<n2; ++j) {
            p(i,j) = 100.0;
        }
    }

    return p;
}

}
