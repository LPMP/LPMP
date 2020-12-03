#include "visitors/standard_visitor.hxx"
#include "test.h"

namespace LPMP {

// UAI test input. Note: not all unaries are present, hence zero unaries must be added.
std::string chain_uai_input_small = 
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
30.000000 9.000000 10.000000 19.000000 

4
9.000000 17.000000 38.000000 8.000000 

MAX-POTENTIALS
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
5.000000 20.000000 2.000000 20.000000 

4
30.000000 33.000000 3.000000 50.000000 
)";

std::vector<std::string> solver_options_small = {
   {"chain 3x1 test"},
   {"--maxIter"}, {"5"},
   {"--timeout"}, {"10"}, // one minute
   {"-v"}, {"2"},
   {"--inputFile"}, chain_uai_input_small
};
}
