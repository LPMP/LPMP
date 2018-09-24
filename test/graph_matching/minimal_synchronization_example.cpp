#include "graph_matching/multigraph_matching.hxx"
#include "visitors/standard_visitor.hxx"
#include "solver.hxx"
#include <string>
#include <vector>
#include "test.h"

using namespace LPMP;

const std::string minimal_synchronization_example = 
R"(gm 0 1
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 -10
a 2 1 0 -10
a 3 1 1 -1

gm 0 2
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 -10
a 2 1 0 -10
a 3 1 1 -1

gm 1 2
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 -10
a 2 1 0 -10
a 3 1 1 -1
)";

const std::vector<std::string> options = 
{
"",
"--standardReparametrization", "anisotropic",
"--roundingReparametrization", "anisotropic",
"--tightenIteration", "10",
"--tightenInterval", "5",
"--tightenReparametrization", "uniform:0.5",
"--tightenConstraintsMax", "5",
"--maxIter", "100",
"--tighten",
"-v", "2"
};

int main(int argc, char** argv)
{
    MpRoundingSolver<Solver<LP<FMC_MGM>,StandardTighteningVisitor>> solver(options);
    auto input = multigraph_matching_input::parse_string(minimal_synchronization_example);
    solver.template GetProblemConstructor<0>().construct(input);
    solver.Solve();
    std::cout << solver.GetLP().number_of_messages() << "\n";
    test( std::abs(-42 - solver.lower_bound()) <= 1e-6 ); 
    std::cout << solver.get_primal() << "\n";
}
