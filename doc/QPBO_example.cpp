#include "QPBO_example.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;

int main(int argc, char** argv)
{
    Solver<LP<QPBO_FMC>,StandardVisitor> solver(argc,argv);
    auto input = parse_QPBO_file(solver.get_input_file());
    solver.GetProblemConstructor().construct(input);
    return solver.Solve(); 
}
