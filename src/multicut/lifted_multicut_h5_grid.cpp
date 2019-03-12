
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LPMP;
int main(int argc, char* argv[])

{
ProblemConstructorRoundingSolver<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutH5Input::ParseLiftedProblem<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>,true>);
return solver.Solve();

}
