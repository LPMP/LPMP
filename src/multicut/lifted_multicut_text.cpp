
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LPMP;
int main(int argc, char* argv[])

{
ProblemConstructorRoundingSolver<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutTextInput::ParseLiftedProblem<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
