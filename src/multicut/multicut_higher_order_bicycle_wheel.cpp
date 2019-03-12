
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LPMP;
int main(int argc, char* argv[])

{
ProblemConstructorRoundingSolver<Solver<FMC_ODD_BICYCLE_WHEEL_MULTICUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutTextInput::parse_higher_order<Solver<FMC_ODD_BICYCLE_WHEEL_MULTICUT,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
