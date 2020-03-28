
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"

#include "multicut/multicut_opengm_input.h"
using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_ODD_BICYCLE_WHEEL_MULTICUT>,StandardTighteningVisitor>> solver(argc,argv);
auto input = LPMP::multicut_opengm_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}