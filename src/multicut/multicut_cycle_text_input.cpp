
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"

#include "multicut/multicut_text_input.h"
using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_MULTICUT>,StandardTighteningVisitor>> solver(argc,argv);
auto input = LPMP::multicut_text_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}