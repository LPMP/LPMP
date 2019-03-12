
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"

#include "multicut/multicut_andres_input.h"
using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_MULTICUT<MessageSendingType::SRMP>>,StandardTighteningVisitor>> solver(argc,argv);
auto input = LPMP::multicut_andres_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}