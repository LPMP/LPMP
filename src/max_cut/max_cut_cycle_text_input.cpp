
#include "max_cut/max_cut.h"
#include "visitors/standard_visitor.hxx"

#include "max_cut/max_cut_text_input.h"
using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_MAX_CUT>,StandardTighteningVisitor>> solver(argc,argv);
auto input = LPMP::max_cut_text_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}