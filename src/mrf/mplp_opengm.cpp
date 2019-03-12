
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_opengm_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
MpRoundingSolver<Solver<LP<FMC_MPLP>,StandardVisitor>> solver(argc,argv);
auto input = mrf_opengm_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}