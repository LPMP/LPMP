
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_uai_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
MpRoundingSolver<Solver<LP<FMC_SRMP_T>,StandardTighteningVisitor>> solver(argc,argv);
auto input = mrf_uai_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}