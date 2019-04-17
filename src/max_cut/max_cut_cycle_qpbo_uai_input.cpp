
#include "max_cut/max_cut.h"
#include "visitors/standard_visitor.hxx"

#include "mrf/mrf_uai_input.h"
using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_MAX_CUT>,StandardTighteningVisitor>> solver(argc,argv);
auto input = LPMP::binary_MRF_uai_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}