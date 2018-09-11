
#include "graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;
int main(int argc, char** argv) {
MpRoundingSolver<Solver<LP<FMC_MP_T<PairwiseConstruction::Left>>,StandardTighteningVisitor>>solver(argc,argv);
auto input = LPMP::TorresaniEtAlInput::parse_file(solver.get_input_file());
solver.template GetProblemConstructor<0>().read_input(input);
solver.template GetProblemConstructor<0>().construct();
return solver.Solve();
}