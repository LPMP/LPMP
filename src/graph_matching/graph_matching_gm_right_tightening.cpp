
#include "graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_GM_T<PairwiseConstruction::Right>>,StandardTighteningVisitor>>solver(argc,argv);
auto input = LPMP::TorresaniEtAlInput::parse_file(solver.get_input_file());
solver.template GetProblemConstructor<0>().construct(input);
return solver.Solve();
}