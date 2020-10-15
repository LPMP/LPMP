
#include "multigraph_matching/multigraph_matching.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_MGM_T<true>>,StandardTighteningVisitor>>solver(argc,argv);
auto input = Torresani_et_al_multigraph_matching_input::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}