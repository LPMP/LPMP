#include "graph_matching/multigraph_matching.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;
int main(int argc, char** argv) {
    MpRoundingSolver<Solver<LP<FMC_MGM>,StandardTighteningVisitor>> solver(argc,argv);
    auto input = multigraph_matching_input::parse_file(solver.get_input_file());
    solver.template GetProblemConstructor<0>().construct(input);
    return solver.Solve();
}
