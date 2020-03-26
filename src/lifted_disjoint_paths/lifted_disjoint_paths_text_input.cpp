#include "lifted_disjoint_paths/lifted_disjoint_paths_input.h"
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;

int main(int argc, char** argv) {
    Solver<LP<FMC_LIFTED_DISJOINT_PATHS,StandardTighteningVisitor>> solver(argc,argv);
    auto input = lifted_disjoint_paths_input::parse_file(solver.get_input_file());
    solver.GetProblemConstructor().construct(input);
    return solver.Solve();
}
