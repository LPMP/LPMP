
#include "horizon_tracking/horizon_tracking.h"
#include "horizon_tracking/horizon_tracking_primal_rounding.hxx"
#include "visitors/standard_visitor.hxx"
#include "LP.h"
#include "tree_decomposition.hxx"
using namespace LPMP;
int main(int argc, char** argv) {
    Solver<LP_subgradient_ascent<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS>,StandardVisitor> solver(argc,argv);
    auto input = horizon_tracking_uai_input::parse_file(solver.get_input_file());
    construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.template GetProblemConstructor<0>(), true);
    solver.Solve();
    round_primal_solution(solver, false, false);
    solver.WritePrimal();
}
