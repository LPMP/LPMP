
#include "horizon_tracking/horizon_tracking.h"
#include "horizon_tracking/horizon_tracking_primal_rounding.hxx"
#include "visitors/standard_visitor.hxx"
#include "LP.h"
#include "LP_conic_bundle.hxx"
using namespace LPMP;
int main(int argc, char** argv) {
Solver<LP_conic_bundle<FMC_HORIZON_TRACKING_CHAINS>,StandardVisitor> solver(argc,argv);
auto input = horizon_tracking_uai_input::parse_file(solver.get_input_file());
construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.template GetProblemConstructor<0>());
order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
solver.Solve();
solver.GetLP().write_back_reparametrization();
round_primal_solution(solver);
solver.WritePrimal();
std::cout<<"\n\n Primal Cost: "<<solver.primal_cost();
std::cout<<"\n Percentage duality gap: "<<100.0 * (solver.primal_cost() - solver.lower_bound()) / solver.lower_bound() <<"\%\n\n";
}