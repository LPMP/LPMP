#include "horizon_tracking/horizon_tracking.h"
#include "horizon_tracking/horizon_tracking_primal_rounding.hxx"
#include "visitors/standard_visitor.hxx"
#include "LP.h"
#include "LP_FWMAP.hxx"
using namespace LPMP;

template<typename SOLVER>
void PrintObjectives(int argc, char** argv, SOLVER& solver)
{
    Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS>,StandardVisitor> solver_reference(argc,argv);
    auto input = horizon_tracking_uai_input::parse_file(solver_reference.get_input_file());
    construct_horizon_tracking_problem_on_grid_to_chains(input, solver_reference, solver_reference.template GetProblemConstructor<0>());
    auto constructor = solver.template GetProblemConstructor<0>();
    auto constructor_reference = solver_reference.template GetProblemConstructor<0>();
    auto numberPairwise = constructor_reference.get_number_of_pairwise_factors();
    assert(numberPairwise == constructor.get_number_of_pairwise_factors());
    for (std::size_t p = 0; p < numberPairwise; p++) {
        auto* pairwise_reference = constructor_reference.get_pairwise_factor(p);
        auto* pairwise = constructor.get_pairwise_factor(p);
        pairwise_reference->get_factor()->primal() = pairwise->get_factor()->primal();
        pairwise_reference->propagate_primal_through_messages();
    }
    REAL mrf_costs = 0;
    for (std::size_t u = 0; u < constructor_reference.get_number_of_variables(); u++) {
        mrf_costs += constructor_reference.get_unary_factor(u)->EvaluatePrimal();
    }
    for (std::size_t p = 0; p < constructor_reference.get_number_of_pairwise_factors(); p++) {
        mrf_costs += constructor_reference.get_pairwise_factor(p)->EvaluatePrimal();
    }
    solver_reference.RegisterPrimal();
    std::cout<<"\n\nMRF Cost: "<<mrf_costs<<std::endl;
    REAL b_costs = 0;
    for (const auto& f : constructor_reference.max_multiple_chains_factors()) {
        b_costs += f->EvaluatePrimal();
    }
    std::cout<<"Bottleneck Cost: "<<b_costs<<std::endl<<std::endl;
}

int main(int argc, char** argv) {
Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS>,StandardVisitor> solver(argc,argv);
auto input = horizon_tracking_uai_input::parse_file(solver.get_input_file());
construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.template GetProblemConstructor<0>());
/*
{
   lp_reparametrization repam {lp_reparametrization_mode::Anisotropic, 0.0};
   solver.GetLP().set_reparametrization(repam); 
   for(std::size_t iter=0; iter<100; ++iter) {
      std::cout << "MRF lower bound =  " << solver.GetLP().get_original_lp().LowerBound() << "\n";
      solver.GetLP().get_original_lp().ComputePass(); 
   }
}
*/
solver.Solve();
round_primal_solution(solver,true);
solver.WritePrimal();
std::cout<<"\n\n Primal Cost: "<<solver.primal_cost();
std::cout<<"\n Percentage duality gap: "<<100.0 * (solver.primal_cost() - solver.lower_bound()) / solver.lower_bound() <<"\%\n\n";
PrintObjectives(argc, argv, solver);
}
