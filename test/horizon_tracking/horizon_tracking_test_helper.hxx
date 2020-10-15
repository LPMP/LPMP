#ifndef LPMP_HORIZON_TRACKING_TEST_HELPER_HXX
#define LPMP_HORIZON_TRACKING_TEST_HELPER_HXX

#include "test.h"
#include "horizon_tracking/horizon_tracking.h"
#include "solver.hxx"
#include "LP_FWMAP.hxx"
#include "LP_conic_bundle.hxx"
#include <type_traits>
#include "horizon_tracking/horizon_tracking_primal_rounding.hxx"

using namespace LPMP;

template<typename SOLVER>
void compute_lower_bound_chains(SOLVER& solver, std::string uai_file, REAL expected_lb, bool check_lb = true)
{
    auto input = horizon_tracking_uai_input::parse_string(uai_file);
    construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.GetProblemConstructor());
    solver.Solve();
    if (check_lb) test(solver.lower_bound(), expected_lb);
}

template<typename SOLVER>
void validate_objective_from_solution(SOLVER& solver, std::vector<std::string> solverOptions, std::string uai_file)
{
    auto constructor = solver.GetProblemConstructor();
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS>, StandardVisitor>;
    solver_type solver_reference(solverOptions);
    auto input = horizon_tracking_uai_input::parse_string(uai_file);
    construct_horizon_tracking_problem_on_grid_to_chains(input, solver_reference, solver_reference.GetProblemConstructor());
    auto constructor_reference = solver_reference.GetProblemConstructor();
    auto numberPairwise = constructor_reference.get_number_of_pairwise_factors();
    assert(numberPairwise == constructor.get_number_of_pairwise_factors());
    for (std::size_t p = 0; p < numberPairwise; p++) {
        auto* pairwise_reference = constructor_reference.get_pairwise_factor(p);
        auto* pairwise = constructor.get_pairwise_factor(p);
        pairwise_reference->get_factor()->primal() = pairwise->get_factor()->primal();
        pairwise_reference->propagate_primal_through_messages();
    }
    for (std::size_t u = 0; u < constructor_reference.get_number_of_variables(); u++) {
        auto* uf = constructor_reference.get_unary_factor(u);
    }
    auto max_potential_on_chain_factors_original = constructor_reference.max_multiple_chains_factors();
    for (auto* current_chain : max_potential_on_chain_factors_original) {
        current_chain->propagate_primal_through_messages();
    }
    solver_reference.RegisterPrimal();
    test(solver_reference.primal_cost(), solver.primal_cost());
}

bool TestUAIChains(std::vector<std::string> solverOptions, std::string uaiFile, double expectedLb) 
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS>, StandardVisitor>;
    solver_type solver(solverOptions);
    compute_lower_bound_chains(solver, uaiFile, expectedLb, false);
    solver.GetLP().write_back_reparametrization();
    test(solver.GetLP().original_factors_lower_bound(), expectedLb);
    test(solver.primal_cost(), std::numeric_limits<REAL>::infinity());

    round_primal_solution(solver, true, true);
    test(solver.primal_cost(), expectedLb);

    // Check if the computed solution actually has the correct objective value:
    validate_objective_from_solution(solver, solverOptions, uaiFile);
}

#endif // HORIZON_TRACKING_TEST_HELPER
