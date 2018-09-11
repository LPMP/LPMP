#include "test.h"
#include "horizon_tracking/horizon_tracking.h"
#include "solver.hxx"
#include "LP_FWMAP.hxx"
#include <type_traits>

using namespace LPMP;

bool TestUAI(std::vector<std::string> solverOptions, std::string uaiFile, double expectedLb) 
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING>, StandardVisitor>;
    solver_type solver(solverOptions);
    auto input = horizon_tracking_uai_input::parse_string(uaiFile);
    construct_horizon_tracking_problem_on_grid(input, solver, solver.template GetProblemConstructor<0>());
    order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());
    solver.Solve();
    test(std::abs(solver.lower_bound() -  expectedLb) <= eps);
    solver.GetLP().write_back_reparametrization();
    // send messages from botleneck potentials down to mrf
    auto chain_constructor = solver.template GetProblemConstructor<0>();
    auto max_potential_factors = chain_constructor.max_graph_factors();
    for(auto* f : max_potential_factors) {
        f->send_messages_to_left<FMC_HORIZON_TRACKING::max_chain_to_max_potential_message_container>(1.0);
    }

    auto max_potential_on_chain_factors = chain_constructor.max_chain_factors();
    for(auto* f : max_potential_on_chain_factors) {
        f->send_messages_to_left<FMC_HORIZON_TRACKING::pairwise_max_message_container>(1.0);
    }
    test(std::abs(solver.GetLP().original_factors_lower_bound() - expectedLb) <= eps);
    auto numF = solver.GetLP().number_of_factors();
    std::vector<FactorTypeAdapter*> factors;
    for (auto i = 0; i < numF; i++)
    {
        auto currentF = solver.GetLP().get_factor(i);
        
        if (dynamic_cast<FMC_HORIZON_TRACKING::UnaryFactor*>(currentF) || 
            dynamic_cast<FMC_HORIZON_TRACKING::PairwiseFactor*>(currentF))
            factors.push_back(currentF);
    }
    test(solver.primal_cost() == std::numeric_limits<REAL>::infinity());
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(factors.begin(), factors.end(), 1);
    solver.RegisterPrimal();
    solver.GetLP().ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(factors.begin(), factors.end(), 2);
    solver.RegisterPrimal();
    test(solver.primal_cost() < std::numeric_limits<REAL>::max());
    test((solver.primal_cost() - expectedLb) <= eps);

    // Check if the computed solution actually has the correct objective value:
    solver_type solver_reference(solverOptions);
    construct_horizon_tracking_problem_on_grid(input, solver_reference, solver_reference.template GetProblemConstructor<0>());
    order_nodes_by_label_space_cadinality(solver_reference.template GetProblemConstructor<0>());
    auto chain_constructor_reference = solver_reference.template GetProblemConstructor<0>();
    auto numberPairwise = chain_constructor_reference.get_number_of_pairwise_factors();
    assert(numberPairwise == chain_constructor.get_number_of_pairwise_factors());
    for (std::size_t p = 0; p < numberPairwise; p++) {
        auto* pairwise_reference = chain_constructor_reference.get_pairwise_factor(p);
        auto* pairwise = chain_constructor.get_pairwise_factor(p);
        pairwise_reference->get_factor()->primal() = pairwise->get_factor()->primal();
        pairwise_reference->propagate_primal_through_messages();
    }
    for (std::size_t u = 0; u < chain_constructor_reference.get_number_of_variables(); u++) {
        auto* uf = chain_constructor_reference.get_unary_factor(u);
    }
    auto max_potential_on_chain_factors_original = chain_constructor_reference.max_chain_factors();
    for (auto* current_chain : max_potential_on_chain_factors_original) {
        current_chain->propagate_primal_through_messages();
    }
    solver_reference.RegisterPrimal();
    test(std::abs(solver_reference.primal_cost() - solver.primal_cost()) <= eps);
}
