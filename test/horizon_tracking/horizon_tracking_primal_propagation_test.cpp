#include "horizon_tracking_test_helper.hxx"
#include "test.h"
#include "data_max_potential_grids_test.hxx"
#include "visitors/standard_visitor.hxx"
#include "horizon_tracking/horizon_tracking.h"
#include "LP_FWMAP.hxx"

using namespace LPMP;

std::vector<std::string> solver_options_primal_prop = {
   {"max potential grid test"},
   {"--maxIter"}, {"10"},
   {"--timeout"}, {"60"},
   {"--roundingReparametrization"}, {"anisotropic"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"-v"}, {"2"}
};

int main(int argc, char** argv)
{
    using solver_type = Solver<LP_tree_FWMAP<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS>, StandardVisitor>;
    solver_type solver(solver_options_primal_prop);
    compute_lower_bound_chains(solver, grid_uai_input_medium, grid_uai_input_medium_lb, false);

    solver.GetLP().ComputePass();
    auto& constructor = solver.GetProblemConstructor();
    std::vector<std::size_t> labeling;
    for(std::size_t i=0; i<constructor.get_number_of_variables(); ++i) {
	    auto* f = constructor.get_unary_factor(i);
	    std::size_t minimal_label = std::min_element(f->get_factor()->begin(), f->get_factor()->end()) - f->get_factor()->begin();
	    labeling.push_back(minimal_label);
    }

	solver.GetLP().write_back_reparametrization();
    for (auto i = 0; i < solver.GetLP().number_of_factors(); i++) {
		solver.GetLP().get_factor(i)->init_primal();
	}

    for(std::size_t i=0; i<constructor.get_number_of_variables(); ++i) {
	    auto* f = constructor.get_unary_factor(i);
	    f->get_factor()->primal() = labeling[i];
	    f->propagate_primal_through_messages();
	    // check consistency of pairwise factors
	    for(auto* m : f->template get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryPairwiseMessageLeftContainer>()) {
		    assert(m->CheckPrimalConsistency());
	    }
	    for(auto* m : f->template get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryPairwiseMessageRightContainer>()) {
		    assert(m->CheckPrimalConsistency());
	    }
	    for(std::size_t p=0; p<constructor.get_number_of_pairwise_factors(); ++p) {
		    const auto idx = constructor.get_pairwise_variables(p);
		    if(idx[0] <= i && idx[1] <= i) {
			    auto* p_f = constructor.get_pairwise_factor(p);
			    assert(p_f->get_factor()->primal()[0] == labeling[idx[0]]);
			    assert(p_f->get_factor()->primal()[1] == labeling[idx[1]]);
			    auto msgs_to_chain = p_f->get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::PairwiseMultipleChainsMessageContainer>();
			    assert(msgs_to_chain.size() == 1);
			    assert(msgs_to_chain[0]->CheckPrimalConsistency()); 
		    }
	    }
    }

    assert(solver.GetLP().CheckPrimalConsistency());
}



