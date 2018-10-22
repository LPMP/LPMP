#ifndef LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX
#define LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX

#include "horizon_tracking.h"
#include "solver.hxx"
#include "LP.h"
#include "LP_FWMAP.hxx"

using namespace LPMP;

template<typename SOLVER>
std::vector<FactorTypeAdapter*> get_mrf_factors(SOLVER& solver)
{
    std::vector<FactorTypeAdapter*> mrf_factors;
    for (auto i = 0; i < solver.GetLP().number_of_factors(); i++) {
        auto f = solver.GetLP().get_factor(i);      
        if (dynamic_cast<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryFactor*>(f) || 
            dynamic_cast<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::PairwiseFactor*>(f))
            mrf_factors.push_back(f);
    }    
    return mrf_factors;
}

template<typename SOLVER>
void round_primal_solution(SOLVER& solver, bool do_rounding_on_mrf = false, bool send_backward = false)
{
    solver.GetLP().write_back_reparametrization();
    auto multiple_chain_constructor = solver.template GetProblemConstructor<0>();
    if (send_backward) {
        auto prevLb = solver.GetLP().LowerBound();
        auto olb1 = solver.GetLP().original_factors_lower_bound();
        std::cout<<"Lower bound before send message left: "<<prevLb<<std::endl;
        // Send messages from Multiple Chain Linear pairwise potentials to MRF pairwise potentials:
        for (INDEX passes = 0; passes < 1; passes++) {
            for (INDEX c = 0; c < multiple_chain_constructor.num_chains(); c++) {
                for(auto* m : multiple_chain_constructor.pairwise_to_multiple_chain_messages()) {
        #ifndef NDEBUG
                    const REAL before_left_lb = m->GetLeftFactor()->LowerBound();
                    const REAL before_right_lb = m->GetRightFactor()->LowerBound();
        #endif
                    if (c != m->GetMessageOp().GetChainIndex()) {
                        continue; // utilize the computation on one chain first.
                    }
                    m->send_message_to_left();
        #ifndef NDEBUG
                    const REAL after_left_lb = m->GetLeftFactor()->LowerBound();
                    const REAL after_right_lb = m->GetRightFactor()->LowerBound();
                    assert(before_left_lb + before_right_lb <= after_left_lb + after_right_lb + eps);
                    std::cout<<"LB Change:"<<-before_left_lb - before_right_lb + after_left_lb + after_right_lb <<std::endl;
        #endif
                }
            }
        }
        auto olb3 = solver.GetLP().original_factors_lower_bound();
        auto newLb = solver.GetLP().LowerBound();
        if (prevLb - newLb > eps) {
            std::cout<<"Previous Lower Bound: "<<prevLb<<std::endl;
            std::cout<<"New Lower Bound: "<<newLb<<std::endl;
			throw std::runtime_error("Lower Bound decreased by send backward!");
        }
        std::cout<<"Lower bound after send message left: "<<newLb<<std::endl;
    }

    auto lb = solver.GetLP().LowerBound();
    if (!do_rounding_on_mrf) {
        for(std::size_t p=0; p<multiple_chain_constructor.get_number_of_pairwise_factors(); ++p) {
            auto [i, j] = multiple_chain_constructor.get_pairwise_variables(p);
            if (j - i > 1)
                continue;

            auto* f = multiple_chain_constructor.get_pairwise_factor(p);
            auto msgs_left = f->template get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryPairwiseMessageLeftContainer>();
            auto msgs_right =  f->template get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryPairwiseMessageRightContainer>();
            for(auto* m : msgs_left) {
                m->send_message_to_right(1.0);
            }
            for(auto* m : msgs_right) {
                m->send_message_to_right(1.0);
            }
        }

        for(auto* m : multiple_chain_constructor.pairwise_to_multiple_chain_messages()) {
            m->send_message_to_right();
        }

        for(std::size_t i=0; i<multiple_chain_constructor.get_number_of_variables(); ++i) {
            auto* f = multiple_chain_constructor.get_unary_factor(i);
            f->init_primal();
        }
        
        for(std::size_t p=0; p<multiple_chain_constructor.get_number_of_pairwise_factors(); ++p) {
            auto* f = multiple_chain_constructor.get_pairwise_factor(p);
            f->init_primal();
        }

        for (const auto& chainsFactor : multiple_chain_constructor.max_multiple_chains_factors()) {
            chainsFactor->init_primal();
            chainsFactor->get_factor()->ComputeAndSetPrimal();
        }
        for(auto* m : multiple_chain_constructor.pairwise_to_multiple_chain_messages()) {
            m->ComputeLeftFromRightPrimal();
        }
        for(std::size_t i=0; i<multiple_chain_constructor.get_number_of_variables(); ++i) {
            auto* f = multiple_chain_constructor.get_unary_factor(i);
            for(auto* m : f->template get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryPairwiseMessageLeftContainer>()) {
                m->ComputeLeftFromRightPrimal();
            }
            for(auto* m : f->template get_messages<typename FMC_HORIZON_TRACKING_MULTIPLE_CHAINS::UnaryPairwiseMessageRightContainer>()) {
                m->ComputeLeftFromRightPrimal();
            }
        }
        solver.RegisterPrimal();
    }
    
    else {
        solver.GetLP().set_reparametrization(lp_reparametrization(lp_reparametrization_mode::Anisotropic, 0.0));
        for(std::size_t i=0; i<1; ++i) {
        solver.GetLP().ComputeForwardPassAndPrimal();
        solver.RegisterPrimal();
        solver.GetLP().ComputeBackwardPassAndPrimal();
        solver.RegisterPrimal();
        }
    }

    std::cout<<"Primal Cost: "<<solver.primal_cost()<<std::endl;
    std::cout<<"Perentage Gap: "<<100*(solver.primal_cost() - lb)/std::abs(lb)<<std::endl;
}

#endif //LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_NEW_HXX
