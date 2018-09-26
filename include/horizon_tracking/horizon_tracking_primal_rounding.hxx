#ifndef LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX
#define LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX

#include "horizon_tracking.h"
#include "solver.hxx"
#include "LP.h"

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
void round_primal_solution(SOLVER& solver, bool send_backward = true)
{
    solver.GetLP().write_back_reparametrization();
    auto multiple_chain_constructor = solver.template GetProblemConstructor<0>();
    if (send_backward) {
        auto prevLb = solver.GetLP().LowerBound();
        auto olb1 = solver.GetLP().original_factors_lower_bound();
        std::cout<<"Lower bound before send message left: "<<prevLb<<std::endl;
        // Send messages from Multiple Chain Linear pairwise potentials to MRF pairwise potentials:
        for(auto* m : multiple_chain_constructor.pairwise_to_multiple_chain_messages()) {
#ifndef NDEBUG
            const REAL before_left_lb = m->GetLeftFactor()->LowerBound();
            const REAL before_right_lb = m->GetRightFactor()->LowerBound();
#endif
            m->send_message_to_left();
#ifndef NDEBUG
            const REAL after_left_lb = m->GetLeftFactor()->LowerBound();
            const REAL after_right_lb = m->GetRightFactor()->LowerBound();
            assert(before_left_lb + before_right_lb <= after_left_lb + after_right_lb + eps);
            std::cout<<"LB Change:"<<-before_left_lb - before_right_lb + after_left_lb + after_right_lb <<std::endl;
#endif
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

    solver.GetLP().set_reparametrization(lp_reparametrization(lp_reparametrization_mode::Anisotropic, 0.0));
    for(std::size_t i=0; i<30; ++i) {
       solver.GetLP().ComputeForwardPassAndPrimal();
       solver.RegisterPrimal();
       solver.GetLP().ComputeBackwardPassAndPrimal();
       solver.RegisterPrimal();
    }
}

#endif //LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_NEW_HXX
