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
        if (dynamic_cast<FMC_HORIZON_TRACKING_CHAINS::UnaryFactor*>(f) || 
            dynamic_cast<FMC_HORIZON_TRACKING_CHAINS::PairwiseFactor*>(f))
            mrf_factors.push_back(f);
    }    
    return mrf_factors;
}

template<typename SOLVER>
void round_primal_solution(SOLVER& solver, bool send_backward = true)
{
    solver.GetLP().write_back_reparametrization();
    auto chain_constructor = solver.template GetProblemConstructor<0>();
    if (send_backward) {
        auto prevLb = solver.GetLP().LowerBound();
        std::cout<<"Lower bound before send message left: "<<prevLb<<std::endl;
        for(auto* m : chain_constructor.max_chain_to_graph_messages()) {
            m->send_message_to_left(); 
        }
        assert(std::abs(solver.GetLP().LowerBound() - prevLb) <= eps);
        for(auto* m : chain_constructor.pairwise_to_chain_messages()) {
            m->send_message_to_left(); 
            assert(std::abs(solver.GetLP().LowerBound() - prevLb) <= eps);
        }
        auto newLb = solver.GetLP().LowerBound();
        if (prevLb - newLb > eps) {
            std::cout<<"Previous Lower Bound: "<<prevLb<<std::endl;
            std::cout<<"New Lower Bound: "<<newLb<<std::endl;
			throw std::runtime_error("Lower Bound decreased by send backward!");
        }
        std::cout<<"Lower bound after send message left: "<<newLb<<std::endl;
    }

    //TO ADDRESS: Changing this primal pass make primal solution on 5x5 grid test very weak.
    for(std::size_t i=0; i<3; ++i) {
       solver.GetLP().ComputeForwardPassAndPrimal();
       solver.RegisterPrimal();
       solver.GetLP().ComputeBackwardPassAndPrimal();
       solver.RegisterPrimal();
    }

    //auto mrf_factors = get_mrf_factors(solver);
    //solver.GetLP().template ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::forward>(mrf_factors.begin(), mrf_factors.end(), std::numeric_limits<INDEX>::max()-2);
    //solver.RegisterPrimal();
    //solver.GetLP().template ComputePassAndPrimal<std::vector<FactorTypeAdapter*>::iterator, Direction::backward>(mrf_factors.begin(), mrf_factors.end(), std::numeric_limits<INDEX>::max()-1);
    //solver.RegisterPrimal();
}

#endif //LPMP_HORIZON_TRACKING_PRIMAL_ROUNDING_HXX
