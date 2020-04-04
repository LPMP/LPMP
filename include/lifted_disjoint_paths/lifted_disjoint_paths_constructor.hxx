#pragma once

#include "LP.h"
#include "solver.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"

namespace LPMP {
//namespace lifted_disjoint_paths{

    template<class FACTOR_MESSAGE_CONNECTION, std::size_t MCF_FACTOR, std::size_t SINGLE_NODE_CUT_FACTOR, std::size_t MCF_SINGLE_NODE_CUT_MESSAGE> 
    class lifted_disjoint_paths_constructor
    {
        public:
            using FMC = FACTOR_MESSAGE_CONNECTION;
            template<typename SOLVER>
                lifted_disjoint_paths_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

            //void construct(const lifted_disjoint_paths_instance& i);
            void construct(const lifted_disjoint_paths::LdpInstance& i);

            void ComputePrimal();

            void Tighten(const std::size_t nr_constraints_to_add);

        private:
            LP<FMC>* lp_; 
    };

//}
}
