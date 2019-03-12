#pragma once

#include "factors_messages.hxx"
#include "LP.h"
#include "mrf/simplex_factor.hxx"
#include "mrf/simplex_marginalization_message.hxx"
#include "horizon_tracking_factors.hxx"
#include "mrf/mrf_problem_construction.hxx"
#include "horizon_tracking_chain_constructor.hxx"
#include "visitors/standard_visitor.hxx"

namespace LPMP {

struct FMC_HORIZON_TRACKING_MULTIPLE_CHAINS {
   constexpr static const char* name = "Horizon tracking problem on multiple chains";

   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 1>;
   using MultipleChainsContainer = FactorContainer<max_potential_on_multiple_chains, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 2>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, MultipleChainsContainer>;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,true>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,true>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 1 >;
   using PairwiseMultipleChainsMessageContainer = MessageContainer<pairwise_max_potential_on_multiple_chains_message, 1, 2, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 2>;

   using MessageList = meta::list<UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, PairwiseMultipleChainsMessageContainer>;

   using mrf_c = mrf_constructor<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS,0,1,0,1>;
   using constructor = max_multiple_chains_constructor<mrf_c, MultipleChainsContainer, PairwiseMultipleChainsMessageContainer>;
   using problem_constructor = constructor;
};
}
