#ifndef LPMP_HORIZON_TRACKING_H
#define LPMP_HORIZON_TRACKING_H

#include "factors_messages.hxx"
#include "LP.h"
#include "mrf/simplex_factor.hxx"
#include "mrf/simplex_marginalization_message.hxx"
#include "horizon_tracking_factors.hxx"
#include "mrf/mrf_problem_construction.hxx"
#include "horizon_tracking_chain_constructor.hxx"
#include "horizon_tracking_tree_constructor.hxx"
#include "visitors/standard_visitor.hxx"

namespace LPMP {

struct FMC_HORIZON_TRACKING_MULTIPLE_CHAINS {
   constexpr static const char* name = "Horizon tracking problem on multiple chains";

   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 1>;
   using MultipleChainsContainer = FactorContainer<max_potential_on_multiple_chains, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 2>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, MultipleChainsContainer>;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 1 >;
   using PairwiseMultipleChainsMessageContainer = MessageContainer<pairwise_max_potential_on_multiple_chains_message, 1, 2, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_HORIZON_TRACKING_MULTIPLE_CHAINS, 2>;

   using MessageList = meta::list<UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, PairwiseMultipleChainsMessageContainer>;

   using mrf_c = mrf_constructor<FMC_HORIZON_TRACKING_MULTIPLE_CHAINS,0,1,0,1>;
   using constructor = max_multiple_chains_constructor<mrf_c, MultipleChainsContainer, PairwiseMultipleChainsMessageContainer>;
   using ProblemDecompositionList = meta::list<constructor>;
};

// struct FMC_HORIZON_TRACKING_TREES {
//    constexpr static const char* name = "Horizon tracking problem trees";

//    using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_HORIZON_TRACKING_TREES, 0, true>;
//    using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_HORIZON_TRACKING_TREES, 1>;
//    using max_tree_container = FactorContainer<max_potential_on_tree_iterative, FMC_HORIZON_TRACKING_TREES, 2>;
//    using max_potential_container = FactorContainer<max_potential_on_graph, FMC_HORIZON_TRACKING_TREES, 3>;

//    using FactorList = meta::list< UnaryFactor, PairwiseFactor, max_tree_container, max_potential_container >;

//    using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING_TREES, 0 >;
//    using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_HORIZON_TRACKING_TREES, 1 >;
//    using pairwise_max_message_container = MessageContainer<pairwise_max_factor_tree_message, 1, 2, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, FMC_HORIZON_TRACKING_TREES, 2>;
//    using max_tree_to_max_potential_message_container = MessageContainer<max_factor_tree_graph_message, 2, 3, message_passing_schedule::left, 1, variableMessageNumber, FMC_HORIZON_TRACKING_TREES, 3>;

//    using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, pairwise_max_message_container, max_tree_to_max_potential_message_container >;

//    using mrf_c = mrf_constructor<FMC_HORIZON_TRACKING_TREES,0,1,0,1>;
//    using constructor = max_tree_constructor<mrf_c, max_tree_container, max_potential_container, pairwise_max_message_container, max_tree_to_max_potential_message_container>;
//    using ProblemDecompositionList = meta::list<constructor>;
// };

}

#endif // LPMP_HORIZON_TRACKING_H
