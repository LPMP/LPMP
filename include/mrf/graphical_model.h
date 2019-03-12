#ifndef LPMP_GRAPH_MATCHING_H
#define LPMP_GRAPH_MATCHING_H

#include "factors_messages.hxx"
#include "LP.h"
#include "simplex_factor.hxx"
//#include "messages/unary_pairwise_mcf_message.hxx"
#include "simplex_marginalization_message.hxx"
#include "mrf_problem_construction.hxx"

// this file contains solvers for graphical models, among them SRMP and MPLP

namespace LPMP {

struct FMC_SRMP { // equivalent to SRMP or TRWS
   constexpr static const char* name = "SRMP for pairwise case = TRWS";

   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_SRMP, 0, true >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_SRMP, 1, false >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP, 1 >;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer >;

   using mrf = mrf_constructor<FMC_SRMP,0,1,0,1>;
   using problem_constructor = mrf;
};

struct FMC_SRMP_T { // equivalent to SRMP or TRWS
   constexpr static const char* name = "SRMP for pairwise case with tightening triplets";

   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_SRMP_T, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_SRMP_T, 1, false>;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP_T, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP_T, 1 >;
   // tightening
   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_SRMP_T, 2, false>;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP_T, 2>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP_T, 3>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_SRMP_T, 4>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer,
         PairwiseTriplet12MessageContainer, PairwiseTriplet13MessageContainer, PairwiseTriplet23MessageContainer
         >;

   using mrf = mrf_constructor<FMC_SRMP_T,0,1,0,1>;
   using tighteningMrf = tightening_mrf_constructor<mrf,2,2,3,4>;
   using problem_constructor = tighteningMrf;
};



struct FMC_MPLP {
   constexpr static const char* name = "MPLP for pairwise case";

   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MPLP, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MPLP, 1, false>;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::right, variableMessageNumber, 1, FMC_MPLP, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::right, variableMessageNumber, 1, FMC_MPLP, 1 >;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer >;

   using mrf = mrf_constructor<FMC_MPLP,0,1,0,1>;
   using problem_constructor = mrf;
};

} // namespace LPMP

#endif // LPMP_GRAPH_MATCHING_H

