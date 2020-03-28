#pragma once

#include "factors_messages.hxx"
#include "LP.h"
#include "solver.hxx"
#include "mrf/simplex_factor.hxx"
#include "mrf/simplex_marginalization_message.hxx"
#include "mrf/mrf_problem_construction.hxx"
#include "mrf/cycle_inequalities.hxx"
#include "graph_matching_constructor.hxx"
#include "graph_matching_input.h"

#include "equality_message.hxx"
#include "graph_matching_inter_quadratic_message.hxx"
#include "min_cost_flow_factor_ssp.hxx"
#include "graph_matching_local_problem.hxx"
#include "tree_decomposition.hxx"

#include <vector>
#include <fstream>

// this file contains definitions of various graph matching solvers and grammars.
//
// solvers:
// FMC_MP implements graph matching with the uniqueness constraints implemented via messages.
// FMC_MCF implements graph matching with a global min cost flow factor.
// FMC_GM amounts to TRWS with infinity on diagonals.
// FMC_${MODEL}_T implements tightening version of all three solvers using violated cycle tightening of Sontag et al.
//
// input grammars:
// TorresaniEtAlInput contains the grammar used by the dual decomposition algorithm of Torresani, Kolmogorov and Rother.
// UaiGraphMatchingInput contains the grammar in uai MRF format plus constraints section.

// do zrobienia: remove mcf constructors from FMCs and from include

namespace LPMP {

// graph matching with assignment via message passing
struct FMC_MP {
   constexpr static const char* name = "AMP";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MP, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP, 2 >;

   using FactorList = meta::list<UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer>;

   using mrf = mrf_constructor<FMC_MP,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using problem_constructor = gm_constructor;
};

// + inter quadratic messages
struct FMC_MP_Q {
   constexpr static const char* name = "AMP-Q";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_Q, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_Q, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_Q, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q, 2 >;
   using inter_quadratic_message_container = MessageContainer<graph_matching_inter_quadratic_message, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MP_Q, 3 >;

   using FactorList = meta::list<UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, inter_quadratic_message_container>;

   using mrf = mrf_constructor<FMC_MP_Q,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using gm_constructor_q = graph_matching_inter_quadratic_message_constructor<gm_constructor, inter_quadratic_message_container>;
   using problem_constructor = gm_constructor_q;
};

// graph matching with assignment via message passing + tightening triplets
struct FMC_MP_T {
   constexpr static const char* name = "AMP-T";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_T, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_T, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_T, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_T, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_T, 2 >;

   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MP_T, 2 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_T, 3>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_T, 4>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_T, 5>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list< 
      AssignmentConstraintMessage,
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = mrf_constructor<FMC_MP_T,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,2,3,4,5>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage >;
   using problem_constructor = gm_constructor;
};

// graph matching with assignment via message passing + tightening triplets
struct FMC_MP_Q_T {
   constexpr static const char* name = "AMP-Q-T";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_Q_T, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_Q_T, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_Q_T, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q_T, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q_T, 2 >;
   using inter_quadratic_message_container = MessageContainer<graph_matching_inter_quadratic_message, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MP_Q_T, 3 >;

   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MP_Q_T, 2 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q_T, 4>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q_T, 5>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_Q_T, 6>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list< 
      AssignmentConstraintMessage,
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      inter_quadratic_message_container, 
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = mrf_constructor<FMC_MP_Q_T,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,2,4,5,6>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage >;
   using gm_constructor_q = graph_matching_inter_quadratic_message_constructor<gm_constructor, inter_quadratic_message_container>;
   using problem_constructor = gm_constructor_q;
};

// graph matching with assignment via minimum cost flow solver
// only used for subgradient based method, since messages from mcf factor to unary factors are not implemented fully.

struct FMC_MCF {
   constexpr static const char* name = "AMCF";
      
   constexpr static INDEX McfCoveringFactor = 2; // TODO: covering factor not needed anymore

   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_MCF, 0, false >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MCF, 1, true >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MCF, 2 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2,  message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0,  message_passing_schedule::none, 1, variableMessageNumber, FMC_MCF, 2>;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MCF, 3 >;  // message passign schedule should be none, as messages should not be sent at all, only receive restricted message should be called.

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer,
       AssignmentConstraintMessage
      >;

   using mrf = mrf_constructor<FMC_MCF,1,2,0,1>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using mcf_gm_constructor = graph_matching_mcf_constructor<gm_constructor, MinCostFlowAssignmentFactor, UnaryToAssignmentMessageContainer>;
   using problem_constructor = mcf_gm_constructor;
};

// naive version where the assignment is enforced through inf on pairwise diagonals. One has to insert all possible diagonals then.
// this results in a dense standard pairwise graphical model
struct FMC_GM {
   constexpr static const char* name = "GM";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_GM, 0, true >; // make true, if primal rounding similar to TRW-S is required
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_GM, 1, false >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM, 1 >;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer>;

   using mrf = mrf_constructor<FMC_GM,0,1,0,1>;
   using gm_constructor = graph_matching_mrf_constructor<mrf>;
   using problem_constructor = gm_constructor;
};

// + tightening triplets
struct FMC_GM_T {
   constexpr static const char* name = "GM-T";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_GM_T, 0, true >; // make true, if primal rounding similar to TRW-S is required
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_GM_T, 1, false >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_T, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_T, 1 >;

   // tightening
   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_GM_T, 2 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_T, 2>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_T, 3>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_T, 4>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list<
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = mrf_constructor<FMC_GM_T,0,1,0,1>;
   using tightening_mrf = tightening_mrf_constructor<mrf,2,2,3,4>;
   using gm_constructor = graph_matching_mrf_constructor<tightening_mrf>;
   using problem_constructor = gm_constructor;
};

struct FMC_HUNGARIAN_BP {

   constexpr static const char* name = "HBP";
      
   constexpr static INDEX McfCoveringFactor = 2; // TODO: this is not needed anymore

   // rounding is done by the mcf factor
   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_HUNGARIAN_BP, 0, true >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_HUNGARIAN_BP, 1, false >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_HUNGARIAN_BP, 2 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0, message_passing_schedule::right, 1, variableMessageNumber, FMC_HUNGARIAN_BP, 2>;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer
      >;

   using mrf = mrf_constructor<FMC_HUNGARIAN_BP,1,2,0,1>;
   using gm_constructor = graph_matching_mrf_constructor<mrf>;
   using mcf_gm_constructor = graph_matching_mcf_constructor<gm_constructor, MinCostFlowAssignmentFactor, UnaryToAssignmentMessageContainer>;
   using problem_constructor = mcf_gm_constructor;
};

struct FMC_HUNGARIAN_BP_T {

   constexpr static const char* name = "HBP-T";
      
   constexpr static INDEX McfCoveringFactor = 2; // TODO: not needed anymore

   // rounding is done by the mcf factor
   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_HUNGARIAN_BP_T, 0, true >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_HUNGARIAN_BP_T, 1, false >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_HUNGARIAN_BP_T, 2 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP_T, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP_T, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0, message_passing_schedule::right, 1, variableMessageNumber, FMC_HUNGARIAN_BP_T, 2>;

   // tightening
   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_HUNGARIAN_BP_T, 3 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 2, 3, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP_T, 3>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 2, 3, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP_T, 4>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 2, 3, message_passing_schedule::right, variableMessageNumber, 1, FMC_HUNGARIAN_BP_T, 5>;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer,
       PairwiseTriplet12MessageContainer, 
       PairwiseTriplet13MessageContainer, 
       PairwiseTriplet23MessageContainer 
      >;

   using mrf = mrf_constructor<FMC_HUNGARIAN_BP_T,1,2,0,1>;
   using tightening_mrf = tightening_mrf_constructor<mrf,3,3,4,5>;
   using gm_constructor = graph_matching_mrf_constructor<tightening_mrf>;
   using mcf_gm_constructor = graph_matching_mcf_constructor<gm_constructor, MinCostFlowAssignmentFactor, UnaryToAssignmentMessageContainer>;
   using problem_constructor = mcf_gm_constructor;
};

/*
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> 
struct FMC_LOCAL_SUBPROBLEM {
   using FMC_MCF_PARAM = FMC_LOCAL_SUBPROBLEM<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "DD-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "DD-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "DD-B"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_MCF_PARAM, 0, false >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MCF_PARAM, 1, true >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MCF_PARAM, 2 >;
   using local_subproblem_container = FactorContainer<graph_matching_local_subproblem, FMC_MCF_PARAM, 3>; 

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2,  message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 1 >;
   
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0,  message_passing_schedule::right, 1, variableMessageNumber, FMC_MCF_PARAM, 2>;

   using unary_local_subproblem_message_container = MessageContainer<unary_local_problem_message, 1, 3, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, FMC_MCF_PARAM, 3>; 
   using pairwise_local_subproblem_message_container = MessageContainer<pairwise_local_problem_message, 2, 3, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, FMC_MCF_PARAM, 4>; 

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor, local_subproblem_container>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer,
       unary_local_subproblem_message_container,
       pairwise_local_subproblem_message_container
      >;

   //using assignment = AssignmentViaMinCostFlowConstructor<FMC_MCF_PARAM,0>;
   //using mcf = AssignmentConstructor<MinCostFlowConstructorCS2<FMC_MCF_PARAM,0>>;
   //using mrf = StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>;
   using mrf = mrf_constructor<FMC_MCF_PARAM,1,2,0,1>;

   using local_subproblem_constructor_type = local_subproblem_constructor<local_subproblem_container, unary_local_subproblem_message_container, pairwise_local_subproblem_message_container>;
   using local_subproblem_constructor_left = local_subproblem_constructor_type;
   using local_subproblem_constructor_right = local_subproblem_constructor_type;

   using problem_constructor = meta::list<>; //mrf_left,mrf_right,local_subproblem_constructor_left,local_subproblem_constructor_right>;
};
*/

} // namespace LPMP
