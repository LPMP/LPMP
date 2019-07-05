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

enum class PairwiseConstruction {Left,Right,BothSides}; // Indicates whether pairwise potentials should be built on {left|right|both} side(s) of assignment graph.

// graph matching with assignment via message passing
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP {
   using FMC_MP_PARAM = FMC_MP<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B"
      : "unknown variant"));
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_PARAM, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_PARAM, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 2 >;

   using FactorList = meta::list<UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer>;

   using mrf = mrf_constructor<FMC_MP_PARAM,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using problem_constructor = gm_constructor;
};

// + inter quadratic messages
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP_Q {
   using FMC_MP_PARAM = FMC_MP_Q<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B"
      : "unknown variant"));
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_PARAM, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_PARAM, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 2 >;
   using inter_quadratic_message_container = MessageContainer<graph_matching_inter_quadratic_message, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 3 >;

   using FactorList = meta::list<UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, inter_quadratic_message_container>;

   using mrf = mrf_constructor<FMC_MP_PARAM,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using gm_constructor_q = graph_matching_inter_quadratic_message_constructor<gm_constructor, inter_quadratic_message_container>;
   using problem_constructor = gm_constructor_q;
};

// graph matching with assignment via message passing + tightening triplets
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP_T {
   using FMC_MP_PARAM = FMC_MP_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B-T"
      : "unknown variant"));
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_PARAM, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_PARAM, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 2 >;

   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MP_PARAM, 2 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 3>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 4>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 5>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list< 
      AssignmentConstraintMessage,
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = mrf_constructor<FMC_MP_PARAM,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,2,3,4,5>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage >;
   using problem_constructor = gm_constructor;
};

// graph matching with assignment via message passing + tightening triplets
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP_Q_T {
   using FMC_MP_PARAM = FMC_MP_Q_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B-T"
      : "unknown variant"));
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MP_PARAM, 0, true >; // set to true if labeling by unaries is desired
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MP_PARAM, 1, false >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::none, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 2 >;
   using inter_quadratic_message_container = MessageContainer<graph_matching_inter_quadratic_message, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 3 >;

   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MP_PARAM, 2 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 4>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 5>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MP_PARAM, 6>;

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

   using mrf = mrf_constructor<FMC_MP_PARAM,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,2,4,5,6>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage >;
   using gm_constructor_q = graph_matching_inter_quadratic_message_constructor<gm_constructor, inter_quadratic_message_container>;
   using problem_constructor = gm_constructor_q;
};

// graph matching with assignment via minimum cost flow solver

// first good option: construct pairwise potentials on both sides, only send messages from unary to assignment (no receiving) and adjust all factors.
// another good algorithm can be obtained by just using the left side, only receiving messages from unaries (no sending) and adjusting all factors.
// For both algorithms: In maximal perturbation problem capacities have to be set [-1,-inf] for edges = 1 and [0,inf] for edges = 0

template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MCF {
   using FMC_MCF_PARAM = FMC_MCF<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMCF-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMCF-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMCF-B"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_MCF_PARAM, 0, false >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MCF_PARAM, 1, true >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MCF_PARAM, 2 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2,  message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0,  message_passing_schedule::right, 1, variableMessageNumber, FMC_MCF_PARAM, 2>;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MCF_PARAM, 3 >;  // message passign schedule should be none, as messages should not be sent at all, only receive restricted message should be called.

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer,
       AssignmentConstraintMessage
      >;

   using mrf = mrf_constructor<FMC_MCF_PARAM,1,2,0,1>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using mcf_gm_constructor = graph_matching_mcf_constructor<gm_constructor, MinCostFlowAssignmentFactor, UnaryToAssignmentMessageContainer>;
   using problem_constructor = mcf_gm_constructor;
};

// + tightening
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MCF_T {
   using FMC_MCF_PARAM = FMC_MCF_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMCF-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMCF-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMCF-B-T"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_MCF_PARAM, 0, false >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MCF_PARAM, 1, true >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MCF_PARAM, 2 >;
   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MCF_PARAM, 3 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0, message_passing_schedule::left, 1, variableMessageNumber, FMC_MCF_PARAM, 2>;

   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 3>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 4>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MCF_PARAM, 5>;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 1, 1, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, FMC_MCF_PARAM, 6 >;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer, 
       PairwiseTriplet12MessageContainer, 
       PairwiseTriplet13MessageContainer, 
       PairwiseTriplet23MessageContainer,
       AssignmentConstraintMessage
      >;

   using mrf = mrf_constructor<FMC_MCF_PARAM,1,2,0,1>;
   using tightening_mrf = tightening_mrf_constructor<mrf,3,3,4,5>;
   using gm_constructor = graph_matching_constructor< graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage>;
   using mcf_gm_constructor = graph_matching_mcf_constructor<gm_constructor, MinCostFlowAssignmentFactor, UnaryToAssignmentMessageContainer>;
   using problem_constructor = mcf_gm_constructor;
};


// naive version where the assignment is enforced through inf on pairwise diagonals. One has to insert all possible diagonals then.
// this results in a dense standard pairwise graphical model
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_GM {
   using FMC_GM_PARAM = FMC_GM<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "GM-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "GM-I"
      : "unknown variant");
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_GM_PARAM, 0, true >; // make true, if primal rounding similar to TRW-S is required
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_GM_PARAM, 1, false >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 1 >;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;//, McfLabelingFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer>;//, UnaryMcfLabelingMessage >;

   using mrf = mrf_constructor<FMC_GM_PARAM,0,1,0,1>;
   using gm_constructor = graph_matching_mrf_constructor<mrf>;
   using problem_constructor = gm_constructor;
};

// + tightening triplets
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_GM_T {
   using FMC_GM_PARAM = FMC_GM_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "GM-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "GM-I-T"
      : "unknown variant");
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_GM_PARAM, 0, true >; // make true, if primal rounding similar to TRW-S is required
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_GM_PARAM, 1, false >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 1 >;

   // tightening
   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_GM_PARAM, 2 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 2>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 3>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_GM_PARAM, 4>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list<
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = mrf_constructor<FMC_GM_PARAM,0,1,0,1>;
   using tightening_mrf = tightening_mrf_constructor<mrf,2,2,3,4>;
   using gm_constructor = graph_matching_mrf_constructor<tightening_mrf>;
   using problem_constructor = gm_constructor;
};

template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_HUNGARIAN_BP {
   using FMC_PARAM = FMC_HUNGARIAN_BP<PAIRWISE_CONSTRUCTION>;

   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "HBP-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "HBP-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "HBP-B"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   // rounding is done by the mcf factor
   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_PARAM, 0, true >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_PARAM, 1, false >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_PARAM, 2 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0, message_passing_schedule::right, 1, variableMessageNumber, FMC_PARAM, 2>;

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

template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_HUNGARIAN_BP_T {
   using FMC_PARAM = FMC_HUNGARIAN_BP_T<PAIRWISE_CONSTRUCTION>;

   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "HBP-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "HBP-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "HBP-B-T"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   // rounding is done by the mcf factor
   using MinCostFlowAssignmentFactor = FactorContainer<min_cost_flow_factor, FMC_PARAM, 0, true >;
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_PARAM, 1, false >;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_PARAM, 2 >;

   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 0 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 1, 2, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 1 >;
   using UnaryToAssignmentMessageType = unary_min_cost_flow_message<McfCoveringFactor>;
   using UnaryToAssignmentMessageContainer = MessageContainer<UnaryToAssignmentMessageType, 1, 0, message_passing_schedule::right, 1, variableMessageNumber, FMC_PARAM, 2>;

   // tightening
   using EmptyTripletFactor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_PARAM, 3 >;
   using PairwiseTriplet12MessageContainer = MessageContainer<PairwiseTripletMessage<0,1,true>, 2, 3, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 3>;
   using PairwiseTriplet13MessageContainer = MessageContainer<PairwiseTripletMessage<0,2,true>, 2, 3, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 4>;
   using PairwiseTriplet23MessageContainer = MessageContainer<PairwiseTripletMessage<1,2,true>, 2, 3, message_passing_schedule::right, variableMessageNumber, 1, FMC_PARAM, 5>;

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


// helper function for extracting types from FMC
template<template <PairwiseConstruction> class FMC, PairwiseConstruction PC>
constexpr PairwiseConstruction FmcConstruction(FMC<PC>)
{
   return PC;
}

template<template <PairwiseConstruction> class FMC1, template <PairwiseConstruction> class FMC2, PairwiseConstruction PC>
constexpr bool FmcTypeCheck(FMC2<PC>)
{
   if(std::is_same<FMC1<PC>,FMC2<PC>>::value) return true;
   else return false;
}



// the UAI mrf format followed by custom constraints describing an underlying assignment problem.
/*
namespace UaiGraphMatchingInput {

   // first use uai input of mrf constructor, afterwards continue with special constraints section
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;

   struct matching_init_line : pegtl::seq< opt_whitespace, pegtl::string<'m','a','t','c','h','i','n','g'>, opt_whitespace, pegtl::eol> {};
   struct variable : positive_integer {};
   struct label : pegtl::sor< positive_integer, pegtl::string<'s','l','a','c','k'> > {};
   struct matching_line : pegtl::seq< variable, mand_whitespace, pegtl::until< pegtl::eolf, pegtl::seq< label, opt_whitespace> > > {};


   struct grammar : pegtl::seq<
                    pegtl::until<matching_init_line>, 
                    pegtl::plus<matching_line> > {};

   using matching = std::vector<std::vector<INDEX> >;
   using input = std::tuple<UaiMrfInput::MrfInput, matching>;

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< variable > {
      template<typename INPUT>
      static void apply(const INPUT & in, matching & m)
      { 
         const INDEX var = std::stoul(in.string());
         assert(m.size() == var);
         m.push_back(std::vector<INDEX>(0));
      }
   };

   template<> struct action< label > {
      template<typename INPUT>
      static void apply(const INPUT & in, matching & m)
      { 
         const std::string slack = "slack";
         if(slack == in.string()) {
            m.back().push_back(std::numeric_limits<INDEX>::max());
         } else {
            const INDEX label = std::stoul(in.string());
            assert(label < 100); // do zrobienia: remove this
            m.back().push_back(label);
         }
      }
   };

   input ParseFile(const std::string& filename)
   {
      UaiMrfInput::MrfInput gm_input;

      pegtl::file_parser problem(filename);
      bool ret = problem.parse< UaiMrfInput::grammar, UaiMrfInput::action >( gm_input );
      if(!ret) {
         throw std::runtime_error("could not read mrf input"); 
      }

      matching m;
      ret = problem.parse< grammar, action >( m );
      if(!ret) {
         throw std::runtime_error("could not read matching input"); 
      }
      return std::move(std::make_tuple(std::move(gm_input), std::move(m)));
   }

   // build datastructure indexed by labels and whose elements are the nodes in the gm
   using constraints = std::vector<std::vector<std::tuple<INDEX,INDEX>>>; // variable label pairs which match to specific label
   constraints invert_matching(const matching& m)
   {
      constraints m_inv;
      for(INDEX i=0; i<m.size(); ++i) {
         for(INDEX j=0; j<m[i].size(); ++j) { // first number is the variable
            if(m[i][j] != std::numeric_limits<INDEX>::max() && m[i][j] >= m_inv.size()) {
               m_inv.resize(m[i][j]+1);
            }
            if(m[i][j] != std::numeric_limits<INDEX>::max()) {
               m_inv[m[i][j]].push_back(std::make_tuple(i,j));
            }
         }
      }
      return std::move(m_inv);
   }

   template<typename SOLVER>
   bool ParseProblemGM(const std::string& filename, SOLVER& s)
   {
      using FMC = typename SOLVER::FMC;
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      auto i = ParseFile(filename);
      auto& mrf = s.template GetProblemConstructor<0>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf, std::get<0>(i));

      // add additional empty pairwise potentials with infty on diagonals for each constraint, if not already there.

      auto& m = std::get<1>(i);
      auto constraints = invert_matching(m);
      std::map<std::tuple<INDEX,INDEX>, matrix<REAL>> pairwisePot;
      for(INDEX c=0; c<constraints.size(); ++c) {
         for(INDEX c1=0; c1<constraints[c].size(); ++c1) {
            for(INDEX c2=0; c2<c1; ++c2) {
               INDEX var1 = std::get<0>(constraints[c][c1]);
               INDEX label1 = std::get<1>(constraints[c][c1]);
               INDEX var2 = std::get<0>(constraints[c][c2]);
               INDEX label2 = std::get<1>(constraints[c][c2]);
               assert(var1 != var2);
               if(var1 > var2) {
                  std::swap(var1,var2);
                  std::swap(label1,label2);
               }
               assert(var2 < mrf_input.number_of_variables_);
               assert(label1 < mrf_input.cardinality_[var1]);
               assert(label2 < mrf_input.cardinality_[var2]);

               if(!mrf.HasPairwiseFactor(var1,var2)) {
                  if(pairwisePot.find(std::make_tuple(var1,var2)) == pairwisePot.end()) {
                     pairwisePot.insert( std::make_pair(std::make_pair(var1,var2), matrix<REAL>(mrf_input.cardinality_[var1], mrf_input.cardinality_[var2], 0.0)) );
                  } 
                  auto it = pairwisePot.find(std::make_pair(var1,var2));
                  assert(it != pairwisePot.end());
                  it->second.operator()(label1, label2) = std::numeric_limits<REAL>::infinity();
               } else {
                  const INDEX factorId = mrf.GetPairwiseFactorId(var1,var2);
                  const REAL val = mrf.GetPairwiseValue(factorId,label1,label2);
                  assert(val > 1000000);
               }
            }
         }
      }
      for(auto it = pairwisePot.cbegin(); it!=pairwisePot.cend(); ++it) {
         const INDEX var1 = std::get<0>(it->first);
         const INDEX var2 = std::get<1>(it->first);
         const matrix<REAL> pot = it->second;
         mrf.AddPairwiseFactor(var1,var2,pot);
      }
      return true;
   }

   template<typename SOLVER>
   void construct_mp(SOLVER& s, const input& i)
   {
      using FMC = typename SOLVER::FMC;
      auto& mrf_left = s.template GetProblemConstructor<0>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // now build unaries for mrf_right. There will be as many unaries (=labels) on the right as constraints
      auto& mrf_right = s.template GetProblemConstructor<1>();
      auto& matching = std::get<1>(i);
      auto constraints = invert_matching(matching);
      for(auto& c : constraints) {
         auto* u_r = mrf_right.AddUnaryFactor(std::vector<REAL>(c.size()+1,0.0)); // extra label is for non-assignment of label
         for(INDEX var_label=0; var_label<c.size(); ++var_label) {
            const INDEX var = std::get<0>(c[var_label]);
            const INDEX label = std::get<1>(c[var_label]);
            auto* u_l = mrf_left.GetUnaryFactor(var);
            auto* m = new typename FMC::AssignmentConstraintMessage( typename FMC::EqualityMessageType(label, var_label), u_l, u_r);
            s.GetLP().AddMessage(m);
         }
      }
      
      std::cout << "Constructed gm with " << mrf_left.GetNumberOfVariables() << " unary factors and " << mrf_left.GetNumberOfPairwiseFactors() << " pairwise factors\n";
   }

   template<typename SOLVER>
   void construct_mcf(SOLVER& s, const input& in)
   {
      using FMC = typename SOLVER::FMC;
      auto& mrf_left = s.template GetProblemConstructor<0>();
      const auto& mrf_input = std::get<0>(in);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // build assignment problem
      // We have two types of nodes: matching nodes and slack nodes, the latter taking care whenever a label says slack. Because CS2 orders edges, we must insert slack nodes after the matching nodes, when the need arises. Hence we may have to shift matchign node numbers after constructing slack nodes and inserting them between the matching nodes.
      const auto& matching = std::get<1>(in);
      for(const auto& m : matching) {
         // do zrobienia: this should not be done when assert is not called
         std::vector<INDEX> m_filtered;
         std::copy_if(m.begin(), m.end(), std::back_inserter(m_filtered), [](auto x) { return x != std::numeric_limits<INDEX>::max(); });
         assert(std::is_sorted( m_filtered.begin(), m_filtered.end()));
      }
      const INDEX no_left_nodes = mrf_left.GetNumberOfVariables();
      INDEX no_right_nodes = 0;
      for(auto& m : matching) {
         for(auto c : m) {
            if(c < std::numeric_limits<INDEX>::max()) { // is not a slack node
               no_right_nodes = std::max(no_right_nodes, c+1);
            }
         }
      }
      std::vector<INDEX> slack_node_count(no_right_nodes+1,0); // denotes number of slack nodes that come before normal matching node
      for(INDEX i=0; i<matching.size(); ++i) {
         SIGNED_INDEX last_matching_node = -1;
         for(INDEX j=0; j<matching[i].size(); ++j) {
            if(matching[i][j] == std::numeric_limits<INDEX>::max()) { // is a slack node
               // increase number of slack nodes coming right after right node matching[i][j-1]
               slack_node_count[last_matching_node+1]++;
            } else {
               last_matching_node = matching[i][j];
            }
         }
      }
         
      std::vector<INDEX> matching_node_idx(no_right_nodes+1);
      std::vector<INDEX> cum_slack_node_count(no_right_nodes+1,0);
      std::partial_sum(slack_node_count.begin(), slack_node_count.end(), cum_slack_node_count.begin()); 
      for(INDEX i=0; i<matching_node_idx.size(); ++i) {
         matching_node_idx[i] = i + cum_slack_node_count[i];
      }
      const INDEX total_no_right_nodes = no_right_nodes + std::accumulate(slack_node_count.begin(), slack_node_count.end(), 0);
      std::vector<INDEX> slack_nodes_used(no_right_nodes+1,0);
      std::vector<typename min_cost_flow_factor::Edge> edges;
      for(INDEX i=0; i<matching.size(); ++i) {
         SIGNED_INDEX last_matching_node = -1;
         for(INDEX j=0; j<matching[i].size(); ++j) {
            if(matching[i][j] == std::numeric_limits<INDEX>::max()) { // is a slack node
               const INDEX right_node_no = matching_node_idx[ last_matching_node+1 ] + slack_nodes_used[ last_matching_node+1 ] - slack_node_count[ last_matching_node+1 ];
               ++(slack_nodes_used[ last_matching_node+1 ]);
               edges.push_back({i, no_left_nodes + right_node_no, 0, 1, 0.0});
            } else {
               edges.push_back({i, no_left_nodes + matching_node_idx[ matching[i][j] ], 0, 1, 0.0});
               last_matching_node = matching[i][j];
            }
         }
      }
      for(INDEX i=0; i<total_no_right_nodes; ++i) {
         edges.push_back({no_left_nodes + i, no_left_nodes + total_no_right_nodes, 0, 1, 0.0});
      }

      std::vector<std::vector<INDEX>> edgeId(no_left_nodes);
      for(INDEX i=0; i<edges.size(); ++i) {
              const INDEX left_node = edges[i].start_node;
              if(left_node < no_left_nodes) {
                      edgeId[left_node].push_back(i);
              }
      }

      std::vector<SIGNED_INDEX> demands(no_left_nodes + total_no_right_nodes + 1);
      std::fill(demands.begin(), demands.begin() + no_left_nodes, 1);
      std::fill(demands.begin() + no_left_nodes, demands.end(), 0);
      demands.back() = -no_left_nodes;

      auto* f = new typename FMC::MinCostFlowAssignmentFactor( min_cost_flow_factor(edges, demands, edges.size()) );
      s.GetLP().AddFactor(f);
      auto* mcf = f->get_factor()->GetMinCostFlowSolver();

      for(INDEX i=0; i<no_left_nodes; ++i) {
         auto *u = mrf_left.GetUnaryFactor(i);
         using MessageType = typename FMC::UnaryToAssignmentMessageType;
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(edgeId[i]), u, f, mrf_left.GetNumberOfLabels(i));
         s.GetLP().AddMessage(m);
      }
   }

   template<typename SOLVER>
   bool ParseProblemMP(const std::string& filename, SOLVER& s)
   {
      // do zrobienia: FMC must be left type -> static_assert this
      using FMC = typename SOLVER::FMC;
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      const auto input = ParseFile(filename);
      construct_mp(s, input);
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemMCF(const std::string& filename, SOLVER& s)
   {
      using FMC = typename SOLVER::FMC;
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      const auto input = ParseFile(filename);
      construct_mcf(s, input);
      return true;
   }
}
*/


} // namespace LPMP
