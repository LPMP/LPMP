#pragma once

#include "LP.h"
#include "config.hxx"
#include "solver.hxx"
#include "factors_messages.hxx"
#include "multicut_factors_messages.h"
//#include "multiway_cut_factors_messages.hxx"
#include "cut_base/lifted_factor.h"
#include "cut_base/lifted_factor_messages.h"
#include "multicut_triplet_constructor.hxx"
#include "multicut_quadruplet_constructor.hxx"
#include "multicut_quintuplet_constructor.hxx"
#include "cut_base/cut_base_lifted_constructor.hxx"

#include <iostream>
#include <vector>


namespace LPMP {

// do zrobienia: possibly rename unary to edge factor

template<MessageSendingType MESSAGE_SENDING>
struct FMC_MULTICUT {
   constexpr static const char* name = "Multicut with cycle constraints";
   //constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_MULTICUT, 0, true>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_MULTICUT, 1>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTICUT, 2 >;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container>;
   using MessageList = meta::list<edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;

   using multicut = multicut_triplet_constructor<FMC_MULTICUT,edge_factor_container,triplet_factor_container,edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;
   using problem_constructor = multicut;
};

// It would be nice to be able to derive from FMC_MULTICUT. This is not possible due to deviating FMCs. Possibly parametrize above FMC with template
template<MessageSendingType MESSAGE_SENDING>
struct FMC_ODD_WHEEL_MULTICUT {
   constexpr static const char* name = "Multicut with cycle and odd wheel constraints";

   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_ODD_WHEEL_MULTICUT, 0, true>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_ODD_WHEEL_MULTICUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_ODD_WHEEL_MULTICUT, 2>;
      
   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 6>;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container, odd_3_wheel_factor_container>;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,  
      triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123 
      >;

   using multicut_c = multicut_triplet_constructor<FMC_ODD_WHEEL_MULTICUT,edge_factor_container,triplet_factor_container,edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;
   using multicut_cow = multicut_quadruplet_constructor<multicut_c, odd_3_wheel_factor_container, triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123>;
   using problem_constructor = multicut_cow;
};

struct FMC_ODD_BICYCLE_WHEEL_MULTICUT {
   constexpr static const char* name = "Multicut with cycle, odd wheel and odd bicycle wheel constraints";

   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 0>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 2>;
   using odd_bicycle_3_wheel_factor_container = FactorContainer<multicut_odd_bicycle_3_wheel_factor, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 3>;
      
   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 6>;

   using odd_3_wheel_odd_bicycle_wheel_message_0123 = MessageContainer<multicut_odd_3_wheel_odd_bicycle_message_0123, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 7>;
   using odd_3_wheel_odd_bicycle_wheel_message_0124 = MessageContainer<multicut_odd_3_wheel_odd_bicycle_message_0124, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 8>;
   using odd_3_wheel_odd_bicycle_wheel_message_0134 = MessageContainer<multicut_odd_3_wheel_odd_bicycle_message_0134, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 9>;
   using odd_3_wheel_odd_bicycle_wheel_message_0234 = MessageContainer<multicut_odd_3_wheel_odd_bicycle_message_0234, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 10>;
   using odd_3_wheel_odd_bicycle_wheel_message_1234 = MessageContainer<multicut_odd_3_wheel_odd_bicycle_message_1234, 2, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MULTICUT, 11>;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container, odd_3_wheel_factor_container, odd_bicycle_3_wheel_factor_container>;
   using MessageList = meta::list<
      edge_triplet_message_0_container,
      edge_triplet_message_1_container,
      edge_triplet_message_2_container,  

      triplet_odd_wheel_message_012, 
      triplet_odd_wheel_message_013,
      triplet_odd_wheel_message_023,
      triplet_odd_wheel_message_123,
      
      odd_3_wheel_odd_bicycle_wheel_message_0123,
      odd_3_wheel_odd_bicycle_wheel_message_0124,
      odd_3_wheel_odd_bicycle_wheel_message_0134,
      odd_3_wheel_odd_bicycle_wheel_message_0234,
      odd_3_wheel_odd_bicycle_wheel_message_1234
      >;

   using multicut_c = multicut_triplet_constructor<FMC_ODD_BICYCLE_WHEEL_MULTICUT,edge_factor_container,triplet_factor_container,edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;
   using multicut_cow = multicut_quadruplet_constructor<multicut_c, odd_3_wheel_factor_container, triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123>;
   using multicut_cowobw = multicut_quintuplet_constructor<multicut_cow, odd_bicycle_3_wheel_factor_container, odd_3_wheel_odd_bicycle_wheel_message_0123, odd_3_wheel_odd_bicycle_wheel_message_0124, odd_3_wheel_odd_bicycle_wheel_message_0134, odd_3_wheel_odd_bicycle_wheel_message_0234, odd_3_wheel_odd_bicycle_wheel_message_1234>;
   using problem_constructor = multicut_cowobw;
};

struct FMC_LIFTED_MULTICUT {
   constexpr static const char* name = "Lifted Multicut with cycle constraints";
   constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   // no rounding performed: do it via GAEC and K&L, called from problem constructor
   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_LIFTED_MULTICUT, 0, true>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_LIFTED_MULTICUT, 1>;
   using cut_factor_container = FactorContainer<LiftedMulticutCutFactor, FMC_LIFTED_MULTICUT, 2>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_LIFTED_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_LIFTED_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_LIFTED_MULTICUT, 2 >;

   using CutEdgeLiftedMulticutFactorMessageContainer = MessageContainer<CutEdgeLiftedMulticutFactorMessage, 0, 2, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, FMC_LIFTED_MULTICUT, 3 >;
   using LiftedEdgeLiftedMulticutFactorMessageContainer = MessageContainer<LiftedEdgeLiftedMulticutFactorMessage, 0, 2, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, FMC_LIFTED_MULTICUT, 4 >;

   using FactorList = meta::list<
      edge_factor_container, 
      triplet_factor_container, 
      cut_factor_container
         >;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,
      CutEdgeLiftedMulticutFactorMessageContainer, 
      LiftedEdgeLiftedMulticutFactorMessageContainer
         >;

   using base_multicut_c = multicut_triplet_constructor<FMC_LIFTED_MULTICUT,edge_factor_container,triplet_factor_container,edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;
   using lifted_multicut_c = class lifted_constructor<base_multicut_c, base_multicut_c, cut_factor_container, CutEdgeLiftedMulticutFactorMessageContainer, LiftedEdgeLiftedMulticutFactorMessageContainer>;
   using problem_constructor = lifted_multicut_c;
};

// make own project out of this. This project would depend on mrf and cut project
/*
// also only separate with violated cycles only in multiway cut
struct FMC_MULTIWAY_CUT {
   constexpr static const char* name = "Multiway cut with cycle and odd wheel constraints";

   // multicut
   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_MULTIWAY_CUT, 0>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_MULTIWAY_CUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_MULTIWAY_CUT, 2>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 6>;

   // mrf
   using unary_factor_container = FactorContainer<UnarySimplexFactor, FMC_MULTIWAY_CUT, 4, true>;
   using potts_factor_container = FactorContainer<pairwise_potts_factor, FMC_MULTIWAY_CUT, 5>;

   using unary_pairwise_message_0_container = MessageContainer<UnaryPairwiseMessage<Chirality::left>, 4, 5, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 7>;
   using unary_pairwise_message_1_container = MessageContainer<UnaryPairwiseMessage<Chirality::right>, 4, 5, message_passing_schedule::left, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 8>;

   // join multicut edge and Potts factor
   using multicut_edge_potts_message_container = MessageContainer<multicut_edge_potts_message, 0, 5, message_passing_schedule::full, atMostOneMessage, atMostOneMessage, FMC_MULTIWAY_CUT, 9>; 
   // when we tighten, additional edges may not be connected to any MRF factor. Also, before we tighten we actually

   using FactorList = meta::list< 
      edge_factor_container,
      triplet_factor_container,
      odd_3_wheel_factor_container,

      unary_factor_container,
      potts_factor_container 
         >;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,  
      triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123,

      unary_pairwise_message_0_container, unary_pairwise_message_1_container,
      multicut_edge_potts_message_container 
      >;

   using multicut_c = multicut_triplet_constructor<FMC_MULTIWAY_CUT,0,1, 0,1,2, 3>;
   using multicut_cow = multicut_quadruplet_constructor<multicut_c,2, 3,4,5,6>;
   using mrf = StandardMrfConstructor<FMC_MULTIWAY_CUT, 4, 5, 7, 8>;
   using multiway_cut_c = multiway_cut_constructor<FMC_MULTIWAY_CUT,0,1,9>;
   using problem_constructor = meta::list<multicut_cow, mrf, multiway_cut_c>; 
};
*/

/*
struct FMC_ASYMMETRIC_MULTIWAY_CUT {
   constexpr static const char* name = "Asymmetric multiway cut with cycle and odd wheel constraints";

   // multicut
   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 0>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 2>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 6>;

   // mrf
   using unary_factor_container = FactorContainer<UnarySimplexFactor, FMC_ASYMMETRIC_MULTIWAY_CUT, 4, true>;
   using potts_factor_container = FactorContainer<amwc_pairwise_potts_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 5>;

   using unary_pairwise_message_0_container = MessageContainer<UnaryPairwiseMessage<Chirality::left>, 4, 5, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 7>;
   using unary_pairwise_message_1_container = MessageContainer<UnaryPairwiseMessage<Chirality::right>, 4, 5, message_passing_schedule::left, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 8>;

   // join multicut edge and Potts factor
   using multicut_edge_potts_message_container = MessageContainer<multicut_edge_potts_message, 0, 5, message_passing_schedule::full, atMostOneMessage, atMostOneMessage, FMC_ASYMMETRIC_MULTIWAY_CUT, 9>; 
   // when we tighten, additional edges may not be connected to any MRF factor. Also, before we tighten we actually

   using FactorList = meta::list< 
      edge_factor_container,
      triplet_factor_container,
      odd_3_wheel_factor_container,

      unary_factor_container,
      potts_factor_container 
         >;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,  
      triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123,

      unary_pairwise_message_0_container, unary_pairwise_message_1_container,
      multicut_edge_potts_message_container 
      >;

   using multicut_c = multicut_triplet_constructor<FMC_ASYMMETRIC_MULTIWAY_CUT,0,1, 0,1,2, 3>;
   using multicut_cow = multicut_quadruplet_constructor<multicut_c,2, 3,4,5,6>;
   using mrf = StandardMrfConstructor<FMC_ASYMMETRIC_MULTIWAY_CUT, 4, 5, 7, 8>;
   using multiway_cut_c = multiway_cut_constructor<FMC_ASYMMETRIC_MULTIWAY_CUT,0,1,9>;
   using problem_constructor = meta::list<multicut_cow, mrf, multiway_cut_c>; 
};
*/



} // end namespace LPMP
