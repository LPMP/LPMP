#ifndef LPMP_MULTIGRAPH_MATCHING_H
#define LPMP_MULTIGRAPH_MATCHING_H

#include "graph_matching.h"
#include "multigraph_matching_triplet_consistency_factor.h"
#include "multigraph_matching_simplex_triplet_consistency_messages.h"
#include "multigraph_matching_constructor.hxx"
#include "multigraph_matching_input.h"

namespace LPMP {

template<bool SEND_ASSIGNMENT_MESSAGES = true>
struct FMC_MGM { // factor message connection
   constexpr static const char* name = "multigraph matching";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM, 3>;
   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero>;

   constexpr static auto assignment_message_schedule = SEND_ASSIGNMENT_MESSAGES ? message_passing_schedule::full : message_passing_schedule::none;
   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, assignment_message_schedule, variableMessageNumber, variableMessageNumber, FMC_MGM, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 2 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 3>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 4>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 2, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM, 5>;

   using pq_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 5>; 
   using qr_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 6>; 

   using MessageList = meta::list< AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, 
          pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
          pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container >;

   using mrf = mrf_constructor<FMC_MGM,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor, 
         triplet_consistency_factor, triplet_consistency_factor_zero,
         pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
         pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container >;
   using problem_constructor = mgm_constructor;
};

template<bool SEND_ASSIGNMENT_MESSAGES = true>
struct FMC_MGM_Q { // factor message connection
   constexpr static const char* name = "multigraph matching";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM_Q, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM_Q, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM_Q, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM_Q, 3>;
   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero>;

   constexpr static auto assignment_message_schedule = SEND_ASSIGNMENT_MESSAGES ? message_passing_schedule::full : message_passing_schedule::none;
   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, assignment_message_schedule, variableMessageNumber, variableMessageNumber, FMC_MGM_Q, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q, 2 >;
   using inter_quadratic_message_container = MessageContainer<graph_matching_inter_quadratic_message, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MGM_Q, 3 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q, 4>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q, 5>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 2, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM_Q, 6>;

   using pq_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q, 7>; 
   using qr_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q, 8>; 

   using MessageList = meta::list< AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, inter_quadratic_message_container, 
          pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
          pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container >;

   using mrf = mrf_constructor<FMC_MGM_Q,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using gm_constructor_q = graph_matching_inter_quadratic_message_constructor<gm_constructor, inter_quadratic_message_container>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor_q, 
         triplet_consistency_factor, triplet_consistency_factor_zero,
         pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
         pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container >;
   using problem_constructor = mgm_constructor;
};

template<bool SEND_ASSIGNMENT_MESSAGES = true>
struct FMC_MGM_T { // factor message connection with tightening triplets for underlying MRFs
   constexpr static const char* name = "multigraph matching with mrf tightening";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM_T, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM_T, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM_T, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM_T, 3>;
   using empty_mrf_triplet_factor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MGM_T, 4>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero, empty_mrf_triplet_factor >;

   constexpr static auto assignment_message_schedule = SEND_ASSIGNMENT_MESSAGES ? message_passing_schedule::full : message_passing_schedule::none;
   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, assignment_message_schedule, variableMessageNumber, variableMessageNumber, FMC_MGM_T, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,true>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,true>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 2 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 3>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 4>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 2, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM_T, 5>;

   using pq_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 6>; 
   using qr_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 7>; 

   using pairwise_triplet_12_message = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 8>;
   using pairwise_triplet_13_message = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 9>;
   using pairwise_triplet_23_message = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 10>; 

   using MessageList = meta::list<
      AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, 
      pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
      pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container,
      pairwise_triplet_12_message, pairwise_triplet_13_message, pairwise_triplet_23_message
         >;

   using mrf = mrf_constructor<FMC_MGM_T,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,4,8,9,10>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor,
         triplet_consistency_factor, triplet_consistency_factor_zero,
         pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
         pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container>;
   using problem_constructor = mgm_constructor;
};

template<bool SEND_ASSIGNMENT_MESSAGES = true>
struct FMC_MGM_Q_T { // factor message connection with tightening triplets for underlying MRFs
   constexpr static const char* name = "multigraph matching with mrf tightening";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM_Q_T, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM_Q_T, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM_Q_T, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM_Q_T, 3>;
   using empty_mrf_triplet_factor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MGM_Q_T, 4>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero, empty_mrf_triplet_factor >;

   constexpr static auto assignment_message_schedule = SEND_ASSIGNMENT_MESSAGES ? message_passing_schedule::full : message_passing_schedule::none;
   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, assignment_message_schedule, variableMessageNumber, variableMessageNumber, FMC_MGM_Q_T, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,true>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,true>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 2 >;
   using inter_quadratic_message_container = MessageContainer<graph_matching_inter_quadratic_message, 1, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MGM_Q_T, 3 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 4>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 5>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 2, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM_Q_T, 6>;

   using pq_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 7>; 
   using qr_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 8>; 

   using pairwise_triplet_12_message = MessageContainer<PairwiseTripletMessage<0,1,true>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 9>;
   using pairwise_triplet_13_message = MessageContainer<PairwiseTripletMessage<0,2,true>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 10>;
   using pairwise_triplet_23_message = MessageContainer<PairwiseTripletMessage<1,2,true>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_Q_T, 11>; 

   using MessageList = meta::list<
      AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, inter_quadratic_message_container,
      pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
      pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container,
      pairwise_triplet_12_message, pairwise_triplet_13_message, pairwise_triplet_23_message
         >;

   using mrf = mrf_constructor<FMC_MGM_Q_T,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,4,9,10,11>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage>;
   using gm_constructor_q = graph_matching_inter_quadratic_message_constructor<gm_constructor, inter_quadratic_message_container>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor_q,
         triplet_consistency_factor, triplet_consistency_factor_zero,
         pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
         pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container>;
   using problem_constructor = mgm_constructor;
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_H
