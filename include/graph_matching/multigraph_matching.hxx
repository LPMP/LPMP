#ifndef LPMP_MULTIGRAPH_MATCHING_H
#define LPMP_MULTIGRAPH_MATCHING_H

#include "graph_matching.h"
#include "multigraph_matching_consistency_constraint.hxx"
#include "multigraph_matching_constructor.hxx"
#include "multigraph_matching_input.h"

namespace LPMP {

struct FMC_MGM { // factor message connection
   constexpr static const char* name = "multigraph matching";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM, 3>;
   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero>;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MGM, 0 >;
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
   using ProblemDecompositionList = meta::list<mgm_constructor>;
};

struct FMC_MGM_LINEAR { // factor message connection
   constexpr static const char* name = "multigraph matching";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM_LINEAR, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM_LINEAR, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM_LINEAR, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM_LINEAR, 3>;
   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero>;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MGM_LINEAR, 0 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_LINEAR, 1>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_LINEAR, 2>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 1, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM_LINEAR, 3>;

   using pq_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_LINEAR, 4>; 
   using qr_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_LINEAR, 5>; 

   using MessageList = meta::list< AssignmentConstraintMessage, pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container, pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container >;

   // TODO: make mrf constructor correct for case of no pairwise potentials
   using mrf = mrf_constructor<FMC_MGM_LINEAR,0,0,0,0>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor, triplet_consistency_factor, triplet_consistency_factor_zero, pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container, pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container >;
   using ProblemDecompositionList = meta::list<mgm_constructor>;
}; 

struct FMC_MGM_T { // factor message connection with tightening triplets for underlying MRFs
   constexpr static const char* name = "multigraph matching with mrf tightening";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM_T, 0, true>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM_T, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM_T, 2>; 
   using triplet_consistency_factor_zero = FactorContainer<multigraph_matching_triplet_consistency_factor_zero, FMC_MGM_T, 3>;
   using empty_mrf_triplet_factor = FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MGM_T, 4>;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor, triplet_consistency_factor_zero, empty_mrf_triplet_factor >;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MGM_T, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 2 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 3>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 4>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 2, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM_T, 5>;

   using pq_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 4>; 
   using qr_vector_triplet_consistency_zero_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 3, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 5>; 

   using pairwise_triplet_12_message = MessageContainer<PairwiseTripletMessage<0,1>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 6>;
   using pairwise_triplet_13_message = MessageContainer<PairwiseTripletMessage<0,2>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 7>;
   using pairwise_triplet_23_message = MessageContainer<PairwiseTripletMessage<1,2>, 1, 4, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM_T, 8>; 

   using MessageList = meta::list<
      AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, 
      pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
      pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container,
      pairwise_triplet_12_message, pairwise_triplet_13_message, pairwise_triplet_23_message
         >;

   using mrf = mrf_constructor<FMC_MGM_T,0,1,1,2>;
   using tightening_mrf = tightening_mrf_constructor<mrf,3,6,7,8>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<tightening_mrf>, AssignmentConstraintMessage>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor,
         triplet_consistency_factor, triplet_consistency_factor_zero,
         pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container,
         pq_vector_triplet_consistency_zero_message_container, qr_vector_triplet_consistency_zero_message_container>;
   using ProblemDecompositionList = meta::list<mgm_constructor>;
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_H
