#ifndef LPMP_MULTIGRAPH_MATCHING_H
#define LPMP_MULTIGRAPH_MATCHING_H

#include "graph_matching.h"
#include "multigraph_matching_consistency_constraint.hxx"
#include "multigraph_matching_constructor.hxx"
#include "multigraph_matching_input.h"

namespace LPMP {

struct FMC_MGM { // factor message connection
   constexpr static const char* name = "multigraph matching";
      
   using UnaryFactor = FactorContainer<UnarySimplexFactor, FMC_MGM, 0>;
   using PairwiseFactor = FactorContainer<PairwiseSimplexFactor, FMC_MGM, 1>;
   using triplet_consistency_factor = FactorContainer<multigraph_matching_triplet_consistency_factor, FMC_MGM, 2>; 
   using FactorList = meta::list< UnaryFactor, PairwiseFactor, triplet_consistency_factor>;

   using AssignmentConstraintMessage = MessageContainer<EqualityMessage, 0, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, FMC_MGM, 0 >;
   using UnaryPairwiseMessageLeftContainer = MessageContainer<UnaryPairwiseMessage<Chirality::left,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 1 >;
   using UnaryPairwiseMessageRightContainer = MessageContainer<UnaryPairwiseMessage<Chirality::right,false>, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 2 >;

   using pq_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 3>; 
   using qr_vector_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right>, 0, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_MGM, 4>; 
   using pr_scalar_triplet_consistency_message_container = MessageContainer<simplex_multigraph_matching_triplet_scalar_consistency_message, 0, 2, message_passing_schedule::left, variableMessageNumber, 2, FMC_MGM, 5>;


   using MessageList = meta::list< AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer, 
          pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container>;

   using mrf = mrf_constructor<FMC_MGM,0,1,1,2>;
   using gm_constructor = graph_matching_constructor<graph_matching_mrf_constructor<mrf>, AssignmentConstraintMessage>;
   using mgm_constructor = multigraph_matching_constructor<gm_constructor, triplet_consistency_factor, pq_vector_triplet_consistency_message_container, qr_vector_triplet_consistency_message_container, pr_scalar_triplet_consistency_message_container>;
   using ProblemDecompositionList = meta::list<mgm_constructor>;
};

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_H
