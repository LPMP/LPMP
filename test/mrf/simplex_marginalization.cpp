#include "test.h"
#include "mrf/simplex_factor.hxx"
#include "mrf/simplex_marginalization_message.hxx"
#include "../test_message.hxx"
#include <vector>
#include <random>

using namespace LPMP;

int main()
{
  // "unary/pairwise message", "[unary/pairwise message between simplex factors]" 
   std::vector<double> costUnary {0.1, 0.2, 0.05, 1};
   std::vector<double> costPairwise {  0.1,    0.2,  0.05,
                                  0.3,  0.001,  0.2 ,
                                 -0.3, -0.001, -0.2 ,
                                  0.3,  0.001,  0.2 };
   UnarySimplexFactor simplexUnary(costUnary);
   PairwiseSimplexFactor simplexPairwise(4,3);
   for(INDEX i=0; i<4; ++i) {
      for(INDEX j=0; j<3; ++j) {
         simplexPairwise.cost(i,j) = costPairwise[i*3 + j];
      }
   }
   const auto pairwise_min = *std::min_element(costPairwise.begin(), costPairwise.end());

   UnaryPairwiseMessage<Chirality::left,true> leftMessage(4,3);
   UnaryPairwiseMessage<Chirality::right,true> rightMessage(4,3);

   // must add operators -= and += to vector to support the below things
   // "marginalize pairwise right"
   {
     LPMP::vector<REAL> marg(4,0.0);
     leftMessage.send_message_to_left(simplexPairwise, marg);
     test(marg[0] == -0.05 + pairwise_min);
     test(marg[1] == -0.001 + pairwise_min);
     test(marg[2] ==  0.3 + pairwise_min);
     test(marg[3] == -0.001 + pairwise_min);
   }

   // "marginalize pairwise left"
   {
     LPMP::vector<REAL> marg(3,0.0);
     rightMessage.send_message_to_left(simplexPairwise, marg);
     test(marg[0] == 0.3 + pairwise_min);
     test(marg[1] == 0.001 + pairwise_min);
     test(marg[2] == 0.2 + pairwise_min);
   }

   // do zrobienia: test ComputeLeftFromRightPrimal and reverse, repamLeft,repamRight, Send{Left|Right}Message and also interchange unary and right/left loop

   // randomized message test
   std::random_device rd;
   { 
      for(std::size_t i=1; i<100; ++i) {
         UnarySimplexFactor s(i);
         PairwiseSimplexFactor p(i,i);
         UnaryPairwiseMessage<Chirality::left> msg_left;

         vector<double> msg_vec(i);
         test_repam_message(s, p, msg_left, msg_vec, rd);
      }
   }

   // pairwise triplet message
   {
      for(INDEX i=1; i<30; ++i) {
         for(INDEX j=1; j<30; ++j) {
            PairwiseSimplexFactor p(i,j);
            matrix<double> msg_vec(i,j);

            // 0,1
            {
               SimpleTighteningTernarySimplexFactor t(i,j, (i+j)/2);
               PairwiseTripletMessage<0,1> m;
               test_message(p, t, m, msg_vec, rd);
            }
            // 0,2
            {
               SimpleTighteningTernarySimplexFactor t(i, (i+j)/2, j);
               PairwiseTripletMessage<0,2> m;
               test_message(p, t, m, msg_vec, rd);
            }
            // 1,2
            {
               SimpleTighteningTernarySimplexFactor t((i+j)/2, i, j);
               PairwiseTripletMessage<1,2> m;
               test_message(p, t, m, msg_vec, rd);
            }
         }
      }
   }
}

/*
TEST_CASE("pairwise/triplet message", "[pairwise/triplet message between simplex factors]") {
   PairwiseSimplexFactor simplexPairwise(4,3);
   for(INDEX i=0; i<3; ++i) {
      for(INDEX j=0; j<3; ++j) {
         simplexPairwise(i,j) = costPairwise[i*3 + j];
      }
   }
   SimpleTighteningTernarySimplexFactor simplexTriplet(3,3,3);

   PairwiseTripletMessage12<MessageSendingType::SRMP> message(3,3,3);

   SECTION( "marginalize triplet" ) {
      matrix marg(3,3,0.0);
      message.ReceiveMessageFromRight(simplexTriplet, marg);
   }

   SECTION( "marginalize pairwise" ) {
      matrix marg(3,3,0.0);
      message.ReceiveMessageFromLeft(simplexPairwise, marg);
   }
}
*/
