#include <random>
#include "graph_matching/multigraph_matching_consistency_constraint.hxx"
#include "mrf/simplex_factor.hxx"
#include "../test_message.hxx"

using namespace LPMP;

int main(int argc, char** argv)
{
   std::random_device rd;

   for(std::size_t no_simplex_labels_matched = 1; no_simplex_labels_matched < 50; ++no_simplex_labels_matched) {

      multigraph_matching_triplet_consistency_factor t(no_simplex_labels_matched);

      t.for_each_labeling([&](const std::size_t _x, const std::size_t _y, const std::size_t _z) {
         t.x = _x;
         t.y = _y;
         t.z = _z;
         test(t.primal_feasible());
      });

      test_factor(t,rd);

      // test_send_message_to<Chirality::right> does not work currently, because receive_restricted_message_from_right would need to be called, which test_message does not currently do.
      for(std::size_t no_simplex_labels_not_matched = 1; no_simplex_labels_not_matched < 50; ++no_simplex_labels_not_matched) {
         UnarySimplexFactor s(std::vector<double>(no_simplex_labels_matched + no_simplex_labels_not_matched));

         for(std::size_t i=0; i<s.size(); ++i) {
            simplex_multigraph_matching_triplet_scalar_consistency_message m(i);

            vector<double> msg_vec(1);
            test_repam_message(s,t,m, msg_vec, rd);
            test_send_message_to<Chirality::left>(s,t,m, msg_vec, rd); 
         }

         std::vector<std::size_t> idx(no_simplex_labels_matched);
         std::iota(idx.begin(), idx.end(), 0);
         vector<double> msg_vec(no_simplex_labels_matched);

         {
            simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left> m(idx.begin(), idx.end());
            test_repam_message(s,t,m, msg_vec, rd);
            test_send_message_to<Chirality::left>(s,t,m, msg_vec, rd);
         }

         {
            simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right> m(idx.begin(), idx.end());
            test_repam_message(s,t,m, msg_vec, rd);
            test_send_message_to<Chirality::left>(s,t,m, msg_vec, rd); 
         } 
      } 
   } 
}
