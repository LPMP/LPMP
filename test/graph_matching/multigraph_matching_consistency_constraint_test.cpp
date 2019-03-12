#include <random>
#include "graph_matching/multigraph_matching_triplet_consistency_factor.h"
#include "graph_matching/multigraph_matching_simplex_triplet_consistency_messages.h"
#include "mrf/simplex_factor.hxx"
#include "generate_random_label_set.hxx"
#include "../test_message.hxx"

using namespace LPMP; 

int main(int argc, char** argv)
{
   std::random_device rd;

   const std::size_t no_max_labels = 10;

   for(std::size_t no_x_labels_total = 1; no_x_labels_total<no_max_labels; ++no_x_labels_total) {
      for(std::size_t no_y_labels_total = 1; no_y_labels_total<no_max_labels; ++no_y_labels_total) {
         for(std::size_t no_x_labels_matched = 1; no_x_labels_matched<=no_x_labels_total; ++no_x_labels_matched) {
            for(std::size_t no_y_labels_matched = 1; no_y_labels_matched<=no_y_labels_total; ++no_y_labels_matched) {

               const auto labels_x = generate_random_label_set(no_x_labels_matched, no_x_labels_total);
               const auto labels_y = generate_random_label_set(no_y_labels_matched, no_y_labels_total);
               multigraph_matching_triplet_consistency_factor t(labels_x, labels_y);

               double lb = std::numeric_limits<double>::infinity();
               t.for_each_labeling([&](const std::size_t _x, const std::size_t _y, const std::size_t _z) {
                     t.x = _x;
                     t.y = _y;
                     t.z = _z;
                     test(t.primal_feasible());
                     lb = std::min(lb, t.evaluate(_x,_y,_z));
                     });

               test(lb == t.LowerBound());
               test_factor(t,rd);

               multigraph_matching_triplet_consistency_factor_zero t_zero(labels_x, labels_y);
               t_zero.for_each_labeling([&](const std::size_t _x, const std::size_t _y) {
                     t_zero.x = _x;
                     t_zero.y = _y;
                     test(t.primal_feasible());
                     });

               test_factor(t_zero,rd); 

               // test_send_message_to<Chirality::right> does not work currently, because receive_restricted_message_from_right would need to be called, which test_message does not currently do.

               for(std::size_t i=0; i<no_x_labels_matched; ++i) {
                  UnarySimplexFactor s(std::vector<double>(no_x_labels_matched+1));
                  simplex_multigraph_matching_triplet_scalar_consistency_message m(i);

                  vector<double> msg_vec(1);
                  test_repam_message(s,t,m, msg_vec, rd);
                  test_send_message_to<Chirality::left>(s,t,m, msg_vec, rd); 
               }

               {
                  UnarySimplexFactor s(std::vector<double>(no_x_labels_matched+1));
                  vector<double> msg_vec(no_x_labels_matched);
                  simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::left> m;
                  test_repam_message(s,t,m, msg_vec, rd);
                  test_send_message_to<Chirality::left>(s,t,m, msg_vec, rd);

                  test_repam_message(s,t_zero,m, msg_vec, rd);
                  test_send_message_to<Chirality::left>(s,t_zero,m, msg_vec, rd);
               }


               {
                  UnarySimplexFactor s(std::vector<double>(no_y_labels_matched+1));
                  vector<double> msg_vec(no_y_labels_matched);
                  simplex_multigraph_matching_triplet_vector_consistency_message<Chirality::right> m;
                  test_repam_message(s,t,m, msg_vec, rd);
                  test_send_message_to<Chirality::left>(s,t,m, msg_vec, rd); 

                  test_repam_message(s,t_zero,m, msg_vec, rd);
                  test_send_message_to<Chirality::left>(s,t_zero,m, msg_vec, rd); 
               } 
            } 
         } 
      }
   }
}
