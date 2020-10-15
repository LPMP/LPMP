#include "../test_message.hxx"
#include "horizon_tracking/horizon_tracking_factors.hxx"
#include "three_dimensional_variable_array.hxx"
#include "mrf/simplex_factor.hxx"
#include <array>
#include <random>

using namespace LPMP;

void randomly_initialize_potentials(std::vector<three_dimensional_variable_array<REAL>> pot, std::random_device& rd)
{
   std::normal_distribution<> nd(0.0,1.0);
   for(std::size_t i=0; i<pot.size(); ++i) {
      for(std::size_t j=0; j<pot[i].size(); ++j) {
         pot[i].data()[j] = nd(rd);
      }
   } 
}

int main(int argc, char** argv)
{
   std::random_device rd;
   std::uniform_int_distribution ud(5,15);

   for(std::size_t n=2; n<10; ++n) {
      std::vector<three_dimensional_variable_array<REAL>> linear_potentials(n);
      std::vector<three_dimensional_variable_array<REAL>> max_potentials(n);
      std::vector<std::vector<INDEX>> num_labels(n);

      for(std::size_t i=0; i<n; ++i) {
         const std::size_t chain_length = ud(rd);
         assert(chain_length > 1);
         for(std::size_t j=0; j<chain_length; ++j) {
            num_labels[i].push_back(ud(rd));
         } 
         std::vector<std::array<std::size_t,2>> potentials_size;
         for(std::size_t j=0; j<chain_length-1; ++j) {
            potentials_size.push_back( std::array<std::size_t,2>{num_labels[i][j], num_labels[i][j+1]} );
         }
         linear_potentials[i].resize(potentials_size.begin(), potentials_size.end());
         max_potentials[i].resize(potentials_size.begin(), potentials_size.end());
      } 

      randomly_initialize_potentials(linear_potentials, rd);
      randomly_initialize_potentials(max_potentials, rd); 

      max_potential_on_multiple_chains p(linear_potentials, max_potentials, num_labels, );

      test_factor(p,rd);

      for(std::size_t chain_index=0; chain_index<n; ++chain_index) {
         for(std::size_t pairwise_index = 0; pairwise_index < num_labels[chain_index].size()-1; ++pairwise_index) {
           PairwiseSimplexFactor s(num_labels[chain_index][pairwise_index], num_labels[chain_index][pairwise_index+1]);
           pairwise_max_potential_on_multiple_chains_message m(chain_index, pairwise_index, pairwise_index, pairwise_index + 1);
           const std::size_t msg_size = num_labels[chain_index][pairwise_index] * num_labels[chain_index][pairwise_index+1];
           vector<REAL> msg_vec(msg_size);

           test_repam_message<Chirality::right>(s, p, m, msg_vec, rd);
           test_send_message_to<Chirality::right>(s, p, m, msg_vec, rd);
         }
      }
   }
}
