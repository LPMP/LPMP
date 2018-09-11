#include "test.h"
#include "factors_messages.hxx"
#include "factors_storage.hxx"
#include "messages_storage.hxx"
#include "message_passing_weight_computation.hxx"
#include <iostream>

using namespace LPMP;

struct dummy_factor {
   REAL EvaluatePrimal() const { return 0.0; }
   REAL LowerBound() const { return 0.0; }

   void init_primal() {}

   template<typename ARCHIVE>
   void serialize_dual(ARCHIVE& ar) {}

   template<typename ARCHIVE>
   void serialize_primal(ARCHIVE& ar) {}

   auto export_variables() { return std::tie(); }
};

struct dummy_message {

template<typename FACTOR, typename MSG>
void send_message_to_left(const FACTOR& f, MSG& msg, const REAL omega) {}

template<typename FACTOR, typename MSG>
void send_message_to_right(const FACTOR& f, MSG& msg, const REAL omega) {}

};

struct dummy_mrf_FMC {
   using unary_dummy = FactorContainer<dummy_factor, dummy_mrf_FMC, 0>;
   using pairwise_dummy = FactorContainer<dummy_factor, dummy_mrf_FMC, 1>;
   using triplet_dummy = FactorContainer<dummy_factor, dummy_mrf_FMC, 2>;

   using unary_pairwise_dummy = MessageContainer<dummy_message, 0, 1, message_passing_schedule::left, variableMessageNumber, 2, dummy_mrf_FMC, 0>;
   using pairwise_triplet_dummy = MessageContainer<dummy_message, 1, 2, message_passing_schedule::left, variableMessageNumber, 3, dummy_mrf_FMC, 1>;

   using FactorList = meta::list<unary_dummy, pairwise_dummy, triplet_dummy>;
   using MessageList = meta::list<unary_pairwise_dummy, pairwise_triplet_dummy>;
   using ProblemDecompositionList = meta::list<>;
};

int main(int argc, char** argv)
{
   factors_storage<dummy_mrf_FMC> fs;
   messages_storage<dummy_mrf_FMC> ms;

   using unary_dummy = typename dummy_mrf_FMC::unary_dummy;
   using pairwise_dummy = typename dummy_mrf_FMC::pairwise_dummy;
   using triplet_dummy = typename dummy_mrf_FMC::triplet_dummy;
   using unary_pairwise_dummy = typename dummy_mrf_FMC::unary_pairwise_dummy;
   using pairwise_triplet_dummy = typename dummy_mrf_FMC::pairwise_triplet_dummy;

   // build test model with four unaries and pairwise forming a 2x2 grid and two triplets on {0,1,2} and {1,2,3}.
   auto* f1 = fs.template add_factor<unary_dummy>();
   auto* f2 = fs.template add_factor<unary_dummy>();
   auto* f3 = fs.template add_factor<unary_dummy>();
   auto* f4 = fs.template add_factor<unary_dummy>();
   fs.add_factor_relation(f1,f2);
   fs.add_factor_relation(f2,f3);
   fs.add_factor_relation(f3,f4);

   auto add_unary_pairwise_message = [&fs,&ms](auto* f1, auto* f2) 
   {
      auto* f12 = fs.template add_factor<pairwise_dummy>();
      ms.template add_message<unary_pairwise_dummy>(f1, f12);
      fs.add_factor_relation(f1,f12);
      ms.template add_message<unary_pairwise_dummy>(f2, f12);
      fs.add_factor_relation(f12,f2);
      return f12;
   };

   auto* f12 = add_unary_pairwise_message(f1,f2);
   auto* f23 = add_unary_pairwise_message(f2,f3);
   auto* f34 = add_unary_pairwise_message(f3,f4);
   auto* f14 = add_unary_pairwise_message(f1,f4);

   // test for standard pairwise mrf
   {
      auto [forward_sorting, forward_update_sorting] = fs.get_sorted_factors(Direction::forward);
      test(forward_sorting.size() == 8 && forward_update_sorting.size() == 4);
      auto forward_factor_indices = get_factor_indices(forward_sorting.begin(), forward_sorting.end());
      test(forward_factor_indices[f1] < forward_factor_indices[f12] && forward_factor_indices[f12] < forward_factor_indices[f2]);
      test(forward_factor_indices[f2] < forward_factor_indices[f23] && forward_factor_indices[f23] < forward_factor_indices[f3]);
      test(forward_factor_indices[f3] < forward_factor_indices[f34] && forward_factor_indices[f34] < forward_factor_indices[f4]);
      test(forward_factor_indices[f1] < forward_factor_indices[f14] && forward_factor_indices[f14] < forward_factor_indices[f4]);

      auto forward_update_factor_indices = get_factor_indices(forward_update_sorting.begin(), forward_update_sorting.end());
      test(0 == forward_update_factor_indices[f1]);
      test(1 == forward_update_factor_indices[f2]);
      test(2 == forward_update_factor_indices[f3]);
      test(3 == forward_update_factor_indices[f4]);

      auto [omega_forward, receive_mask_forward] = compute_anisotropic_weights(forward_sorting.begin(), forward_sorting.end(), 0.0);
      test(omega_forward.size() == receive_mask_forward.size() && omega_forward.size() == forward_update_sorting.size());

      auto greater_zero = [](const auto x) { return x > 0.0; };

      for(std::size_t i=0; i<4; ++i) {
         test(omega_forward[i].size() == 2 && receive_mask_forward[i].size() == omega_forward[i].size());
      }

      std::cout << omega_forward[0][0] << "," << omega_forward[0][1] << "\n";
      std::cout << int(receive_mask_forward[0][0]) << "," << int(receive_mask_forward[0][1]) << "\n";
      test(std::count_if(omega_forward[0].begin(), omega_forward[0].end(), greater_zero) == 2);
      test(std::count(receive_mask_forward[0].begin(), receive_mask_forward[0].end(), 0) == 2);

      std::cout << omega_forward[1][0] << "," << omega_forward[1][1] << "\n";
      std::cout << int(receive_mask_forward[1][0]) << "," << int(receive_mask_forward[1][1]) << "\n";
      test(std::count_if(omega_forward[1].begin(), omega_forward[1].end(), greater_zero) == 1);
      test(std::count(receive_mask_forward[1].begin(), receive_mask_forward[1].end(), 0) == 1);

      std::cout << omega_forward[2][0] << "," << omega_forward[2][1] << "\n";
      std::cout << int(receive_mask_forward[2][0]) << "," << int(receive_mask_forward[2][1]) << "\n";
      test(std::count_if(omega_forward[2].begin(), omega_forward[2].end(), greater_zero) == 1);
      test(std::count(receive_mask_forward[2].begin(), receive_mask_forward[2].end(), 0) == 1);

      std::cout << omega_forward[3][0] << "," << omega_forward[3][1] << "\n";
      std::cout << int(receive_mask_forward[3][0]) << "," << int(receive_mask_forward[3][1]) << "\n";
      test(std::count_if(omega_forward[3].begin(), omega_forward[3].end(), greater_zero) == 0);
      test(std::count(receive_mask_forward[3].begin(), receive_mask_forward[3].end(), 0) == 0);
   }

   auto add_pairwise_triplet_message = [&fs,&ms](auto* f12, auto* f13, auto* f23)
   {
      auto* f123 = fs.template add_factor<triplet_dummy>();
      ms.template add_message<pairwise_triplet_dummy>(f12, f123);
      fs.add_factor_relation(f12, f123);
      fs.add_factor_relation(f12,f13);
      ms.template add_message<pairwise_triplet_dummy>(f13, f123);
      ms.template add_message<pairwise_triplet_dummy>(f23, f123);
      fs.add_factor_relation(f123, f23);
      fs.add_factor_relation(f13,f23);
      return f123;
   };

   auto* f13 = add_unary_pairwise_message(f1,f3);
   auto* f123 = add_pairwise_triplet_message(f12,f13,f23);
   auto* f134 = add_pairwise_triplet_message(f13,f14,f34);

   // test for mrf with triplet potentials

   // sorting
   auto [forward_sorting, forward_update_sorting] = fs.get_sorted_factors(Direction::forward);
   test(forward_sorting.size() == 11 && forward_update_sorting.size() == 9);
   auto forward_factor_indices = get_factor_indices(forward_sorting.begin(), forward_sorting.end());
   test(forward_factor_indices[f1] < forward_factor_indices[f12] && forward_factor_indices[f12] < forward_factor_indices[f2]);
   test(forward_factor_indices[f2] < forward_factor_indices[f23] && forward_factor_indices[f23] < forward_factor_indices[f3]);
   test(forward_factor_indices[f3] < forward_factor_indices[f34] && forward_factor_indices[f34] < forward_factor_indices[f4]);
   test(forward_factor_indices[f1] < forward_factor_indices[f13] && forward_factor_indices[f13] < forward_factor_indices[f3]);
   test(forward_factor_indices[f12] < forward_factor_indices[f123] && forward_factor_indices[f123] < forward_factor_indices[f23]);
   test(forward_factor_indices[f13] < forward_factor_indices[f134] && forward_factor_indices[f134] < forward_factor_indices[f34]);

   auto [backward_sorting, backward_update_sorting] = fs.get_sorted_factors(Direction::backward);
   test(backward_sorting.size() == 11 && backward_update_sorting.size() == 9);
   auto backward_factor_indices = get_factor_indices(backward_sorting.begin(), backward_sorting.end());
   test(backward_factor_indices[f1] > backward_factor_indices[f12] && backward_factor_indices[f12] > backward_factor_indices[f2]);
   test(backward_factor_indices[f2] > backward_factor_indices[f23] && backward_factor_indices[f23] > backward_factor_indices[f3]);
   test(backward_factor_indices[f3] > backward_factor_indices[f34] && backward_factor_indices[f34] > backward_factor_indices[f4]);
   test(backward_factor_indices[f1] > backward_factor_indices[f13] && backward_factor_indices[f13] > backward_factor_indices[f3]);
   test(backward_factor_indices[f12] > backward_factor_indices[f123] && backward_factor_indices[f123] > backward_factor_indices[f23]);
   test(backward_factor_indices[f13] > backward_factor_indices[f134] && backward_factor_indices[f134] > backward_factor_indices[f34]);


   // anisotropic weight generation
   {
      auto [omega_forward, receive_mask_forward] = compute_anisotropic_weights(forward_sorting.begin(), forward_sorting.end(), 0.0);
      auto forward_update_factor_indices = get_factor_indices(forward_update_sorting.begin(), forward_update_sorting.end());
      test(omega_forward.size() == receive_mask_forward.size() && omega_forward.size() == forward_update_sorting.size());

      auto greater_zero = [](const auto x) { return x > 0.0; };

      test(0 == forward_update_factor_indices[f1]);
      test(omega_forward[0].size() == 3 && receive_mask_forward[0].size() == omega_forward[0].size());
      test(std::count_if(omega_forward[0].begin(), omega_forward[0].end(), greater_zero) == 3);
      test(std::count(receive_mask_forward[0].begin(), receive_mask_forward[0].end(), 0) == 3);

      test(omega_forward[forward_update_factor_indices[f12]].size() == 1 && receive_mask_forward[forward_update_factor_indices[f12]].size() == omega_forward[forward_update_factor_indices[f12]].size());
      test(std::count_if(omega_forward[forward_update_factor_indices[f12]].begin(), omega_forward[forward_update_factor_indices[f12]].end(), greater_zero) == 1);
      test(std::count(receive_mask_forward[forward_update_factor_indices[f12]].begin(), receive_mask_forward[forward_update_factor_indices[f12]].end(), 0) == 1);

      test(omega_forward[forward_update_factor_indices[f13]].size() == 2 && receive_mask_forward[forward_update_factor_indices[f13]].size() == omega_forward[forward_update_factor_indices[f13]].size());
      test(std::count_if(omega_forward[forward_update_factor_indices[f13]].begin(), omega_forward[forward_update_factor_indices[f13]].end(), greater_zero) == 2);
      test(std::count(receive_mask_forward[forward_update_factor_indices[f13]].begin(), receive_mask_forward[forward_update_factor_indices[f13]].end(), 0) == 1);

      test(omega_forward[forward_update_factor_indices[f2]].size() == 2 && receive_mask_forward[forward_update_factor_indices[f2]].size() == omega_forward[forward_update_factor_indices[f2]].size());
      test(std::count_if(omega_forward[forward_update_factor_indices[f2]].begin(), omega_forward[forward_update_factor_indices[f2]].end(), greater_zero) == 1);
      test(std::count(receive_mask_forward[forward_update_factor_indices[f2]].begin(), receive_mask_forward[forward_update_factor_indices[f2]].end(), 0) == 1);


      auto [omega_backward, receive_mask_backward] = compute_anisotropic_weights(backward_sorting.begin(), backward_sorting.end(), 0.0);
      auto backward_update_factor_indices = get_factor_indices(backward_update_sorting.begin(), backward_update_sorting.end());
      test(omega_backward.size() == receive_mask_backward.size() && omega_backward.size() == backward_update_sorting.size());

      test(omega_backward[backward_update_factor_indices[f12]].size() == 1 && receive_mask_backward[backward_update_factor_indices[f12]].size() == omega_backward[backward_update_factor_indices[f12]].size());
      test(std::count_if(omega_forward[backward_update_factor_indices[f12]].begin(), omega_forward[backward_update_factor_indices[f12]].end(), greater_zero) == 0);
      test(std::count(receive_mask_forward[backward_update_factor_indices[f12]].begin(), receive_mask_forward[backward_update_factor_indices[f12]].end(), 0) == 0); 
   }

   // isotropic weight generation
   {
      auto omega_forward = compute_isotropic_weights(forward_sorting.begin(), forward_sorting.end(), 0.0);
      auto receive_mask_forward = compute_full_receive_mask(forward_sorting.begin(), forward_sorting.end());
      test(omega_forward.size() == receive_mask_forward.size());
      for(std::size_t i=0; i<omega_forward.size(); ++i) {
         test(omega_forward[i].size() == forward_update_sorting[i]->no_send_messages());
         test( std::abs(1.0 - std::accumulate(omega_forward[i].begin(), omega_forward[i].end(), 0.0)) <= 1e-8);
         test(receive_mask_forward[i].size() == forward_update_sorting[i]->no_receive_messages());
         for(std::size_t j=0; j<receive_mask_forward[i].size(); ++j) {
            test(receive_mask_forward(i,j) == 1);
         }
      }

   }

   {
      auto omega_forward = compute_isotropic_weights(forward_sorting.begin(), forward_sorting.end(), 0.5);
      for(std::size_t i=0; i<omega_forward.size(); ++i) {
         test(omega_forward[i].size() == forward_update_sorting[i]->no_send_messages());
         test( std::abs(0.5 - std::accumulate(omega_forward[i].begin(), omega_forward[i].end(), 0.0)) <= 1e-8);
      }
   }
}
