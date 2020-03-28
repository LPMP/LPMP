#ifndef LPMP_TEST_MESSAGE_HXX
#define LPMP_TEST_MESSAGE_HXX

#include <random>
#include "function_existence.hxx"
#include "config.hxx"
#include "template_utilities.hxx"
#include "test.h"

namespace LPMP {

   template<typename T>
   std::enable_if_t<std::is_same_v<T,double>> initialize_data_randomly(T& x, std::random_device& rd)
   {
      std::mt19937 gen(rd());
      std::normal_distribution<> dis(0.0, 1.0);

      x = dis(gen);
   }

   template<typename T>
   std::enable_if_t<!std::is_same_v<T,double>> initialize_data_randomly(T& x, std::random_device& rd)
   {
      std::mt19937 gen(rd());
      std::normal_distribution<> dis(0.0, 1.0);

      for(std::size_t i=0; i<x.size(); ++i) {
         x[i] = dis(gen);
      }
   }

   template<typename FACTOR>
   void randomly_initialize_factor(FACTOR& f, std::random_device& rd)
   {
      auto dual = f.export_variables();

      for_each_tuple(dual, [&](auto& x) {
         initialize_data_randomly(x,rd);
      });
   }

   template<typename MESSAGE_VECTOR>
   void randomly_initialize_message(MESSAGE_VECTOR& msg_vec, std::random_device& rd)
   {
      std::mt19937 gen(rd());
      std::normal_distribution<> dis(0, 1);

      for(std::size_t i=0; i<msg_vec.size(); ++i) {
         msg_vec[i] = dis(gen);
      }
   }

   namespace test_message_detail {
      LPMP_FUNCTION_EXISTENCE_CLASS(has_repam_left,RepamLeft)
      LPMP_FUNCTION_EXISTENCE_CLASS(has_repam_right,RepamRight) 
   }

   template<typename MESSAGE, typename MESSAGE_VECTOR, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void reparametrize(const MESSAGE& m, const MESSAGE_VECTOR& msg_vec, LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      if constexpr(test_message_detail::has_repam_left<MESSAGE,void,LEFT_FACTOR,REAL,INDEX>()) {
         for(std::size_t i=0; i<msg_vec.size(); ++i) {
            m.RepamLeft(l, msg_vec[i], i);
         }
      } else {
         m.RepamLeft(l, msg_vec);
      }


      if constexpr(test_message_detail::has_repam_right<MESSAGE,void,RIGHT_FACTOR&,REAL,INDEX>()) {
         for(std::size_t i=0; i<msg_vec.size(); ++i) {
            m.RepamRight(r, -msg_vec[i], i);
         }
      } else {
         m.RepamRight(r, -1.0*msg_vec);
      }
   }

   template<Chirality C = Chirality::left, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename MESSAGE, typename MESSAGE_VECTOR>
   void test_repam_message(LEFT_FACTOR& l, RIGHT_FACTOR& r, MESSAGE& m, MESSAGE_VECTOR& msg_vec, std::random_device& rd)
   {
      randomly_initialize_factor(l, rd);
      randomly_initialize_factor(r, rd);

      l.init_primal();
      r.init_primal();
      if constexpr(C == Chirality::left) {
         l.MaximizePotentialAndComputePrimal();
         m.ComputeRightFromLeftPrimal(l,r);
         r.MaximizePotentialAndComputePrimal();
      } else {
         r.MaximizePotentialAndComputePrimal();
         m.ComputeLeftFromRightPrimal(l,r);
         l.MaximizePotentialAndComputePrimal();
      }

      // test primal on left and right: cost must stay invariant
      const double prev_cost = l.EvaluatePrimal() + r.EvaluatePrimal();
      randomly_initialize_message(msg_vec, rd);
      reparametrize(m, msg_vec, l, r);
      const double after_cost = l.EvaluatePrimal() + r.EvaluatePrimal();

      test(m.CheckPrimalConsistency(l,r));
      test(std::abs(prev_cost - after_cost) <= 1e-8); 
   }

   template<Chirality C, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename MESSAGE, typename MESSAGE_VECTOR>
   void test_send_message_to(LEFT_FACTOR& l, RIGHT_FACTOR& r, MESSAGE& m, MESSAGE_VECTOR& msg_vec, std::random_device& rd)
   {
      randomly_initialize_factor(l, rd);
      randomly_initialize_factor(r, rd);

      l.init_primal();
      r.init_primal();

      const double prev_lb = l.LowerBound() + r.LowerBound();

      std::fill(msg_vec.begin(), msg_vec.end(), 0.0);
      if constexpr(C == Chirality::left) {
         m.send_message_to_left(r, msg_vec, 1.0);
         reparametrize(m, -msg_vec, l, r);
      } else {
         m.send_message_to_right(l, msg_vec, 1.0);
         reparametrize(m, msg_vec, l, r);
      }


      if constexpr(C == Chirality::left) {
         l.MaximizePotentialAndComputePrimal();
         m.ComputeRightFromLeftPrimal(l,r);
         test(m.CheckPrimalConsistency(l,r)); 
         r.MaximizePotentialAndComputePrimal();
      } else {
         r.MaximizePotentialAndComputePrimal();
         m.ComputeLeftFromRightPrimal(l,r);
         test(m.CheckPrimalConsistency(l,r)); 
         l.MaximizePotentialAndComputePrimal();
      }

      const double after_lb = l.LowerBound() + r.LowerBound();

      test(prev_lb <= after_lb + 1e-8);
      test(std::abs(l.LowerBound() - l.EvaluatePrimal()) <= 1e-8);
      test(std::abs(r.LowerBound() - r.EvaluatePrimal()) <= 1e-8);
      test(m.CheckPrimalConsistency(l,r)); 
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename MESSAGE, typename MESSAGE_VECTOR>
   void test_send_message(LEFT_FACTOR& l, RIGHT_FACTOR& r, MESSAGE& m, MESSAGE_VECTOR& msg_vec, std::random_device& rd)
   {
      test_send_message_to<Chirality::left>(l,r,m,msg_vec,rd);
      test_send_message_to<Chirality::right>(l,r,m,msg_vec,rd);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename MESSAGE, typename MESSAGE_VECTOR>
   void test_message(LEFT_FACTOR& l, RIGHT_FACTOR& r, MESSAGE& m, MESSAGE_VECTOR& msg_vec, std::random_device& rd)
   {
      test_repam_message(l,r,m,msg_vec,rd);
      test_send_message_to<Chirality::left>(l,r,m,msg_vec,rd);
      test_send_message_to<Chirality::right>(l,r,m,msg_vec,rd);
   }


   template<typename FACTOR>
   void test_factor(FACTOR& f, std::random_device& rd)
   {
      randomly_initialize_factor(f, rd);
      f.init_primal();
      f.MaximizePotentialAndComputePrimal();
      test(std::abs(f.EvaluatePrimal() - f.LowerBound()) < 1e-8);
   }

} // namespace LPMP

#endif // LPMP_TEST_MESSAGE_HXX
