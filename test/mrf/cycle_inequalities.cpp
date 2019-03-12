#include "test.h"
#include <vector>
#include "mrf/graphical_model.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "mrf/cycle_inequalities.hxx"

using namespace LPMP;

// cycle inequalities tightening for MAP-MRF"
int main()
{

   // conditional minima
  {
      matrix<REAL> m(4,4);
      INDEX c=0;
      for(INDEX x1=0; x1<4; ++x1) {
         for(INDEX x2=0; x2<4; ++x2) {
            m(x1,x2) = c;
            ++c;
         }
      }
      const auto row_minima = k_ary_cycle_inequalities_search<int>::row_minima(m);
      const auto column_minima = k_ary_cycle_inequalities_search<int>::column_minima(m);
      const auto principal_minima = k_ary_cycle_inequalities_search<int>::principal_minima(m, column_minima);

      test(row_minima(0,0) == 0);
      test(row_minima(0,1) == 1);
      test(row_minima(1,0) == 4);
      test(row_minima(1,1) == 5);
      test(row_minima(2,0) == 8);
      test(row_minima(2,1) == 9);
      test(row_minima(3,0) == 12);
      test(row_minima(3,1) == 13);

      test(column_minima(0,0) == 0);
      test(column_minima(0,1) == 4);
      test(column_minima(1,0) == 1);
      test(column_minima(1,1) == 5);
      test(column_minima(2,0) == 2);
      test(column_minima(2,1) == 6);
      test(column_minima(3,0) == 3);
      test(column_minima(3,1) == 7);

      test(principal_minima(0,0) == 5);
      test(principal_minima(0,1) == 4);
      test(principal_minima(0,2) == 4);
      test(principal_minima(0,3) == 4);

      test(principal_minima(1,0) == 1);
      test(principal_minima(1,1) == 0);
      test(principal_minima(1,2) == 0);
      test(principal_minima(1,3) == 0);

      test(principal_minima(2,0) == 1);
      test(principal_minima(2,1) == 0);
      test(principal_minima(2,2) == 0);
      test(principal_minima(2,3) == 0);
   
      test(principal_minima(3,0) == 1);
      test(principal_minima(3,1) == 0);
      test(principal_minima(3,2) == 0);
      test(principal_minima(3,3) == 0);
   }

   std::vector<std::string> i = {
      "",
      "-i", "",
      "-v", "2"
   };

   matrix<REAL> negPotts2(2,2);
   negPotts2(0,0) = 1.0;
   negPotts2(1,1) = 1.0;
   negPotts2(0,1) = 0.0;
   negPotts2(1,0) = 0.0;

   matrix<REAL> posPotts2(2,2); 
   posPotts2(0,0) = 0.0;
   posPotts2(1,1) = 0.0;
   posPotts2(0,1) = 1.0;
   posPotts2(1,0) = 1.0;

   // binary violated 4 cycle
   {
      Solver<LP<FMC_SRMP_T>,StandardVisitor> s(i);
      auto& mrf = s.GetProblemConstructor();
      s.GetLP().set_reparametrization(lp_reparametrization(lp_reparametrization_mode::Anisotropic, 0.1));

      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});
      
      mrf.add_pairwise_factor(0,1,negPotts2);
      mrf.add_pairwise_factor(1,2,posPotts2);
      mrf.add_pairwise_factor(2,3,posPotts2);
      mrf.add_pairwise_factor(0,3,posPotts2);

      mrf.Tighten(100);
      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type> cycle_search(mrf);

      auto triplets = cycle_search.search();
      test(triplets.size() >= 2);
      mrf.add_triplets(triplets);

      for(std::size_t i=0; i<100; ++i) {
         s.GetLP().ComputePass();
      }
      std::cout << "omega forward\n";
      std::cout << s.GetLP().LowerBound() << std::endl;
      const auto w = s.GetLP().get_message_passing_weight(lp_reparametrization(lp_reparametrization_mode::Uniform, 0.1));
      for(std::size_t i=0; i<w.omega_forward.size(); ++i) {
         for(std::size_t j=0; j<w.omega_forward[i].size(); ++j) {
            std::cout << w.omega_forward[i][j] << ", ";
         }
         std::cout << "\n";
      }
      std::cout << "receive mask forward\n";
      for(std::size_t i=0; i<w.receive_mask_forward.size(); ++i) {
         for(std::size_t j=0; j<w.receive_mask_forward[i].size(); ++j) {
            std::cout << int(w.receive_mask_forward[i][j]) << ", ";
         }
         std::cout << "\n";
      }
      test(s.GetLP().LowerBound() > 1.0-eps);
   }

   // binary violated 5 cycle
   {
     Solver<LP<FMC_SRMP_T>,StandardVisitor> s(i);
     auto& mrf = s.GetProblemConstructor();
     s.GetLP().set_reparametrization(lp_reparametrization(lp_reparametrization_mode::Uniform, 0.5)); // setting reparametrization mode to anisotropic leads to suboptimal fixed point

      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});
      mrf.add_unary_factor({0,0});

      
      mrf.add_pairwise_factor(0,1,negPotts2);
      mrf.add_pairwise_factor(1,2,posPotts2);
      mrf.add_pairwise_factor(2,3,posPotts2);
      mrf.add_pairwise_factor(3,4,posPotts2);
      mrf.add_pairwise_factor(0,4,posPotts2);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type> cycle_search(mrf);

      auto triplets = cycle_search.search();

      test(triplets.size() >= 3);
      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass();
      }
      std::cout << s.GetLP().LowerBound() << std::endl;
      test(s.GetLP().LowerBound() > 1.0-eps);
      
   }

   matrix<REAL> negPotts4(4,4);
   negPotts4(0,0) = 1.0; negPotts4(0,1) = 1.0; negPotts4(0,2) = 0.0; negPotts4(0,3) = 0.0;
   negPotts4(1,0) = 1.0; negPotts4(1,1) = 1.0; negPotts4(1,2) = 0.0; negPotts4(1,3) = 0.0;
   negPotts4(2,0) = 0.0; negPotts4(2,1) = 0.0; negPotts4(2,2) = 1.0; negPotts4(2,3) = 1.0;
   negPotts4(3,0) = 0.0; negPotts4(3,1) = 0.0; negPotts4(3,2) = 1.0; negPotts4(3,3) = 1.0;

   matrix<REAL> posPotts4(4,4);
   posPotts4(0,0) = 0.0; posPotts4(0,1) = 0.0; posPotts4(0,2) = 1.0; posPotts4(0,3) = 1.0;
   posPotts4(1,0) = 0.0; posPotts4(1,1) = 0.0; posPotts4(1,2) = 1.0; posPotts4(1,3) = 1.0;
   posPotts4(2,0) = 1.0; posPotts4(2,1) = 1.0; posPotts4(2,2) = 0.0; posPotts4(2,3) = 0.0;
   posPotts4(3,0) = 1.0; posPotts4(3,1) = 1.0; posPotts4(3,2) = 0.0; posPotts4(3,3) = 0.0;

   // expanded k-ary cycle search 4-cycle
   {
     Solver<LP<FMC_SRMP_T>,StandardVisitor> s(i);
     auto& mrf = s.GetProblemConstructor();
     s.GetLP().set_reparametrization(lp_reparametrization(lp_reparametrization_mode::Uniform, 0.5)); // setting reparametrization mode to anisotropic leads to suboptimal fixed point

      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});

      
      mrf.add_pairwise_factor(0,1,negPotts4);
      mrf.add_pairwise_factor(1,2,posPotts4);
      mrf.add_pairwise_factor(2,3,posPotts4);
      mrf.add_pairwise_factor(0,3,posPotts4);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,false> cycle_search(mrf);
      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,true> cycle_search2(mrf);

      auto triplets = cycle_search.search();
      test(triplets.size() == 0);
      triplets = cycle_search2.search();
      test(triplets.size() >= 2);

      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass();
      }
      test(s.GetLP().LowerBound() > 1.0-eps);
   }

   // expanded k-ary cycle search 5-cycle
   {
     Solver<LP<FMC_SRMP_T>,StandardVisitor> s(i);
     auto& mrf = s.GetProblemConstructor();
     s.GetLP().set_reparametrization(lp_reparametrization(lp_reparametrization_mode::Uniform, 0.5)); // setting reparametrization mode to anisotropic leads to suboptimal fixed point

      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});
      mrf.add_unary_factor({0,0,0,0});

      
      mrf.add_pairwise_factor(0,1,negPotts4);
      mrf.add_pairwise_factor(1,2,posPotts4);
      mrf.add_pairwise_factor(2,3,posPotts4);
      mrf.add_pairwise_factor(3,4,posPotts4);
      mrf.add_pairwise_factor(0,4,posPotts4);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,false> cycle_search(mrf);
      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,true> cycle_search2(mrf);

      auto triplets = cycle_search.search();
      test(triplets.size() == 0);
      triplets = cycle_search2.search();
      test(triplets.size() >= 3);

      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass();
      }
      test(s.GetLP().LowerBound() > 1.0-eps);
   }
}

