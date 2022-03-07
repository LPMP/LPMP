#include "test_mrf.hxx"
#include "LP_external_interface.hxx"

using namespace LPMP;

int main()
{
   auto sat_solver_options = solver_options;
   sat_solver_options[12] = std::string("uniform:0.5");

   using VisitorType = StandardVisitor;

   auto pos_potts = construct_potts(2,2, 0.0, 1.0);
   auto neg_potts = construct_potts(2,2, 1.0, 0.0);

   std::cout << "sat based reduction does not work yet\n";
   return 0;

   // SAT rounding with duality gap
   {
       using FMC = FMC_SRMP;
       using SolverType = MpRoundingSolver<Solver<LP_external_solver<DD_ILP::sat_solver,LP<FMC>>,VisitorType>>;

       SolverType s(sat_solver_options);
       auto& mrf = s.GetProblemConstructor();

       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));

       // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts

       mrf.add_pairwise_factor(0,1,neg_potts);
       mrf.add_pairwise_factor(1,2,pos_potts);
       mrf.add_pairwise_factor(2,3,pos_potts);
       mrf.add_pairwise_factor(3,4,pos_potts);
       mrf.add_pairwise_factor(0,4,pos_potts);

       s.Solve();
       test(std::abs(s.lower_bound() - 0.0) <= eps);
       test(std::abs(s.primal_cost() - 1.0) <= eps);
   }

   // SAT rounding without duality gap
   {
       using FMC = FMC_SRMP;
       using SolverType = MpRoundingSolver<Solver<LP_external_solver<DD_ILP::sat_solver,LP<FMC>>,VisitorType>>;

       SolverType s(sat_solver_options);
       auto& mrf = s.GetProblemConstructor();

       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));
       mrf.add_unary_factor(std::vector<REAL>(2,0.0));

       // make a chain of length 4 visiting each label once -> one negative Potts and three positive Potts

       mrf.add_pairwise_factor(0,1,neg_potts);
       mrf.add_pairwise_factor(1,2,pos_potts);
       mrf.add_pairwise_factor(2,3,pos_potts);
       mrf.add_pairwise_factor(3,4,pos_potts);

       s.Solve();
       test(std::abs(s.lower_bound() - 0.0) <= eps);
       test(std::abs(s.primal_cost() - 0.0) <= eps);
   }

    // SAT rounding with triplets
    {
       using FMC = FMC_SRMP_T;
        using SatSolverType = MpRoundingSolver<Solver<LP_external_solver<DD_ILP::sat_solver,LP<FMC>>,VisitorType>>;

        SatSolverType s_sat(solver_options);
        auto& mrf_sat = s_sat.GetProblemConstructor();

        mrf_sat.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf_sat.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf_sat.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf_sat.add_unary_factor(std::vector<REAL>(2,0.0));

        // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts
        mrf_sat.add_pairwise_factor(0,1,neg_potts);
        mrf_sat.add_pairwise_factor(1,2,pos_potts);
        mrf_sat.add_pairwise_factor(2,3,pos_potts);
        mrf_sat.add_pairwise_factor(0,3,pos_potts);

        mrf_sat.add_tightening_triplet(0,1,2);
        mrf_sat.add_tightening_triplet(0,1,3);
        mrf_sat.add_tightening_triplet(0,2,3);
        mrf_sat.add_tightening_triplet(1,2,3);

        s_sat.Solve();
        test(std::abs(s_sat.lower_bound() - 1.0) <= eps);
        test(std::abs(s_sat.primal_cost() - 1.0) <= eps);
    }
}
