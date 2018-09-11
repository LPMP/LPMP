#include "config.hxx"
#include "factors_messages.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "test.h"
#include "LP_external_interface.hxx"
#include "test_model.hxx"
#include <random>

using namespace LPMP; 

// build simple model and check properties of factors/messages
int main()
{

  //std::vector<std::string> options({{"test solver"}, {"--maxIter"}, {"5"}});
  //MpRoundingSolver<Solver<test_FMC, LP, StandardVisitor>> s(options);
    {
        Solver<LP_external_solver<DD_ILP::problem_export, LP<test_FMC>>, StandardVisitor> s;
        auto& lp = s.GetLP();
        auto* f1 = lp.template add_factor<typename test_FMC::factor>(0,1);
        auto* f2 = lp.template add_factor<typename test_FMC::factor>(1,0);
        auto* f3 = lp.template add_factor<typename test_FMC::factor>(0,0);

        auto* m12 = lp.template add_message<typename test_FMC::message>(f1,f2);
        auto* m13 = lp.template add_message<typename test_FMC::message>(f1,f3);

        test(lp.number_of_factors() == 3);
        test(lp.number_of_messages() == 2);
        test(lp.get_factor(0) == f1);
        test(lp.get_factor(1) == f2);
        test(lp.get_factor(2) == f3);

        test(f1->no_messages() == 2);
        test(f1->no_send_messages() == 2);

        test(f2->no_messages() == 1);
        test(f2->no_send_messages() == 0);

        test(f3->no_messages() == 1);
        test(f3->no_send_messages() == 0);

        std::cout << "lower bound before optimization = " << s.GetLP().LowerBound() << "\n";
        s.Solve();
        std::cout << "lower bound after optimization = " << s.GetLP().LowerBound() << "\n";
        test(std::abs(s.GetLP().LowerBound() - 1.0) <= eps);
        s.GetLP().get_external_solver().write_to_file("test_problem.lp");
    }

   {
       //Solver<LP<test_FMC>, StandardVisitor> s;
       //auto& lp = s.GetLP();

       //build_test_model(lp);

       //s.Solve();

   }
}
