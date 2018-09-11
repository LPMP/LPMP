#include "test_model.hxx"
#include "LP_FWMAP.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "test.h"

using namespace LPMP;

int main(int argc, char** argv)
{
  Solver<LP_tree_FWMAP<test_FMC>, StandardVisitor> s;
  auto& lp = s.GetLP();

  build_test_model(lp);

  s.Solve();

  test( std::abs(s.GetLP().decomposition_lower_bound() - 1.0) <= eps );

  s.GetLP().write_back_reparametrization();
  test(std::abs(s.GetLP().original_factors_lower_bound() - 1.0) <= eps);
}


