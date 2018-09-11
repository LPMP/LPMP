#include "graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "LP_external_interface.hxx"
#include "solver.hxx"
#include "combiLP.hxx"
#include "hdf5_routines.hxx"
#include "gurobi_interface.hxx"

using namespace LPMP;

int main(int argc, char** argv)
{
  using solver_type = Solver<FMC_SRMP,LP_external_solver<DD_ILP::gurobi_interface,LP>,StandardTighteningVisitor>;
  std::vector<std::string> options = { 
    {""},
    {"-i"}, {"combiLP.uai"}
    //{"-i"}, {"test.h5"}
    //{"-i"}, {"6000032.h5"}
  };
  solver_type solver(options);
  solver.ReadProblem(UaiMrfInput::ParseProblem<solver_type>);
  //solver.ReadProblem(LPMP::ParseOpenGM<solver_type>);
  solver.GetLP().write_to_file("dd_ilp_gurobi.lp");
  solver.GetLP().solve();
}
