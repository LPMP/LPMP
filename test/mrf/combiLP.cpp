#include "graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "solver.hxx"
#include "combiLP.hxx"
#include "hdf5_routines.hxx"
#include "gurobi_interface.hxx"

int main(int argc, char** argv)
{
  using solver_type = LPMP::MpRoundingSolver<LPMP::Solver<LPMP::FMC_SRMP_T,LPMP::combiLP<DD_ILP::gurobi_interface, LPMP::LP>,LPMP::StandardTighteningVisitor>>;
  std::vector<std::string> options = { 
    {""},
    {"-i"}, {"combiLP.uai"},
    //{"-i"}, {"combiLP_test.uai"},
    //{"-i"}, {"test.h5"},
    {"--maxIter"}, {"250"},
    {"--roundingReparametrization"}, {"anisotropic"},
    {"--standardReparametrization"}, {"anisotropic"},
    //{"--tighten"},
    {"--tightenReparametrization"}, {"uniform:0.5"},
    {"--tightenIteration"}, {"100"},
    {"--tightenInterval"}, {"50"},
    {"--tightenConstraintsMax"}, {"10"}
  };
  solver_type solver(options);
  solver.ReadProblem(LPMP::UaiMrfInput::ParseProblem<LPMP::Solver<LPMP::FMC_SRMP_T,LPMP::combiLP<DD_ILP::gurobi_interface, LPMP::LP>,LPMP::StandardTighteningVisitor>>);
  //solver.ReadProblem(LPMP::ParseOpenGM<LPMP::Solver<LPMP::FMC_SRMP_T,LPMP::combiLP<DD_ILP::gurobi_interface, LPMP::LP>,LPMP::StandardTighteningVisitor>>);

  solver.Solve(); 
}
