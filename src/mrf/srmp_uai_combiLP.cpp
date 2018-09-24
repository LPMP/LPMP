
#include "mrf/graphical_model.h"
#include "mrf/combiLP.hxx"
#include "gurobi_interface.hxx"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_uai_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
MpRoundingSolver<Solver<combiLP<DD_ILP::gurobi_interface, LP<FMC_SRMP>>,StandardVisitor>> solver(argc,argv);
auto input = mrf_uai_input::parse_file(solver.get_input_file());
solver.template GetProblemConstructor<0>().construct(input);
return solver.Solve();
}