#include "lifted_disjoint_paths/lifted_disjoint_paths_input.h"
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "visitors/standard_visitor.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "solver.hxx"
#include "LP.h"

using namespace LPMP;

int main(int argc, char** argv) {
	//Solver<LP<FMC_LIFTED_DISJOINT_PATHS,StandardTighteningVisitor>> solver(argc,argv); //orig

	Solver<LP<lifted_disjoint_paths_FMC>,StandardTighteningVisitor> solver(argc,argv);
		std::string inputFileName=solver.get_input_file();
    const lifted_disjoint_paths::LdpInstance input = lifted_disjoint_paths::parse_file(inputFileName);
    solver.GetProblemConstructor().construct(input);
    return solver.Solve();
}
