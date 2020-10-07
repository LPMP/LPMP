#include "lifted_disjoint_paths/lifted_disjoint_paths_input.h"
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "visitors/standard_visitor.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "solver.hxx"
#include "LP.h"
#include "andres/graph/digraph.hxx"
#include "lifted_disjoint_paths/ldp_config.hxx"

using namespace LPMP;

int main(int argc, char** argv) {

//    std::vector<std::string> args;
//    for (int i = 0; i < argc; i++){
//        args.push_back(argv[i]);
//    }
//    for (int i = 0; i < args.size(); ++i) {
//        std::cout<<args[i]<<":";
//    }


	ProblemConstructorRoundingSolver<Solver<LP<lifted_disjoint_paths_FMC>,StandardTighteningVisitor>> solver(argc,argv);
	std::string inputFileName=solver.get_input_file();

	LPMP::lifted_disjoint_paths::ConfigDisjoint<> configParams(inputFileName);
	LPMP::lifted_disjoint_paths::LdpInstance ldpInstance(configParams);

	//const lifted_disjoint_paths::LdpInstance input = lifted_disjoint_paths::parse_file(inputFileName);
	solver.GetProblemConstructor().construct(ldpInstance);

	return solver.Solve();

//    const andres::graph::Digraph<>& graph=input.getGraph();
//    std::cout<<"base graph"<<std::endl;
//    for (int i = 0; i < 50; ++i) {
//    	std::cout<<graph.vertexOfEdge(i,0)<<", "<<graph.vertexOfEdge(i,1)<<": "<<input.getEdgeScore(i)<<std::endl;
//    }
//
//    const andres::graph::Digraph<>& liftedGraph=input.getGraphLifted();
//    std::cout<<"lifted graph"<<std::endl;
//    for (int i = 0; i < 50; ++i) {
//    	std::cout<<liftedGraph.vertexOfEdge(i,0)<<", "<<liftedGraph.vertexOfEdge(i,1)<<": "<<input.getLiftedEdgeScore(i)<<std::endl;
//    }



}
