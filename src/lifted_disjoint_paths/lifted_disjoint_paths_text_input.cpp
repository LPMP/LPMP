#include "lifted_disjoint_paths/lifted_disjoint_paths_input.h"
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "visitors/standard_visitor.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "solver.hxx"
#include "LP.h"
#include "andres/graph/digraph.hxx"

using namespace LPMP;

int main(int argc, char** argv) {
	ProblemConstructorRoundingSolver<Solver<LP<lifted_disjoint_paths_FMC>,StandardTighteningVisitor>> solver(argc,argv);
	std::string inputFileName=solver.get_input_file();

	lifted_disjoint_paths::LdpInstance input = lifted_disjoint_paths::parse_file(inputFileName);
	solver.GetProblemConstructor().construct(input);
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
