
#include "lifted_disjoint_paths/lifted_disjoint_paths_fmc.h"
#include "visitors/standard_visitor.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "solver.hxx"
#include "LP.h"
#include "andres/graph/digraph.hxx"
#include "lifted_disjoint_paths/ldp_parameters.hxx"
#include "lifted_disjoint_paths/ldp_complete_structure.hxx"

using namespace LPMP;

int main(int argc, char** argv) {
//    std::map<size_t,std::map<size_t,double>> edges;
//    edges[0][50]=-5;
//    edges[0][70]=-1;
//    edges[0][80]=2;
//    edges[10][50]=-1;
//    edges[10][60]=-3;
//    edges[10][70]=2;
//    edges[20][60]=-2;
//    edges[20][70]=-3;
//    edges[20][80]=2;
//    edges[30][50]=1;
//    edges[30][70]=3;
//    edges[30][80]=-4;
//    edges[40][50]=2;
//    edges[40][60]=3;
//    edges[40][80]=1;

//    size_t v=30;
//    size_t w=70;
//    double liftedCost=-10;
//    ldp_cut_factor<lifted_disjoint_paths::LdpInstance> cutFactor(v,w,liftedCost,edges);
//    double lb=cutFactor.LowerBound();
//    std::cout<<"lower bound "<<lb<<std::endl;
//    for(size_t i=0;i<cutFactor.cutGraph.getNumberOfInputs();i++){
//        auto * it=cutFactor.cutGraph.forwardNeighborsBegin(i);
//        auto * end=cutFactor.cutGraph.forwardNeighborsEnd(i);
//        for(;it!=end;it++){
//            size_t j=it->first;
//            std::cout<<"min marginal "<<i<<","<<j<<std::endl;
//            double mm=cutFactor.getOneEdgeMinMarginal(i,j);
//            std::cout<<": "<<mm<<std::endl;
//        }
//    }




    ProblemConstructorRoundingSolver<Solver<LP<lifted_disjoint_paths_FMC>,StandardTighteningVisitor>> solver(argc,argv);
    std::string inputFileName=solver.get_input_file();

    LPMP::lifted_disjoint_paths::LdpParameters<> configParams(inputFileName);

    LPMP::CompleteStructure<> completeStructure(configParams);

    LPMP::lifted_disjoint_paths::LdpInstance ldpInstance(configParams,completeStructure);


    solver.GetProblemConstructor().construct(ldpInstance);

    return solver.Solve();



}
