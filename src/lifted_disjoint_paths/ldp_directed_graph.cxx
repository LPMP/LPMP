#include "lifted_disjoint_paths/ldp_directed_graph.hxx"

namespace LPMP {
//LdpDirectedGraph::LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts){
//    //TODO the same as in the previous constructor but addd s and t nodes and edges with given cost
//    LdpDirectedGraph(edges,inputEdgeCosts,std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
//}


LdpDirectedGraph::LdpDirectedGraph(const andres::graph::Digraph<>& inputGraph,const std::vector<double>& inputEdgeCosts){
    numberOfVertices=inputGraph.numberOfVertices();
    std::vector<std::size_t> adjacencyForward(numberOfVertices);
    std::vector<std::size_t> adjacencyBackward(numberOfVertices);

    assert(inputGraph.numberOfEdges()==inputEdgeCosts.size());
    for(size_t i=0;i<numberOfVertices;i++){
        size_t neForward=inputGraph.numberOfEdgesFromVertex(i);
        adjacencyForward[i]=neForward;
        size_t neBackward=inputGraph.numberOfEdgesToVertex(i);
        adjacencyBackward[i]=neBackward;
    }

    forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());

    backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());

    for(size_t i=0;i<numberOfVertices;i++){
        size_t neForward=inputGraph.numberOfEdgesFromVertex(i);
        for(size_t j=0;j<neForward;j++){
            size_t e=inputGraph.edgeFromVertex(i,j);
            size_t w=inputGraph.vertexFromVertex(i,j);
            forwardEdges[i][j]={w,inputEdgeCosts[e]};

        }
        size_t neBackward=inputGraph.numberOfEdgesToVertex(i);
        for(size_t j=0;j<neBackward;j++){
            size_t w=inputGraph.vertexToVertex(i,j);
            size_t e=inputGraph.edgeToVertex(i,j);
            backwardEdges[i][j]={w,inputEdgeCosts[e]};
        }
    }

    for (size_t i=0;i<numberOfVertices;i++) {
        std::sort(forwardEdges[i].begin(),forwardEdges[i].end());
        std::sort(backwardEdges[i].begin(),backwardEdges[i].end());
    }
    //Need to sort within edges?
}


//LdpDirectedGraph::LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts,double inCost,double outCost){

//    bool addST=inCost<=std::numeric_limits<double>::max();

//    std::vector<std::size_t> adjacencyForward;
//    std::vector<std::size_t> adjacencyBackward;


//    assert(edges.size()==inputEdgeCosts.size());

//    // first determine size for adjacency_list
//    for(const std::array<size_t,2>& e:edges){
//        size_t i=e[0];
//        size_t j=e[1];
//        size_t size=std::max({i+1,j+1,adjacencyForward.size()});
//        adjacencyBackward.resize(size);
//        adjacencyForward.resize(size);
//        adjacencyForward[i]++;
//        adjacencyBackward[j]++;
//    }
//    numberOfVertices=adjacencyForward.size();

//    if(addST){
//        for(size_t i=0;i<numberOfVertices;i++){
//            adjacencyForward[i]++;
//            adjacencyBackward[i]++;
//        }

//        adjacencyForward.push_back(numberOfVertices); //for s
//        adjacencyBackward.push_back(0);
//        adjacencyForward.push_back(0);  //for t
//        adjacencyBackward.push_back(numberOfVertices);

//    }
//    forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());

//    backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());



//     std::fill(adjacencyForward.begin(), adjacencyForward.end(), 0);
//     std::fill(adjacencyBackward.begin(), adjacencyBackward.end(), 0);


//    for(size_t i=0;i<edges.size();i++){
//        size_t v=edges[i][0];
//        size_t w=edges[i][1];
//        forwardEdges[v][adjacencyForward[v]]={w,inputEdgeCosts[i]};
//        adjacencyForward[v]++;
//        backwardEdges[w][adjacencyBackward[w]]={v,inputEdgeCosts[i]};
//        adjacencyBackward[w]++;
//    }
//    if(addST){
//        size_t s=numberOfVertices;
//        size_t t=numberOfVertices+1;

//        for(size_t i=0;i<numberOfVertices;i++){
//            forwardEdges[i][adjacencyForward[i]]={t,outCost};
//            backwardEdges[i][adjacencyBackward[i]]={s,inCost};

//        }
//        numberOfVertices+=2;
//    }

//    for (size_t i=0;i<numberOfVertices;i++) {
//        std::sort(forwardEdges[i].begin(),forwardEdges[i].end());
//        std::sort(backwardEdges[i].begin(),backwardEdges[i].end());
//    }


//    //Need to sort within edges?
//}
}
