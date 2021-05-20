#include "lifted_disjoint_paths/ldp_directed_graph.hxx"

namespace LPMP {
//LdpDirectedGraph::LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts){
//    //TODO the same as in the previous constructor but addd s and t nodes and edges with given cost
//    LdpDirectedGraph(edges,inputEdgeCosts,std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
//}


LdpDirectedGraph::LdpDirectedGraph(const LdpDirectedGraph& inputGraph,double inputEdgeCost,double outputEdgeCost){


    //std::cout<<"graph constructor "<<std::endl;
    numberOfVertices=inputGraph.getNumberOfVertices()+2;
    numberOfEdges=0;
    size_t numberOfBackwardEdges=0;

    //std::cout<<"vertices "<<numberOfVertices<<std::endl;
    size_t s=numberOfVertices-2;
    size_t t=numberOfVertices-1;
    std::vector<std::size_t> adjacencyForward(numberOfVertices);
    std::vector<std::size_t> adjacencyBackward(numberOfVertices);
    for(size_t i=0;i<numberOfVertices-2;i++){
        //std::cout<<"for "<<i<<std::endl;
        adjacencyForward[i]=inputGraph.getNumberOfEdgesFromVertex(i)+1;
        //std::cout<<"adj forward "<<i<<": "<<adjacencyForward[i]<<std::endl;
        adjacencyBackward[i]=inputGraph.getNumberOfEdgesToVertex(i)+1;
        //std::cout<<"adj backward "<<i<<": "<<adjacencyBackward[i]<<std::endl;
    }
    //std::cout<<"after for "<<std::endl;
    adjacencyForward[s]=numberOfVertices-2;  //vertices from s
    adjacencyBackward[t]=numberOfVertices-2;              //vertices to t
    forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());
    backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());

   // std::cout<<"adjacency ok"<<std::endl;
    for (size_t i = 0; i < numberOfVertices-2; ++i) {
        forwardEdges[s][i]={i,inputEdgeCost,numberOfVertices};

        backwardEdges[t][i]={i,outputEdgeCost,numberOfVertices};
        size_t j = 0;
        for (; j < adjacencyForward[i]-1; ++j) {
            forwardEdges[i][j]=inputGraph.getForwardEdges()[i][j];
            numberOfEdges++;
        }
        forwardEdges[i][j]={t,outputEdgeCost,numberOfVertices};
        numberOfEdges++;
        j = 0;
        for (; j < adjacencyBackward[i]-1; ++j) {
            backwardEdges[i][j]=inputGraph.getBackwardEdges()[i][j];
            numberOfBackwardEdges++;
        }
        backwardEdges[i][j]={s,inputEdgeCost,numberOfVertices};
        numberOfBackwardEdges++;

    }
   // std::cout<<"edges added"<<std::endl;

    assert(numberOfBackwardEdges==numberOfEdges);

    for (size_t i=0;i<numberOfVertices;i++) {
        std::sort(forwardEdges[i].begin(),forwardEdges[i].end());
        std::sort(backwardEdges[i].begin(),backwardEdges[i].end());
    }

    setNeighborPointers();
}


void LdpDirectedGraph::setNeighborPointers(){
    std::vector<size_t> backCounters(numberOfVertices,0);
    for (size_t i = 0; i < numberOfVertices; ++i) {
        edge* iterForward=forwardEdges[i].begin();
        size_t counter=0;
        while(iterForward!=forwardEdges[i].end()){
            size_t neighborID=iterForward->first;
            edge& backwardEdge=backwardEdges[neighborID][backCounters[neighborID]];
            assert(backwardEdge.first==i);
            iterForward->reverse_neighbor_index=backCounters[neighborID];
            backwardEdge.reverse_neighbor_index=counter;
            backCounters[neighborID]++;
            counter++;
            iterForward++;
        }
    }
}

}
