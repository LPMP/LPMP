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
        forwardEdges[s][i]={i,inputEdgeCost};

        backwardEdges[t][i]={i,outputEdgeCost};
        size_t j = 0;
        for (; j < adjacencyForward[i]-1; ++j) {
            forwardEdges[i][j]=inputGraph.getForwardEdges()[i][j];
            numberOfEdges++;
        }
        forwardEdges[i][j]={t,outputEdgeCost};
        numberOfEdges++;
        j = 0;
        for (; j < adjacencyBackward[i]-1; ++j) {
            backwardEdges[i][j]=inputGraph.getBackwardEdges()[i][j];
            numberOfBackwardEdges++;
        }
        backwardEdges[i][j]={s,inputEdgeCost};
        numberOfBackwardEdges++;

    }
   // std::cout<<"edges added"<<std::endl;

    assert(numberOfBackwardEdges==numberOfEdges);

    for (size_t i=0;i<numberOfVertices;i++) {
        std::sort(forwardEdges[i].begin(),forwardEdges[i].end());
        std::sort(backwardEdges[i].begin(),backwardEdges[i].end());
    }


}

//LdpDirectedGraph::LdpDirectedGraph(const andres::graph::Digraph<>& inputGraph,const std::vector<double>& inputEdgeCosts){
//    numberOfVertices=inputGraph.numberOfVertices();
//    numberOfEdges=0;
//    std::vector<std::size_t> adjacencyForward(numberOfVertices);
//    std::vector<std::size_t> adjacencyBackward(numberOfVertices);

//    assert(inputGraph.numberOfEdges()==inputEdgeCosts.size());
//    for(size_t i=0;i<numberOfVertices;i++){
//        size_t neForward=inputGraph.numberOfEdgesFromVertex(i);
//        adjacencyForward[i]=neForward;
//        size_t neBackward=inputGraph.numberOfEdgesToVertex(i);
//        adjacencyBackward[i]=neBackward;
//    }

//    forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());

//    backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());
//    size_t numberOfBackwardEdges=0;

//    for(size_t i=0;i<numberOfVertices;i++){
//        size_t neForward=inputGraph.numberOfEdgesFromVertex(i);
//        for(size_t j=0;j<neForward;j++){
//            size_t e=inputGraph.edgeFromVertex(i,j);
//            size_t w=inputGraph.vertexFromVertex(i,j);
//            forwardEdges[i][j]={w,inputEdgeCosts[e]};
//            numberOfEdges++;

//        }
//        size_t neBackward=inputGraph.numberOfEdgesToVertex(i);
//        for(size_t j=0;j<neBackward;j++){
//            size_t w=inputGraph.vertexToVertex(i,j);
//            size_t e=inputGraph.edgeToVertex(i,j);
//            backwardEdges[i][j]={w,inputEdgeCosts[e]};
//            numberOfBackwardEdges++;
//        }
//    }
//    assert(numberOfEdges==numberOfBackwardEdges);

//    for (size_t i=0;i<numberOfVertices;i++) {
//        std::sort(forwardEdges[i].begin(),forwardEdges[i].end());
//        std::sort(backwardEdges[i].begin(),backwardEdges[i].end());
//    }
//    //Need to sort within edges?
//}


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
