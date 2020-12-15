#include "lifted_disjoint_paths/ldp_two_layer_graph.hxx"

namespace LPMP {

LdpTwoLayerGraph::LdpTwoLayerGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts){

    std::vector<std::size_t> adjacencyForward;
    std::vector<std::size_t> adjacencyBackward;

    assert(edges.size()==inputEdgeCosts.size());


    // first determine size for adjacency_list
    numberOfOutputs=0;
    numberOfInputs=0;
    for(const std::array<size_t,2>& e:edges){
        size_t i=e[0];
        size_t j=e[1];
        numberOfInputs=std::max(i+1,adjacencyForward.size());
        numberOfOutputs=std::max(j+1,adjacencyBackward.size());
        adjacencyForward.resize(numberOfInputs);
        adjacencyForward[i]++;
        adjacencyBackward.resize(numberOfOutputs);
        adjacencyBackward[j]++;
    }
    numberOfInputs=adjacencyForward.size();
    numberOfOutputs=adjacencyBackward.size();


    forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());
    backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());


     std::fill(adjacencyForward.begin(), adjacencyForward.end(), 0);
     std::fill(adjacencyBackward.begin(), adjacencyBackward.end(), 0);

    for(size_t i=0;i<edges.size();i++){
        size_t v=edges[i][0];
        size_t w=edges[i][1];
        forwardEdges[v][adjacencyForward[v]]={w,inputEdgeCosts[i],adjacencyBackward[w]};
        backwardEdges[w][adjacencyBackward[w]]={v,inputEdgeCosts[i],adjacencyForward[v]};
        //std::cout<<"add edge to layer graph "<<v<<", "<<w<<std::endl;
        adjacencyForward[v]++;
        adjacencyBackward[w]++;
    }

//        for (size_t i=0;i<numberOfInputs;i++) {
//            std::sort(forwardEdges[i].begin(),forwardEdges[i].end());

//        }
    //Need to sort within edges?
}


}
