#ifndef LDP_DIRECTED_GRAPH_HXX
#define LDP_DIRECTED_GRAPH_HXX

#include "two_dimensional_variable_array.hxx"
#include "andres/graph/digraph.hxx"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace LPMP {

class LdpDirectedGraph{
public:
    LdpDirectedGraph(){}


  //  LdpDirectedGraph(const andres::graph::Digraph<>& inputGraph,const std::vector<double>& inputEdgeCosts);
//    LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts,double inCost,double outCost);

//    LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts);

    template<class EDGES, class COSTS>
    LdpDirectedGraph(const EDGES& edges,const COSTS& inputEdgeCosts,double inCost,double outCost,size_t _numberOfVertices=0);

    template<class EDGES, class COSTS>
    LdpDirectedGraph(const EDGES& edges,const COSTS& inputEdgeCosts,size_t _numberOfVertices=0);

    LdpDirectedGraph(const LdpDirectedGraph& inputGraph,double inputEdgeCost,double outputEdgeCost);

    template<class EDGES, class COSTS>
    void initFromEdgesAndCosts(const EDGES& edges,const COSTS& inputEdgeCosts,double inCost,double outCost,size_t _numberOfVertices);


    double getForwardEdgeCost(size_t vertex,size_t neighborIndex) const{
        assert(vertex<numberOfVertices);
        assert(neighborIndex<forwardEdges[vertex].size());
        return forwardEdges[vertex][neighborIndex].second;
    }


    double getBackwardEdgeCost(size_t vertex,size_t neighborIndex) const{
        assert(vertex<numberOfVertices);
        assert(neighborIndex<backwardEdges[vertex].size());
        return backwardEdges[vertex][neighborIndex].second;
    }

    size_t getForwardEdgeVertex(size_t vertex,size_t neighborIndex) const{
        assert(vertex<numberOfVertices);
        assert(neighborIndex<forwardEdges[vertex].size());
        return forwardEdges[vertex][neighborIndex].first;
    }


    size_t getBackwardEdgeVertex(size_t vertex,size_t neighborIndex) const{
         assert(vertex<numberOfVertices);
         assert(neighborIndex<backwardEdges[vertex].size());
        return backwardEdges[vertex][neighborIndex].first;
    }

    const std::pair<size_t,double> * forwardNeighborsBegin(size_t i)const {
         assert(i<numberOfVertices);
        return forwardEdges[i].begin();
    }

    const std::pair<size_t,double> * forwardNeighborsEnd(size_t i)const {
         assert(i<numberOfVertices);
        return forwardEdges[i].end();
    }

    const std::pair<size_t,double> * backwardNeighborsBegin(size_t i)const {
         assert(i<numberOfVertices);
        return backwardEdges[i].begin();
    }

    const std::pair<size_t,double> * backwardNeighborsEnd(size_t i)const {
         assert(i<numberOfVertices);
        return backwardEdges[i].end();
    }

    std::pair<size_t,double> * forwardNeighborsBegin(size_t i) {
        if(i>=numberOfVertices) std::cout<<"i is "<<i<<", number of vertices "<<numberOfVertices<<std::endl;
         assert(i<numberOfVertices);
        return forwardEdges[i].begin();
    }

     std::pair<size_t,double> * forwardNeighborsEnd(size_t i) {
          assert(i<numberOfVertices);
        return forwardEdges[i].end();
    }

     std::pair<size_t,double> * backwardNeighborsBegin(size_t i) {
         assert(i<numberOfVertices);
        return backwardEdges[i].begin();
    }

     std::pair<size_t,double> * backwardNeighborsEnd(size_t i) {
        assert(i<numberOfVertices);
        return backwardEdges[i].end();
    }

    const size_t & getNumberOfVertices()const{
        return numberOfVertices;
    }

    const size_t & getNumberOfEdges()const{
        return numberOfEdges;
    }

    size_t getNumberOfEdgesFromVertex(const size_t& i)const{
        assert(i<numberOfVertices);
        return forwardEdges[i].size();
    }

    size_t getNumberOfEdgesToVertex(const size_t& i)const{
        assert(i<numberOfVertices);
        return backwardEdges[i].size();
    }


    const two_dim_variable_array<std::pair<size_t,double>>& getForwardEdges()const{
        return forwardEdges;
    }

    const two_dim_variable_array<std::pair<size_t,double>>& getBackwardEdges()const{
        return backwardEdges;
    }


private:
    two_dim_variable_array<std::pair<size_t,double>> forwardEdges;
    two_dim_variable_array<std::pair<size_t,double>> backwardEdges;
    size_t numberOfVertices;
    size_t numberOfEdges;
    //two_dim_variable_array<double> forwardCost;
    //two_dim_variable_array<double> backwardCost;

};


template<class EDGES, class COSTS>
LdpDirectedGraph::LdpDirectedGraph(const EDGES& edges,const COSTS& inputEdgeCosts,size_t _numberOfVertices){
    //TODO the same as in the previous constructor but addd s and t nodes and edges with given cost
    initFromEdgesAndCosts(edges,inputEdgeCosts,std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),_numberOfVertices);
   // std::cout<<"small constructor vertices "<<numberOfVertices<<std::endl;
}



template<class EDGES, class COSTS>
LdpDirectedGraph::LdpDirectedGraph(const EDGES& edges, const COSTS& inputEdgeCosts, double inCost, double outCost, size_t _numberOfVertices){
    initFromEdgesAndCosts(edges,inputEdgeCosts,inCost,outCost,_numberOfVertices);
}

template<class EDGES, class COSTS>
void LdpDirectedGraph::initFromEdgesAndCosts(const EDGES& edges,const COSTS& inputEdgeCosts,double inCost,double outCost,size_t _numberOfVertices){

    bool addST=inCost<std::numeric_limits<double>::max();

    numberOfEdges=0;

    std::vector<std::size_t> adjacencyForward;
    std::vector<std::size_t> adjacencyBackward;


    if(_numberOfVertices>0){
        adjacencyForward=std::vector<size_t>(_numberOfVertices);
        adjacencyBackward=std::vector<size_t>(_numberOfVertices);
    }


    assert(edges.shape(0)==inputEdgeCosts.shape(0));

    // first determine size for adjacency_list
    for(size_t i=0;i<edges.shape(0);i++){
        size_t v=edges(i,0);
        size_t w=edges(i,1);
        size_t size=std::max({v+1,w+1,adjacencyForward.size()});
        assert(_numberOfVertices==0||size==adjacencyForward.size());
        adjacencyBackward.resize(size);
        adjacencyForward.resize(size);
        adjacencyForward[v]++;
        adjacencyBackward[w]++;
//        if(w==0){
//            std::cout<<"back to w "<<adjacencyBackward[w]<<std::endl;
//        }
    }
    numberOfVertices=adjacencyForward.size();

    if(addST){
        for(size_t i=0;i<numberOfVertices;i++){
            adjacencyForward[i]++;
            adjacencyBackward[i]++;
        }

        adjacencyForward.push_back(numberOfVertices); //for s
        adjacencyBackward.push_back(0);
        adjacencyForward.push_back(0);  //for t
        adjacencyBackward.push_back(numberOfVertices);

    }
    forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());

    backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());



     std::fill(adjacencyForward.begin(), adjacencyForward.end(), 0);
     std::fill(adjacencyBackward.begin(), adjacencyBackward.end(), 0);


    for(size_t i=0;i<edges.shape(0);i++){
        size_t v=edges(i,0);
        size_t w=edges(i,1);
        forwardEdges[v][adjacencyForward[v]]={w,inputEdgeCosts(i)};
        adjacencyForward[v]++;
        numberOfEdges++;
        backwardEdges[w][adjacencyBackward[w]]={v,inputEdgeCosts(i)};
        adjacencyBackward[w]++;
    }
    if(addST){
        size_t s=numberOfVertices;
        size_t t=numberOfVertices+1;

        for(size_t i=0;i<numberOfVertices;i++){
            forwardEdges[s][i]={i,inCost};
            backwardEdges[t][i]={i,outCost};
            forwardEdges[i][adjacencyForward[i]]={t,outCost};
            backwardEdges[i][adjacencyBackward[i]]={s,inCost};
            numberOfEdges=numberOfEdges+2;

        }

        if(numberOfEdges!=(edges.shape(0)+2*numberOfVertices)){
            std::cout<<"number of edges "<<numberOfEdges<<std::endl;
            std::cout<<"expected "<<(edges.shape(0)+2*numberOfVertices)<<std::endl;
        }
        assert(numberOfEdges==(edges.shape(0)+2*numberOfVertices));
        numberOfVertices+=2;
    }
    else{
        assert(numberOfEdges==edges.shape(0));
    }

    for (size_t i=0;i<numberOfVertices;i++) {
        std::sort(forwardEdges[i].begin(),forwardEdges[i].end());
        std::sort(backwardEdges[i].begin(),backwardEdges[i].end());
    }

    //std::cout<<"ldp graph, vertices: "<<numberOfVertices<<std::endl;

    //Need to sort within edges?
}


}

#endif // LDP_DIRECTED_GRAPH_HXX
