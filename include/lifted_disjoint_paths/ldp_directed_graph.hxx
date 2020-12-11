#ifndef LDP_DIRECTED_GRAPH_HXX
#define LDP_DIRECTED_GRAPH_HXX

#include "two_dimensional_variable_array.hxx"
#include "andres/graph/digraph.hxx"

namespace LPMP {

class LdpDirectedGraph{
public:
    LdpDirectedGraph(){}


    LdpDirectedGraph(const andres::graph::Digraph<>& inputGraph,const std::vector<double>& inputEdgeCosts);
    LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts,double inCost,double outCost);

    LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts);


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

    size_t getNumberOfEdgesFromVertex(const size_t& i)const{
        assert(i<numberOfVertices);
        return forwardEdges[i].size();
    }




private:
    two_dim_variable_array<std::pair<size_t,double>> forwardEdges;
    two_dim_variable_array<std::pair<size_t,double>> backwardEdges;
    size_t numberOfVertices;
    //two_dim_variable_array<double> forwardCost;
    //two_dim_variable_array<double> backwardCost;

};




}

#endif // LDP_DIRECTED_GRAPH_HXX
