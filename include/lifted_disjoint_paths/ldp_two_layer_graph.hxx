#ifndef LDP_TWO_LAYER_GRAPH_HXX
#define LDP_TWO_LAYER_GRAPH_HXX

#include "two_dimensional_variable_array.hxx"
namespace LPMP {

class LdpTwoLayerGraph{
public:
    LdpTwoLayerGraph(){}
    LdpTwoLayerGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts);


    void setForwardEdgeCost(size_t vertex,size_t neighborIndex,double value) {
        //todo asserts
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        forwardEdges[vertex][neighborIndex].second=value;
    }


    void updateForwardEdgeCost(size_t vertex,size_t neighborIndex,double value) {
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        forwardEdges[vertex][neighborIndex].second+=value;
    }

    double getForwardEdgeCost(size_t vertex,size_t neighborIndex) const{
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        return forwardEdges[vertex][neighborIndex].second;
    }

    size_t getForwardEdgeVertex(size_t vertex,size_t neighborIndex) const{
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        return forwardEdges[vertex][neighborIndex].first;
    }



    const std::pair<size_t,double> * forwardNeighborsBegin(size_t i)const {
        assert(i<numberOfInputs);
        return forwardEdges[i].begin();
    }

    const std::pair<size_t,double> * forwardNeighborsEnd(size_t i)const {
        assert(i<numberOfInputs);
        return forwardEdges[i].end();
    }



    std::pair<size_t,double> * forwardNeighborsBegin(size_t i){
        assert(i<numberOfInputs);
        return forwardEdges[i].begin();
    }

    std::pair<size_t,double> * forwardNeighborsEnd(size_t i){
        assert(i<numberOfInputs);
        return forwardEdges[i].end();
    }


    const size_t & getNumberOfInputs()const{
        return numberOfInputs;
    }


private:

    two_dim_variable_array<std::pair<size_t,double>> forwardEdges;
    size_t numberOfInputs;
    size_t numberOfOutputs;

};

}



#endif // LDP_TWO_LAYER_GRAPH_HXX
