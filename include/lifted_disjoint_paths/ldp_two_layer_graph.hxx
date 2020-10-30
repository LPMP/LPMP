#ifndef LDP_TWO_LAYER_GRAPH_HXX
#define LDP_TWO_LAYER_GRAPH_HXX

#include "two_dimensional_variable_array.hxx"
namespace LPMP {

class LdpTwoLayerGraph{
public:
    LdpTwoLayerGraph(){}
    LdpTwoLayerGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts){

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
            numberOfOutputs=std::max(j+1,adjacencyForward.size());
            adjacencyForward.resize(numberOfInputs);
            adjacencyForward[i]++;
        }
        numberOfInputs=adjacencyForward.size();


        forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());


         std::fill(adjacencyForward.begin(), adjacencyForward.end(), 0);

        for(size_t i=0;i<edges.size();i++){
            size_t v=edges[i][0];
            size_t w=edges[i][1];
            forwardEdges[v][adjacencyForward[v]]={w,inputEdgeCosts[i]};
            std::cout<<"add edge to layer graph "<<v<<", "<<w<<std::endl;
            adjacencyForward[v]++;
        }

//        for (size_t i=0;i<numberOfInputs;i++) {
//            std::sort(forwardEdges[i].begin(),forwardEdges[i].end());

//        }
        //Need to sort within edges?
    }


    void setForwardEdgeCost(size_t vertex,size_t neighborIndex,double value) {
        forwardEdges[vertex][neighborIndex].second=value;
    }



    double getForwardEdgeCost(size_t vertex,size_t neighborIndex) const{
        return forwardEdges[vertex][neighborIndex].second;
    }

    size_t getForwardEdgeVertex(size_t vertex,size_t neighborIndex) const{
        return forwardEdges[vertex][neighborIndex].first;
    }



    const std::pair<size_t,double> * forwardNeighborsBegin(size_t i)const {
        return forwardEdges[i].begin();
    }

    const std::pair<size_t,double> * forwardNeighborsEnd(size_t i)const {
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
