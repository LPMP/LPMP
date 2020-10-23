#ifndef LDP_DIRECTED_GRAPH_HXX
#define LDP_DIRECTED_GRAPH_HXX

#include "two_dimensional_variable_array.hxx"
#include "andres/graph/digraph.hxx"

namespace LPMP {

class LdpDirectedGraph{
public:
    LdpDirectedGraph(const andres::graph::Digraph<>& inputGraph,const std::vector<double>& inputEdgeCosts){
        size_t numberOfVertices=inputGraph.numberOfVertices();
        std::vector<std::size_t> adjacencyForward(numberOfVertices);
        std::vector<std::size_t> adjacencyBackward(numberOfVertices);

        for(size_t i=0;i<numberOfVertices;i++){
            size_t neForward=inputGraph.numberOfEdgesFromVertex(i);
            adjacencyForward[i]=neForward;
            size_t neBackward=inputGraph.numberOfEdgesToVertex(i);
            adjacencyBackward[i]=neBackward;
        }

        forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());
        forwardCost.resize(adjacencyForward.begin(),adjacencyForward.end());
        backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());
        backwardCost.resize(adjacencyBackward.begin(),adjacencyBackward.end());

        for(size_t i=0;i<numberOfVertices;i++){
            size_t neForward=inputGraph.numberOfEdgesFromVertex(i);
            for(size_t j=0;j<neForward;j++){
                size_t e=inputGraph.edgeFromVertex(i,j);
                size_t w=inputGraph.vertexFromVertex(i,j);
                forwardEdges[i][j]=w;
                forwardCost[i][j]=inputEdgeCosts[e];
            }
            size_t neBackward=inputGraph.numberOfEdgesToVertex(i);
            for(size_t j=0;j<neBackward;j++){
                size_t w=inputGraph.vertexToVertex(i,j);
                size_t e=inputGraph.edgeToVertex(i,j);
                backwardEdges[i][j]=w;
                backwardCost[i][j]=inputEdgeCosts[e];
            }
        }
        //Need to sort within edges?
    }

    LdpDirectedGraph(const std::vector<std::array<size_t,2>>& edges,const std::vector<double>& inputEdgeCosts){

        std::vector<std::size_t> adjacencyForward;
        std::vector<std::size_t> adjacencyBackward;


        // first determine size for adjacency_list
        for(const std::array<size_t,2>& e:edges){
            size_t i=e[0];
            size_t j=e[j];
            size_t size=std::max({i+1,j+1,adjacencyForward.size()});
            adjacencyBackward.resize(size);
            adjacencyForward.resize(size);
            adjacencyForward[i]++;
            adjacencyBackward[j]++;
        }

        forwardEdges.resize(adjacencyForward.begin(),adjacencyForward.end());
        forwardCost.resize(adjacencyForward.begin(),adjacencyForward.end());
        backwardEdges.resize(adjacencyBackward.begin(),adjacencyBackward.end());
        backwardCost.resize(adjacencyBackward.begin(),adjacencyBackward.end());


         std::fill(adjacencyForward.begin(), adjacencyForward.end(), 0);
         std::fill(adjacencyBackward.begin(), adjacencyBackward.end(), 0);


        size_t numberOfVertices=adjacencyForward.size();
        for(size_t i=0;i<edges.size();i++){
            size_t v=edges[i][0];
            size_t w=edges[i][1];
            forwardEdges[v][adjacencyForward[v]]=w;
            forwardCost[v][adjacencyForward[v]]=inputEdgeCosts[i];
            adjacencyForward[v]++;
            backwardEdges[w][adjacencyBackward[w]]=v;
            backwardCost[w][adjacencyBackward[w]]=inputEdgeCosts[i];
            adjacencyBackward[w]++;
        }

        //Need to sort within edges?
    }

    LdpDirectedGraph(std::vector<std::array<size_t,2>> edges, std::vector<double> inputEdgeCosts,double inCost,double outCost){
        //TODO the same as in the previous constructor but addd s and t nodes and edges with given cost
    }



    double getForwardEdgeCost(size_t vertex,size_t neighborIndex) const{
        return forwardCost[vertex][neighborIndex];
    }


    double getBackwardEdgeCost(size_t vertex,size_t neighborIndex) const{
        return backwardCost[vertex][neighborIndex];
    }

    const size_t * forwardNeighborsBegin(size_t i)const {
        return forwardEdges[i].begin();
    }

    const size_t * forwardNeighborsEnd(size_t i)const {
        return forwardEdges[i].end();
    }

    const size_t * backwardNeighborsBegin(size_t i)const {
        return backwardEdges[i].begin();
    }

    const size_t * backwardNeighborsEnd(size_t i)const {
        return backwardEdges[i].end();
    }

    const double * forwardCostBegin(size_t i)const {
        return forwardCost[i].begin();
    }

    const double * forwardCostEnd(size_t i)const {
        return forwardCost[i].end();
    }

    const double * backwardCostBegin(size_t i)const {
        return backwardCost[i].begin();
    }

    const double * backwardCostEnd(size_t i)const {
        return backwardCost[i].end();
    }

    //iterator for(double *it=forwardCost[i].begin();it!=forwardCost[i].end();it++)




private:
    two_dim_variable_array<size_t> forwardEdges;
    two_dim_variable_array<size_t> backwardEdges;
    two_dim_variable_array<double> forwardCost;
    two_dim_variable_array<double> backwardCost;

};

}

#endif // LDP_DIRECTED_GRAPH_HXX
