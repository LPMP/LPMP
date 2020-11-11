#ifndef LDP_COMPLETE_STRUCTURE_HXX
#define LDP_COMPLETE_STRUCTURE_HXX



#include <stdexcept>
//#include "disjoint-paths/disjointPathsMethods.hxx"
#include "ldp_methods.hxx"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <iterator>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include "lifted_disjoint_paths/ldp_vertex_groups.hxx"
#include "ldp_file_processing_methods.hxx"
#include <array>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>


namespace py = pybind11;
namespace LPMP {

template<class T = size_t>
struct CompleteStructure {
public:
    template<class PAR>
    CompleteStructure(PAR & configParameters)
    {

        pVertexGroups=new VertexGroups<>(configParameters);
        VertexGroups<>& vg=*pVertexGroups;
        maxTime=vg.getMaxTime();
        assert(vg.getMinTime()==1);

        configParameters.getControlOutput()<<"max time "<<maxTime<<std::endl;
        completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);
        verticesScore=std::vector<double>(vg.getMaxVertex()+1);

        configParameters.getControlOutput()<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
        addEdgesFromFile(configParameters.getGraphFileName(),configParameters);
        deleteVG=true;

}


    CompleteStructure(VertexGroups<>& vertexGroups):
        pVertexGroups(&vertexGroups)
{
        VertexGroups<>& vg=*pVertexGroups;

        maxTime=vg.getMaxTime();
        verticesScore=std::vector<double>(vg.getMaxVertex()+1);
       // std::cout<<"max time "<<maxTime<<std::endl;
        completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);
       // std::cout<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
        deleteVG=false;


}
    ~CompleteStructure(){
        if(deleteVG){
            delete pVertexGroups;
        }
    }


    template<class PAR>
    void addEdgesFromFile(const std::string& fileName,PAR& params);
    void addEdgesFromMatrix(size_t time1,size_t time2,const py::array_t<double> inputMatrix);
    template<class PAR>
    void addEdgesFromVectors(const py::array_t<size_t>& edges,const py::array_t<double>& costs,PAR& params);
    void addEdgesFromVectorsAll(const py::array_t<size_t>& edges,const py::array_t<double>& costs);
    void setVerticesCostsVector(const std::vector<double>& costs);
     void setVerticesCosts(const py::array_t<double>& costs);
    std::vector<bool> getGraphEdgeLabels(const std::vector<std::vector<size_t>>& paths) const;
    std::vector<std::array<size_t,2>> getEdgeList() const;
    const VertexGroups<>& getVertexGroups(){
        return *pVertexGroups;
    }


        andres::graph::Digraph<> completeGraph;
        std::vector<double> completeScore;
        std::vector<double> verticesScore;

        size_t maxTime;


private:
    bool deleteVG;
    VertexGroups<>* pVertexGroups;

};


template<class T>
inline void CompleteStructure<T>::setVerticesCostsVector(const std::vector<double>& costs){
    assert(costs.size()==verticesScore.size());
    verticesScore.assign(costs.begin(),costs.end());
}


template<class T>
inline void CompleteStructure<T>::setVerticesCosts(const py::array_t<double>& costs){
    const auto costVector=costs.unchecked<1>();
    const size_t dimOfCosts=costVector.shape(0);
    assert(dimOfCosts==verticesScore.size());
    for (size_t i=0;i<dimOfCosts;i++){
        verticesScore[i]=costVector(i);
    }
}

template<class T>
inline void CompleteStructure<T>::addEdgesFromVectorsAll(const py::array_t<size_t>& edges,const py::array_t<double>& costs){
    char delim=',';
    VertexGroups<>& vg=*pVertexGroups;

    const auto edgeVector=edges.unchecked<2>();
    const std::size_t dim1=edgeVector.shape(0);
    const std::size_t dim2=edgeVector.shape(1);
    const auto costVector=costs.unchecked<1>();
    const size_t dimOfCosts=costVector.shape(0);


    if(dim2!=2){
        std::string message="Wrong dimension of edge array, second dimension 2 expected";
        throw std::invalid_argument(message);
    }
    if(dim1!=dimOfCosts){
        std::string message="Dimension of edge array and edge costs do not match.";
        throw std::invalid_argument(message);
    }


    for (size_t i=0;i<dim1;i++) {
        size_t v=edgeVector(i,0);
        size_t w=edgeVector(i,1);
        double edgeCost=costVector(i);

        if(v>vg.getMaxVertex()||w>vg.getMaxVertex()){
            throw std::invalid_argument("Input edges contain vertices out of the range.");
        }

        completeGraph.insertEdge(v,w);
        completeScore.push_back(edgeCost);

    }
}



template<class T>
template<class PAR>
inline void CompleteStructure<T>::addEdgesFromVectors(const py::array_t<size_t>& edges,const py::array_t<double>& costs,PAR& params){
    char delim=',';
    VertexGroups<>& vg=*pVertexGroups;

    const auto edgeVector=edges.unchecked<2>();
    const std::size_t dim1=edgeVector.shape(0);
    const std::size_t dim2=edgeVector.shape(1);
    const auto costVector=costs.unchecked<1>();
    const size_t dimOfCosts=costVector.shape(0);


    if(dim2!=2){
        std::string message="Wrong dimension of edge array, second dimension 2 expected";
        throw std::invalid_argument(message);
    }
    if(dim1!=dimOfCosts){
        std::string message="Dimension of edge array and edge costs do not match.";
        throw std::invalid_argument(message);
    }

    params.getControlOutput()<<"Reading base edges from vector. "<<std::endl;
    params.writeControlOutput();
    for (size_t i=0;i<dim1;i++) {
        size_t v=edgeVector(i,0);
        size_t w=edgeVector(i,1);
        double edgeCost=costVector(i);


        if(v>vg.getMaxVertex()) break;

        if(w>vg.getMaxVertex()) continue;
        size_t l0=vg.getGroupIndex(v);
        size_t l1=vg.getGroupIndex(w);

        if(l1-l0<=params.getMaxTimeGapComplete()){
            completeGraph.insertEdge(v,w);
            completeScore.push_back(edgeCost);
        }
    }
}




template<class T>
inline void CompleteStructure<T>::addEdgesFromMatrix(size_t time1,size_t time2,const py::array_t<double> inputMatrix){

    VertexGroups<>& vg=*pVertexGroups;
    const auto matrix=inputMatrix.unchecked<2>();
    const std::size_t dim1=matrix.shape(0);
    const std::size_t dim2=matrix.shape(1);

    size_t transformIndex1=vg.getGroupVertices(time1)[0];
    size_t numberOfVertices1=vg.getGroupVertices(time1).size();
    size_t transformIndex2=vg.getGroupVertices(time2)[0];
    size_t numberOfVertices2=vg.getGroupVertices(time2).size();

    if(dim1!=numberOfVertices1||dim2!=numberOfVertices2){
        std::string message="Dimension mismatch, expected dimensions: "+std::to_string(numberOfVertices1)+", "+std::to_string(numberOfVertices2)+", got: "+std::to_string(dim1)+", "+std::to_string(dim2);
        throw std::invalid_argument(message);
    }
    for(std::size_t i=0; i<dim1; ++i) {
        size_t vertex1=i+transformIndex1;
        for(std::size_t j=0; j<dim2; ++j) {
            double score=matrix(i,j);
            if(!isinf(score)){
                size_t vertex2=j+transformIndex2;
                completeGraph.insertEdge(vertex1,vertex2);
                completeScore.push_back(score);
            }
        }
    }

}




template<class T>
template<class PAR>
inline void CompleteStructure<T>::addEdgesFromFile(const std::string& fileName,PAR& params){
    char delim=',';
    VertexGroups<>& vg=*pVertexGroups;

    std::string line;
    std::ifstream data;
    try{
        data.open(fileName);
        if(!data){
            throw std::system_error(errno, std::system_category(), "failed to open graph file "+fileName);
        }

        std::getline(data, line);
        double objValue=0;


        params.getControlOutput()<<"Read big graph" << std::endl;
        std::vector<std::string> strings;

        params.getControlOutput()<<"Reading vertices from file. "<<std::endl;
        params.writeControlOutput();
        //Vertices that are not found have score=0. Appearance and disappearance cost are read here.
        while (std::getline(data, line) && !line.empty()) {
            strings = split(line, delim);

            unsigned int v = std::stoul(strings[0]);

            double c = std::stod(strings[1]);
            assert(v<verticesScore.size());
            verticesScore[v]=c;

        }

        params.getControlOutput()<<"Reading base edges from file. "<<std::endl;
        params.writeControlOutput();
        size_t maxLabel=0;
        while (std::getline(data, line) && !line.empty()) {


            strings = split(line, delim);

            unsigned int v = std::stoul(strings[0]);

            unsigned int w = std::stoul(strings[1]);

            if(v>vg.getMaxVertex()) break;

            if(w>vg.getMaxVertex()) continue;

            size_t l0=vg.getGroupIndex(v);
            size_t l1=vg.getGroupIndex(w);

            //if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;

            if(l1-l0<=params.getMaxTimeGapComplete()){
                double score = std::stod(strings[2]);
                completeGraph.insertEdge(v,w);
                completeScore.push_back(score);
            }

        }


        //objValue+=maxLabel*parameters.getInputCost()+maxLabel*parameters.getOutputCost();

        data.close();
    }
    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }
}




template<class T>
inline std::vector<bool> CompleteStructure<T>::getGraphEdgeLabels(const std::vector<std::vector<size_t>>& paths) const{
    std::vector<bool> labels=getEdgeLabels(completeGraph,paths);
    //std::vector<bool> labels;
    return labels;
}

template<class T>
inline std::vector<std::array<size_t,2>> CompleteStructure<T>::getEdgeList() const{
    std::vector<std::array<size_t,2>> edges(completeGraph.numberOfEdges());
    for (size_t e = 0; e < completeGraph.numberOfEdges(); ++e) {
        size_t v=completeGraph.vertexOfEdge(e,0);
        size_t w=completeGraph.vertexOfEdge(e,1);
        edges[e]={v,w};
    }
    return edges;
}



}//End of namespace

#endif // LDP_COMPLETE_STRUCTURE_HXX
