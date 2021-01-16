#ifndef LDP_COMPLETE_STRUCTURE_HXX
#define LDP_COMPLETE_STRUCTURE_HXX



#include <stdexcept>
//#include "disjoint-paths/disjointPathsMethods.hxx"
#include "ldp_methods.hxx"
#include "ldp_functions.hxx"
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
#include<chrono>
#include"ldp_directed_graph.hxx"


namespace py = pybind11;
namespace LPMP {

template<class T = size_t>
struct CompleteStructure {
public:
    template<class PAR>
    CompleteStructure(PAR & configParameters)
    {
        //TODO min time in parameters
        constructorBegin=std::chrono::steady_clock::now();

        pVertexGroups=new VertexGroups<>(configParameters);
        VertexGroups<>& vg=*pVertexGroups;
        maxTime=vg.getMaxTime();
        assert(vg.getMinTime()==1);

        vertexShiftBack=vg.getVertexShiftBack();

        configParameters.getControlOutput()<<"max time "<<maxTime<<std::endl;

        verticesScore=std::vector<double>(vg.getMaxVertex()+1);

        configParameters.getControlOutput()<<"vg vertices "<<(vg.getMaxVertex()+1)<<std::endl;
        addEdgesFromFile(configParameters.getGraphFileName(),configParameters,vertexShiftBack);
        edgeFileName=configParameters.getGraphFileName();
        deleteVG=true;

}


    CompleteStructure(VertexGroups<>& vertexGroups):
        pVertexGroups(&vertexGroups)
{
        constructorBegin=std::chrono::steady_clock::now();
        VertexGroups<>& vg=*pVertexGroups;

        vertexShiftBack=vg.getVertexShiftBack();  //Carefully with this. There can be vertex shift even if we get zero from VG!

        maxTime=vg.getMaxTime();
        verticesScore=std::vector<double>(vg.getMaxVertex()+1);

        deleteVG=false;
        edgeFileName="";


}
    ~CompleteStructure(){
        if(deleteVG){
            delete pVertexGroups;
        }
    }


    template<class PAR>
    void addEdgesFromFile(const std::string& fileName,PAR& params,size_t minVertexToUse=0);

    template<class PAR>
    void addEdgesFromVectors(const py::array_t<size_t>& edges,const py::array_t<double>& costs,PAR& params,size_t minVertexToUse=0);

    void addEdgesFromVectorsAll(const py::array_t<size_t>& edges,const py::array_t<double>& costs);

    void setVerticesCostsVector(const std::vector<double>& costs);

    void setVerticesCosts(const py::array_t<double>& costs);

    std::vector<bool> getGraphEdgeLabels(const std::vector<std::vector<size_t>>& paths) const;

    std::vector<std::array<size_t,2>> getEdgeList() const;

    const VertexGroups<>& getVertexGroups(){
        return *pVertexGroups;
    }

    const std::chrono::steady_clock::time_point& getContructorBegin()const {
        return constructorBegin;
    }

    const std::string& getEdgeFileName()const {
       return edgeFileName;
    }


      //  andres::graph::Digraph<> completeGraph;
     //   std::vector<double> completeScore; //TODO make local variable
        std::vector<double> verticesScore;

        LdpDirectedGraph myCompleteGraph;
        size_t maxTime;


private:
    bool deleteVG;
    VertexGroups<>* pVertexGroups;
    std::chrono::steady_clock::time_point constructorBegin;
    std::string edgeFileName;
    size_t vertexShiftBack;

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


//    for (size_t i=0;i<dim1;i++) {
//        size_t v=edgeVector(i,0);
//        size_t w=edgeVector(i,1);
//        double edgeCost=costVector(i);

//        if(v>vg.getMaxVertex()||w>vg.getMaxVertex()){
//            throw std::invalid_argument("Input edges contain vertices out of the range.");
//        }

//    //    completeGraph.insertEdge(v,w);
//        completeScore.push_back(edgeCost);

//    }

    myCompleteGraph=LdpDirectedGraph(edgeVector,costVector);
}



template<class T>
template<class PAR>
inline void CompleteStructure<T>::addEdgesFromVectors(const py::array_t<size_t>& edges, const py::array_t<double>& costs, PAR& params, size_t minVertexToUse){
    char delim=',';
    VertexGroups<>& vg=*pVertexGroups;
    vertexShiftBack=std::max(minVertexToUse,vg.getVertexShiftBack());

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

    std::vector<std::array<size_t,2>> edgesToUse;
    std::vector<double> costsToUse;

    size_t i=0;


    for (;i<dim1;i++) {
        size_t v0=edgeVector(i,0);
        size_t w0=edgeVector(i,1);

        assert(w0>v0);

        if(v0>=minVertexToUse){
            size_t v=v0-minVertexToUse;
            size_t w=w0-minVertexToUse;

            if(w<=vg.getMaxVertex()){

                size_t l0=vg.getGroupIndex(v);
                size_t l1=vg.getGroupIndex(w);

                assert(l1>l0);

                if(l1-l0<=params.getMaxTimeGapComplete()){
                    double edgeCost=costVector(i);
                    edgesToUse.push_back({v,w});
                    costsToUse.push_back(edgeCost);

                }
            }
        }
    }


    EdgeVector ev(edgesToUse);
    InfoVector iv(costsToUse);
   // std::cout<<"before my complete graph"<<std::endl;
    myCompleteGraph=LdpDirectedGraph(ev,iv,vg.getMaxVertex()+1);

}






template<class T>
template<class PAR>
inline void CompleteStructure<T>::addEdgesFromFile(const std::string& fileName, PAR& params, size_t minVertexToUse){  //min VertexToUse cannot be always based on VG
    char delim=',';
    VertexGroups<>& vg=*pVertexGroups;
    edgeFileName=fileName;

    vertexShiftBack=std::max(minVertexToUse,vg.getVertexShiftBack());


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

            unsigned int v0 = std::stoul(strings[0]);
            if(v0>=minVertexToUse){
                size_t v=v0-minVertexToUse;
                double c = std::stod(strings[1]);
                assert(v<verticesScore.size());
                verticesScore[v]=c;
            }

        }

        params.getControlOutput()<<"Reading base edges from file. "<<std::endl;
        params.writeControlOutput();
        size_t maxLabel=0;

        std::vector<std::array<size_t,2>> listOfEdges;
        std::vector<double> completeScore;

        if(minVertexToUse>0){
            bool minVertexFound=false;
            while (std::getline(data, line) && !line.empty()) {


                strings = split(line, delim);

                unsigned int v0 = std::stoul(strings[0]);

                if(v0>=minVertexToUse){
                    minVertexFound=true;
                    unsigned int w0 = std::stoul(strings[1]);
                    size_t w=w0-minVertexToUse;
                    size_t v=v0-minVertexToUse;

                    assert(v<=vg.getMaxVertex());

                    if(w<=vg.getMaxVertex()){

                        size_t l0=vg.getGroupIndex(v);
                        size_t l1=vg.getGroupIndex(w);

                        //if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;

                        if(l1-l0<=params.getMaxTimeGapComplete()){
                            double score = std::stod(strings[2]);

                            completeScore.push_back(score);
                            listOfEdges.push_back({v,w});
                        }
                    }

                    break;
                }

            }
            assert(minVertexFound);

        }

        while (std::getline(data, line) && !line.empty()) {


            strings = split(line, delim);

            unsigned int v0 = std::stoul(strings[0]);

            unsigned int w0 = std::stoul(strings[1]);
            assert(w0>v0);
            assert(v0>=minVertexToUse);

            size_t v=v0-minVertexToUse;
            size_t w=w0-minVertexToUse;

            if(v>vg.getMaxVertex()) break;

            if(w>vg.getMaxVertex()) continue;

            size_t l0=vg.getGroupIndex(v);
            size_t l1=vg.getGroupIndex(w);

            //if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;

            if(l1-l0<=params.getMaxTimeGapComplete()){
                double score = std::stod(strings[2]);

                completeScore.push_back(score);
                listOfEdges.push_back({v,w});
            }

        }


        data.close();

        EdgeVector ev(listOfEdges);
        InfoVector iv(completeScore);
       // std::cout<<"before my complete graph"<<std::endl;
        myCompleteGraph=LdpDirectedGraph(ev,iv,vg.getMaxVertex()+1);
        //std::cout<<"after my complete graph "<<myCompleteGraph.getNumberOfVertices()<<std::endl;
    }
    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }
}





template<class T>
inline std::vector<std::array<size_t,2>> CompleteStructure<T>::getEdgeList() const{
    //std::vector<std::array<size_t,2>> edges(myCompleteGraph.getNumberOfEdges());
    std::vector<std::array<size_t,2>> edges;
    size_t edgeCounter=0;
    for (size_t i = 0; i < myCompleteGraph.getNumberOfVertices(); ++i) {
        for (auto iter=myCompleteGraph.forwardNeighborsBegin(i);iter!=myCompleteGraph.forwardNeighborsEnd(i);iter++) {
            edges.push_back({i,iter->first});
        }
    }
    assert(edges.size()==myCompleteGraph.getNumberOfEdges());

    return edges;
}



}//End of namespace

#endif // LDP_COMPLETE_STRUCTURE_HXX
