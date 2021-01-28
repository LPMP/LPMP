#ifndef LDP_METHODS_HXX
#define LDP_METHODS_HXX

/*
 * disjointPathsMethods.hxx
 *
 *  Created on: Jul 3, 2020
 *      Author: fuksova
 */


#include <stdexcept>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <stack>
#include <unordered_set>
#include <iterator>
#include <unordered_map>
#include <string>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
//#include "levinkov/timer.hxx"
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <list>
//#include "disjoint-paths/disjointParams.hxx"
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include "ldp_vertex_groups.hxx"

namespace py = pybind11;
namespace LPMP {



template<class T,class PAR>
    std::vector<std::unordered_set<size_t>> initReachableSet(T & graph,PAR& parameters,VertexGroups<size_t>* vg=0){


//        levinkov::Timer tfw;
//                        tfw.start();


        parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
                const size_t n=graph.numberOfVertices();

                // todo: use some matrix structure
                std::vector<std::unordered_set<size_t>> desc(n);
                std::vector<std::vector<std::bitset<10000>>> descBit(n);
                size_t columns=n/10000;
                size_t mod=n%10000;
                if(mod>0) columns++;


                for (int i = 0; i < n; ++i) {
                        descBit[i]=std::vector<std::bitset<10000>>(columns);
                }

                for (size_t v = 0; v < n; ++v) {

                        descBit[v][v/10000][v%10000]=1; //make this reflexive
                        size_t edges=graph.numberOfEdgesFromVertex(v);
                        for (int j = 0; j < edges; ++j) {
                                size_t w=graph.vertexFromVertex(v,j);
                                descBit[v][w/10000][w%10000]=1;


                        }
                }



                if(vg==0){
                        for (int k1 = 0; k1 <columns; ++k1) {
                                for (int k2 = 0; k2 < 10000; ++k2) {
                                        if(k1*10000+k2<n){
                                                for (int i = 0; i < n; ++i) {
                                                        if(descBit[i][k1][k2]){
                                                                for (int j = 0; j < columns; ++j) {
                                                                        descBit[i][j]|=descBit[k1*10000+k2][j];
                                                                }

                                                        }
                                                }
                                        }
                                        else{
                                                break;
                                        }
                                }

                        }
                }
                else{
                        for (int k1 = 0; k1 <columns; ++k1) {
                                for (int k2 = 0; k2 < 10000; ++k2) {
                                        if(k1*10000+k2<n){
                                                size_t maxTime=vg->getGroupIndex(k1*10000+k2);
                                                for (int t = 0; t < maxTime; ++t) {//TODO use time gap
                                                        const std::vector<size_t>& vertices=vg->getGroupVertices(t);
                                                        for (size_t i:vertices) {
                                                                if(descBit[i][k1][k2]){
                                                                        for (int j = 0; j < columns; ++j) {
                                                                                descBit[i][j]|=descBit[k1*10000+k2][j];
                                                                        }

                                                                }
                                                        }
                                                }
                                        }
                                        else{
                                                break;
                                        }
                                }

                        }

                }

                for (int i = 0; i < n; ++i) {
                        for (int k1 = 0; k1 <columns; ++k1) {
                                for (int k2 = 0; k2 < 10000; ++k2) {
                                        if(k1*10000+k2<n&&descBit[i][k1][k2]){
                                                desc[i].insert(k1*10000+k2);
                                        }
                                }
                        }

                }

             //   tfw.stop();


        //parameters.getControlOutput()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
                return desc;

}







//    std::vector<std::unordered_set<size_t>> initReachableSet(T & graph,PAR& parameters,VertexGroups<size_t>* vg=0){


////        levinkov::Timer tfw;
////                        tfw.start();


//        parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
//                const size_t n=graph.numberOfVertices();

//                // todo: use some matrix structure
//                std::vector<std::unordered_set<size_t>> desc(n);
//                std::vector<std::vector<std::bitset<10000>>> descBit(n);
//                size_t columns=n/10000;
//                size_t mod=n%10000;
//                if(mod>0) columns++;


//                for (int i = 0; i < n; ++i) {
//                        descBit[i]=std::vector<std::bitset<10000>>(columns);
//                }

//                for (size_t v = 0; v < n; ++v) {

//                        descBit[v][v/10000][v%10000]=1; //make this reflexive
//                        size_t edges=graph.numberOfEdgesFromVertex(v);
//                        for (int j = 0; j < edges; ++j) {
//                                size_t w=graph.vertexFromVertex(v,j);
//                                descBit[v][w/10000][w%10000]=1;


//                        }
//                }



//                if(vg==0){
//                        for (int k1 = 0; k1 <columns; ++k1) {
//                                for (int k2 = 0; k2 < 10000; ++k2) {
//                                        if(k1*10000+k2<n){
//                                                for (int i = 0; i < n; ++i) {
//                                                        if(descBit[i][k1][k2]){
//                                                                for (int j = 0; j < columns; ++j) {
//                                                                        descBit[i][j]|=descBit[k1*10000+k2][j];
//                                                                }

//                                                        }
//                                                }
//                                        }
//                                        else{
//                                                break;
//                                        }
//                                }

//                        }
//                }
//                else{
//                        for (int k1 = 0; k1 <columns; ++k1) {
//                                for (int k2 = 0; k2 < 10000; ++k2) {
//                                        if(k1*10000+k2<n){
//                                                size_t maxTime=vg->getGroupIndex(k1*10000+k2);
//                                                for (int t = 0; t < maxTime; ++t) {//TODO use time gap
//                                                        const std::vector<size_t>& vertices=vg->getGroupVertices(t);
//                                                        for (size_t i:vertices) {
//                                                                if(descBit[i][k1][k2]){
//                                                                        for (int j = 0; j < columns; ++j) {
//                                                                                descBit[i][j]|=descBit[k1*10000+k2][j];
//                                                                        }

//                                                                }
//                                                        }
//                                                }
//                                        }
//                                        else{
//                                                break;
//                                        }
//                                }

//                        }

//                }

//                for (int i = 0; i < n; ++i) {
//                        for (int k1 = 0; k1 <columns; ++k1) {
//                                for (int k2 = 0; k2 < 10000; ++k2) {
//                                        if(k1*10000+k2<n&&descBit[i][k1][k2]){
//                                                desc[i].insert(k1*10000+k2);
//                                        }
//                                }
//                        }

//                }

//             //   tfw.stop();


//        //parameters.getControlOutput()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
//                return desc;

//}



template<class T,class PAR>
    std::vector<std::vector<bool>> initReachable(T & graph,PAR& parameters,VertexGroups<size_t>* vg=0){


//        levinkov::Timer tfw;
//                        tfw.start();


        parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
                const size_t n=graph.numberOfVertices();

                // todo: use some matrix structure
                std::vector<std::vector<bool>> desc(n);
                std::vector<std::vector<std::bitset<10000>>> descBit(n);
                size_t columns=n/10000;
                size_t mod=n%10000;
                if(mod>0) columns++;


                for (int i = 0; i < n; ++i) {
                        desc[i]=std::vector<bool>(n,0);
                        descBit[i]=std::vector<std::bitset<10000>>(columns);
                }

                for (size_t v = 0; v < n; ++v) {

                        descBit[v][v/10000][v%10000]=1; //make this reflexive
                        size_t edges=graph.numberOfEdgesFromVertex(v);
                        for (int j = 0; j < edges; ++j) {
                                size_t w=graph.vertexFromVertex(v,j);
                                descBit[v][w/10000][w%10000]=1;


                        }
                }



                if(vg==0){
                        for (int k1 = 0; k1 <columns; ++k1) {
                                for (int k2 = 0; k2 < 10000; ++k2) {
                                        if(k1*10000+k2<n){
                                                for (int i = 0; i < n; ++i) {
                                                        if(descBit[i][k1][k2]){
                                                                for (int j = 0; j < columns; ++j) {
                                                                        descBit[i][j]|=descBit[k1*10000+k2][j];
                                                                }

                                                        }
                                                }
                                        }
                                        else{
                                                break;
                                        }
                                }

                        }
                }
                else{
                        for (int k1 = 0; k1 <columns; ++k1) {
                                for (int k2 = 0; k2 < 10000; ++k2) {
                                        if(k1*10000+k2<n){
                                                size_t maxTime=vg->getGroupIndex(k1*10000+k2);
                                                for (int t = 0; t < maxTime; ++t) {//TODO use time gap
                                                        const std::vector<size_t>& vertices=vg->getGroupVertices(t);
                                                        for (size_t i:vertices) {
                                                                if(descBit[i][k1][k2]){
                                                                        for (int j = 0; j < columns; ++j) {
                                                                                descBit[i][j]|=descBit[k1*10000+k2][j];
                                                                        }

                                                                }
                                                        }
                                                }
                                        }
                                        else{
                                                break;
                                        }
                                }

                        }

                }

                for (int i = 0; i < n; ++i) {
                        for (int k1 = 0; k1 <columns; ++k1) {
                                for (int k2 = 0; k2 < 10000; ++k2) {

                                        if(k1*10000+k2<n){
                                                desc[i][k1*10000+k2]=descBit[i][k1][k2];
                                        }
                                }
                        }

                }


                //tfw.stop();

        //parameters.getControlOutput()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;

                return desc;

}

   template<class T>
   std::vector<bool> getEdgeLabels(const T& graph,const std::vector<std::vector<size_t>>& paths){
        char delim=',';
        std::vector<bool> activeEdges(graph.numberOfEdges());


        for (size_t i=0;i<paths.size();i++) {
            std::vector<bool> activeVertices(graph.numberOfVertices());
            for (size_t j = 0; j < paths.at(i).size(); ++j) {
                activeVertices[paths.at(i).at(j)]=1;
            }
            for (size_t j = 0; j < paths.at(i).size(); ++j) {
                size_t v=paths.at(i).at(j);
                size_t n=graph.numberOfEdgesFromVertex(v);
                for (int k = 0; k < n; ++k) {
                    size_t w=graph.vertexFromVertex(v,k);
                    if(activeVertices[w]){
                        size_t e=graph.edgeFromVertex(v,k);
                        activeEdges[e]=1;
                    }

                }
            }
        }
        return activeEdges;
    }




   template<class T>
   std::vector<bool> getLiftedEdgeLabels(const T& edges,const std::vector<std::vector<size_t>>& paths,size_t numberOfVertices){
        char delim=',';
        std::vector<bool> activeEdges(edges.size());
        std::vector<size_t> vertexLabels(numberOfVertices);



        for (size_t i=0;i<paths.size();i++) {
            for (size_t j = 0; j < paths.at(i).size(); ++j) {
                vertexLabels[paths.at(i).at(j)]=i+1;
            }
       }
        for(size_t i=0;i<activeEdges.size();i++){
            size_t label1=vertexLabels[edges.at(i).at(0)];
            size_t label2=vertexLabels[edges.at(i).at(1)];
            if(label1==label2) activeEdges[i]=1;
        }
        return activeEdges;
    }


   template<class T>
   std::vector<bool> getBaseEdgeLabels(const T& edges,const std::vector<std::vector<size_t>>& paths,size_t numberOfVertices){
        char delim=',';
        std::vector<bool> activeEdges(edges.size());
        std::vector<size_t> vertexDescendants(numberOfVertices);
       // std::vector<bool> needCheck(numberOfVertices);

        //TODO: make check if all inner base edges stored in vertexDescendants, were present in edges contrainer
        size_t t=numberOfVertices+1;


        for (size_t i=0;i<paths.size();i++) {
            for (size_t j = 0; j < paths.at(i).size()-1; ++j) {
                size_t vertex=paths.at(i).at(j);
                vertexDescendants[vertex]=paths.at(i).at(j+1);
              //  needCheck[vertex]=1;
            }
            size_t lastVertex=*(paths.at(i).rbegin());
            vertexDescendants[lastVertex]=t;
       }
        for(size_t i=0;i<activeEdges.size();i++){
            size_t firstVertex=edges.at(i).at(0);
            size_t secondVertex=edges.at(i).at(1);
            size_t descendant=vertexDescendants[firstVertex];
            if(descendant==secondVertex){
                activeEdges[i]=1;
              //  needCheck[firstVertex]=0;
            }
        }
//        if(debug()){
//            for(size_t i=0;i<needCheck.size();i++){
//                if(needCheck[i]) //TODO exception
//            }

//        }

        return activeEdges;
    }


 template<class INSTANCE, class PAR>
   void createKnnBaseGraph(INSTANCE& instance,  PAR& parameters){
        parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
        parameters.writeControlOutput();

        std::vector<double> newBaseCosts;
        //std::vector<size_t> inOutEdges;
        size_t k=parameters.getKnnK();
        //std::vector<size_t> goodLongEdges;


        const andres::graph::Digraph<>& graph_=instance.getGraph();
        const size_t t_=instance.getTerminalNode();
        const size_t s_=instance.getSourceNode();
        const VertexGroups<size_t>& vg=instance.getVertexGroups();
        const std::vector<double>& edgeScore=instance.getEdgesScore();
        andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());



        std::vector<bool> finalEdges(graph_.numberOfEdges(),false);
        for (int v0 = 0; v0 < graph_.numberOfVertices(); ++v0) {
            std::unordered_map<int,std::list<size_t>> edgesToKeep;
            size_t l0=vg.getGroupIndex(v0);
            for (size_t ne = 0; ne < graph_.numberOfEdgesFromVertex(v0); ++ne) {
                size_t e=graph_.edgeFromVertex(v0,ne);
                size_t v1=graph_.vertexFromVertex(v0,ne);
    //			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
                if(v0==s_||v1==t_){
                    //tempGraph.insertEdge(v0,v1);
                    //newBaseCosts.push_back(edgeScore[e]);
                    finalEdges[e]=true;
                }
                else{
                    size_t l1=vg.getGroupIndex(v1);
                    size_t gap=l1-l0;
                    if(gap<=parameters.getMaxTimeBase()){
                    //if(gap<=parameters.getKnnTimeGap()){
                        //gap=std::min(parameters.getKnnTimeGap()+1,gap);
                        double cost=edgeScore[e];
                        if(edgesToKeep.count(gap)>0){
                            std::list<size_t>& smallList=edgesToKeep[gap];
                            auto it=smallList.begin();
                            double bsf=edgeScore[*it];
                            //std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
                            while(bsf>cost&&it!=smallList.end()){
                                it++;
                                size_t index=*it;
                                if(it!=smallList.end()){
                                    bsf=edgeScore[index];
                                    //	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
                                }
                            }
                            if(it!=smallList.begin()){
                                smallList.insert(it,e);
                                if(smallList.size()>k) smallList.pop_front();
                            }
                            else if(smallList.size()<k){
                                smallList.push_front(e);
                            }
                        }
                        else{
                            edgesToKeep[gap].push_front(e);
                        }
                    }
    //				else if(gap<=parameters.getMaxTimeBase()){
    //					if(getEdgeScore(e)<=parameters.getBaseUpperThreshold()){
    //						//tempGraph.insertEdge(v0,v1);
    //						//newBaseCosts.push_back(getEdgeScore(e));
    //						finalEdges[e]=true;
    //					}
    //
    //				}
                }
            }
            //std::cout.precision(4);
            double bsf=0;
            for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
                if(edgesToKeep.count(gap)>0){
                    auto& smallList=edgesToKeep[gap];
                    for(size_t e:smallList){
                        finalEdges[e]=true;
                        if(edgeScore[e]<bsf){
                            bsf=edgeScore[e];
                        }
                    }
                }
            }

            for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
                if(edgesToKeep.count(gap)>0){
                    auto& smallList=edgesToKeep[gap];
                    for(size_t e:smallList){
                        double score=edgeScore[e];
                        if(score<=parameters.getBaseUpperThreshold()){
                            finalEdges[e]=true;
                        }
                    }

                }
            }

        }

        for (int e = 0; e < graph_.numberOfEdges(); ++e) {
            if(finalEdges[e]){
                size_t v0=graph_.vertexOfEdge(e,0);
                size_t v1=graph_.vertexOfEdge(e,1);
                tempGraph.insertEdge(v0,v1);
                newBaseCosts.push_back(edgeScore[e]);
            }
        }

        if(newBaseCosts.size()!=tempGraph.numberOfEdges()){
            throw std::runtime_error("Error in base graph sparsification.");
        }



        instance.setGraph(tempGraph);
        instance.setEdgesScore(newBaseCosts);


        if(graph_.numberOfEdges()!=newBaseCosts.size()){
            parameters.getControlOutput()<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
            parameters.writeControlOutput();

        }
        else{

            parameters.getControlOutput()<<"edge number and graph size match "<<std::endl;
            parameters.writeControlOutput();
        }

        parameters.getControlOutput()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
        parameters.writeControlOutput();

    }






    template<class INSTANCE, class PAR>
    void keepFractionOfLifted(INSTANCE& instance,  PAR& parameters){


        parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
        parameters.writeControlOutput();
        //TODO run automaticLifted to find candidates first

        double negMaxValue=0;
        double posMinValue=0;
        //TODO adaptive lifted threshold in config file
        std::vector<double> newLiftedCosts;

        const andres::graph::Digraph<>& graph_=instance.getGraph();
        const andres::graph::Digraph<>& graphLifted_=instance.getGraphLifted();
        const std::vector<double>& liftedCosts=instance.getLiftedEdgesScore();
        const std::vector<std::unordered_set<size_t>>* pReachable =instance.getPReachable();
        const std::vector<std::unordered_set<size_t>>& reachable=*pReachable;
        const size_t t_=instance.getTerminalNode();
        const VertexGroups<size_t>& vg=instance.getVertexGroups();

        std::vector<double> newBaseEdgeScore=instance.getEdgesScore();


        andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


        negMaxValue=parameters.getNegativeThresholdLifted();
        posMinValue=parameters.getPositiveThresholdLifted();



        std::unordered_map<size_t,std::set<size_t>> liftedEdges;
        for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
            std::unordered_set<size_t> alternativePath;
            for (int i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
                size_t w=graph_.vertexFromVertex(v,i);
                for(size_t u:reachable[w]){
                    if(u!=w) alternativePath.insert(u);
                }
            }
            for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
                size_t w=graphLifted_.vertexFromVertex(v,i);
                if(w!=t_){
                    if(alternativePath.count(w)>0) liftedEdges[v].insert(w);

                }

            }
        }


        parameters.getControlOutput()<<"done"<<std::endl;
        parameters.writeControlOutput();


        for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
            size_t v0=graphLifted_.vertexOfEdge(i,0);
            size_t v1=graphLifted_.vertexOfEdge(i,1);
            int l0=vg.getGroupIndex(v0);
            int l1=vg.getGroupIndex(v1);
            double cost=liftedCosts.at(i);
            bool goodCost=(cost<=negMaxValue)||(cost>=posMinValue);
            //if(isReachable(v0,v1)){
            if(instance.isReachable(v0,v1)&&goodCost){

                int timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
                bool timeConstraint=l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
                //int modulo=timeGapDiff%parameters.getLongerIntervalLifted();
                //if(l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0)){
                if(timeConstraint){
                    if(liftedEdges[v0].count(v1)>0){
                        tempGraphLifted.insertEdge(v0,v1);
                        newLiftedCosts.push_back(cost);
                    }
                    else{
                        auto edgeTest=graph_.findEdge(v0,v1);
                        if(edgeTest.first){
                            newBaseEdgeScore[edgeTest.second]+=cost;  //Compensate that the lifted edge has been removed
                        }

                    }

                }
            }

        }



        instance.setGraphLifted(tempGraphLifted);
        instance.setEdgesScore(newBaseEdgeScore);
        instance.setLiftedEdgesScore(newLiftedCosts);

        parameters.getControlOutput()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
        parameters.writeControlOutput();

        if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){

            parameters.getControlOutput()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
        }
        else{

            parameters.getControlOutput()<<"lifted edge number and lifted graph size match "<<std::endl;

        }
        parameters.writeControlOutput();

    }




}





#endif // LDP_METHODS_HXX
