/*
 * disjoint-paths-init.hxx
 *
 *  Created on: Sep 10, 2018
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_
#define INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_

#include <stdexcept>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <lifted_disjoint_paths/ldp_parameters.hxx>
#include <stack>
#include <unordered_set>
#include <iterator>
#include <unordered_map>
#include <string>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <list>
#include "ldp_functions.hxx"
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <disjoint-paths/disjointPathsMethods.hxx>
#include <disjoint-paths/completeStructure.hxx>
#include <disjoint-paths/twoGraphsInputStructure.hxx>
#include "lifted_disjoint_paths/ldp_directed_graph.hxx"
#include "config.hxx"

namespace py = pybind11;
namespace LPMP{
namespace lifted_disjoint_paths {



class LdpInstance {
public:


    LdpInstance(LdpParameters<>& configParameters);
     LdpInstance(LdpParameters<>& configParameters,disjointPaths::CompleteStructure<>& cs);
//     LdpInstance(LdpParameters<>& configParameters,const disjointPaths::TwoGraphsInputStructure& twoGraphsIS);
     LdpInstance(LdpParameters<>& configParameters,const py::array_t<size_t>& baseEdges,const py::array_t<size_t>& liftedEdges,const  py::array_t<double>& baseCosts,const  py::array_t<double>& liftedCosts,disjointPaths::VertexGroups<>& pvg);
   // LdpInstance(const ConfigDisjoint<>& configParameters,char delim=',',disjointPaths::CompleteStructure<>* cs=0,size_t minTime=0,size_t maxTime=0);

	bool isReachable(size_t i,size_t j) const{
		if(i==t_||j==s_) return false;  //Assume no path from the terminal node
		if(i==s_||j==t_) return true;
		if(reachable.size()==0) return true;
		return reachable[i].count(j)>0;
	}

	const std::unordered_set<size_t>& reachableFromVertex(size_t v)const{
		return reachable.at(v);
	}


	const andres::graph::Digraph<>& getGraph() const  {
		return graph_;
	}

    void setGraph(const andres::graph::Digraph<>& newGraph){
        graph_=newGraph;
    }


	const andres::graph::Digraph<>& getGraphLifted() const  {
		return graphLifted_;
	}

    void setGraphLifted(const andres::graph::Digraph<>& newGraphLifted){
        graphLifted_=newGraphLifted;
    }

	const size_t getGapLifted() const {
		return parameters.getMaxTimeLifted();
	}


    const std::vector<std::unordered_set<size_t>>* getPReachable(){
        return &reachable;
    }

	size_t getSourceNode() const {
		return s_;
	}

	size_t getTerminalNode() const {
		return t_;
	}

	double getEdgeScore(size_t e) const {
		return edgeScore[e];
	}

	double getEdgeScore(size_t v0,size_t v1) const {
		auto findEdge=graph_.findEdge(v0,v1);
		assert(findEdge.first);
		return edgeScore[findEdge.second];
	}

	const std::vector<double>& getEdgesScore() {
		return edgeScore;
	}

    void setEdgesScore(const std::vector<double>& newEdgesScore){
        edgeScore=newEdgesScore;
    }

	const std::vector<double>& getLiftedEdgesScore() {
			return liftedEdgeScore;
	}

    void setLiftedEdgesScore(const std::vector<double>& newLiftedScore){
        liftedEdgeScore=newLiftedScore;
    }



	const std::vector<double>& getVerticesScore() const{
		return vertexScore;
	}

	double getLiftedEdgeScore(size_t e) const {
		return liftedEdgeScore[e];
	}

	double getLiftedEdgeScore(size_t v0,size_t v1) const {
		auto findEdge=graphLifted_.findEdge(v0,v1);
		assert(findEdge.first);
		return liftedEdgeScore[findEdge.second];
	}

	double getVertexScore(size_t v) const {
		return vertexScore[v];
	}

    const disjointPaths::VertexGroups<size_t>& getVertexGroups()const {
        return vertexGroups;
	}

	size_t getGroupIndex(size_t v)const {
		return vertexGroups.getGroupIndex(v);
	}


	size_t getEdgeVarIndex(size_t edgeIndex)const {
		return edgeIndex+numberOfVertices;
	}


	size_t getLiftedEdgeVarIndex(size_t liftedEdgeIndex)const {
		return liftedEdgeIndex+numberOfEdges+numberOfVertices;
	}

	size_t getVertexVarIndex(size_t vertexIndex)const{
		return vertexIndex;
	}

    bool existLiftedEdge(const size_t v,const size_t w)const{
        bool value=liftedStructure.at(v).isWithinBounds(w)&&liftedStructure.at(v)[w]>0;
        return value;
    }

//	std::vector<bool>& getBaseEdgeLabels()  {
//		return baseEdgeLabels;
//	}

    bool isStrongBase(size_t v,size_t w) const;

    void sparsifyLiftedGraph();

    LdpParameters<>& parameters;
    disjointPaths::VertexGroups<size_t> vertexGroups;
	size_t minV=0;
	size_t maxV=0;

    mutable std::vector<size_t> sncNeighborStructure;
    mutable std::vector<size_t> sncBUNeighborStructure;

    mutable std::vector<double> sncTDStructure;
    mutable std::vector<double> sncBUStructure;
    mutable std::vector<char> sncClosedVertices;
    mutable std::vector<double> sncLiftedMessages;


    size_t getNumberOfVertices() const{
        return numberOfVertices;
    }

    const LdpDirectedGraph& getMyGraph() const{
        return myGraph;
    }

    const LdpDirectedGraph& getMyGraphLifted() const{
        return myGraphLifted;
    }

private:

    //LdpInstance(const LdpInstance& ldpI);
    void init();

	size_t s_;
	size_t t_;

	std::vector<double> vertexScore;
	std::vector<double> edgeScore;
	std::vector<double> liftedEdgeScore;
	//std::vector<std::vector<bool>> desc;
	std::vector<std::unordered_set<size_t>> reachable;

	andres::graph::Digraph<> graph_;
	andres::graph::Digraph<> graphLifted_;

    LdpDirectedGraph myGraph;
    LdpDirectedGraph myGraphLifted;

    std::vector<std::unordered_set<size_t>> strongBaseEdges;

	std::vector<bool> baseEdgeLabels;

	size_t numberOfVertices;
	size_t numberOfEdges;
	size_t numberOfLiftedEdges;


	void readGraph(std::ifstream& data,size_t maxVertex,char delim);
    void readGraphWithTime(size_t minTime,size_t maxTime,disjointPaths::CompleteStructure<>* cs);
   // void sparsifyBaseGraph();
   // void sparsifyLiftedGraph();
    void initLiftedStructure();

    std::vector<ShiftedVector<char>> liftedStructure;


};





}
}
#endif /* INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_ */
