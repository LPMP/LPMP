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
#include "lifted_disjoint_paths/ldp_directed_graph.hxx"
#include "config.hxx"
#include "lifted_disjoint_paths/ldp_complete_structure.hxx"
#include "lifted_disjoint_paths/ldp_vertex_groups.hxx"
#include "ldp_batch_process.hxx"
#include <chrono>


namespace py = pybind11;
namespace LPMP{
namespace lifted_disjoint_paths {



class LdpInstance {
public:



     LdpInstance(LdpParameters<>& configParameters,CompleteStructure<>& cs);
     LdpInstance(LdpParameters<>& configParameters,LdpBatchProcess& BP);
     LdpInstance(LdpParameters<>& configParameters, const py::array_t<size_t>& baseEdges, const py::array_t<size_t>& liftedEdges, const  py::array_t<double>& baseCosts, const  py::array_t<double>& liftedCosts, const py::array_t<double> &verticesCosts, VertexGroups<>& pvg);

	bool isReachable(size_t i,size_t j) const{
		if(i==t_||j==s_) return false;  //Assume no path from the terminal node
		if(i==s_||j==t_) return true;
		if(reachable.size()==0) return true;
		return reachable[i].count(j)>0;
	}

	const std::unordered_set<size_t>& reachableFromVertex(size_t v)const{
		return reachable.at(v);
	}



	const size_t getGapLifted() const {
		return parameters.getMaxTimeLifted();
	}

    const size_t getGapBase() const {
        return parameters.getMaxTimeBase();
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


	const std::vector<double>& getVerticesScore() const{
		return vertexScore;
	}


	double getVertexScore(size_t v) const {
                assert(v<vertexScore.size());
		return vertexScore[v];
	}

    const VertexGroups<size_t>& getVertexGroups()const {
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


    bool isStrongBase(size_t v,size_t w) const;

    double evaluateClustering(const std::vector<size_t>& labels) const;

    std::vector<std::unordered_set<size_t>> initReachableLdp(const LdpDirectedGraph &graph, LdpParameters<> &parameters, const VertexGroups<size_t> *vg=nullptr);




    LdpParameters<>& parameters;
    LPMP::VertexGroups<> vertexGroups;
	size_t minV=0;
	size_t maxV=0;

    mutable std::vector<size_t> sncNeighborStructure;
    mutable std::vector<size_t> sncBUNeighborStructure;

    mutable std::vector<double> sncTDStructure;
    mutable std::vector<double> sncBUStructure;
    mutable std::vector<char> sncClosedVertices;
    mutable std::vector<double> sncLiftedMessages;
    mutable std::vector<char> sncVerticesInScope;


    const size_t& getNumberOfVertices() const{
        return numberOfVertices;
    }

    const LdpDirectedGraph& getMyGraph() const{
        return myGraph;
    }

    const LdpDirectedGraph& getMyGraphLifted() const{
        return myGraphLifted;
    }

    bool checkStrongBase(const size_t& v,const size_t& w)const;

    bool canJoin(const size_t& v,const size_t& w)const{
        if(parameters.isMustCutMissing()){
            bool value=canJoinStructure.at(v).isWithinBounds(w)&&canJoinStructure.at(v)[w]>0;
            return value;
        }
        else{
            return true;
        }

    }

private:


    void initAdaptiveThresholds(const LdpDirectedGraph *pBaseGraph, const LdpDirectedGraph *pLiftedGraph);
    void init();

    void sparsifyBaseGraphNew(const LdpDirectedGraph& inputGraph,bool zeroCost=false);


    void sparsifyLiftedGraphNew(const LdpDirectedGraph& inputLiftedGraph);


    void initLiftedStructure();

    void initCanJoinStructure(const LdpDirectedGraph &completeGraph);

	size_t s_;
	size_t t_;

	std::vector<double> vertexScore;
   std::vector<std::unordered_set<size_t>> reachable;


    LdpDirectedGraph myGraph;
    LdpDirectedGraph myGraphLifted;
    const LdpDirectedGraph* pCompleteGraph;

    std::vector<std::unordered_set<size_t>> strongBaseEdges;

	std::vector<bool> baseEdgeLabels;

	size_t numberOfVertices;
	size_t numberOfEdges;
	size_t numberOfLiftedEdges;


    double negativeLiftedThreshold;
    double positiveLiftedThreshold;
    double baseThreshold;



    std::vector<ShiftedVector<char>> liftedStructure;
    std::vector<ShiftedVector<char>> canJoinStructure;


};





}
}
#endif /* INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_ */
