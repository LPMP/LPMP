/*
 * ldp_min_assignment_factor.hxx
 *
 *  Created on: Apr 7, 2020
 *      Author: fuksova
 */

#ifndef INCLUDE_LIFTED_DISJOINT_PATHS_LDP_MIN_ASSIGNMENT_FACTOR_HXX_
#define INCLUDE_LIFTED_DISJOINT_PATHS_LDP_MIN_ASSIGNMENT_FACTOR_HXX_



#include <cstdlib>
#include "graph_matching/graph_matching_input.h"
#include "graph_matching/min_cost_flow_factor_ssp.hxx"
#include "graph_matching/matching_problem_input.h"
#include <vector>
#include <unordered_map>

namespace LPMP {

template<class GRAPH_STRUCT,class EDGE_ITERATOR>
class ldp_min_assignment_factor
{
public:

	//By default, all edges are lifted. However vu or uw can be base too.
	ldp_min_assignment_factor(GRAPH_STRUCT inputGraph,EDGE_ITERATOR begin,EDGE_ITERATOR end): //Maybe remember vertex indices instead of edge indices?
		graph(inputGraph)
    {
		EDGE_ITERATOR it=begin;
		size_t leftCounter=0;
		size_t rightCounter=0;
		for(;it!=end;it++){
			size_t v0=inputGraph.vertexOfEdge(*it,0);
			size_t v1=inputGraph.vertexOfEdge(*it,1);
			if(leftVertices.count(v0)==0){
				leftVertices[v0]=leftCounter;
				leftCounter++;
			}
			if(rightVertices.count(v1)==0){
				leftVertices[v1]=rightCounter;
				rightCounter++;
			}
			edges.push_back(std::pair<size_t,size_t>(v0,v1));
		}
		edgeCosts=std::vector<double>(edges.size());
    }

	double LowerBound() const{
		double minValue=0;
		return minValue;

	}

	double EvaluatePrimal() const{
		double value=0;
		for (int i = 0; i < edgeCosts.size(); ++i) {
			if(primal_[i]) value+=edgeCosts[i];
		}
		return value;
	}

	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }

	void init_primal() { primal_ = std::vector<bool>(edges.size()); }

	double delta(short edgeId){  //difference: label[edgeId]=1-label[edgeId]=0
        return 0;
	}

	void updateCost(size_t edgeId,double updateValue){  //update cost of one edge, assumed indices 0-2
		edgeCosts[edgeId]+=updateValue;
	}

private:
    std::vector<bool> primal_;
    GRAPH_STRUCT graph;
    std::unordered_map<size_t,size_t> leftVertices;
    std::unordered_map<size_t,size_t> rightVertices;
    std::vector<std::pair<size_t,size_t>> edges;
    std::vector<double> edgeCosts;

};





}





#endif /* INCLUDE_LIFTED_DISJOINT_PATHS_LDP_MIN_ASSIGNMENT_FACTOR_HXX_ */
