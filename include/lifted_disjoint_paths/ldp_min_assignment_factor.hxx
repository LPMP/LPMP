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
#include <array>
#include <unordered_map>

namespace LPMP {

template<class GRAPH_STRUCT,class EDGE_ITERATOR>
class ldp_min_assignment_factor
{
public:

	//By default, all edges are lifted. However vu or uw can be base too.
	ldp_min_assignment_factor(GRAPH_STRUCT inputGraph,EDGE_ITERATOR begin,EDGE_ITERATOR end); //Maybe remember vertex indices instead of edge indices?

	double LowerBound();

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


	void updateCost(size_t edgeId,double updateValue){  //update cost of one edge, assumed indices 0-2
		edgeCosts[edgeId]+=updateValue;
	}

private:
    std::vector<bool> primal_;
    GRAPH_STRUCT graph;
    std::unordered_map<size_t,size_t> leftVertices;
    std::unordered_map<size_t,size_t> rightVertices;
    //std::vector<std::pair<size_t,size_t>> edges;
    std::vector<std::array<size_t,2>> edges;
    std::vector<double> edgeCosts;

};

template<class GRAPH_STRUCT,class EDGE_ITERATOR>
inline ldp_min_assignment_factor<GRAPH_STRUCT,EDGE_ITERATOR>::ldp_min_assignment_factor(GRAPH_STRUCT inputGraph,EDGE_ITERATOR begin,EDGE_ITERATOR end): //Maybe remember vertex indices instead of edge indices?
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
		//edges.push_back(std::pair<size_t,size_t>(v0,v1));
		edges.push_back({v0,v1});
	}
	edgeCosts=std::vector<double>(edges.size());
}



template<class GRAPH_STRUCT,class EDGE_ITERATOR>
inline double ldp_min_assignment_factor<GRAPH_STRUCT,EDGE_ITERATOR>::LowerBound(){
	double minValue=0;
	LPMP::linear_assignment_problem_input lapInput;
	for (int i = 0; i < edges.size(); ++i) {
		if(edgeCosts[i]<0){
			//lapInput.add_assignment(edges[i].first,edges[i].second,edgeCosts[i]);
			lapInput.add_assignment(edges[i][0],edges[i][1],edgeCosts[i]);
		}
	}
	MCF::SSP<long,double> mcf(lapInput.no_mcf_nodes(),lapInput.no_mcf_edges());
	lapInput.initialize_mcf(mcf);
	double result=mcf.solve();

	std::vector<size_t> labeling(leftVertices.size());
	for (int i = 0; i < mcf.no_edges(); ++i) {
		if(mcf.flow(i)>0.5){
			int label=mcf.head(i)-leftVertices.size();
			int vertex=mcf.tail(i);
			//	std::cout<<"label "<<label<<"vertex "<<vertex<<std::endl;
			if(vertex<leftVertices.size()&&label<rightVertices.size()){
				labeling[vertex]=label;
			}
		}
	}
	for (int i = 0; i < edges.size(); ++i) {
		size_t v0=leftVertices[edges[i][0]];
		size_t v1=rightVertices[edges[i][1]];
		if(labeling[v0]==v1){
			primal_[i]=1;
			minValue+=edgeCosts[i];
		}
		else{
			primal_[i]=0;
		}
	}
	return minValue;
}



}//namespace LPMP end





#endif /* INCLUDE_LIFTED_DISJOINT_PATHS_LDP_MIN_ASSIGNMENT_FACTOR_HXX_ */
