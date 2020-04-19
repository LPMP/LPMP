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
#include <union_find.hxx>

namespace LPMP {


class ldp_min_assignment_factor
{
public:

	//By default, all edges are lifted. However vu or uw can be base too.
	ldp_min_assignment_factor(size_t minVertex1, size_t vertexNumber1,size_t minVertex2,size_t vertexNumber2):
		minV1(minVertex1),
		minV2(minVertex2),
		vertexN1(vertexNumber1),
		vertexN2(vertexNumber2),
		unassigned(vertexNumber2+vertexNumber1),
		uf(vertexNumber1+vertexNumber2)
    {
		edgeCosts=std::vector<std::unordered_map<size_t,double>>(vertexNumber1);
    }

	//double LowerBound() const;

	double EvaluatePrimal() const{
		double value=0;
		for (int i = 0; i < edgeCosts.size(); ++i) {
			if(primal_[i]) value+=edgeCosts[i];
		}
		return value;
	}

	void updateCost(size_t v1,size_t v2,double value){
		edgeCosts[v1-vertexN1][v2-vertexN2]+=value;
	}


	LPMP::linear_assignment_problem_input createLAStandard() const{
		LPMP::linear_assignment_problem_input lapInput;
		for (size_t i = 0; i < vertexN1; ++i) {
			for(std::pair<size_t,double> e:edgeCosts.at(i)){
				if(e.second<0){
					lapInput.add_assignment(i,e.first,e.second);
				}
			}
		}
		return lapInput;
	}

	LPMP::linear_assignment_problem_input laExcludeEdge(size_t v1,size_t v2) const{
		LPMP::linear_assignment_problem_input lapInput;
		for (size_t i = 0; i < vertexN1; ++i) {
			for(std::pair<size_t,double> e:edgeCosts.at(i)){
				if(e.second<0&&(v1!=i||v2!=e.first)){
					lapInput.add_assignment(i,e.first,e.second);
				}
			}
		}
		return lapInput;
	}

	LPMP::linear_assignment_problem_input laExcludeVertices(size_t v1,size_t v2) const{
		LPMP::linear_assignment_problem_input lapInput;
		for (size_t i = 0; i < vertexN1; ++i) {
			if(i==v1) continue;
			size_t vertex=i;
			if(i>v1) vertex--;
			for(std::pair<size_t,double> e:edgeCosts.at(i)){
				if(e.second<0&&e.second!=v2){
					size_t vertex2=e.second;
					if(vertex2>v2) vertex2--;
					lapInput.add_assignment(vertex,vertex2,e.second);
				}
			}
		}
		return lapInput;
	}

	std::vector<size_t> getLabeling(linear_assignment_problem_input& lapInput,size_t shift1=vertexN1,size_t shift2=vertexN2)const {
		MCF::SSP<long,double> mcf(lapInput.no_mcf_nodes(),lapInput.no_mcf_edges());
		lapInput.initialize_mcf(mcf);
		double result=mcf.solve();

		std::vector<size_t> labeling(vertexN1);
		for (int i = 0; i < mcf.no_edges(); ++i) {
			if(mcf.flow(i)>0.5){
				int label=mcf.head(i);
				if(label>shift2) label++;
				int vertex=mcf.tail(i);
				if(vertex>shift1) vertex++;
				//	std::cout<<"label "<<label<<"vertex "<<vertex<<std::endl;
				if(vertex<vertexN1&&label<vertexN2){
					labeling[vertex]=label;
				}
			}
		}
		return labeling;
	}

	double evaluateLabeling(std::vector<size_t>& labeling)const {
		double cost=0;
		for (int i = 0; i < vertexN1; ++i) {
			cost+=edgeCosts[i].at(labeling[i]);
		}
		return cost;
	}

	double LowerBound() const{
		LPMP::linear_assignment_problem_input lapInput=createLAStandard();
		std::vector<size_t> labeling=getLabeling(lapInput);
		return evaluateLabeling(labeling);
	}


	bool notAllowed(size_t v1,size_t v2){  //TODO make this based on reachability
		return false;
	}

	std::unordered_map<std::pair<size_t,size_t>,double> adjustCostAndGetMessages(){
		LPMP::linear_assignment_problem_input lapInput=createLAStandard();
		std::vector<size_t> labeling=getLabeling(lapInput);
		double optimalValue= evaluateLabeling(labeling);
		std::vector<size_t> isNotOptZero; //enough to have single vertex instead of set.
		std::vector<std::unordered_set<size_t>> isOptOne;
		for (int vertex = 0; vertex < vertexN1; ++vertex) {
			isOptOne[vertex].insert(labeling[vertex]);
			isNotOptZero[vertex]=labeling[vertex];
		}
		for (int v1 = 0; v1 < vertexN1; ++v1) {
			for (int v2 = 0; v2 < vertexN2; ++v2) {
				if(isOptOne[v1].count(v2)){ //has label 1 in optimum
					if(isNotOptZero[v1]==v2){
						//do nothing delta is zero
					}
					else{
						linear_assignment_problem_input lapInput=laExcludeEdge(v1,v2);
						labeling=getLabeling(lapInput);
						double value=evaluateLabeling(labeling);
						double delta=optimalValue-value;
						for (int vertex = 0; vertex < vertexN1; ++vertex) {
							isOptOne[vertex].insert(labeling[vertex]);
						}
						for (int vertex = 0; vertex < vertexN1; ++vertex) {
							size_t vertex2=isNotOptZero[vertex];
							if(labeling[vertex]!=vertex2){
								isNotOptZero=vertexN2;
							}
						}

						//TODO update structures isOptOne and isNotOptZero
					}

				}
			}
		}

	}



	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }

	void init_primal() { primal_ = std::vector<size_t>(vertexN1,unassigned); }


	void updateCost(size_t edgeId,double updateValue){  //update cost of one edge, assumed indices 0-2
		edgeCosts[edgeId]+=updateValue;
	}

private:
    std::vector<size_t> primal_;
    std::vector<std::unordered_map<size_t,double>> edgeCosts;   //Ideally store edges and costs as node->(node,cost) in a map.
    size_t minV1,vertexN1,minV2, vertexN2,unassigned;

    union_find uf;

};








}//namespace LPMP end





#endif /* INCLUDE_LIFTED_DISJOINT_PATHS_LDP_MIN_ASSIGNMENT_FACTOR_HXX_ */
