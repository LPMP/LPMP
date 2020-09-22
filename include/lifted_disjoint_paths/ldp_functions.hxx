/*
 * ldp_functions.hxx
 *
 *  Created on: Mar 30, 2020
 *      Author: fuksova
 */

#ifndef INCLUDE_LIFTED_DISJOINT_PATHS_LDP_FUNCTIONS_HXX_
#define INCLUDE_LIFTED_DISJOINT_PATHS_LDP_FUNCTIONS_HXX_


#include <stdexcept>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <lifted_disjoint_paths/ldp_config.hxx>
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
#include <utility>
#include <config.hxx>






namespace LPMP{





namespace lifted_disjoint_paths {




template<class T=char>
 std::vector<std::string> split(
		std::string inputString, T delim) {
	size_t occurence = 0;
	size_t newOccurence = 0;
	std::vector<std::string> strings;
	while (newOccurence < inputString.size()) {
		newOccurence = std::min(inputString.find_first_of(delim, occurence),
				inputString.size());

		std::string newString(inputString, occurence, newOccurence - occurence);
		strings.push_back(newString);
		newOccurence = newOccurence + 1;
		occurence = newOccurence;
	}

	return strings;
}



//template<class T = size_t>
//struct VertexGroups {
//public:
//	VertexGroups(std::unordered_map<size_t,std::vector<size_t>> groups_,std::vector<size_t> vToGroup_):
//		groups(groups_),vToGroup(vToGroup_)
//	{
//		maxVertex=vToGroup_.size()-3;
//		maxTime=*(vToGroup.rbegin());
//	}
//	VertexGroups(){
//		maxVertex=0;
//		maxTime=0;
//	}

//	VertexGroups(const ConfigDisjoint<>& parameters,char delim){
//		//std::vector<std::vector<size_t>> groups;

//		std::string line;
//		std::vector<std::string> strings;
//		std::ifstream timeData(parameters.getTimeFileName());

//		size_t lineCounter=0;
//		std::vector<size_t> currentGroup;

//		unsigned int previousTime=0;
//		unsigned int time=0;

//		while (std::getline(timeData, line) && !line.empty()) {

//			//std::cout<<"time file line "<<lineCounter<<" delim "<<delim<<std::endl;
//			lineCounter++;

//			strings = split(line,delim);
//			//std::cout<<"split "<<std::endl;
//			if (strings.size() < 2) {
//				throw std::runtime_error(
//						std::string("Vertex and time frame expected"));
//			}

//			unsigned int v = std::stoul(strings[0]);
//			time = std::stoul(strings[1]);
//			if(time>parameters.getMaxTimeFrame()){
//				groups[previousTime]=currentGroup;
//				break;
//			}
//			if(vToGroup.size()!=v){
//				throw std::runtime_error(
//						std::string("Wrong vertex numbering in time file"));
//			}
//			else{

//				vToGroup.push_back(time);
//				//if(time==groups.size()+1){  //Change groups to map?
//				if(time==previousTime||previousTime==0){
//					currentGroup.push_back(v);
//				}
//				//else if(time==groups.size()+2){
//				else{
//					//groups.push_back(currentGroup);
//					groups[previousTime]=currentGroup;
//					currentGroup=std::vector<size_t>();
//					currentGroup.push_back(v);
//				}
//				previousTime=time;
//				//				else{
//				//					throw std::runtime_error("Wrong time frame numbering.");
//				//				}

//			}
//		}


//		maxTime=*(vToGroup.rbegin());

//		//time frame of s is zero
//		vToGroup.push_back(0);
//		currentGroup=std::vector<size_t>();
//		currentGroup.push_back(vToGroup.size()-1);
//		groups[0]=currentGroup;

//		//time frame of t is maxTime
//		vToGroup.push_back(maxTime+1);
//		currentGroup=std::vector<size_t>();
//		currentGroup.push_back(vToGroup.size()-1);
//		groups[maxTime+1]=currentGroup;

//		maxVertex=vToGroup.size()-3;

//	}



//	size_t getGroupIndex(size_t v) const{
//		return vToGroup[v];
//	}

//	const std::vector<size_t>& getGroupVertices(size_t index){
//		return groups[index];
//	}

//	//ID of maximal valid vertex (i.e. without s and t)
//	size_t getMaxVertex() const {
//		return maxVertex;
//	}

//	//Time of the last video frame
//	size_t getMaxTime() const {
//			return maxTime;
//    }

//private:
//	//std::vector<std::vector<size_t>> groups;
//	std::unordered_map<size_t,std::vector<size_t>> groups;
//	std::vector<size_t> vToGroup;
//	size_t maxVertex;
//	size_t maxTime;


//};




//template<class T>
//	std::vector<std::unordered_set<size_t>> initReachableSet(T & graph,VertexGroups<size_t>* vg=0){


//		std::cout<<"Run Floyd Warshall"<<std::endl;

//		const size_t n=graph.numberOfVertices();

//		// todo: use some matrix structure
//		std::vector<std::unordered_set<size_t>> desc(n);
//		std::vector<std::vector<std::bitset<10000>>> descBit(n);
//		size_t columns=n/10000;
//		size_t mod=n%10000;
//		if(mod>0) columns++;


//		for (int i = 0; i < n; ++i) {
//			descBit[i]=std::vector<std::bitset<10000>>(columns);
//		}

//		for (size_t v = 0; v < n; ++v) {

//			descBit[v][v/10000][v%10000]=1; //make this reflexive
//			size_t edges=graph.numberOfEdgesFromVertex(v);
//			for (int j = 0; j < edges; ++j) {
//				size_t w=graph.vertexFromVertex(v,j);
//				descBit[v][w/10000][w%10000]=1;


//			}
//		}



//		if(vg==0){
//			for (int k1 = 0; k1 <columns; ++k1) {
//				for (int k2 = 0; k2 < 10000; ++k2) {
//					if(k1*10000+k2<n){
//						for (int i = 0; i < n; ++i) {
//							if(descBit[i][k1][k2]){
//								for (int j = 0; j < columns; ++j) {
//									descBit[i][j]|=descBit[k1*10000+k2][j];
//								}

//							}
//						}
//					}
//					else{
//						break;
//					}
//				}

//			}
//		}
//		else{
//			for (int k1 = 0; k1 <columns; ++k1) {
//				for (int k2 = 0; k2 < 10000; ++k2) {
//					if(k1*10000+k2<n){
//						size_t maxTime=vg->getGroupIndex(k1*10000+k2);
//						for (int t = 0; t < maxTime; ++t) {//TODO use time gap
//							const std::vector<size_t>& vertices=vg->getGroupVertices(t);
//							for (size_t i:vertices) {
//								if(descBit[i][k1][k2]){
//									for (int j = 0; j < columns; ++j) {
//										descBit[i][j]|=descBit[k1*10000+k2][j];
//									}

//								}
//							}
//						}
//					}
//					else{
//						break;
//					}
//				}

//			}

//		}

//		for (int i = 0; i < n; ++i) {
//			for (int k1 = 0; k1 <columns; ++k1) {
//				for (int k2 = 0; k2 < 10000; ++k2) {
//					if(k1*10000+k2<n&&descBit[i][k1][k2]){
//						desc[i].insert(k1*10000+k2);
//					}
//				}
//			}

//		}


//		return desc;

//}



//template<class T>
//	std::vector<std::vector<bool>> initReachable(T & graph,ConfigDisjoint<>& parameters,VertexGroups<size_t>* vg=0){


//		std::cout<<"Run Floyd Warshall"<<std::endl;
//		parameters.infoFile()<<"Run Floyd Warshall"<<std::endl;
//		const size_t n=graph.numberOfVertices();

//		// todo: use some matrix structure
//		std::vector<std::vector<bool>> desc(n);
//		std::vector<std::vector<std::bitset<10000>>> descBit(n);
//		size_t columns=n/10000;
//		size_t mod=n%10000;
//		if(mod>0) columns++;


//		for (int i = 0; i < n; ++i) {
//			desc[i]=std::vector<bool>(n,0);
//			descBit[i]=std::vector<std::bitset<10000>>(columns);
//		}

//		for (size_t v = 0; v < n; ++v) {

//			descBit[v][v/10000][v%10000]=1; //make this reflexive
//			size_t edges=graph.numberOfEdgesFromVertex(v);
//			for (int j = 0; j < edges; ++j) {
//				size_t w=graph.vertexFromVertex(v,j);
//				descBit[v][w/10000][w%10000]=1;


//			}
//		}



//		if(vg==0){
//			for (int k1 = 0; k1 <columns; ++k1) {
//				for (int k2 = 0; k2 < 10000; ++k2) {
//					if(k1*10000+k2<n){
//						for (int i = 0; i < n; ++i) {
//							if(descBit[i][k1][k2]){
//								for (int j = 0; j < columns; ++j) {
//									descBit[i][j]|=descBit[k1*10000+k2][j];
//								}

//							}
//						}
//					}
//					else{
//						break;
//					}
//				}

//			}
//		}
//		else{
//			for (int k1 = 0; k1 <columns; ++k1) {
//				for (int k2 = 0; k2 < 10000; ++k2) {
//					if(k1*10000+k2<n){
//						size_t maxTime=vg->getGroupIndex(k1*10000+k2);
//						for (int t = 0; t < maxTime; ++t) {//TODO use time gap
//							const std::vector<size_t>& vertices=vg->getGroupVertices(t);
//							for (size_t i:vertices) {
//								if(descBit[i][k1][k2]){
//									for (int j = 0; j < columns; ++j) {
//										descBit[i][j]|=descBit[k1*10000+k2][j];
//									}

//								}
//							}
//						}
//					}
//					else{
//						break;
//					}
//				}

//			}

//		}

//		for (int i = 0; i < n; ++i) {
//			for (int k1 = 0; k1 <columns; ++k1) {
//				for (int k2 = 0; k2 < 10000; ++k2) {

//					if(k1*10000+k2<n){
//						desc[i][k1*10000+k2]=descBit[i][k1][k2];
//					}
//				}
//			}

//		}


//		return desc;

//}


//template<class T = size_t>
//struct CompleteStructure {
//public:

//	CompleteStructure(ConfigDisjoint<T> & configParameters,char delim=','){

//		vg=VertexGroups<>(configParameters,delim);
//		maxTime=vg.getMaxTime();
//		std::cout<<"max time "<<maxTime<<std::endl;
//		configParameters.infoFile()<<"max time "<<maxTime<<std::endl;
//		completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);
//		std::cout<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
//		configParameters.infoFile()<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
//		std::string line;
//		std::ifstream data(configParameters.getGraphFileName());

//		std::getline(data, line);
//		double objValue=0;

//		std::cout << "Read big graph" << std::endl;
//        configParameters.infoFile()<<"Read big graph" << std::endl;
//		std::vector<std::string> strings;


//		std::cout<<"Reading vertices from file. "<<std::endl;
//		configParameters.infoFile()<<"Reading vertices from file. "<<std::endl;
//		configParameters.infoFile().flush();
//		//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
//		while (std::getline(data, line) && !line.empty()) {

//		}

//		std::cout<<"Reading base edges from file. "<<std::endl;
//		configParameters.infoFile()<<"Reading base edges from file. "<<std::endl;
//		configParameters.infoFile().flush();
//		size_t maxLabel=0;
//		while (std::getline(data, line) && !line.empty()) {

//			strings = split(line, delim);

//			unsigned int v = std::stoul(strings[0]);
//			unsigned int w = std::stoul(strings[1]);
//			size_t l0=vg.getGroupIndex(v);
//			size_t l1=vg.getGroupIndex(w);
//			if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;
//			//std::cout<<"edge "<<v<<" "<<w<<std::endl;

//			if(l1-l0<=configParameters.getMaxTimeGapComplete()){
//				double score = std::stod(strings[2]);
//				completeGraph.insertEdge(v,w);
//				completeScore.push_back(score);
//			}

//		}


//		//objValue+=maxLabel*parameters.getInputCost()+maxLabel*parameters.getOutputCost();

//		data.close();
//	}

//	andres::graph::Digraph<> completeGraph;
//	std::vector<double> completeScore;
//	VertexGroups<> vg;
//	size_t maxTime;



//};




}  //end of namespace








template<class T>
struct ShiftedVector{
public:
    ShiftedVector<T>(size_t boundary1,size_t boundary2,const T& value): //inclusive both min and max vertex
    minVertex(std::min(boundary1,boundary2)),maxVertex(std::max(boundary1,boundary2))
    {
        //const size_t min = std::min(boundary1,boundary2);
        //const size_t max = std::max(boundary1,boundary2);
        //const size_t size = max - min + 1;
        //myVector = new T[size];
        myVector=std::vector<T>(maxVertex-minVertex+1,value);
    }
    ShiftedVector<T>(size_t boundary1,size_t boundary2): //inclusive both min and max vertex
    minVertex(std::min(boundary1,boundary2)),maxVertex(std::max(boundary1,boundary2))
    {
        myVector=std::vector<T>(maxVertex-minVertex+1);
    }
    ShiftedVector<T>(){
        minVertex=0;
        maxVertex=0;
        ShiftedVector(minVertex,maxVertex);
    }


    T& operator[](size_t idx){
        if(debug()&&( idx<minVertex||idx>maxVertex)) std::cout<<"out of bounds, index "<<idx<<", interval "<<minVertex<<","<<maxVertex<<std::endl;
        assert(idx>=minVertex&&idx<=maxVertex);
        size_t shiftedIndex=idx-minVertex;
        return myVector[shiftedIndex];
    }

    void setValue(size_t index,T value){
        assert(index>=minVertex&&index<=maxVertex);
        myVector[index-minVertex]=value;

    }

    T getValue(size_t idx)const {
        assert(idx>=minVertex&&idx<=maxVertex);
        size_t shiftedIndex=idx-minVertex;
        return myVector[shiftedIndex];
    }

    const T& operator [](size_t idx) const {
        if(debug()&&( idx<minVertex||idx>maxVertex)) std::cout<<"out of bounds, index "<<idx<<", interval "<<minVertex<<","<<maxVertex<<std::endl;
        assert(idx>=minVertex&&idx<=maxVertex);
        return myVector[idx-minVertex];
    }

    void fillWith(const T& value){
        myVector=std::vector<T>(maxVertex-minVertex+1,value);
    }

    void fillWith(const T& value,size_t index1,size_t index2){ //inclusive both boundary indices
        size_t minIndex=std::min(index1,index2);
        size_t maxIndex=std::max(index1,index2);
        assert(minIndex>=minVertex&&maxIndex<=maxVertex);
        for (size_t i = minIndex; i <= maxIndex; ++i) {
            myVector[i-minVertex]=value;
        }
    }

    bool isWithinBounds(const size_t index)const{
        return index>=minVertex&&index<=maxVertex;

    }

private:
    //T* myVector; // TODO: might be faster depending on whether initializing myVector is expensive
    std::vector<T> myVector;
    size_t minVertex;
    size_t maxVertex;
};

} //end of namespace


#endif /* INCLUDE_LIFTED_DISJOINT_PATHS_LDP_FUNCTIONS_HXX_ */
