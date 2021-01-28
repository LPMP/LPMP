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
#include <utility>
#include <config.hxx>






namespace LPMP{

template<class SOLVER,class INSTANCE>
void constructProblemFromSolver(SOLVER& solver,const INSTANCE& instance){
    solver.GetProblemConstructor().construct(instance);
}



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

//template<class T>
//inline bool edgeCompare(const std::tuple<T,size_t,size_t,bool>& t1,const std::tuple<T,size_t,size_t,bool>& t2) {
//     if(std::get<0>(t1)==std::get<0>(t2)){
//         if(std::get<1>(t1)==std::get<1>(t2)){
//             if(std::get<2>(t1)==std::get<2>(t2)){
//                 if(std::get<3>(t1)==std::get<3>(t2)){
//                     std::cout<<std::get<0>(t1)<<" ";
//                     std::cout<<std::get<1>(t1)<<" ";
//                     std::cout<<std::get<2>(t1)<<" ";
//                     std::cout<<std::get<3>(t1)<<" ";
//                     std::cout<<std::endl;
//                     throw std::runtime_error("duplicated edge in edge list");
//                 }
//                 else {
//                     return std::get<3>(t1)>std::get<3>(t2); //lifted first
//                 }

//             }
//             else{
//                 return std::get<2>(t1)<std::get<2>(t2);
//             }

//         }
//         else{
//             return std::get<1>(t1)<std::get<1>(t2);
//         }
//     }
//     else{
//         return std::get<0>(t1)<std::get<0>(t2);
//     }
// }


// template<class T>
// inline bool queueCompare(const std::pair<double,T>& t1,const std::pair<double,T>& t2) {
//      if(abs(t1.first-t2.first)>=eps){
//          return std::get<0>(t1)<std::get<0>(t2);
//      }
//      else{
//          return false;
//      }
//  }

// template <class T>
// class ForQueueComparator{
// public:
//     bool operator()(const std::tuple<double,size_t,T>& lhs, const  std::tuple<double,size_t,T>& rhs) {
//         if(std::abs(std::get<0>(rhs)-std::get<0>(lhs))<eps){
//             return std::get<1>(rhs) < std::get<1>(lhs);
//         }
//         else{
//             return std::get<0>(lhs)<std::get<0>(rhs);
//         }
//     }
// };




template<class T>
inline bool edgeCompare(const std::tuple<T,size_t,size_t,bool>& t1,const std::tuple<T,size_t,size_t,bool>& t2) {
     if(abs(std::get<0>(t1)-std::get<0>(t2))<eps){
         if(std::get<1>(t1)==std::get<1>(t2)){
             if(std::get<2>(t1)==std::get<2>(t2)){
                 if(std::get<3>(t1)==std::get<3>(t2)){
                     std::cout<<std::get<0>(t1)<<" ";
                     std::cout<<std::get<1>(t1)<<" ";
                     std::cout<<std::get<2>(t1)<<" ";
                     std::cout<<std::get<3>(t1)<<" ";
                     std::cout<<std::endl;
                     throw std::runtime_error("duplicated edge in edge list");
                 }
                 else {
                     return std::get<3>(t1)>std::get<3>(t2); //lifted first
                 }

             }
             else{
                 return std::get<2>(t1)<std::get<2>(t2);
             }

         }
         else{
             return std::get<1>(t1)<std::get<1>(t2);
         }
     }
     else{
         return std::get<0>(t1)<std::get<0>(t2);
     }
 }

template<class T>
inline bool baseEdgeCompare(const std::tuple<T,size_t,size_t>& t1,const std::tuple<T,size_t,size_t>& t2) {
    //if(std::get<0>(t1)==std::get<0>(t2)){
     if(abs(std::get<0>(t1)-std::get<0>(t2))<eps){
        if(std::get<1>(t1)==std::get<1>(t2)){
            if(std::get<2>(t1)==std::get<2>(t2)){
                throw std::runtime_error("duplicated edge in list");
            }
            else{
                return std::get<2>(t1)<std::get<2>(t2);
            }

        }
        else{
            return std::get<1>(t1)<std::get<1>(t2);
        }
    }
    else{
        return std::get<0>(t1)<std::get<0>(t2);
    }
 }

}  //end of namespace




template<class T=size_t>
struct EdgeVector{
public:
    EdgeVector(std::vector<std::array<T,2>>& inputVector):
        edges(inputVector)
    {}

    T& operator()(size_t i,size_t j){
        assert(j<2);
        assert(i<edges.size());
        return edges[i][j];
    }

    const T& operator()(size_t i,size_t j)const{
        assert(j<2);
        assert(i<edges.size());
        return edges[i][j];
    }

    const size_t shape(size_t i)const{
        assert(i<=1);
        if(i==0){
            return edges.size();
        }
        if(i==1){
            return 2;
        }
    }

private:
    std::vector<std::array<T,2>>& edges;

};

template<class T=size_t>
struct InfoVector{
public:
    InfoVector(std::vector<T>& inputVector):
        myVector(inputVector)
    {}

    T& operator()(size_t i){
        assert(i<myVector.size());
        return myVector[i];
    }

    const T& operator()(size_t i)const {
        assert(i<myVector.size());
        return myVector[i];
    }

    const size_t shape(size_t i)const{
        assert(i==0);
        return myVector.size();

    }

private:
    std::vector<T>& myVector;

};



template<typename T>
void fillWithValue(std::vector<T>& myVector,size_t first,size_t last, T value){
    for(size_t i=first;i<last;i++){
        myVector[i]=value;
    }
}



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
    ShiftedVector<T>(size_t boundary1,size_t boundary2,const std::vector<T>& sourceVector): //inclusive both min and max vertex
    minVertex(std::min(boundary1,boundary2)),maxVertex(std::max(boundary1,boundary2))
    {
        assert(sourceVector.size()>maxVertex);
        auto itBegin=sourceVector.begin();
        itBegin+=minVertex;
        auto itEnd=sourceVector.begin();
        itEnd+=maxVertex+1;
        myVector=std::vector<T>(itBegin,itEnd);
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
       // if(debug()&&( idx<minVertex||idx>maxVertex)) std::cout<<"out of bounds, index "<<idx<<", interval "<<minVertex<<","<<maxVertex<<std::endl;
        assert(idx>=minVertex&&idx<=maxVertex);
        size_t shiftedIndex=idx-minVertex;
        return myVector[shiftedIndex];
    }

//    void setValue(size_t index,T value){
//        assert(index>=minVertex&&index<=maxVertex);
//        myVector[index-minVertex]=value;

//    }

//    T getValue(size_t idx)const {
//        assert(idx>=minVertex&&idx<=maxVertex);
//        size_t shiftedIndex=idx-minVertex;
//        return myVector[shiftedIndex];
//    }

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
