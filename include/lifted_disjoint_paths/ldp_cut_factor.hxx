#ifndef LDP_CUT_FACTOR_HXX
#define LDP_CUT_FACTOR_HXX

#endif // LDP_CUT_FACTOR_HXX
#include <utility>
#include<cstdlib>
#include<vector>
#include<array>
#include<set>
#include<map>
#include "ldp_directed_graph.hxx"
#include <config.hxx>
#include "graph_matching/graph_matching_input.h"
#include "graph_matching/min_cost_flow_factor_ssp.hxx"
#include "graph_matching/matching_problem_input.h"
#include "ldp_two_layer_graph.hxx"

namespace LPMP {

template<class INSTANCE>
class ldp_cut_factor
{
public:

    ldp_cut_factor(size_t v_,size_t w_, double liftedCost_,std::map<size_t,std::map<size_t,double>> inputEdges);
   // ldp_cut_factor(size_t v_, size_t w_, const std::vector<std::array<size_t,2>>& edgesInDirectedGraph,const std::vector<double>& edgeCosts,const INSTANCE& ldpInstance);
    double LowerBound() const;


    double EvaluatePrimal() const;

    void setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t> &vertexLabels);

    const std::vector<size_t>& getPrimal();

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

    auto export_variables() { return std::tie(); }

    void init_primal();

   // std::array<double,3> getAllMinMarginals();

    double getOneMinMarginal(size_t edgeId) const;

    void updateCostBase(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value);

    double advancedMinimizer(const size_t& index1, const size_t& neighborIndex, bool restrictToOne)const;

    LPMP::linear_assignment_problem_input createLAStandard()const;
    LPMP::linear_assignment_problem_input laExcludeEdge(const size_t& v1,const size_t& v2)const ;
    LPMP::linear_assignment_problem_input laExcludeVertices(const size_t& v1, const size_t& v2)const;

    double getOneEdgeMinMarginal(const size_t & index1,const size_t & neighborIndex) const;

    LdpTwoLayerGraph cutGraph;

private:
//std::vector<double> costs;

std::vector<size_t> inputVertices;
std::vector<size_t> outputVertices;
size_t v;
size_t w;
//double liftedEdgeCost;
std::vector<size_t> primalSolution; //inputNodeIndex ->index of neighbor within cutGraph
std::array<size_t,2> baseCoveringLifted; //{input node index,index in output nodes+numberOfInput}
size_t numberOfInput;
size_t numberOfOutput;
bool baseCoverLiftedExists;
size_t unassignedLabel;
mutable std::vector<size_t> storeLabeling;
double liftedCost;
mutable bool liftedActive;






};

//template<class INSTANCE>
//inline ldp_cut_factor<INSTANCE>::ldp_cut_factor(size_t v_, size_t w_, const std::vector<std::array<size_t,2>>& edgesInDirectedGraph,const std::vector<double>& edgeCosts,const INSTANCE& ldpInstance):
//    v(v_),
//    w(w_)
//{ //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
//   // std::set<size_t> outVertices;
//    size_t edgeCounter=0;
//    const LdpDirectedGraph& baseGraph=ldpInstance.getMyGraph();
//    //TODO leave out later if lexicographical order of edges is ensured
//    std::set<size_t> inputVertices;
//    std::set<size_t> outputVertices;
//    for(size_t i=0;i<edgesInDirectedGraph.size();i++){
//        const std::array<size_t,2>& edge=edgesInDirectedGraph[i];
//        inputVertices.insert(edge[0]);
//        outputVertices.insert(edge[1]);
//    }



//}

template<class INSTANCE>
inline ldp_cut_factor<INSTANCE>::ldp_cut_factor(size_t v_, size_t w_, double liftedCost_, std::map<size_t,std::map<size_t,double>> inputEdges):
    v(v_),
    w(w_),
    liftedCost(liftedCost_)
{ //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
    std::set<size_t> outVertices;
    numberOfInput=inputEdges.size();
    unassignedLabel=std::numeric_limits<size_t>::max();
    baseCoveringLifted={unassignedLabel,unassignedLabel};
    baseCoverLiftedExists=0;
    liftedActive=false;
    size_t edgeCounter=0;
    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
        size_t v1=iter->first;
        inputVertices.push_back(v1);
        std::map<size_t,double>& neighbors=iter->second;
        for(auto iter2=neighbors.begin();iter2!=neighbors.end();iter2++){
            outVertices.insert(iter2->first);
            edgeCounter++;
        }
    }
    for(size_t vert:outVertices){
        outputVertices.push_back(vert);
    }
    numberOfOutput=outVertices.size();
    size_t counterInput=0;
    std::vector<double> edgeCosts;

    std::vector<std::array<size_t,2>> localEdges(edgeCounter);
    edgeCounter=0;

    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
        size_t v1=iter->first;
        std::map<size_t,double>& neighbors=iter->second;
        for(size_t i=0;i<outputVertices.size();i++){
            auto it=neighbors.find(outputVertices[i]);
            if(it!=neighbors.end()){
                if(v1==v&&it->first==w){
                    baseCoveringLifted={counterInput,i};
                    baseCoverLiftedExists=1;
                }
                assert(edgeCounter<localEdges.size());
                localEdges[edgeCounter]={counterInput,i};

                edgeCounter++;
                std::cout<<"adding edge "<<counterInput<<", "<<i<<std::endl;
                edgeCosts.push_back(it->second);
            }
        }
        counterInput++;
    }

    cutGraph=LdpTwoLayerGraph(localEdges,edgeCosts);
    primalSolution=std::vector<size_t>(inputVertices.size()+1,unassignedLabel); //plus one for lifted edge

//    if(diagnostics()){
//        size_t edgeCounter=0;
//        for (size_t v=0;v<cutGraph.getNumberOfVertices();v++) {
//            assert(localEdges[edgeCounter][0]==v);
//            const auto * it=cutGraph.forwardNeighborsBegin(v);
//            const auto * end=cutGraph.forwardNeighborsEnd(v);
//            for(;it!=end;it++){
//                assert(edgeCounter<localEdges.size());
//                assert(localEdges[edgeCounter][1]==it->first);
//                edgeCounter++;
//            }
//            if(edgeCounter>=localEdges.size()) break;
//        }

//    }
    storeLabeling=std::vector<size_t>(numberOfInput+1,unassignedLabel);

}

//ldp_cut_factor::ldp_cut_factor(size_t v_, size_t w_, double liftedCost, std::map<size_t,std::map<size_t,double>> inputEdges):
//    v(v_),
//    w(w_)
//{ //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
//    std::set<size_t> outVertices;
//    size_t edgeCounter=0;
//    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
//        size_t v1=iter->first;
//        inputVertices.push_back(v1);
//        std::map<size_t,double>& neighbors=iter->second;
//        for(auto iter2=neighbors.begin();iter2!=neighbors.end();iter2++){
//            outVertices.insert(iter2->first);
//            edgeCounter++;
//        }
//    }
//    for(size_t vert:outVertices){
//        outputVertices.push_back(vert);
//    }
//    size_t counterInput=0;

//    baseCoveringLifted=edgeCounter; //meaning there is none

//    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
//        size_t v1=iter->first;
//        std::map<size_t,double>& neighbors=iter->second;
//        for(size_t i=0;i<outputVertices.size();i++){
//            auto it=neighbors.find(outputVertices[i]);
//            if(it!=neighbors.end()){
//                if(v1==v&&it->first==w){
//                    baseCoveringLifted=edges.size();
//                }
//                edges.push_back({counterInput,i});
//                costs.push_back(it->second);
//            }
//        }
//        counterInput++;
//    }
//    costs.push_back(liftedCost);
//    primalSolution=std::vector<bool>(edges.size()+1);
//    assert(edgeCounter==edges.size());

//}


template<class INSTANCE>
inline void ldp_cut_factor<INSTANCE>::setPrimal(const std::vector<size_t>& primalDescendants,const std::vector<size_t>& vertexLabels) {
    //TODO: missing lifted edge label!
    for(size_t i=0;i<inputVertices.size();i++){
        const size_t& vertexID=inputVertices[i];
        primalSolution[i]=unassignedLabel;
        const auto* it=cutGraph.forwardNeighborsBegin(i);
        const auto* end=cutGraph.forwardNeighborsEnd(i);
        size_t counter=0;
        for(;it!=end;it++){
            size_t index=it->first;

            assert(index>=0&&index<outputVertices.size());
            size_t neighborID=outputVertices[index];
            if(neighborID==primalDescendants[vertexID]){
                assert(i<primalSolution.size());
                primalSolution[i]=counter;
            }
            counter++;
        }
    }
    assert(v<vertexLabels.size()&&w<vertexLabels.size());
    if(vertexLabels[v]==vertexLabels[w]){
        primalSolution.back()=w;
    }
    else{
        primalSolution.back()=unassignedLabel;
    }

}

template <class INSTANCE>
inline void ldp_cut_factor<INSTANCE>::init_primal(){
    for(size_t& v :primalSolution){
        v=unassignedLabel;
    }
}


template<class INSTANCE>
inline const std::vector<size_t>& ldp_cut_factor<INSTANCE>::getPrimal() {
    return primalSolution;
}

template<class INSTANCE>
inline void ldp_cut_factor<INSTANCE>::updateCostBase(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value){
    double oldValue=cutGraph.getForwardEdgeCost(inputVertexIndex,neighborIndex);
    cutGraph.setForwardEdgeCost(inputVertexIndex,neighborIndex,oldValue+value);
}


template<class INSTANCE>
inline double ldp_cut_factor<INSTANCE>::EvaluatePrimal() const{
    double value=0;
    for(size_t i=0;i<primalSolution.size()-1;i++){
        if(primalSolution[i]!=unassignedLabel){
            value+=cutGraph.getForwardEdgeCost(i,primalSolution[i]);
        }
    }
    return value;
}


template<class INSTANCE>
inline double ldp_cut_factor<INSTANCE>::advancedMinimizer(const size_t& index1,const size_t& index2,bool restrictToOne)const {
    // not restrictToOne .. block only edge
    // restrictToOne ... block both vertices
    //index1 .. index in inputVertices, index2.. index in outputVertices + numberOfInput
    LPMP::linear_assignment_problem_input lapInput;
    double minValue=0;
   std::fill(storeLabeling.begin(),storeLabeling.end(),unassignedLabel);
   bool activeExists=false;
    //TODO: Assert of neighborIndex?
   //TODO add lifted edge to evaluation, based on
    if(index1==unassignedLabel){
        lapInput=createLAStandard();
    }
    else if(restrictToOne){
        assert(index1<inputVertices.size());
        lapInput=laExcludeVertices(index1,index2);
        std::cout<<"lap input for one created "<<std::endl;
        //minValue=cutGraph.getForwardEdgeCost(index1,neighborIndex);
        assert(index2<outputVertices.size());
        bool found=false;
        for (auto* it=cutGraph.forwardNeighborsBegin(index1);it!=cutGraph.forwardNeighborsEnd(index1);it++) {
            if(it->first==index2){
                minValue=it->second;
                found=true;
                break;
            }
        }
        assert(found);
        storeLabeling[index1]=index2;
        activeExists=true;

    }
    else{//edge is restricted to zero
        assert(index1<inputVertices.size());
        lapInput=laExcludeEdge(index1,index2);
        std::cout<<"lap input for zero created "<<std::endl;
    }
    MCF::SSP<long,double> mcf(lapInput.no_mcf_nodes(),lapInput.no_mcf_edges());
    lapInput.initialize_mcf(mcf);

    std::cout<<"mcf init"<<std::endl;
    minValue+=mcf.solve();


    for (int i = 0; i < mcf.no_edges(); ++i) {
        if(mcf.flow(i)>0.99){
           // int label=mcf.head(i)-numberOfInput;
            assert(lapInput.no_left_nodes==numberOfInput);
             int label=mcf.head(i)-numberOfInput;
            int vertex=mcf.tail(i);
            if(vertex<numberOfInput&&label<numberOfOutput){
                std::cout<<"vertex "<<vertex<<", label "<<label<<std::endl;
                storeLabeling[vertex]=label;
                activeExists=true;
            }
        }
    }

    liftedActive=false;
    if(baseCoverLiftedExists&&storeLabeling[baseCoveringLifted[0]]==baseCoveringLifted[1]){
        liftedActive=true; //cost already present in LAP
        std::cout<<"lifted forced one"<<std::endl;
    }
    else if(activeExists&&liftedCost<0){
        minValue+=liftedCost;
        liftedActive=true;
    }

    std::cout<<"lifted active "<<liftedActive<<std::endl;

    return minValue;

}





template<class INSTANCE>
inline LPMP::linear_assignment_problem_input   ldp_cut_factor<INSTANCE>::createLAStandard() const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
            double value=it->second;
            if(baseCoverLiftedExists&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=liftedCost;
            }
            if(value<0){
                assert(outputIndex<outputVertices.size());
                std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
                lapInput.add_assignment(i,outputIndex,value);
            }
        }

    }
    return lapInput;
}


template<class INSTANCE>
inline LPMP::linear_assignment_problem_input   ldp_cut_factor<INSTANCE>::laExcludeEdge(const size_t& v1,const size_t& v2) const{
    LPMP::linear_assignment_problem_input lapInput;
    bool edgeFound=false;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);        
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
             double value=it->second;
             if(baseCoverLiftedExists&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                 value+=liftedCost;
             }
             if(i==v1&&outputIndex==v2){
                 edgeFound=true;
                 continue;
             }
            if(value<0){
                assert(outputIndex<outputVertices.size());
                lapInput.add_assignment(i,outputIndex,value);
                std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
            }
        }

    }
    assert(edgeFound);
    return lapInput;
}

template<class INSTANCE>
inline LPMP::linear_assignment_problem_input   ldp_cut_factor<INSTANCE>::laExcludeVertices(const size_t& v1, const size_t& v2) const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        if(i==v1){
            continue;
        }
        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
            double value=it->second;
            if(baseCoverLiftedExists&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=liftedCost;
            }
            if(outputIndex==v2) continue;
            if(value<0){
                assert(outputIndex<outputVertices.size());
                lapInput.add_assignment(i,outputIndex,value);
                std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
            }
        }


    }
    return lapInput;
}

template<class INSTANCE>
inline double ldp_cut_factor<INSTANCE>::LowerBound() const{
    double value=advancedMinimizer(unassignedLabel,unassignedLabel,false);
    return value;
}

template<class INSTANCE>
inline double ldp_cut_factor<INSTANCE>::getOneEdgeMinMarginal(const size_t & index1,const size_t & index2) const{
    double restrictOne= advancedMinimizer(index1,index2,true);
    double restrictZero=advancedMinimizer(index1,index2,false);
    return restrictOne-restrictZero;
}



//double ldp_cut_factor::LowerBound() const{
//    size_t liftedIndex=edges.size();
//    double lb=0;
//    if(costs.at(liftedIndex)>=0){
//        std::vector<bool> usedInput(inputVertices.size(),0);
//        std::vector<bool> usedOutput(outputVertices.size(),0);
//        bool canUseSimple=true;
//        for(size_t i=0;i<costs.size()-1;i++){
//            if(costs[i]<0){
//                if(i==baseCoveringLifted){
//                    if(costs[i]<-costs[liftedIndex]){
//                        lb+=costs[i]+costs[liftedIndex];
//                        if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
//                            canUseSimple=false;
//                            break;
//                        }
//                        else{
//                            usedInput.at(edges[i][0])=1;
//                            usedOutput.at(edges[i][1])=1;
//                        }
//                    }

//                }
//                else{
//                    lb+=costs[i];
//                    if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
//                        canUseSimple=false;
//                        break;
//                    }
//                    else{
//                        usedInput.at(edges[i][0])=1;
//                        usedOutput.at(edges[i][1])=1;
//                    }

//                }
//            }
//        }
//        if(canUseSimple){
//            return lb;
//        }
//        else{
//            //TODO solve min cost assignment
//        }
//    }
//    else{
//        std::vector<bool> usedInput(inputVertices.size(),0);
//        std::vector<bool> usedOutput(outputVertices.size(),0);
//        bool canUseSimple=true;
//        double minValue=std::numeric_limits<double>::max();
//        for(size_t i=0;i<costs.size()-1;i++){
//            minValue=std::min(minValue,costs[i]);
//            if(costs[i]<0){
//                lb+=costs[i];
//                if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
//                    canUseSimple=false;
//                    break;
//                }
//                else{
//                    usedInput.at(edges[i][0])=1;
//                    usedOutput.at(edges[i][1])=1;
//                }
//            }
//        }
//        if(canUseSimple){
//            if(minValue<0){
//                lb+=costs[liftedIndex];
//            }
//            else if(costs[liftedIndex]<-minValue){
//                lb=costs[liftedIndex]+minValue;
//            }
//            return lb;
//        }
//        else{
//            //TODO solve min cost assignment
//        }

//    }

//}


//double ldp_cut_factor::getOneMinMarginal(size_t edgeId) const{
//    assert(edgeId<=costs.size());
//    size_t liftedIndex=edges.size();
//    std::vector<bool> usedInput(inputVertices.size(),0);
//    std::vector<bool> usedOutput(outputVertices.size(),0);
//    //std::vector<bool> activeEdges(edges.size(),0); //maybe change to set of indices
//    std::set<size_t> activeEdges;
//    double lb=0;
//    bool canUseSimple=true;
//    for(size_t i=0;i<costs.size()-1;i++){  //TODO maybe make one method to be used in LowerBound too!
//        if(costs[i]<0){
//            lb+=costs[i];
//            //activeEdges[i]=1;
//            activeEdges.insert(i);
//            if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
//                canUseSimple=false;
//                break;
//            }
//            else{
//                usedInput.at(edges[i][0])=1;
//                usedOutput.at(edges[i][1])=1;
//            }
//        }
//    }
//    if(canUseSimple){

//        double inactiveLiftedCost=0;
//        double activeLiftedCost=0;
//        if(costs[baseCoveringLifted]>=0){
//            inactiveLiftedCost=lb;
//        }
//        else{
//            inactiveLiftedCost=lb-costs[baseCoveringLifted]; //maybe remove vertices from used input and used output
//        }
//        if(lb<0){
//            activeLiftedCost=lb+costs[liftedIndex];
//        }
//        else{
//            double minValue=costs[0];
//            for (size_t i = 1; i < costs.size()-1; ++i) {
//                minValue=std::min(minValue,costs[i]);
//            }
//            activeLiftedCost=minValue+costs[liftedIndex];
//        }
//        //TODO up to here use for lower bound, return lower from the two
//        if(edgeId==liftedIndex){
//            return activeLiftedCost-inactiveLiftedCost;
//        }
//        else if(edgeId!=baseCoveringLifted){
//            double optValue=inactiveLiftedCost;
//            bool liftedActiveInOpt=false;
//            if(activeLiftedCost<inactiveLiftedCost){
//                optValue=activeLiftedCost;
//                liftedActiveInOpt=true;
//            }
//            double inactiveCost=0;
//            double activeCost=0;
//            if(costs[edgeId]>=0){
//                inactiveCost=optValue;
//                activeCost=costs[edgeId];
//                size_t blockV1=edges[edgeId][0];
//                size_t blockV2=edges[edgeId][1];
//                for(size_t i=0;i<costs.size()-1;i++){
//                    if(costs[i]<0&&edges[i][0]!=blockV1&&edges[i][1]!=blockV2){
//                        activeCost+=costs[i];
//                    }
//                }
//                if(costs[liftedIndex]<0){
//                    activeCost+=costs[liftedIndex];
//                }




//            }
//            else{
//                activeCost=optValue;
//            }


//        }

//    }
//    else{
//        //TODO use min cost assignment
//    }

//    if(edgeId==liftedIndex){

//    }
//}


//    ldp_cut_factor(size_t v_,size_t w_, std::vector<std::array<size_t,2>> inputEdges,std::vector<double> inputCosts) { //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
//        std::set<size_t> inpVertices;
//        std::set<size_t> outVertices;
//        for (size_t i = 0; i < inputEdges.size(); ++i) {
//            size_t v1=inputEdges[i][0];
//            size_t v2=inputEdges[i][1];
//            if(inpVertices.count(v1)==0){
//                inpVertices.insert(v1);
//            }

//            if(outVertices.count(v2)==0){
//                outVertices.insert(v2);
//            }

//        }
//        //inputVertices=std::vector<size_t>(inpVertices.size());
//        for(size_t vert:inpVertices){
//            inputVertices.push_back(vert);
//        }
//        for(size_t vert:outVertices){
//            outputVertices.push_back(vert);
//        }

//    }



}
