#ifndef LDP_CUT_FACTOR_HXX
#define LDP_CUT_FACTOR_HXX

#endif // LDP_CUT_FACTOR_HXX
#include <utility>
#include<cstdlib>
#include<vector>
#include<array>
#include<set>
#include<map>
//#include "ldp_directed_graph.hxx"
#include <config.hxx>
#include "graph_matching/graph_matching_input.h"
#include "graph_matching/min_cost_flow_factor_ssp.hxx"
#include "graph_matching/matching_problem_input.h"
#include "ldp_two_layer_graph.hxx"

namespace LPMP {


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


    void updateCostBase(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value);

    void updateCostLifted(const double& value);

    std::tuple<double, size_t, size_t, char> minCutEdge(size_t index1, size_t index2,const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const;

    double getLiftedMinMarginal(const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const;
    double getOneEdgeMinMarginal(const size_t & index1, const size_t & neighborIndex, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const;


    const LdpTwoLayerGraph& getCutGraph()const {
        return cutGraph;
    }

    const std::vector<size_t>& getInputVertices()const{
        return inputVertices;
    }

    const std::vector<size_t>& getOutputVertices()const{
        return outputVertices;
    }

    const size_t& getNumberOfInputs()const{
        return numberOfInput;
    }

    const size_t& getNumberOfOutputs()const{
        return numberOfOutput;
    }

    const size_t& getLiftedInputVertex() const{
        return v;
    }
    const size_t& getLiftedOutputVertex() const{
        return w;
    }

    const double& getLiftedCost()const{
        return liftedCost;
    }

    const double& getPrimalLiftedCost()const{
        return primalLiftedCost;
    }
    const double& getPrimalBaseCost()const{
        return primalBaseCost;
    }




    void print()const ;
private:
    double advancedMinimizer(const size_t& index1, const size_t& neighborIndex, bool restrictToOne, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost)const;

    LPMP::linear_assignment_problem_input createLAStandard(bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost)const;
    LPMP::linear_assignment_problem_input laExcludeEdge(const size_t& v1,const size_t& v2,bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost)const ;
    LPMP::linear_assignment_problem_input laExcludeVertices(const size_t& v1, const size_t& v2,bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost)const;


//std::vector<double> costs;

std::vector<size_t> inputVertices;
std::vector<size_t> outputVertices;
 LdpTwoLayerGraph cutGraph;
size_t v;
size_t w;
//double liftedEdgeCost;
std::vector<size_t> primalSolution; //inputNodeIndex ->index of neighbor within cutGraph
bool liftedActiveInPrimal;
std::array<size_t,2> baseCoveringLifted; //{input node index,index in output nodes+numberOfInput}
size_t numberOfInput;
size_t numberOfOutput;
bool baseCoverLiftedExists;
size_t unassignedLabel;
mutable std::vector<size_t> storeLabeling;  //index to index
double liftedCost;
mutable bool liftedActive;
mutable double primalBaseCost;
mutable double primalLiftedCost;







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


ldp_cut_factor::ldp_cut_factor(size_t v_, size_t w_, double liftedCost_, std::map<size_t,std::map<size_t,double>> inputEdges):
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
    liftedActiveInPrimal=false;
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
               // std::cout<<"adding edge "<<counterInput<<", "<<i<<std::endl;
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

//Returning neighbors indices, not output indices!
//Last entry for lifted edge, different encoding
void ldp_cut_factor::setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t>& vertexLabels) {
   bool setLiftedActive=(vertexLabels[v]!=0&&vertexLabels[v]==vertexLabels[w]);
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

    //Is it a proper encoding?
    if(setLiftedActive){
        liftedActiveInPrimal=true;
        primalSolution.back()=w;
    }
    else{
        liftedActiveInPrimal=false;
        primalSolution.back()=unassignedLabel;
    }
}


double ldp_cut_factor::EvaluatePrimal() const{
    double value=0;
    for(size_t i=0;i<primalSolution.size()-1;i++){
        if(primalSolution[i]!=unassignedLabel){
            value+=cutGraph.getForwardEdgeCost(i,primalSolution[i]);
        }
    }
    primalBaseCost=value;
    primalLiftedCost=0;

    if(primalSolution.back()==w){
        primalLiftedCost=liftedCost;
        value+=liftedCost;  //meaning lifted edge is active
    }
    return value;
}

void ldp_cut_factor::init_primal(){
    for(size_t& v :primalSolution){
        v=unassignedLabel;
    }
}


const std::vector<size_t>& ldp_cut_factor::getPrimal() {
    return primalSolution;
}


void ldp_cut_factor::print()const {
    std::cout<<"CUT "<<v<<", "<<w<<":"<<liftedCost<<std::endl;
    for (size_t i = 0; i < numberOfInput; ++i) {

        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
             double value=it->second;
            std::cout<<inputVertices[i]<<", "<<outputVertices[outputIndex]<<": "<<value<<std::endl;
        }

    }
}


void ldp_cut_factor::updateCostLifted(const double& value){
    liftedCost+=value;
}



void ldp_cut_factor::updateCostBase(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value){
    //double oldValue=cutGraph.getForwardEdgeCost(inputVertexIndex,neighborIndex);
    cutGraph.updateForwardEdgeCost(inputVertexIndex,neighborIndex,value);
}





double ldp_cut_factor::advancedMinimizer(const size_t& index1, const size_t& index2, bool restrictToOne, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost)const {
    // not restrictToOne .. block only edge
    // restrictToOne ... block both vertices
    //index1 .. index in inputVertices, index2.. index in outputVertices + numberOfInput
    LPMP::linear_assignment_problem_input lapInput;
    double minValue=0;
    std::fill(storeLabeling.begin(),storeLabeling.end(),unassignedLabel);
    liftedActive=false;
    bool activeExists=false;
    const double& localLiftedCost=*pLiftedCost;
    const LdpTwoLayerGraph& localCutGraph=*pCutGraph;
    //TODO: Assert of neighborIndex?
    //TODO add lifted edge to evaluation, based on


    std::tuple<double,size_t,size_t,char> myTuple;
    if(index1!=unassignedLabel&&!restrictToOne){
       if(debug()) std::cout<<"restric to zero"<<std::endl;
        myTuple=minCutEdge(index1,index2,pCutGraph,pLiftedCost);
    }
    else{
        myTuple=minCutEdge(unassignedLabel,unassignedLabel,pCutGraph,pLiftedCost);
    }
    double minCutValue=std::get<0>(myTuple);
    if(minCutValue>=0){
    if(debug())    std::cout<<"simple method"<<std::endl;
        if(restrictToOne){
            if(debug())std::cout<<"restric to one"<<std::endl;
            double valueToReturn=0;
            auto* iter=localCutGraph.forwardNeighborsBegin(index1);
            auto* end=localCutGraph.forwardNeighborsEnd(index1);
            bool found=false;
            for(;iter!=end;iter++){
                if(iter->first==index2){
                    valueToReturn=iter->second;
                    found=true;
                    break;
                }

            }
            storeLabeling[index1]=index2;
            if(localLiftedCost<0||(index1==baseCoveringLifted[0]&&index2==baseCoveringLifted[1])) valueToReturn+=localLiftedCost;
            return valueToReturn;
        }
        else{
            if(debug())std::cout<<"not restric to one"<<std::endl;

            double activeCost=localLiftedCost+minCutValue;
            if(activeCost<0){
                size_t vertex=std::get<1>(myTuple);
                size_t label=std::get<2>(myTuple);
                storeLabeling[vertex]=label;
                liftedActive=true;
                return activeCost;
            }
            else{
                return 0;
            }
        }

    }
    else{
       // std::cout<<"ADVANCED METHOD FOR CUT FACTOR"<<std::endl;
        if(debug())std::cout<<"advanced method"<<std::endl;
        if(index1==unassignedLabel){
            lapInput=createLAStandard(addLiftedCost,pCutGraph,pLiftedCost);
        }
        else if(restrictToOne){
            assert(index1<inputVertices.size());
            lapInput=laExcludeVertices(index1,index2,addLiftedCost,pCutGraph,pLiftedCost);
            //std::cout<<"lap input for one created "<<std::endl;
            //minValue=cutGraph.getForwardEdgeCost(index1,neighborIndex);
            assert(index2<outputVertices.size());
            bool found=false;
            for (auto* it=localCutGraph.forwardNeighborsBegin(index1);it!=localCutGraph.forwardNeighborsEnd(index1);it++) {
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
            assert(index2<outputVertices.size());
            lapInput=laExcludeEdge(index1,index2,addLiftedCost,pCutGraph,pLiftedCost);
            // std::cout<<"lap input for zero created "<<std::endl;
        }
        MCF::SSP<long,double> mcf(lapInput.no_mcf_nodes(),lapInput.no_mcf_edges());
        lapInput.initialize_mcf(mcf);

        //  std::cout<<"mcf init"<<std::endl;
        minValue+=mcf.solve();


        for (int i = 0; i < mcf.no_edges(); ++i) {
            if(mcf.flow(i)>0.99){
                // int label=mcf.head(i)-numberOfInput;
                assert(lapInput.no_left_nodes==numberOfInput);
                int label=mcf.head(i)-numberOfInput;
                int vertex=mcf.tail(i);
                if(vertex<numberOfInput&&label<numberOfOutput){
                    //  std::cout<<"vertex "<<vertex<<", label "<<label<<std::endl;
                    storeLabeling[vertex]=label;
                    activeExists=true;
                }
            }
        }


        if(baseCoverLiftedExists&&storeLabeling[baseCoveringLifted[0]]==baseCoveringLifted[1]){
            liftedActive=true;
            if(!addLiftedCost) minValue+=localLiftedCost;
        }
        else if(activeExists&&localLiftedCost<0){
            minValue+=localLiftedCost;
            liftedActive=true;
        }

        // std::cout<<"lifted active "<<liftedActive<<std::endl;

        return minValue;
    }

}



LPMP::linear_assignment_problem_input   ldp_cut_factor::createLAStandard(bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
            double value=it->second;
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=*pLiftedCost;
            }
            if(value<0){
                assert(outputIndex<outputVertices.size());
               // std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
                lapInput.add_assignment(i,outputIndex,value);
            }
        }

    }
    return lapInput;
}


LPMP::linear_assignment_problem_input   ldp_cut_factor::laExcludeEdge(const size_t& v1, const size_t& v2, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    bool edgeFound=false;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
             double value=it->second;
             if(i==v1&&outputIndex==v2){
                 edgeFound=true;
                 continue;
             }
             if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                 value+=*pLiftedCost;
             }

            if(value<0){
                assert(outputIndex<outputVertices.size());
                lapInput.add_assignment(i,outputIndex,value);
                if(debug())std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
            }
        }

    }
    assert(edgeFound);
    return lapInput;
}


LPMP::linear_assignment_problem_input   ldp_cut_factor::laExcludeVertices(const size_t& v1, const size_t& v2, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        if(i==v1){
            continue;
        }
        auto * it=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
            double value=it->second;
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=*pLiftedCost;
            }
            if(outputIndex==v2) continue;
            if(value<0){
                assert(outputIndex<outputVertices.size());
                lapInput.add_assignment(i,outputIndex,value);
              //  std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
            }
        }


    }
    return lapInput;
}


double ldp_cut_factor::LowerBound() const{
    bool addLifted=(baseCoverLiftedExists&&liftedCost>0);
    double value=advancedMinimizer(unassignedLabel,unassignedLabel,false,addLifted,&cutGraph,&liftedCost);
    return value;
}


double ldp_cut_factor::getOneEdgeMinMarginal(const size_t & index1,const size_t & index2,const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const{
    if(debug())std::cout<<"base edge min marginal in cut"<<std::endl;
    bool addLiftedCost=baseCoverLiftedExists&&liftedCost>0;
    double restrictOne= advancedMinimizer(index1,index2,true,addLiftedCost,pCutGraph,pLiftedCost);
    double restrictZero=advancedMinimizer(index1,index2,false,addLiftedCost,pCutGraph,pLiftedCost);
    return restrictOne-restrictZero;
}


std::tuple<double,size_t,size_t,char> ldp_cut_factor::minCutEdge(size_t index1,size_t index2,const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const{
    double minValue=std::numeric_limits<double>::max();
    size_t v1=0;
    size_t v2=0;
    char cutCoverLiftedNegative=false;
    for(size_t i=0;i<inputVertices.size();i++){
        auto * iter=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;iter!=end;iter++){
            if(iter->second<minValue){
                if(index1!=i||index2!=iter->first){
                    minValue=iter->second;
                    v1=i;
                    v2=iter->first;
                }
            }
            if(baseCoverLiftedExists&&i==baseCoveringLifted[0]&&iter->first==baseCoveringLifted[1]){
                cutCoverLiftedNegative=(iter->second<0);
            }
        }

    }

    std::tuple<double,size_t,size_t,char>t (minValue,v1,v2,cutCoverLiftedNegative);
    return t;
}


double ldp_cut_factor::getLiftedMinMarginal(const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const{
    if(debug())std::cout<<"cut lifted min marginal"<<std::endl;
//    bool nonPositiveEdgeExists=false;
//    bool baseCoverCost=std::numeric_limits<double>::min();
//    double minValue=std::numeric_limits<double>::max();
//    bool minValueSet=false;
//    for(size_t i=0;i<inputVertices.size();i++){
//        auto * iter=cutGraph.forwardNeighborsBegin(i);
//        auto *end=cutGraph.forwardNeighborsEnd(i);
//        for(;iter!=end;iter++){
//            if(iter->second<minValue){
//                minValue=iter->second;
//                minValueSet=true;
//            }
//            if(i==baseCoveringLifted[0]&&iter->first==baseCoveringLifted[1]){
//                baseCoverCost=iter->second;
//            }
//        }
//    }
    double localLiftedCost=*pLiftedCost;
    std::tuple<double,size_t,size_t,char>t=minCutEdge(unassignedLabel,unassignedLabel,pCutGraph,pLiftedCost);
    bool baseCoverNegative=std::get<3>(t)>0;
    double minValue=std::get<0>(t);

    if(!baseCoverLiftedExists||!baseCoverNegative){
        if(minValue>=0){
            //assert(minValueSet);
            if(debug())std::cout<<"min value big"<<std::endl;
            double value=localLiftedCost+minValue;
            if(debug())std::cout<<value<<std::endl;
            return value;
        }
        else{
            if(debug())std::cout<<"min value negative"<<std::endl;
            return localLiftedCost;
        }
    }
    else{
        if(debug())std::cout<<" base cover lifted"<<std::endl;
        double lowerBound=advancedMinimizer(unassignedLabel,unassignedLabel,false,false,pCutGraph,pLiftedCost); //some cut edge must be active
        double restrictedOne=0;
        double restrictedZero=0;
        if(storeLabeling[baseCoveringLifted[0]]==baseCoveringLifted[1]){
            if(debug())std::cout<<"base cover lifted "<<std::endl;
            assert(liftedActive);
            restrictedOne=lowerBound;
            restrictedZero=advancedMinimizer(baseCoveringLifted[0],baseCoveringLifted[1],false,false,pCutGraph,pLiftedCost);
            if(liftedActive) restrictedZero-=localLiftedCost;
            return restrictedOne-restrictedZero;
        }
        else{
            if(debug())std::cout<<"base cover lifted not active"<<std::endl;
            assert(lowerBound<0);
            return localLiftedCost;
        }
    }


 //   return 0;




}


class ldp_snc_cut_message
{
public:
    ldp_snc_cut_message(  std::vector<size_t> _nodeIndicesInCut,  //empty if it is a message only for lifted edge
                          std::vector<size_t> _nodeIndicesInSnc,
                          size_t _sncNodeIDindexInCut,
                          bool _sncIsOut,
                          bool _containsLiftedEdge,
                          size_t _nodeIndexOfLiftedEdge):
        nodeIndicesInCut(_nodeIndicesInCut),  //empty if it is a message only for lifted edge
        nodeIndicesInSnc(_nodeIndicesInSnc),
        sncNodeIDindexInCut(_sncNodeIDindexInCut),
        sncIsOut(_sncIsOut),   //if true, central node is in inputs, other nodes in outputs
        containsLiftedEdge(_containsLiftedEdge),
        nodeIndexOfLiftedEdge(_nodeIndexOfLiftedEdge){

          }

    template<typename CUT_FACTOR>
    void RepamLeft(CUT_FACTOR& l, const double msg, const std::size_t msg_dim) const
    {
        if(debug())std::cout<<"cut repam left ";
        double before=l.LowerBound();
        if(debug())std::cout<<"before lb "<<before<<std::endl;
       if(debug()){

          l.print();
       }
//       printIndices();

      // std::cout<<std::endl;
        assert(msg_dim <=nodeIndicesInCut.size());
       if(msg_dim==nodeIndicesInCut.size()){
           lastV1=l.getLiftedInputVertex();
           lastV2=l.getLiftedOutputVertex();
           lastValue=msg;

           l.updateCostLifted(msg);
       }
       else{
           if(sncIsOut){
               l.updateCostBase(sncNodeIDindexInCut,nodeIndicesInCut.at(msg_dim),msg);
               lastV1=l.getInputVertices().at(sncNodeIDindexInCut);
               size_t otherVertex=l.getCutGraph().getForwardEdgeVertex(sncNodeIDindexInCut,nodeIndicesInCut[msg_dim]);
               lastV2=l.getOutputVertices().at(otherVertex);
               lastValue=msg;
           }
           else{
               l.updateCostBase(nodeIndicesInCut.at(msg_dim),sncNodeIDindexInCut,msg);
           }
       }//std::cout<<"index to update "<<indicesInTriangle.at(msg_dim);
       double after=l.LowerBound();
       if(debug())std::cout<<"after lb "<<after<<std::endl;

     }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {
        if(debug())std::cout<<"snc repam right "<<std::endl;
        double before=r.LowerBound();
        if(debug())std::cout<<"before lb "<<before<<std::endl;
        if(debug()){
            std::cout<<"snc repam right "<<std::endl;
            r.print();
        }
        assert(msg_dim <=nodeIndicesInSnc.size());
        size_t secondVertex;
        if(msg_dim==nodeIndicesInSnc.size()){
            {
                assert(containsLiftedEdge);
                if(debug())std::cout<<"node id of lifted edge "<<r.getLiftedIDs()[nodeIndexOfLiftedEdge]<<std::endl;
                r.updateEdgeCost(msg,nodeIndexOfLiftedEdge,true);
                secondVertex=r.getLiftedIDs().at(nodeIndexOfLiftedEdge);
            }
        }
        else{
                r.updateEdgeCost(msg,nodeIndicesInSnc.at(msg_dim),false);
                secondVertex=r.getBaseIDs().at(nodeIndicesInSnc.at(msg_dim));
        }
        assert(lastV1==r.nodeID);
        if(lastV2!=secondVertex){
            std::cout<<"cut mismatch, first vertex "<<lastV1<<", vertex in cut "<<lastV2<<",  vertex in snc "<<secondVertex<<", index in cut"<<nodeIndicesInCut[msg_dim]<<", is lifted "<<(msg_dim==nodeIndicesInSnc.size())<<std::endl;

        }
//        else{
//            std::cout<<"cut ok, first vertex "<<lastV1<<", vertex in cut "<<lastV2<<",  vertex in snc "<<secondVertex<<", is lifted "<<(msg_dim==nodeIndicesInSnc.size())<<std::endl;
//        }
        assert(lastV2==secondVertex);
        assert(std::abs(lastValue+msg)<eps);
        double after=r.LowerBound();
        if(debug())std::cout<<"abter lb "<<after<<std::endl;

    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
         if(debug())std::cout<<"snc send one to left "<<std::endl;
        if(debug()) std::cout<<"triangle send one to left ";
        std::vector<double> baseCosts=r.getBaseCosts();
        const std::vector<double>& liftedCosts=r.getLiftedCosts();
        size_t i=0;
        double delta=0;
        for (;i<nodeIndicesInSnc.size();i++) {

            delta = r.getOneBaseEdgeMinMarginal(nodeIndicesInSnc[i],&baseCosts,&liftedCosts);
            baseCosts[nodeIndicesInSnc[i]]-=delta;

            msg[i] -= omega * delta;
        }
        if(containsLiftedEdge){
            delta=r.getOneLiftedMinMarginal(nodeIndexOfLiftedEdge,&baseCosts,&liftedCosts);
            msg[i]-=omega * delta;

        }

        //  std::cout<<" done "<<std::endl;
    }

    template<typename CUT_FACTOR, typename MSG>
    void send_message_to_right(const CUT_FACTOR& l, MSG& msg, const double omega)
    {
        if(debug())std::cout<<"cut send one to right"<<std::endl;
        if(debug()) std::cout<<"triangle send one to right"<<std::endl;
        //for (size_t i=0;i<1;i++) {
        double delta;
        size_t i=0;
        double liftedCost=l.getLiftedCost();
        LdpTwoLayerGraph cutGraph=l.getCutGraph();
        for (;i<nodeIndicesInCut.size();i++) {
            if(sncIsOut){
                size_t otherVertex=l.getCutGraph().getForwardEdgeVertex(sncNodeIDindexInCut,nodeIndicesInCut[i]);
                delta = l.getOneEdgeMinMarginal(sncNodeIDindexInCut,otherVertex,&cutGraph,&liftedCost);
                cutGraph.updateForwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i],-delta);
            }
            else{  //just for output nodes now
                assert(false);
                delta = l.getOneEdgeMinMarginal(nodeIndicesInCut[i],sncNodeIDindexInCut,&cutGraph,&liftedCost);
                cutGraph.updateForwardEdgeCost(nodeIndicesInCut[i],sncNodeIDindexInCut,-delta);
            }
            if(debug())std::cout<<"min marginal value "<<delta<<std::endl;
            msg[i] -= omega * delta;
        }

        if(containsLiftedEdge){
            if(debug())std::cout<<"send to lifted "<<l.getLiftedOutputVertex()<<std::endl;
            delta=l.getLiftedMinMarginal(&cutGraph,&liftedCost);
            msg[i]-=omega*delta;
        }
    }





    template<typename SINGLE_NODE_CUT_FACTOR,typename CUT_FACTOR>
    bool check_primal_consistency(const CUT_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
    {
        return true;

    }

private:


    const std::vector<size_t> nodeIndicesInCut;  //empty if it is a message only for lifted edge
    const std::vector<size_t> nodeIndicesInSnc;
    const size_t sncNodeIDindexInCut;
    const bool sncIsOut;   //if true, central node is in inputs, other nodes in outputs
    const bool containsLiftedEdge;
    const size_t nodeIndexOfLiftedEdge;


    mutable size_t lastV1;
    mutable size_t lastV2;
    mutable double lastValue;
};










}
