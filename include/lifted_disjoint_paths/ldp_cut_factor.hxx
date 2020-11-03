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

    std::tuple<double, size_t, size_t, char> minCutEdge(size_t index1, size_t index2) const;

    double getLiftedMinMarginal() const;
    double getOneEdgeMinMarginal(const size_t & index1,const size_t & neighborIndex) const;


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

    void print()const ;
private:
    double advancedMinimizer(const size_t& index1, const size_t& neighborIndex, bool restrictToOne,bool addLiftedCost=false)const;

    LPMP::linear_assignment_problem_input createLAStandard(bool addLiftedCost=false)const;
    LPMP::linear_assignment_problem_input laExcludeEdge(const size_t& v1,const size_t& v2,bool addLiftedCost=false)const ;
    LPMP::linear_assignment_problem_input laExcludeVertices(const size_t& v1, const size_t& v2,bool addLiftedCost=false)const;


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
    double oldValue=cutGraph.getForwardEdgeCost(inputVertexIndex,neighborIndex);
    cutGraph.setForwardEdgeCost(inputVertexIndex,neighborIndex,oldValue+value);
}


double ldp_cut_factor::EvaluatePrimal() const{
    double value=0;
    for(size_t i=0;i<primalSolution.size()-1;i++){
        if(primalSolution[i]!=unassignedLabel){
            value+=cutGraph.getForwardEdgeCost(i,primalSolution[i]);
        }
    }
    if(primalSolution.back()==w) value+=liftedCost;  //meaning lifted edge is active
    return value;
}


double ldp_cut_factor::advancedMinimizer(const size_t& index1, const size_t& index2, bool restrictToOne, bool addLiftedCost)const {
    // not restrictToOne .. block only edge
    // restrictToOne ... block both vertices
    //index1 .. index in inputVertices, index2.. index in outputVertices + numberOfInput
    LPMP::linear_assignment_problem_input lapInput;
    double minValue=0;
   std::fill(storeLabeling.begin(),storeLabeling.end(),unassignedLabel);
   liftedActive=false;
   bool activeExists=false;
    //TODO: Assert of neighborIndex?
   //TODO add lifted edge to evaluation, based on


   std::tuple<double,size_t,size_t,char> myTuple(0,0,0,0);
   if(index1!=unassignedLabel&&!restrictToOne){
       myTuple=minCutEdge(index1,index2);
   }
   else{
       myTuple=minCutEdge(unassignedLabel,unassignedLabel);
   }
   double minCutValue=std::get<0>(myTuple);
   if(minCutValue>=0){
       if(restrictToOne){
           double valueToReturn=0;
           auto* iter=cutGraph.forwardNeighborsBegin(index1);
           auto* end=cutGraph.forwardNeighborsEnd(index1);
           bool found=false;
           for(;iter!=end;iter++){
               if(iter->first==index2){
                   valueToReturn=iter->second;
                   found=true;
                   break;
               }

           }
           storeLabeling[index1]=index2;
           if(liftedCost<0) valueToReturn+=liftedCost;
           return valueToReturn;
       }
       else{

           double activeCost=liftedCost+minCutValue;
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

       if(index1==unassignedLabel){
           lapInput=createLAStandard(addLiftedCost);
    }
    else if(restrictToOne){
        assert(index1<inputVertices.size());
        lapInput=laExcludeVertices(index1,index2,addLiftedCost);
        //std::cout<<"lap input for one created "<<std::endl;
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
        assert(index2<outputVertices.size());
        lapInput=laExcludeEdge(index1,index2,addLiftedCost);
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
        if(!addLiftedCost) minValue+=liftedCost;
    }
    else if(activeExists&&liftedCost<0){
        minValue+=liftedCost;
        liftedActive=true;
    }

   // std::cout<<"lifted active "<<liftedActive<<std::endl;

    return minValue;
   }

}



LPMP::linear_assignment_problem_input   ldp_cut_factor::createLAStandard(bool addLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
            double value=it->second;
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=liftedCost;
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


LPMP::linear_assignment_problem_input   ldp_cut_factor::laExcludeEdge(const size_t& v1, const size_t& v2, bool addLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    bool edgeFound=false;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);        
        for(;it!=end;it++){
            const size_t& outputIndex=it->first;
             double value=it->second;
             if(i==v1&&outputIndex==v2){
                 edgeFound=true;
                 continue;
             }
             if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                 value+=liftedCost;
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


LPMP::linear_assignment_problem_input   ldp_cut_factor::laExcludeVertices(const size_t& v1, const size_t& v2, bool addLiftedCost) const{
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
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=liftedCost;
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
    double value=advancedMinimizer(unassignedLabel,unassignedLabel,false,addLifted);
    return value;
}


double ldp_cut_factor::getOneEdgeMinMarginal(const size_t & index1,const size_t & index2) const{
    bool addLiftedCost=baseCoverLiftedExists&&liftedCost>0;
    double restrictOne= advancedMinimizer(index1,index2,true,addLiftedCost);
    double restrictZero=advancedMinimizer(index1,index2,false,addLiftedCost);
    return restrictOne-restrictZero;
}


std::tuple<double,size_t,size_t,char> ldp_cut_factor::minCutEdge(size_t index1,size_t index2) const{
    double minValue=std::numeric_limits<double>::max();
    size_t v1=0;
    size_t v2=0;
    char cutCoverLiftedNegative=false;
    for(size_t i=0;i<inputVertices.size();i++){
        auto * iter=cutGraph.forwardNeighborsBegin(i);
        auto *end=cutGraph.forwardNeighborsEnd(i);
        for(;iter!=end;iter++){
            if(iter->second<minValue){
                if(index1!=i||index2!=iter->first){
                    minValue=iter->second;
                    v1=i;
                    v2=iter->first;
                }
            }
        }
        if(baseCoverLiftedExists&&i==baseCoveringLifted[0]&&iter->first==baseCoveringLifted[1]){
            cutCoverLiftedNegative=iter->second<0;
        }
    }

    std::tuple<double,size_t,size_t,char>t (minValue,v1,v2,cutCoverLiftedNegative);
    return t;
}


double ldp_cut_factor::getLiftedMinMarginal() const{
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
    std::tuple<double,size_t,size_t,char>t=minCutEdge(unassignedLabel,unassignedLabel);
    bool baseCoverNegative=std::get<3>(t)>0;
    double minValue=std::get<0>(t);

    if(!baseCoverLiftedExists||!baseCoverNegative){
        if(minValue>=0){
            //assert(minValueSet);
            std::cout<<"min value big"<<std::endl;
            double value=liftedCost+minValue;
            std::cout<<value<<std::endl;
            return value;
        }
        else{
            return liftedCost;
        }
    }
    else{
        double lowerBound=advancedMinimizer(unassignedLabel,unassignedLabel,false,false); //some cut edge must be active
        double restrictedOne=0;
        double restrictedZero=0;
        if(storeLabeling[baseCoveringLifted[0]]==baseCoveringLifted[1]){
            std::cout<<"base cover lifted "<<std::endl;
            assert(liftedActive);
            restrictedOne=lowerBound;
            restrictedZero=advancedMinimizer(baseCoveringLifted[0],baseCoveringLifted[1],false,false);
            if(liftedActive) restrictedZero-=liftedCost;
            return restrictedOne-restrictedZero;
        }
        else{
            std::cout<<"base cover lifted not active"<<std::endl;
            assert(lowerBound<0);
            return liftedCost;
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
        std::cout<<"snc repam left ";
        double before=l.LowerBound();
        std::cout<<"before lb "<<before<<std::endl;
       if(debug()){

          l.print();
       }
//       printIndices();

      // std::cout<<std::endl;
        assert(msg_dim <=nodeIndicesInCut.size());
       if(msg_dim==nodeIndicesInCut.size()){

           l.updateCostLifted(msg);
       }
       else{
           if(sncIsOut){
               l.updateCostBase(sncNodeIDindexInCut,nodeIndicesInCut.at(msg_dim),msg);
           }
           else{
               l.updateCostBase(nodeIndicesInCut.at(msg_dim),sncNodeIDindexInCut,msg);
           }
       }//std::cout<<"index to update "<<indicesInTriangle.at(msg_dim);
       double after=l.LowerBound();
       std::cout<<"after lb "<<after<<std::endl;

     }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {
        std::cout<<"cut repam right "<<std::endl;
        double before=r.LowerBound();
        std::cout<<"before lb "<<before<<std::endl;
        if(debug()){
            std::cout<<"cut repam right "<<std::endl;
            r.print();
        }
        assert(msg_dim <=nodeIndicesInSnc.size());
        if(msg_dim==nodeIndicesInSnc.size()){
            {
                assert(containsLiftedEdge);
                std::cout<<"node id of lifted edge "<<r.getLiftedIDs()[nodeIndexOfLiftedEdge]<<std::endl;
                r.updateEdgeCost(msg,nodeIndexOfLiftedEdge,true);
            }
        }
        else{
                r.updateEdgeCost(msg,nodeIndicesInSnc.at(msg_dim),false);
        }
        double after=r.LowerBound();
        std::cout<<"abter lb "<<after<<std::endl;

    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
         std::cout<<"snc send one to left "<<std::endl;
        if(debug()) std::cout<<"triangle send one to left ";
        //printIndices();
        //for (size_t i=0;i<1;i++) {
        size_t i=0;
        double delta=0;
        for (;i<nodeIndicesInSnc.size();i++) {

            delta = r.getOneBaseEdgeMinMarginal(nodeIndicesInSnc[i]);

            msg[i] -= omega * delta;
        }
        if(containsLiftedEdge){
            delta=r.getOneLiftedMinMarginal(nodeIndexOfLiftedEdge);
            msg[i]-=omega * delta;

        }

        //  std::cout<<" done "<<std::endl;
    }

    template<typename CUT_FACTOR, typename MSG>
    void send_message_to_right(const CUT_FACTOR& l, MSG& msg, const double omega)
    {
        std::cout<<"cut send one to right"<<std::endl;
        if(debug()) std::cout<<"triangle send one to right"<<std::endl;
        //for (size_t i=0;i<1;i++) {
        double delta;
        size_t i=0;
        for (;i<nodeIndicesInCut.size();i++) {
            if(sncIsOut){
                delta = l.getOneEdgeMinMarginal(sncNodeIDindexInCut,nodeIndicesInCut[i]);
            }
            else{
                delta = l.getOneEdgeMinMarginal(nodeIndicesInCut[i],sncNodeIDindexInCut);
            }
            msg[i] -= omega * delta;
        }
        if(containsLiftedEdge){
            std::cout<<"send to lifted "<<l.getLiftedOutputVertex()<<std::endl;
            delta=l.getLiftedMinMarginal();
            msg[i]-=omega*delta;
        }
    }



//    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
//    static void SendMessagesToLeft(const SINGLE_NODE_CUT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
//    {
//        if(debug()) std::cout<<"triangle send all to left"<<std::endl;

//        //TODO take only one half from the base min marginals!
//        std::vector<double> msg_vec_base = r.getAllBaseMinMarginals();
//        std::vector<double> localBaseCost=r.getBaseCosts();
//        for (size_t i = 0; i < localBaseCost.size(); ++i) {
//            msg_vec_base[i]=0.5*msg_vec_base[i];
//            localBaseCost[i]-=msg_vec_base[i];
//        }
//        const std::vector<double> msg_vec_lifted = r.getAllLiftedMinMarginals(&localBaseCost);

//        for(auto it=msg_begin; it!=msg_end; ++it)
//        {
//            auto& msg = (*it).GetMessageOp();
//            for (int i = 0; i <msg.verticesInSnc.size(); ++i) {
//                const size_t vertex = msg.verticesInSnc.at(i);
//                bool lifted=msg.isLifted.at(i);
//                if(lifted){
//                    (*it)[i] -=  omega * msg_vec_lifted.at(vertex);
//                }
//                else{
//                    (*it)[i] -= omega * msg_vec_base.at(vertex);
//                }

//            }
//        }
//    }





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
};





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
