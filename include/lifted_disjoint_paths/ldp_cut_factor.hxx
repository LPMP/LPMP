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


    void updateCostBaseForward(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value);


    void updateCostBaseBackward(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value);

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



    std::pair<LdpTwoLayerGraph,double> getAllMinMarginals()const;

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
        sncNodeIDindexInCut(_sncNodeIDindexInCut), //Warning if only lifted edge is incident, this is equal to number of inputs
        sncIsOut(_sncIsOut),   //if true, central node is in inputs, other nodes in outputs
        containsLiftedEdge(_containsLiftedEdge),
        nodeIndexOfLiftedEdge(_nodeIndexOfLiftedEdge){


    }

    template<typename CUT_FACTOR>
    void RepamLeft(CUT_FACTOR& l, const double msg, const std::size_t msg_dim) const
    {

       // l.LowerBound();
        assert(msg_dim <=nodeIndicesInCut.size());
        if(msg_dim==nodeIndicesInCut.size()){
            if(debug()){
                if(sncIsOut){
                    lastV1=l.getLiftedInputVertex();
                    lastV2=l.getLiftedOutputVertex();
                }
                else{
                    lastV2=l.getLiftedInputVertex();
                    lastV1=l.getLiftedOutputVertex();
                }
                lastValue=msg;
            }

            l.updateCostLifted(msg);
        }
        else{
            if(sncIsOut){
                l.updateCostBaseForward(sncNodeIDindexInCut,nodeIndicesInCut.at(msg_dim),msg);
                if(debug()){
                    lastV1=l.getInputVertices().at(sncNodeIDindexInCut);
                    size_t otherVertex=l.getCutGraph().getForwardEdgeVertex(sncNodeIDindexInCut,nodeIndicesInCut[msg_dim]);
                    lastV2=l.getOutputVertices().at(otherVertex);
                    lastValue=msg;
                }
            }
            else{
                l.updateCostBaseBackward(sncNodeIDindexInCut,nodeIndicesInCut.at(msg_dim),msg);
                if(debug()){
                    lastV1=l.getOutputVertices().at(sncNodeIDindexInCut);
                    size_t otherVertex=l.getCutGraph().getBackwardEdgeVertex(sncNodeIDindexInCut,nodeIndicesInCut[msg_dim]);
                    lastV2=l.getInputVertices().at(otherVertex);
                    lastValue=msg;
                }
            }

        }
      //  l.LowerBound();

    }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {

        assert(msg_dim <=nodeIndicesInSnc.size());
      // if(debug()) r.LowerBound();
        size_t secondVertex;
        if(msg_dim==nodeIndicesInSnc.size()){

            assert(containsLiftedEdge);
            r.updateEdgeCost(msg,nodeIndexOfLiftedEdge,true);
            if(debug()) secondVertex=r.getLiftedIDs().at(nodeIndexOfLiftedEdge);

        }
        else{
            r.updateEdgeCost(msg,nodeIndicesInSnc.at(msg_dim),false);
            if(debug()) secondVertex=r.getBaseIDs().at(nodeIndicesInSnc.at(msg_dim));
        }
        if(debug()){
            assert(lastV1==r.nodeID);
            if(lastV2!=secondVertex){
                std::cout<<"cut mismatch, first vertex "<<lastV1<<", vertex in cut "<<lastV2<<",  vertex in snc "<<secondVertex<<", index in cut"<<nodeIndicesInCut[msg_dim]<<", is lifted "<<(msg_dim==nodeIndicesInSnc.size())<<std::endl;

            }
            assert(lastV2==secondVertex);
            assert(std::abs(lastValue+msg)<eps);
        }
      //  if(debug()) r.LowerBound();

    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
        //std::cout<<"send cut to left"<<std::endl;

        std::vector<double> baseCosts=r.getBaseCosts();
        std::vector<double> liftedCosts=r.getLiftedCosts(); //TODO change to const reference
        size_t i=0;
        double delta=0;
        for (;i<nodeIndicesInSnc.size();i++) {

            delta = r.getOneBaseEdgeMinMarginal(nodeIndicesInSnc[i],&baseCosts,&liftedCosts);
            baseCosts[nodeIndicesInSnc[i]]-=delta;
            //if(debug()){
                double controlDelta=r.getOneBaseEdgeMinMarginal(nodeIndicesInSnc[i],&baseCosts,&liftedCosts);
                if(abs(controlDelta)>eps){
                    std::cout<<"base control value "<<controlDelta<<", original delta: "<<delta<<std::endl;
                    throw std::runtime_error("wrong min marginal check");
                }
            //}

            msg[i] -= omega * delta;
        }
        if(containsLiftedEdge){
            delta=r.getOneLiftedMinMarginal(nodeIndexOfLiftedEdge,&baseCosts,&liftedCosts);
            //if(debug()){
                liftedCosts[nodeIndexOfLiftedEdge]-=delta;
                double controlDelta=r.getOneLiftedMinMarginal(nodeIndexOfLiftedEdge,&baseCosts,&liftedCosts);
                if(abs(controlDelta)>eps){
                    std::cout<<"lifted control value "<<controlDelta<<", original delta: "<<delta<<std::endl;
                    throw std::runtime_error("wrong lifted min marginal check");
                }
            //}
            msg[i]-=omega * delta;

        }


    }

    template<typename CUT_FACTOR, typename MSG>
    void send_message_to_right(const CUT_FACTOR& l, MSG& msg, const double omega)
    {
        //std::cout<<"send cut to right"<<std::endl;

        double delta;
        size_t i=0;
        double liftedCost=l.getLiftedCost();
        LdpTwoLayerGraph cutGraph=l.getCutGraph();
        double controlDelta=0;
        for (;i<nodeIndicesInCut.size();i++) {

            if(sncIsOut){
                size_t otherVertex=l.getCutGraph().getForwardEdgeVertex(sncNodeIDindexInCut,nodeIndicesInCut[i]);
                //std::cout<<"base edge mm, neighbor node "<<l.getOutputVertices()[otherVertex]<<std::endl;
                delta = l.getOneEdgeMinMarginal(sncNodeIDindexInCut,otherVertex,&cutGraph,&liftedCost);
                //std::cout<<"orig edge cost "<<cutGraph.getForwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i])<<std::endl;;
                cutGraph.updateForwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i],-delta);
                //std::cout<<"changed edge cost "<<cutGraph.getForwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i])<<std::endl;
                if(debug()){
                    controlDelta=l.getOneEdgeMinMarginal(sncNodeIDindexInCut,otherVertex,&cutGraph,&liftedCost);
                }
            }
            else{  //just for output nodes now

                size_t otherVertex=l.getCutGraph().getBackwardEdgeVertex(sncNodeIDindexInCut,nodeIndicesInCut[i]);
                //std::cout<<"base edge mm, neighbor node "<<l.getOutputVertices()[otherVertex]<<std::endl;
                delta = l.getOneEdgeMinMarginal(otherVertex,sncNodeIDindexInCut,&cutGraph,&liftedCost);
                //std::cout<<"orig edge cost "<<cutGraph.getForwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i])<<std::endl;;
                cutGraph.updateBackwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i],-delta);
                //std::cout<<"changed edge cost "<<cutGraph.getForwardEdgeCost(sncNodeIDindexInCut,nodeIndicesInCut[i])<<std::endl;
                if(debug()){
                    controlDelta=l.getOneEdgeMinMarginal(otherVertex,sncNodeIDindexInCut,&cutGraph,&liftedCost);
                }
            }
            if(debug()){
                if(abs(controlDelta)>eps){
                    std::cout<<"control value from cut factor"<<controlDelta<<", original delta: "<<delta<<std::endl;
                    throw std::runtime_error("wrong cut min marginal - check base edge");
                }
            }
           // std::cout<<"delta "<<delta<<", sending "<<(omega*delta)<<std::endl;

            msg[i] -= omega * delta;
        }

        if(containsLiftedEdge){

            // std::cout<<"lifted edge mm, neighbor node "<<l.getLiftedOutputVertex()<<std::endl;
            delta=l.getLiftedMinMarginal(&cutGraph,&liftedCost);
            if(debug()){
            //std::cout<<"orig lifted edge cost "<<liftedCost<<std::endl;
                liftedCost-=delta;
              //  std::cout<<"new lifted edge cost "<<liftedCost<<std::endl;
                controlDelta=l.getLiftedMinMarginal(&cutGraph,&liftedCost);
                if(abs(controlDelta)>eps){
                    std::cout<<"control value from cut factor"<<controlDelta<<", original delta: "<<delta<<std::endl;
                    throw std::runtime_error("wrong cut min marginal - check lifted edge");
                }
            }
            //    std::cout<<"delta "<<delta<<", sending "<<(omega*delta)<<std::endl;
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
