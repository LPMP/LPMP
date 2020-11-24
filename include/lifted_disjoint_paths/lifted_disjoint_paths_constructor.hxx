#pragma once

#include "LP.h"
#include "solver.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
//#include "ldp_triangle_factor.hxx"
#include "MCF-SSP/mcf_ssp.hxx"
#include <unordered_map>
#include<andres/graph/components.hxx>
#include<map>
//#include"ldp_cut_factor.hxx"
#include <memory>
#include "ldp_min_marginals_extractor.hxx"
#include "ldp_path_separator.hxx"
#include "ldp_cut_message_creator.hxx"
#include "ldp_cut_factor_separator.hxx"

namespace LPMP {

template<class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
class lifted_disjoint_paths_constructor
{
public:
    using FMC = FACTOR_MESSAGE_CONNECTION;
    using ldp_cut_factor = typename CUT_FACTOR_CONT::FactorType;
       using ldp_path_factor_type = typename PATH_FACTOR::FactorType;

    template<typename SOLVER>
    lifted_disjoint_paths_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

    //void construct(const lifted_disjoint_paths_instance& i);
    void construct(const lifted_disjoint_paths::LdpInstance& instance);

    void ComputePrimal();

    void WritePrimal(std::stringstream& strStream)const;

    size_t Tighten(const std::size_t nr_constraints_to_add);
  //    size_t separateTriangles(const std::size_t nr_constraints_to_add);
        size_t separateCuts(const std::size_t nr_constraints_to_add);
    void pre_iterate() { reparametrize_snc_factors(); }

    std::vector<std::vector<size_t>> getBestPrimal()const { return bestPrimalSolution;}

    double getBestPrimalValue()const { return bestPrimalValue;}
private:
    std::size_t mcf_node_to_graph_node(std::size_t i) const;
    void read_in_mcf_costs(const bool change_marginals = false);
    void write_back_mcf_costs();
    void reparametrize_snc_factors();

    std::size_t nr_nodes() const { assert(single_node_cut_factors_.size() == (mcf_->no_nodes() - 2) / 2); return single_node_cut_factors_.size(); }
    std::size_t incoming_mcf_node(const std::size_t i) const { assert(i < nr_nodes()); return i*2; }
    std::size_t outgoing_mcf_node(const std::size_t i) const { assert(i < nr_nodes()); return i*2+1; }
    std::size_t mcf_source_node() const { return mcf_->no_nodes()-2; }
    std::size_t mcf_terminal_node() const { return mcf_->no_nodes()-1; }
    std::size_t base_graph_node(const std::size_t mcf_node) const;
    std::size_t base_graph_source_node() const { return nr_nodes(); }
    std::size_t base_graph_terminal_node() const { return nr_nodes() + 1; }

    void adjustLiftedLabels();
    void adjustCutLabels(size_t startPointer);
    void adjustPathLabels(size_t startPointer);
//    void adjustTriangleLabels(size_t firstIndex=0);

    bool checkFeasibilityInSnc();
    bool checkFeasibilityLiftedInSnc();
    bool checkFeasibilityBaseInSnc();

    void sncDebug();

    LP<FMC> *lp_;
    using mcf_solver_type = MCF::SSP<long, double>;
    std::unique_ptr<mcf_solver_type> mcf_; // minimum cost flow factor for base edges
    std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>> single_node_cut_factors_;
   // std::vector<CUT_FACTOR_CONT*> triangle_factors_;
    std::vector<CUT_FACTOR_CONT*> cut_factors_;
    std::vector<PATH_FACTOR*> path_factors_;
    std::vector<SINGLE_NODE_CUT_LIFTED_MESSAGE*> snc_lifted_messages_;
   // std::vector<SNC_CUT_MESSAGE*> snc_triangle_messages_;
     std::vector<SNC_CUT_MESSAGE*> snc_cut_messages_;
     std::vector<SNC_PATH_MESSAGE*> snc_path_messages_;
    //std::vector<std::vector<std::unordered_set<size_t>>> usedTriangles;
    const lifted_disjoint_paths::LdpInstance * pInstance;
    double bestPrimalValue;
    std::map<size_t,std::set<size_t>> addedCutFactorLiftedEdges;
    std::vector<std::vector<size_t>> bestPrimalSolution;
    std::vector<size_t> currentPrimalDescendants;
    std::vector<size_t> currentPrimalStartingVertices;
    std::vector<size_t> currentPrimalLabels;
    ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR> minMarginalsExtractor;

};

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE, class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE, SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::base_graph_node(const std::size_t mcf_node) const
{
    assert(mcf_node < mcf_->no_nodes());
    if(mcf_node == mcf_source_node())
        return base_graph_source_node();
    if (mcf_node == mcf_terminal_node())
        return base_graph_terminal_node();
    return mcf_node / 2;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::checkFeasibilityInSnc(){
    bool isFeasible=checkFeasibilityBaseInSnc();
    if(isFeasible){
        isFeasible=checkFeasibilityLiftedInSnc();
    }
    return isFeasible;
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::WritePrimal(std::stringstream& strStream) const{

    for (int i = 0; i < nr_nodes(); ++i) {
        size_t vertex=i;
        auto* snc=single_node_cut_factors_[vertex][0]->get_factor();
        if(snc->isNodeActive()&&(snc->getPrimalBaseVertexID()==base_graph_source_node())){

            while(vertex!=base_graph_terminal_node()){
                strStream<<vertex<<" ";
                auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
                vertex=sncFactorOut->getPrimalBaseVertexID();
            }
            strStream<<"\n";
        }

    }
}



template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::checkFeasibilityLiftedInSnc(){
  //Assumes primal feasible solution w.r.t. base edges and node labels
    bool isFeasible=true;


    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
        size_t vertex=i;
        auto* sncOut=single_node_cut_factors_[vertex][1]->get_factor();
        auto* sncIn=single_node_cut_factors_[vertex][0]->get_factor();
        if(sncOut->isNodeActive()!=sncIn->isNodeActive()){
            isFeasible=false;
            break;
        }
        if(!sncOut->isNodeActive()){
            if(!sncOut->getPrimalLiftedIndices().empty()){
                isFeasible=false;
            }
            if(!sncIn->getPrimalLiftedIndices().empty()){
                isFeasible=false;
            }
        }
        assert(isFeasible);

    }

    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
        size_t vertex=i;
        auto* sncOut=single_node_cut_factors_[vertex][1]->get_factor();
        if(!sncOut->isNodeActive()||(sncOut->getPrimalBaseVertexID()==base_graph_terminal_node())){
            if(!sncOut->getPrimalLiftedIndices().empty()){
                isFeasible=false;
                assert(isFeasible);
            }
        }
        else{
            const std::vector<size_t>& liftedIDs= sncOut->getLiftedIDs();
            size_t centralNodeLabel=currentPrimalLabels[i];
            assert(centralNodeLabel!=0);
            auto iter=sncOut->getPrimalLiftedIndices().begin();
            auto end=sncOut->getPrimalLiftedIndices().end();
            for(size_t j=0;j<liftedIDs.size();j++){
                size_t id=liftedIDs[j];
                if(iter==end||j!=*iter){
                    if(currentPrimalLabels[id]==centralNodeLabel){
                        isFeasible=false;
                        assert(isFeasible);
                    }
                }
                else{
                    if(currentPrimalLabels[id]!=centralNodeLabel){
                        isFeasible=false;
                        assert(isFeasible);
                    }
                    else{
                        iter++;
                    }
                }
            }

        }

        auto* sncIn=single_node_cut_factors_[vertex][0]->get_factor();
        if(!sncIn->isNodeActive()||(sncIn->getPrimalBaseVertexID()==base_graph_source_node())){
            if(!sncIn->getPrimalLiftedIndices().empty()){
                isFeasible=false;
                assert(isFeasible);
            }
        }
        else{
            const std::vector<size_t>& liftedIDs= sncIn->getLiftedIDs();
            size_t centralNodeLabel=currentPrimalLabels[i];
            assert(centralNodeLabel!=0);
            auto iter=sncIn->getPrimalLiftedIndices().begin();
            auto end=sncIn->getPrimalLiftedIndices().end();
            for(size_t j=0;j<liftedIDs.size();j++){
                size_t id=liftedIDs[j];
                if(iter==end||j!=*iter){
                    if(currentPrimalLabels[id]==centralNodeLabel){
                        isFeasible=false;
                        assert(isFeasible);
                    }
                }
                else{
                    if(currentPrimalLabels[id]!=centralNodeLabel){
                        isFeasible=false;
                        assert(isFeasible);
                    }
                    else{
                        iter++;
                    }
                }
            }
        }
    }

 return isFeasible;
}


//template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
//bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::checkFeasibilityLiftedInSnc(){
//  //Assumes primal feasible solution w.r.t. base edges and node labels
//    bool isFeasible=true;


//    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
//        size_t vertex=i;
//        auto* sncOut=single_node_cut_factors_[vertex][1]->get_factor();
//        auto* sncIn=single_node_cut_factors_[vertex][0]->get_factor();
//        if(sncOut->isNodeActive()!=sncIn->isNodeActive()){
//            isFeasible=false;
//            break;
//        }
//        if(!sncOut->isNodeActive()){
//            if(!sncOut->getPrimalLiftedIndices().empty()){
//                isFeasible=false;
//            }
//            if(!sncIn->getPrimalLiftedIndices().empty()){
//                isFeasible=false;
//            }
//        }

//    }

//    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
//        size_t vertex=i;
//        auto* snc=single_node_cut_factors_[vertex][1]->get_factor();
//        if(snc->isNodeActive()&&(snc->getPrimalBaseVertexID()==base_graph_terminal_node())){

//            std::vector<bool> isOnPath(nr_nodes(),0);
//            std::list<size_t> path;

//            while(vertex!=base_graph_source_node()&&isFeasible){
//                auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
//                const std::vector<size_t>& liftedIDs= sncFactorOut->getLiftedIDs();
//                const std::unordered_set<size_t>& activeInFactor=sncFactorOut->getPrimalLiftedIndices();
//                for (int j = 0; j < liftedIDs.size(); ++j) {
//                    if(isOnPath.at(liftedIDs.at(j))){
//                        if(activeInFactor.count(j)==0){
//                            isFeasible=false;
//                            break;
//                        }
//                    }
//                }
//                if(!isFeasible) break;
//                for(size_t activeV:activeInFactor){
//                    size_t activeVertexID=liftedIDs.at(activeV);
//                    if(!isOnPath.at(activeVertexID)){
//                        isFeasible=false;
//                        break;
//                    }
//                }
//                if(!isFeasible) break;
//                isOnPath[vertex]=1;
//                path.push_front(vertex);
//                auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
//                vertex=sncFactorIn->getPrimalBaseVertexID();
//            }
//            if(isFeasible){
//                for(size_t activeVertex:path){
//                    auto* sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();

//                    const std::vector<size_t>& liftedIDs= sncFactorIn->getLiftedIDs();
//                    const std::unordered_set<size_t>& activeInFactor=sncFactorIn->getPrimalLiftedIndices();
//                    for (int j = 0; j < liftedIDs.size(); ++j) {
//                        if(isOnPath.at(liftedIDs.at(j))){
//                            if(activeInFactor.count(j)==0){
//                                isFeasible=false;
//                                break;
//                            }
//                        }
//                    }
//                    if(!isFeasible) break;
//                    for(size_t activeV:activeInFactor){
//                        size_t activeVertexID=liftedIDs.at(activeV);
//                        if(!isOnPath.at(activeVertexID)){
//                            isFeasible=false;
//                            break;
//                        }
//                    }
//                    if(!isFeasible) break;
//                }
//            }
//        }
//    }

//    return isFeasible;

//}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::checkFeasibilityBaseInSnc(){

    //Check flow conservation
    bool isFeasible=true;

    std::unordered_map<size_t,size_t> predecessors; //Used later for lifted edges consistency
    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {

        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        return isFeasible;
        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();

        bool shouldBeActive=currentPrimalLabels[i]!=0;

        if(sncFactorIn->isNodeActive()!=shouldBeActive||sncFactorOut->isNodeActive()!=shouldBeActive){
            isFeasible=false;

        }
        else if(shouldBeActive){
            size_t descI=sncFactorOut->getPrimalBaseVertexID();
            if(currentPrimalDescendants[i]!=descI){
                isFeasible=false;
            }
            if(descI!=base_graph_terminal_node()){
                const auto* pairSncInput=single_node_cut_factors_[descI][0]->get_factor();
                if(pairSncInput->getPrimalBaseVertexID()!=i){
                    isFeasible=false;

                }
                if(currentPrimalLabels[descI]!=currentPrimalLabels[i]) isFeasible=false;
            }
            size_t predI=sncFactorIn->getPrimalBaseVertexID();
            if(predI!=base_graph_source_node()){
                if(currentPrimalDescendants[predI]!=i) isFeasible=false;
                const auto* pairSncOutput=single_node_cut_factors_[predI][1]->get_factor();
                if(pairSncOutput->getPrimalBaseVertexID()!=i){
                    isFeasible=false;

                }
                if(currentPrimalLabels[predI]!=currentPrimalLabels[i]) isFeasible=false;
            }

        }
    }
    for(const size_t& s:currentPrimalStartingVertices){
        const auto* sncFactorIn=single_node_cut_factors_[s][0]->get_factor();
        if(sncFactorIn->getPrimalBaseVertexID()!=base_graph_source_node()) isFeasible=false;
    }


    return isFeasible;
}




//template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
//bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::checkFeasibilityBaseInSnc(){

//    //Check flow conservation
//    bool isFeasible=true;

//    std::unordered_map<size_t,size_t> predecessors; //Used later for lifted edges consistency
//    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {

//        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
//        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();

//        if(sncFactorIn->isNodeActive()!=sncFactorOut->isNodeActive()){
//            isFeasible=false;

//        }
//        else if(sncFactorIn->isNodeActive()&&sncFactorOut->isNodeActive()){
//            size_t inputNode=sncFactorIn->getPrimalBaseVertexID();
//            if(inputNode!=base_graph_source_node()){
//                auto* sncFactorInputNodeOut=single_node_cut_factors_[inputNode][1]->get_factor();
//                if(sncFactorInputNodeOut->getPrimalBaseVertexID()!=i){
//                    isFeasible=false;

//                }
//            }

//            size_t outputNode=sncFactorOut->getPrimalBaseVertexID();
//            if(outputNode!=base_graph_terminal_node()){
//                auto* sncFactorOuputNodeIn=single_node_cut_factors_[outputNode][0]->get_factor();
//                if(sncFactorOuputNodeIn->getPrimalBaseVertexID()!=i){
//                    isFeasible=false;

//                }
//            }
//        }
//    }


//    return isFeasible;
//}



template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::adjustLiftedLabels(){
    //Assumes primal feasible solution w.r.t. base edges and node labels
    bool isFeasible=true;


    for (int i = 0; i < nr_nodes(); ++i) {
        size_t vertex=i;
        auto* sncOut=single_node_cut_factors_[vertex][1]->get_factor();
        auto* sncIn=single_node_cut_factors_[vertex][0]->get_factor();
        size_t centralNodeLabel=currentPrimalLabels[i];
        if(centralNodeLabel==0){
            std::vector<size_t> liftedIndices;
            sncOut->setPrimalLifted(liftedIndices);
            sncIn->setPrimalLifted(liftedIndices);
        }
        else {
            std::vector<size_t> liftedIndicesOut;
            if(sncOut->getPrimalBaseVertexID()!=base_graph_terminal_node()) {
                const std::vector<size_t>& liftedIDs= sncOut->getLiftedIDs();
                for(size_t j=0;j<liftedIDs.size();j++){
                    if(currentPrimalLabels[liftedIDs[j]]==centralNodeLabel){
                        liftedIndicesOut.push_back(j);
                    }
                }

            }
            sncOut->setPrimalLifted(liftedIndicesOut);

            std::vector<size_t> liftedIndicesIn;
            if(sncIn->getPrimalBaseVertexID()!=base_graph_source_node()) {
                const std::vector<size_t>& liftedIDs= sncIn->getLiftedIDs();
                for(size_t j=0;j<liftedIDs.size();j++){
                    if(currentPrimalLabels[liftedIDs[j]]==centralNodeLabel){
                        liftedIndicesIn.push_back(j);
                    }
                }

            }
            sncIn->setPrimalLifted(liftedIndicesIn);


        }

    }


}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::adjustPathLabels(size_t startPointer){
    //Assumes primal feasible solution w.r.t. base edges and node labels
    bool isFeasible=true;


    for (int i = startPointer; i < path_factors_.size(); ++i) {
        auto * pathFactorr=path_factors_[i]->get_factor();
        pathFactorr->setPrimal(currentPrimalDescendants,currentPrimalLabels);
    }

}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::adjustCutLabels(size_t startPointer){
    //Assumes primal feasible solution w.r.t. base edges and node labels
    bool isFeasible=true;


    for (int i = startPointer; i < cut_factors_.size(); ++i) {
        auto * cutFactor=cut_factors_[i]->get_factor();
        cutFactor->setPrimal(currentPrimalDescendants,currentPrimalLabels);
    }

}

//template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
//void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::adjustLiftedLabels(){
////Assumes primal feasible solution w.r.t. base edges and node labels
//    for (int i = 0; i < nr_nodes(); ++i) {
//        size_t vertex=i;
//        auto* snc=single_node_cut_factors_[vertex][1]->get_factor();
//        if(snc->isNodeActive()&&(snc->getPrimalBaseVertexID()==base_graph_terminal_node())){

//            std::vector<bool> isOnPath(nr_nodes(),0);
//            isOnPath[vertex]=1;
//            std::list<size_t> path;
//           // path.push_front(vertex);

//            while(vertex!=base_graph_source_node()){
//                auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
//                std::unordered_set<size_t> activeEndpointIndices;
//                const std::vector<size_t>& liftedIDs= sncFactorOut->getLiftedIDs();
//                for (int j = 0; j < liftedIDs.size(); ++j) {
//                    if(isOnPath[liftedIDs.at(j)]) activeEndpointIndices.insert(j);
//                }

//                sncFactorOut->setPrimalLifted(activeEndpointIndices);
//                isOnPath[vertex]=1;
//                path.push_front(vertex);
//                auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
//                vertex=sncFactorIn->getPrimalBaseVertexID();
//            }
//            for(size_t activeVertex:path){
//                auto* sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();
//                std::unordered_set<size_t> activeEndpointIndices;
//                const std::vector<size_t>& liftedIDs= sncFactorIn->getLiftedIDs();
//                for (int j = 0; j < liftedIDs.size(); ++j) {
//                    if(isOnPath[liftedIDs.at(j)]) activeEndpointIndices.insert(j);
//                }
//                sncFactorIn->setPrimalLifted(activeEndpointIndices);
//            }
//        }
//        else if(!(snc->isNodeActive())){
//            auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
//            std::unordered_set<size_t> verticesOfActiveEdges;
//            sncFactorOut->setPrimalLifted(verticesOfActiveEdges);
//            auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
//            sncFactorIn->setPrimalLifted(verticesOfActiveEdges);

//        }
//    }
//}


/*
template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::adjustTriangleLabels(size_t firstIndex){
//Assumes primal feasible solution w.r.t. base edges and node labels and adjusted lifted labels
    for (size_t i = firstIndex; i < triangle_factors_.size(); ++i) {
        auto * trFactor=triangle_factors_[i]->get_factor();
        std::bitset<3> primalSolution("000");
        size_t v1=trFactor->getV1();
        size_t v2=trFactor->getV2();
        size_t v3=trFactor->getV3();
        auto * sncV1=single_node_cut_factors_[v1][1]->get_factor();
        auto * sncV2=single_node_cut_factors_[v2][1]->get_factor();

        if(sncV1->isNodeActive()){
            if(trFactor->isV1V2Base()){
                size_t primalV1=sncV1->getPrimalBaseVertexID();
                if(primalV1==v2){
                    primalSolution[0]=1;
                }
            }
            else{
                if(sncV1->isActiveInPrimalLifted(sncV1->getLiftedIDToOrder(v2))){
                    primalSolution[0]=1;
                }
            }
            if(sncV1->isActiveInPrimalLifted(sncV1->getLiftedIDToOrder(v3))){
                primalSolution[2]=1;
            }
        }
        if(sncV2->isNodeActive()){
            if(trFactor->isV2V3Base()){
                size_t primalV2=sncV2->getPrimalBaseVertexID();
                if(primalV2==v3){
                    primalSolution[1]=1;
                }
            }
            else{
                if(sncV2->isActiveInPrimalLifted(sncV2->getLiftedIDToOrder(v3))){
                    primalSolution[1]=1;
                }
            }
        }

        trFactor->setPrimal(primalSolution);


    }
}
*/



template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::construct(const lifted_disjoint_paths::LdpInstance &instance)
{
    pInstance=&instance;
    bestPrimalValue=std::numeric_limits<double>::max();




    // first construct minimum cost flow factor for base edges
    const std::size_t nr_base_graph_nodes = instance.getGraph().numberOfVertices() - 2;
    const std::size_t base_graph_source = instance.getSourceNode();
    const std::size_t base_graph_terminal = instance.getTerminalNode();
    //std::cout << "source = " << base_graph_source << ", terminal = " << base_graph_terminal << "\n";

    const std::size_t nr_mcf_nodes = 2*nr_base_graph_nodes + 2; // source/terminal vertex + 2*ordinary vertices to ensure unit capacity vertices
    const std::size_t nr_mcf_edges = 3*nr_base_graph_nodes + instance.getGraph().numberOfEdges() + 1; // appearance/disappearance/uniqueness edge + connection edges + source/terminal edge
    mcf_ = std::make_unique<mcf_solver_type>(nr_mcf_nodes, nr_mcf_edges);

    const std::size_t mcf_source_node = nr_mcf_nodes - 2;
    const std::size_t mcf_terminal_node = nr_mcf_nodes - 1;
    auto incoming_edge = [](const std::size_t i) { return 2*i; };
    auto outgoing_edge = [](const std::size_t i) { return 2*i+1; };

    mcf_->add_edge(mcf_source_node, mcf_terminal_node, 0, nr_base_graph_nodes, 0.0);
    mcf_->add_node_excess(mcf_source_node, nr_base_graph_nodes);
    mcf_->add_node_excess(mcf_terminal_node, -std::ptrdiff_t(nr_base_graph_nodes));

    // ensure unit capacity on vertices, appearance and disappearance edges
    for(std::size_t i=0; i<nr_base_graph_nodes; ++i)
    {
        mcf_->add_edge(incoming_edge(i), outgoing_edge(i), 0, 1, 0.0);
        mcf_->add_edge(mcf_source_node, incoming_edge(i), 0, 1, 0.0);
        mcf_->add_edge(outgoing_edge(i), mcf_terminal_node, 0, 1, 0.0);
    }

    for(std::size_t i=0; i<nr_base_graph_nodes; ++i)
    {
        for (std::size_t ne = 0; ne < instance.getGraph().numberOfEdgesFromVertex(i); ++ne)
        {
            const std::size_t e = instance.getGraph().edgeFromVertex(i, ne);
            const std::size_t j = instance.getGraph().vertexFromVertex(i, ne);
            //std::cout << "base graph edge " << i << "," << j << "\n";
            if (j != base_graph_source && j != base_graph_terminal)
            {
                assert(i < j);
                const double c = instance.getEdgeScore(e);
                mcf_->add_edge(outgoing_edge(i), incoming_edge(j), 0, 1, c);
            }
        }
    }

    mcf_->order();


    // next add all single node cut factors
    single_node_cut_factors_.reserve(instance.getGraph().numberOfVertices());
    for(std::size_t i=0; i<nr_base_graph_nodes; ++i)
    {
        auto* incoming_snc = lp_->template add_factor<SINGLE_NODE_CUT_FACTOR>(instance, i, false);
        incoming_snc->get_factor()->initBaseCosts(0.5);
        incoming_snc->get_factor()->initLiftedCosts(0.5);
        incoming_snc->get_factor()->initNodeCost(0.5);
        auto* outgoing_snc = lp_->template add_factor<SINGLE_NODE_CUT_FACTOR>(instance, i, true);
        outgoing_snc->get_factor()->initBaseCosts(0.5);
        outgoing_snc->get_factor()->initLiftedCosts(0.5);
        outgoing_snc->get_factor()->initNodeCost(0.5);
        single_node_cut_factors_.push_back({incoming_snc, outgoing_snc});
    }

    assert(base_graph_source == this->base_graph_source_node());
    assert(base_graph_terminal == this->base_graph_terminal_node());
    assert(mcf_source_node == this->mcf_source_node());
    assert(mcf_terminal_node == this->mcf_terminal_node());

    if(debug()) std::cout<<"factors created "<<std::endl;

    // add messages between lifted edges of snc factors
    // for all lifted messages
    // get outgoing factor container -> left_snc, left_node
    // get incoming factor container -> right_snc, right_node

    //lp_->add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>(left_snc, right_snc, left_node, right_node);

    for (int vertex1 = 0; vertex1 < single_node_cut_factors_.size(); ++vertex1) {
        auto* left_snc=single_node_cut_factors_[vertex1][1];
        auto* left_snc_factor=left_snc->get_factor();
        //	std::cout<<"vertex 1 "<<vertex1<<std::endl;
        for (int j = 0; j < left_snc_factor->getLiftedIDs().size(); ++j) {
            size_t vertex2=left_snc_factor->getLiftedIDs().at(j);
            //std::cout<<"vertex 2 "<<vertex2<<std::endl;
            assert(instance.getGraphLifted().findEdge(vertex1,vertex2).first);
            if(vertex2==this->base_graph_source_node()||vertex2==this->base_graph_source_node()) continue;
            auto* right_snc=single_node_cut_factors_[vertex2][0];
            size_t i=right_snc->get_factor()->getLiftedIDToOrder(vertex1);
            assert((right_snc->get_factor()->getLiftedIDs().size())>i);
            auto * message=lp_->template add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE>(left_snc, right_snc, i, j);
            snc_lifted_messages_.push_back(message);

        }

    }


    currentPrimalDescendants=std::vector<size_t>(nr_nodes(),base_graph_terminal_node());
    currentPrimalLabels=std::vector<size_t>(nr_nodes(),0);

    minMarginalsExtractor =ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>(&single_node_cut_factors_,pInstance) ;

   /* if(debug()) std::cout<<"messages added"<<std::endl;
    usedTriangles=std::vector<std::vector<std::unordered_set<size_t>>>(nr_nodes());
    for (int i = 0; i < nr_nodes(); ++i) {
        size_t nrIncomingEdges=instance.getGraph().numberOfEdgesToVertex(i);
        nrIncomingEdges+=instance.getGraphLifted().numberOfEdgesToVertex(i);
        usedTriangles[i]=std::vector<std::unordered_set<size_t>>(nrIncomingEdges);
    }
    //Tighten(200);
    */
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::sncDebug()
{


    std::map<size_t,std::map<size_t,double>> baseEdgesFromFactors;
    std::map<size_t,std::map<size_t,double>> liftedEdgesFromFactors;

    for (int i = 0; i < nr_nodes(); ++i) {
        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        const std::vector<double>& bCosts=sncFactorIn->getBaseCosts();
        const std::vector<size_t>& bEdges=sncFactorIn->getBaseIDs();
        double nodeCost=sncFactorIn->getNodeCost();
        for (int j = 0; j < bCosts.size(); ++j) {
            baseEdgesFromFactors[bEdges[j]][i]+=bCosts[j];
           // baseEdgesFromFactors[bEdges[j]][i]+=nodeCost;
        }
        const std::vector<double>& lCosts=sncFactorIn->getLiftedCosts();
        const std::vector<size_t>& lEdges=sncFactorIn->getLiftedIDs();
        for (int j = 0; j < lCosts.size(); ++j) {
            liftedEdgesFromFactors[lEdges[j]][i]+=lCosts[j];
        }
    }

    for (int i = 0; i < nr_nodes(); ++i) {
        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        const std::vector<double>& bCosts=sncFactorOut->getBaseCosts();
        const std::vector<size_t>& bEdges=sncFactorOut->getBaseIDs();
        double nodeCost=sncFactorOut->getNodeCost();
        for (int j = 0; j < bCosts.size(); ++j) {
            baseEdgesFromFactors[i][bEdges[j]]+=bCosts[j];
            //baseEdgesFromFactors[i][bEdges[j]]+=nodeCost;
        }
        const std::vector<double>& lCosts=sncFactorOut->getLiftedCosts();
        const std::vector<size_t>& lEdges=sncFactorOut->getLiftedIDs();
        for (int j = 0; j < lCosts.size(); ++j) {
            liftedEdgesFromFactors[i][lEdges[j]]+=lCosts[j];
        }
    }

    for (int i = 0; i < cut_factors_.size(); ++i) {
        auto * cFactor=cut_factors_[i]->get_factor();
        const std::vector<size_t>& inputs=cFactor->getInputVertices();
        const std::vector<size_t>& outputs=cFactor->getOutputVertices();
        auto & cutGraph=cFactor->getCutGraph();
        for (int j = 0; j < inputs.size(); ++j) {
            size_t inputVertex=inputs[j];
            auto * iter=cutGraph.forwardNeighborsBegin(j);
            auto * end=cutGraph.forwardNeighborsEnd(j);
            for (;iter!=end;iter++) {
                size_t outputVertex=outputs[iter->first];
                baseEdgesFromFactors[inputVertex][outputVertex]+=iter->second;
            }

        }
        size_t v=cFactor->getLiftedInputVertex();
        size_t w=cFactor->getLiftedOutputVertex();
        double liftedCost=cFactor->getLiftedCost();
        liftedEdgesFromFactors[v][w]+=liftedCost;

    }



    for (int i = 0; i < path_factors_.size(); ++i) {
        auto * pFactor=path_factors_[i]->get_factor();
        const std::vector<size_t>& vertices= pFactor->getListOfVertices();
        const std::vector<double>& costs= pFactor->getCosts();
        const std::vector<char>& liftedInfo=pFactor-> getLiftedInfo();
        for (int j = 0; j < vertices.size()-1; ++j) {
            size_t v1=vertices[j];
            size_t v2=vertices[j+1];
            double cost=costs[j];
            if(liftedInfo.at(j)>0){
                liftedEdgesFromFactors[v1][v2]+=cost;
            }
            else{
                baseEdgesFromFactors[v1][v2]+=cost;
            }
        }
        liftedEdgesFromFactors[vertices[0]][vertices.back()]+=costs.back();
        //TODO: Add longest edge!
    }

    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();

    for (int i = 0; i < baseGraph.getNumberOfVertices(); ++i) {
        const auto* iter=baseGraph.forwardNeighborsBegin(i);
        const auto* end=baseGraph.forwardNeighborsEnd(i);
        for (;iter!=end;iter++) {
            size_t vertex=iter->first;
            double cost=iter->second;
            double factorsCost=baseEdgesFromFactors[i][vertex];
            if(abs(factorsCost-cost)>eps){
                std::cout<<"base cost mismatch "<<i<<", "<<vertex<<", orig: "<<cost<<", factor cost: "<<factorsCost<<std::endl;
            }
        }
    }


    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        const auto* iter=liftedGraph.forwardNeighborsBegin(i);
        const auto* end=liftedGraph.forwardNeighborsEnd(i);
        for (;iter!=end;iter++) {
            size_t vertex=iter->first;
            double cost=iter->second;
            double factorsCost=liftedEdgesFromFactors[i][vertex];
            if(abs(factorsCost-cost)>eps){
                std::cout<<"lifted cost mismatch "<<i<<", "<<vertex<<", orig: "<<cost<<", factor cost: "<<factorsCost<<std::endl;
            }
        }
    }




}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::ComputePrimal()
{
    read_in_mcf_costs();
    mcf_->solve();

    std::vector<size_t> startingNodes;
    std::vector<size_t> descendants(nr_nodes(),base_graph_terminal_node());


    for (std::size_t graph_node = 0; graph_node < nr_nodes(); ++graph_node) {
        single_node_cut_factors_[graph_node][0]->get_factor()->setNoBaseEdgeActive();
        single_node_cut_factors_[graph_node][1]->get_factor()->setNoBaseEdgeActive();
        const std::size_t mcf_incoming_node = incoming_mcf_node(graph_node);

        std::size_t vertex_index = 0;
        for (std::size_t j = 0; j < mcf_->no_outgoing_arcs(mcf_incoming_node); ++j) {
            const std::size_t edge_id=j+mcf_->first_outgoing_arc(mcf_incoming_node);
            const std::size_t node = mcf_->head(edge_id);
            if(base_graph_node(node) != graph_node) {
                auto* pSNC=single_node_cut_factors_[graph_node][0]->get_factor();
                if(mcf_->flow(edge_id) == -1) {
                    pSNC->setBaseEdgeActive(vertex_index);
                    size_t nodeID=pSNC->getBaseID(vertex_index);
                    if(nodeID==base_graph_source_node()){
                        startingNodes.push_back(graph_node);
                    }

                }
                ++vertex_index;
            }
        }

        const std::size_t mcf_outgoing_node = outgoing_mcf_node(graph_node);
        vertex_index = 0;

        for (std::size_t j = 0; j < mcf_->no_outgoing_arcs(mcf_outgoing_node); ++j) {
            const std::size_t edge_id=j+mcf_->first_outgoing_arc(mcf_outgoing_node);
            const std::size_t node = mcf_->head(edge_id);
            if(base_graph_node(node) != graph_node) {
                auto* pSNC=single_node_cut_factors_[graph_node][1]->get_factor();
                if(mcf_->flow(edge_id) == 1) {
                    pSNC->setBaseEdgeActive(vertex_index);
                    size_t nodeID=pSNC->getBaseID(vertex_index);
                    assert(descendants[graph_node]==base_graph_terminal_node());
                    descendants[graph_node]=nodeID;
                }
                ++vertex_index;
            }
        }
    }

    currentPrimalDescendants=descendants;
    currentPrimalStartingVertices=startingNodes;
    std::fill(currentPrimalLabels.begin(),currentPrimalLabels.end(),0);

    std::vector<std::vector<size_t>> paths;
    double toAdd=0;
    std::vector<size_t>  labels(nr_nodes());
    for (size_t i = 0; i < startingNodes.size(); ++i) {
        std::vector<size_t> path;
        size_t currentNode=startingNodes[i];

        while(currentNode!=base_graph_terminal_node()){
            path.push_back(currentNode);
            currentPrimalLabels[currentNode]=i+1;
            currentNode=descendants[currentNode];
        }
        paths.push_back(path);
    }

    //TODO check feasibility base: compare set active with the theree above vectors
    //Set lifted based on primal labels . Carefull with zero labels!
    //Set primal for cuts: based on the three vectors

    bool isFeasible=this->checkFeasibilityBaseInSnc();
    if(diagnostics()) std::cout<<"checked feasibility: "<<isFeasible<<std::endl;
    assert(isFeasible);
    adjustLiftedLabels();
    isFeasible=this->checkFeasibilityLiftedInSnc();
  //  adjustTriangleLabels();
    assert(isFeasible);
    adjustCutLabels(0);
    adjustPathLabels(0);
    double primalValue=0;
    double primalBaseValue=0;
    double primalLiftedValue=0;

    for (int i = 0; i < nr_nodes(); ++i) {

        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        primalValue+=sncFactorIn->EvaluatePrimal();
        primalValue+=sncFactorOut->EvaluatePrimal();

        primalBaseValue+=sncFactorIn->getPrimalBaseCost();
        primalBaseValue+=sncFactorOut->getPrimalBaseCost();

        primalLiftedValue+=sncFactorIn->getPrimalLiftedCost();
        primalLiftedValue+=sncFactorOut->getPrimalLiftedCost();


    }




  /*  for (int i = 0; i < triangle_factors_.size(); ++i) {
        auto * trFactor=triangle_factors_[i]->get_factor();
        primalValue+=trFactor->EvaluatePrimal();

    }
    */


    std::map<size_t,std::map<size_t,double>> costsFromOtherFactors;
    std::map<size_t,std::map<size_t,std::set<size_t>>> indicesOfCutFactors;
    std::map<size_t,std::map<size_t,std::set<size_t>>> indicesOfPathFactors;


    for (int i = 0; i < cut_factors_.size(); ++i) {
        auto * cFactor=cut_factors_[i]->get_factor();
        primalValue+=cFactor->EvaluatePrimal();

        primalBaseValue+=cFactor->getPrimalBaseCost();
        primalLiftedValue+=cFactor->getPrimalLiftedCost();
        const auto & inputs=cFactor->getInputVertices();
        const auto & outputs=cFactor->getOutputVertices();
        const auto & cutGraph=cFactor->getCutGraph();
        for (int j = 0; j < inputs.size(); ++j) {
            size_t v1=inputs[j];
            const auto * iter=cutGraph.forwardNeighborsBegin(j);
            const auto * end=cutGraph.forwardNeighborsEnd(j);
            for (;iter!=end;iter++) {
                size_t v2=outputs.at(iter->first);
                double cost=iter->second;
                costsFromOtherFactors[v1][v2]+=cost;
                indicesOfCutFactors[v1][v2].insert(i);
            }
        }

    }

    for (int i = 0; i < path_factors_.size(); ++i) {
        auto * pFactor=path_factors_[i]->get_factor();
        primalValue+=pFactor->EvaluatePrimal();

        primalBaseValue+=pFactor->getPrimalBaseCost();
        primalLiftedValue+=pFactor->getPrimalLiftedCost();

        const auto& pathVertices=pFactor->getListOfVertices();
        for (int j = 0; j < pathVertices.size()-1; ++j) {
            if(pFactor->getLiftedInfo().at(j)==0){
                size_t v1=pathVertices[j];
                size_t v2=pathVertices[j+1];
                costsFromOtherFactors[v1][v2]+=pFactor->getCosts().at(j);
                indicesOfPathFactors[v1][v2].insert(i);

            }
        }

    }


//TODO swap order and use these paths for adjusting lifted and triangle labels



        if(primalValue < bestPrimalValue){

            bestPrimalValue=primalValue;


            bestPrimalSolution=paths;
        }

        if(diagnostics()){
            //sncDebug();
           // std::cout<<"computed primal value "<<primalValue<<std::endl;
            double controlPrimalValue=0;
            double controlLiftedValue=0;
            double controlBaseValue=0;
            //std::vector<size_t> labels(nr_nodes());
            for (size_t i = 0; i < startingNodes.size(); ++i) {
                // std::cout<<"path "<<i<<std::endl;
                double pathCost=0;

                size_t currentNode=startingNodes[i];



                double factorPathCost=0;
                double sncFactorsCost=0;
                controlPrimalValue+=pInstance->getEdgeScore(pInstance->getSourceNode(),currentNode);
                pathCost+=pInstance->getEdgeScore(pInstance->getSourceNode(),currentNode);
                bool hasPathFactor=false;
                bool hasCutFactor=false;
                std::set<size_t> usedPathFactors;
                while(currentNode!=base_graph_terminal_node()){
                    assert(currentNode<labels.size());
                   // labels[currentNode]=i+1;
                    double valueToAdd=pInstance->getVertexScore(currentNode);
                    valueToAdd+=pInstance->getEdgeScore(currentNode,descendants[currentNode]);
                    controlPrimalValue+=valueToAdd;
                    pathCost+=valueToAdd;


                    const auto* sncFactorIn=single_node_cut_factors_[currentNode][0]->get_factor();
                    const auto* sncFactorOut=single_node_cut_factors_[currentNode][1]->get_factor();
                    sncFactorsCost+=sncFactorIn->getPrimalBaseCost();
                    sncFactorsCost+=sncFactorOut->getPrimalBaseCost();
                    factorPathCost+=costsFromOtherFactors[currentNode][descendants[currentNode]];

                    if(!indicesOfPathFactors[currentNode][descendants[currentNode]].empty()){
                        std::set<size_t>& pf=indicesOfPathFactors[currentNode][descendants[currentNode]];
                        hasPathFactor=true;
                        usedPathFactors.insert(pf.begin(),pf.end());

                    }
                    if(!indicesOfCutFactors[currentNode][descendants[currentNode]].empty()) hasCutFactor=true;
                    if(!hasPathFactor&&!hasCutFactor) assert(costsFromOtherFactors[currentNode][descendants[currentNode]]==0);
                    currentNode=descendants[currentNode];
                }
                if(std::abs(pathCost-factorPathCost-sncFactorsCost)>eps){
                    std::cout<<"cost mismatch in path "<<i<<", factor complete "<<(factorPathCost+sncFactorsCost)<<", control "<<pathCost<<"has path "<<hasPathFactor<<", cost snc "<<sncFactorsCost<<", other factors: "<<factorPathCost<<", has cut "<<hasCutFactor<<std::endl;
                    currentNode=startingNodes[i];
                    while(currentNode!=base_graph_terminal_node()){
                        std::cout<<currentNode<<", ";
                        currentNode=descendants[currentNode];
                    }
                    std::cout<<std::endl;
                    for(auto f:usedPathFactors){
                        std::cout<<"pf index "<<f<<std::endl;
                        auto * pathF=path_factors_[f]->get_factor();
                        pathF->print();
                        std::cout<<" base cost "<<pathF->getPrimalBaseCost()<<std::endl;
                    }

                    assert(false);



                }
//                else{
//                    std::cout<<"path "<<i<<" ok"<<std::endl;
//                }

            }
            controlBaseValue=controlPrimalValue;




//            //std::cout<<"lifted edges"<<std::endl;
//            for(size_t e=0;e<pInstance->getGraph().numberOfEdges();e++){
//                size_t v0=pInstance->getGraph().vertexOfEdge(e,0);
//                size_t v1=pInstance->getGraph().vertexOfEdge(e,1);
//                if(v0==base_graph_source_node()){

//                }
//                else if(v1==base_graph_terminal_node()){

//                }
//                else if
//                if(currentPrimalDescendants[v0]==v1){
//                    controlPrimalValue+=pInstance->getEdgeScore(e);
//                }
//            }

            for(size_t e=0;e<pInstance->getGraphLifted().numberOfEdges();e++){
                size_t v0=pInstance->getGraphLifted().vertexOfEdge(e,0);
                size_t v1=pInstance->getGraphLifted().vertexOfEdge(e,1);
                if(currentPrimalLabels[v0]==currentPrimalLabels[v1]&&currentPrimalLabels[v0]!=0){

                    controlPrimalValue+=pInstance->getLiftedEdgeScore(e);
                    controlLiftedValue+=pInstance->getLiftedEdgeScore(e);

                }
            }
            if(std::abs(primalValue-controlPrimalValue)>=eps){
                std::cout<<"computed primal value "<<primalValue<<std::endl;
                std::cout<<"control primal value "<<controlPrimalValue<<std::endl;
                std::cout<<"computed base primal value "<<primalBaseValue<<std::endl;
                std::cout<<"control primal base value "<<controlBaseValue<<std::endl;
                std::cout<<"computed lifted primal value "<<primalLiftedValue<<std::endl;
                std::cout<<"control lifted base value "<<controlLiftedValue<<std::endl;

            }
            assert(std::abs(primalValue-controlPrimalValue)<eps);


            if(debug()){
                for(auto path: bestPrimalSolution){
                    for(auto v:path){
                        std::cout<<v<<" ";
                    }
                    std::cout<<std::endl;
                }
            }
        }

//    best_primal_solution = ...;

    if(diagnostics()) std::cout<<"primal value: "<<primalValue<<std::endl;
}

//std::vector<size_t> get_best_solution() const
//{
//    return best_primal_solution;
//}






template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::Tighten(const std::size_t nr_constraints_to_add)
{
   size_t numberOfCutsToSeparate=nr_constraints_to_add/2;
   size_t numberOfPathsToSeparate=nr_constraints_to_add-numberOfCutsToSeparate;

   // size_t cutsSeparated=  separateCuts(numberOfCutsToSeparate);
   size_t counterAdded=0;
   LdpCutSeparator<ldp_cut_factor,SINGLE_NODE_CUT_FACTOR> cutSeparator(pInstance,minMarginalsExtractor);
   cutSeparator.separateCutInequalities(numberOfCutsToSeparate);
   std::priority_queue<std::pair<double,ldp_cut_factor*>>& queueWithCuts=cutSeparator.getPriorityQueue();
   while(!queueWithCuts.empty()&&counterAdded<numberOfCutsToSeparate){
       ldp_cut_factor* pCutFromQueue=queueWithCuts.top().second;
       auto * newCutFactor=lp_->template add_factor<CUT_FACTOR_CONT>(*pCutFromQueue);
       cut_factors_.push_back(newCutFactor);
       delete pCutFromQueue;
       pCutFromQueue=nullptr;
       queueWithCuts.pop();
       auto * pCutFactor=newCutFactor->get_factor();
       const std::vector<size_t>& inputs=pCutFactor->getInputVertices();
       const std::vector<size_t>& outputs=pCutFactor->getOutputVertices();
       bool liftedAdded=false;
       for(size_t i=0;i<inputs.size();i++){
           size_t inputVertex=inputs[i];
           auto * snc=single_node_cut_factors_[inputVertex][1];
           LdpCutMessageInputs<ldp_cut_factor,SINGLE_NODE_CUT_FACTOR> messageInputs;
           messageInputs.init(pCutFactor,snc,i);


           if(messageInputs.containsLifted) liftedAdded=true;
           auto * message1=lp_->template add_message<SNC_CUT_MESSAGE>(newCutFactor,snc,messageInputs._nodeIndicesInCut,messageInputs._nodeIndicesInSnc,i,true,messageInputs.containsLifted,messageInputs._nodeIndexOfLiftedEdge);
           snc_cut_messages_.push_back(message1);
       }
       if(!liftedAdded){
           std::vector<size_t> _nodeIndicesInCut;
           std::vector<size_t> _nodeIndicesInSnc;
           size_t v=pCutFactor->getLiftedInputVertex();
           size_t w=pCutFactor->getLiftedOutputVertex();
           auto * snc=single_node_cut_factors_[v][1];
           size_t _nodeIndexOfLiftedEdge=snc->get_factor()->getLiftedIDToOrder(w);
           auto * message1=lp_->template add_message<SNC_CUT_MESSAGE>(newCutFactor,snc,_nodeIndicesInCut,_nodeIndicesInSnc,inputs.size(),true,true,_nodeIndexOfLiftedEdge);
           snc_cut_messages_.push_back(message1);
       }
       counterAdded++;
   }
   cutSeparator.clearPriorityQueue();



    //size_t counterAdded=cutsSeparated;

    ldp_path_separator<ldp_path_factor_type,SINGLE_NODE_CUT_FACTOR> pathSeparator(pInstance, minMarginalsExtractor);
    pathSeparator.separatePathInequalities(nr_constraints_to_add);
    std::priority_queue<std::pair<double,ldp_path_factor_type*>>& queueWithPaths=pathSeparator.getPriorityQueue();

    double possibleImprovement=0;
   // std::cout<<"queue filled "<<queueWithPaths.size()<<std::endl;


    size_t pathFactorsOriginalSize=path_factors_.size();
    while(!queueWithPaths.empty()&&counterAdded<nr_constraints_to_add){
       // std::pair<double,ldp_path_factor_type*>& p=queueWithPaths.top();

        ldp_path_factor_type* pPathFactor=queueWithPaths.top().second;
        double improvement=queueWithPaths.top().first;
        possibleImprovement+=improvement;
        //auto* newPathFactor = lp_->template add_factor<PATH_FACTOR>(*pPathFactor);
        auto* newPathFactor = lp_->template add_factor<PATH_FACTOR>(pPathFactor->getListOfVertices(),pPathFactor->getCosts(),pPathFactor->getLiftedInfo());
        path_factors_.push_back(newPathFactor);
       // std::cout<<"factor added, number of vertices "<<pPathFactor->getNumberOfEdges()<<std::endl;
        delete pPathFactor;
        pPathFactor=nullptr;
        queueWithPaths.pop();
        auto *myPathFactor=newPathFactor->get_factor();
        //myPathFactor->print();
        const std::vector<size_t>& pathVertices=myPathFactor->getListOfVertices();
        for (size_t i=0;i<myPathFactor->getNumberOfEdges();i++) {
            size_t pathVertex=pathVertices[i];
            if(i>0){
                auto * pSNC=single_node_cut_factors_[pathVertex][0];
                LdpPathMessageInputs<ldp_path_factor_type,SINGLE_NODE_CUT_FACTOR> messageInputs;
                messageInputs.init(myPathFactor,pSNC,i);
                //LdpPathMessageInputs messageInputs=pathSeparator.getMessageInputsToPathFactor(myPathFactor,pSNC,i);
                bool debugInfo=path_factors_.size()==17||path_factors_.size()==19;
                auto* newMessage = lp_->template add_message<SNC_PATH_MESSAGE>(newPathFactor,pSNC,messageInputs.edgeIndicesInPath,messageInputs.indicesInSnc,messageInputs.isLiftedForMessage,debugInfo);
                snc_path_messages_.push_back(newMessage);
            }
            if(i<myPathFactor->getNumberOfEdges()-1){
                auto * pSNCOut=single_node_cut_factors_[pathVertex][1];
                LdpPathMessageInputs<ldp_path_factor_type,SINGLE_NODE_CUT_FACTOR> messageInputsOut;
                messageInputsOut.init(myPathFactor,pSNCOut,i);
                bool debugInfo=path_factors_.size()==17||path_factors_.size()==19;
              //  LdpPathMessageInputs messageInputsOut=pathSeparator.getMessageInputsToPathFactor(myPathFactor,pSNCOut,i);
                auto* newMessageOut = lp_->template add_message<SNC_PATH_MESSAGE>(newPathFactor,pSNCOut,messageInputsOut.edgeIndicesInPath,messageInputsOut.indicesInSnc,messageInputsOut.isLiftedForMessage,debugInfo);
                snc_path_messages_.push_back(newMessageOut);
            }
            //std::cout<<"vertex "<<i<<" of path solved "<<std::endl;
        }
        counterAdded ++;
        //std::cout<<"added "<<counterAdded<<" path ineq"<<std::endl;
    }

    pathSeparator.clearPriorityQueue();
    adjustPathLabels(pathFactorsOriginalSize);

    if(diagnostics()) std::cout<<"added "<<counterAdded<<" path factors with improvement "<<possibleImprovement<<std::endl;
     return counterAdded;

}





template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::mcf_node_to_graph_node(std::size_t i) const
{
    assert(i < mcf_->no_nodes());
    if(i == mcf_source_node())
        return base_graph_source_node();
    if(i == mcf_terminal_node())
        return base_graph_terminal_node();
    else
        return i/2;
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE, class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE, SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::read_in_mcf_costs(const bool change_marginals)
{
    mcf_->reset_costs();

    const lifted_disjoint_paths::LdpInstance &instance=*pInstance;
    //const andres::graph::Digraph<>& baseGraph=instance.getGraph();
   // std::vector<std::unordered_map<size_t,double>> costsFromCuts(nr_nodes());


 /*   if(debug()) std::cout<<"triangle factors size "<<triangle_factors_.size()<<std::endl;
    for(size_t i=0;i<triangle_factors_.size();i++){
        auto * trFactor=triangle_factors_[i]->get_factor();
        size_t v1=trFactor->getV1();
        size_t v2=trFactor->getV2();
        size_t v3=trFactor->getV3();
        if(trFactor->isV1V2Base()){
            double m=trFactor->getOneMinMarginal(0);
            trFactor->updateCost(0,-m);
            costsFromTriangles[v1][v2]+=m;
        }
        if(trFactor->isV2V3Base()){
            double m=trFactor->getOneMinMarginal(1);
            trFactor->updateCost(1,-m);
            costsFromTriangles[v2][v3]+=m;
        }
    }
    */

//   for(size_t i=0;i<cut_factors_.size();i++){ //Need fix: update costs with already created marginals
//        auto * cFactor=cut_factors_[i]->get_factor();
//        for(size_t j=0;j<cFactor->getNumberOfInputs();j++){
//            for(size_t k=0;k<cFactor->getNumberOfOutputs();k++){
//                double value=cFactor->getOneEdgeMinMarginal(j,k);
//                costsFromCuts[cFactor->getInputVertices()[j]][cFactor->getOutputVertices()[k]]=value;

//            }
//        }
//    }


    for(std::size_t i=0; i<nr_nodes(); ++i)
    {
        {

            auto *incoming_snc = single_node_cut_factors_[i][0]->get_factor();
            const auto incoming_min_marg = incoming_snc->getAllBaseMinMarginalsForMCF();
            const auto incoming_edge_ids = incoming_snc->getBaseIDs();
            assert(incoming_min_marg.size() == incoming_edge_ids.size());
            assert(incoming_min_marg.size() + 1 == mcf_->no_outgoing_arcs(incoming_mcf_node(i))); // one extra mcf node for capacity one arc.
            //std::vector<std::tuple<std::size_t, double>> incoming_min_margs_sorted;
            assert(std::is_sorted(incoming_edge_ids.begin(), incoming_edge_ids.end()));
            std::size_t e = mcf_->first_outgoing_arc(incoming_mcf_node(i));
            assert(incoming_min_marg.size() <= mcf_->no_outgoing_arcs(incoming_mcf_node(i)));
            for (std::size_t l = 0; l < incoming_min_marg.size(); ++l)
            {
                const std::size_t j = incoming_edge_ids[l];
                assert(j != base_graph_terminal_node());
               // std::cout << "incoming min marginal edge " << i << "," << j << "\n";
                while (mcf_node_to_graph_node(mcf_->head(e)) != j)
                {
                    //std::cout << "mcf edge " << e << ": " << mcf_->tail(e) << "," << mcf_->head(e) << "; base edge: " << mcf_node_to_graph_node(mcf_->tail(e)) << "," << mcf_node_to_graph_node(mcf_->head(e)) << "\n";
                    assert(mcf_->tail(e) == incoming_mcf_node(i));
                    ++e;
                }
                const std::size_t start_node = mcf_->tail(e);
                assert(mcf_->tail(e) == incoming_mcf_node(i));
                const double m = incoming_min_marg[l];

                assert(mcf_->lower_bound(e) == 1 && mcf_->upper_bound(e) == 0);
                if (j != base_graph_source_node())
                    mcf_->update_cost(e, -m);
                else
                    mcf_->update_cost(e, -m);
            }

            if(change_marginals)
            {
                for(std::size_t l = 0; l < incoming_min_marg.size(); ++l)
                {
                    incoming_snc->updateEdgeCost(-incoming_min_marg[l], l, false);
                }
            }
        }

        {

            auto *outgoing_snc = single_node_cut_factors_[i][1]->get_factor();
            const auto outgoing_min_marg = outgoing_snc->getAllBaseMinMarginalsForMCF();
            const auto outgoing_edge_ids = outgoing_snc->getBaseIDs();
            assert(outgoing_min_marg.size() + 1 == mcf_->no_outgoing_arcs(outgoing_mcf_node(i))); // one extra mcf node for capacity one arc.
            assert(std::is_sorted(outgoing_edge_ids.begin(), outgoing_edge_ids.end()));
            std::size_t e = mcf_->first_outgoing_arc(outgoing_mcf_node(i));
            for (std::size_t l = 0; l < outgoing_min_marg.size(); ++l)
            {
                const std::size_t j = outgoing_edge_ids[l];
                assert(j != base_graph_source_node());
                while (mcf_node_to_graph_node(mcf_->head(e)) != j)
                {
                    assert(mcf_->tail(e) == outgoing_mcf_node(i));
                    ++e;
                }
                const std::size_t start_node = mcf_->tail(e);
                assert(mcf_->tail(e) == outgoing_mcf_node(i));
                const double m = outgoing_min_marg[l];

                assert(mcf_->lower_bound(e) == 0 && mcf_->upper_bound(e) == 1);
                mcf_->update_cost(e, m);
//                auto iter=costsFromCuts[i].find(j);
//                if(iter!=costsFromCuts[i].end()){
//                    mcf_->update_cost(e,iter->second);
//                }
            }
            if(change_marginals)
            {
                for (std::size_t l = 0; l < outgoing_min_marg.size(); ++l)
                {
                    outgoing_snc->updateEdgeCost(-outgoing_min_marg[l], l, false);
                }
            }
        }
    }
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::write_back_mcf_costs()
{
    assert(std::abs(mcf_->potential(mcf_source_node()) - mcf_->potential(mcf_terminal_node())) <= 1e-8);
    for(std::size_t i=0; i<nr_nodes(); ++i)
    {
        {
            const std::size_t current_node = incoming_mcf_node(i);
            auto *incoming_snc = single_node_cut_factors_[i][0]->get_factor();
            //assert(mcf_->no_outgoing_arcs(incoming_mcf_node(i)) == incoming_snc->getBaseCosts().size() + 1); // TODO: enable again
            std::size_t vertex_index = 0;
            for(std::size_t e = mcf_->first_outgoing_arc(incoming_mcf_node(i)); e<mcf_->first_outgoing_arc(incoming_mcf_node(i)) + mcf_->no_outgoing_arcs(incoming_mcf_node(i)); ++e)
            {
                const int flow = mcf_->flow(e);
                const std::size_t head = mcf_->head(e);
                const double orig_cost = mcf_->cost(e);
                const double reduced_cost = mcf_->reduced_cost(e);
                const std::size_t incoming_node = mcf_->head(e);
                if(base_graph_node(incoming_node) == base_graph_source_node()) // from source
                {
                    assert(flow == 0 || flow == -1);
                    if(flow == 0)
                        assert(reduced_cost <= 1e-8);
                    if (flow == -1)
                        assert(reduced_cost >= -1e-8);
                    //incoming_snc->updateEdgeCost(- reduced_cost, base_graph_node(incoming_node), false);
                    //assert(vertex_index == incoming_snc->getBaseIDs().size()-1);
                    assert(incoming_snc->getBaseIDs()[vertex_index] == base_graph_node(incoming_node));
                   // incoming_snc->updateEdgeCost(-orig_cost, vertex_index , false);
                     incoming_snc->updateEdgeCost(-reduced_cost, vertex_index , false);
                    ++vertex_index;
                }
                else if (base_graph_node(incoming_node) != i) // regular edge
                {
                    assert(flow == 0 || flow == -1);
                    if(flow == 0)
                        assert(reduced_cost <= 1e-8);
                    if (flow == -1)
                        assert(reduced_cost >= -1e-8);
                    assert(incoming_snc->getBaseIDs()[vertex_index] == base_graph_node(incoming_node));
                   // incoming_snc->updateEdgeCost(-0.5 * orig_cost, vertex_index, false);
                     incoming_snc->updateEdgeCost(-0.5 * reduced_cost, vertex_index, false);
                    ++vertex_index;
                }
                else // bottleneck edge
                {
                    assert(flow == 0 || flow == 1);
                    if(flow == 0)
                        assert(reduced_cost >= -1e-8);
                    if (flow == 1)
                        assert(reduced_cost <= 1e-8);
                    assert(base_graph_node(incoming_node) == i);
                    incoming_snc->updateNodeCost(0.5 * reduced_cost);
                }
            }
        }

        {
            const std::size_t current_node = outgoing_mcf_node(i);
            auto *outgoing_snc = single_node_cut_factors_[i][1]->get_factor();
            // assert(mcf_->no_outgoing_arcs(outgoing_mcf_node(i)) == outgoing_snc->getBaseCosts().size() + 1); // TODO: enable again
            std::size_t  vertex_index = 0;
            for(std::size_t e = mcf_->first_outgoing_arc(outgoing_mcf_node(i)); e<mcf_->first_outgoing_arc(outgoing_mcf_node(i)) + mcf_->no_outgoing_arcs(outgoing_mcf_node(i)); ++e)
            {
                const int flow = mcf_->flow(e);
                const std::size_t head = mcf_->head(e);
                const double orig_cost = mcf_->cost(e);
                const double reduced_cost = mcf_->reduced_cost(e);
                const std::size_t outgoing_node = mcf_->head(e);
                if(base_graph_node(outgoing_node) == base_graph_terminal_node()) // terminal edge
                {
                    assert(flow == 0 || flow == 1);
                    if(flow == 0)
                        assert(reduced_cost >= -1e-8);
                    if (flow == 1)
                        assert(reduced_cost <= 1e-8);
                    //outgoing_snc->updateEdgeCost(reduced_cost, base_graph_node(outgoing_node), false);
                    assert(outgoing_snc->getBaseIDs()[vertex_index] == base_graph_node(outgoing_node));
                    //outgoing_snc->updateEdgeCost(orig_cost, vertex_index , false);
                    outgoing_snc->updateEdgeCost(reduced_cost, vertex_index , false);
                    ++vertex_index;
                }
                else if (base_graph_node(outgoing_node) != i) // regular edge
                {
                    assert(flow == 0 || flow == 1);
                    if(flow == 0)
                        assert(reduced_cost >= -1e-8);
                    if (flow == 1)
                        assert(reduced_cost <= 1e-8);
                    //outgoing_snc->updateEdgeCost(0.5 * reduced_cost, base_graph_node(outgoing_node), false);
                    assert(outgoing_snc->getBaseIDs()[vertex_index] == base_graph_node(outgoing_node));
                    //outgoing_snc->updateEdgeCost(0.5 * orig_cost, vertex_index, false);
                    outgoing_snc->updateEdgeCost(0.5 * reduced_cost, vertex_index, false);
                    ++vertex_index;
                }
                else // bottleneck edge
                {
                    assert(flow == 0 || flow == -1);
                    if(flow == 0)
                        assert(reduced_cost <= 1e-8);
                    if (flow == -1)
                        assert(reduced_cost >= -1e-8);
                    assert(base_graph_node(outgoing_node) == i);
                    outgoing_snc->updateNodeCost(-0.5 * reduced_cost);
                }
            }
        }
    }
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class CUT_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_CUT_MESSAGE,class PATH_FACTOR,class SNC_PATH_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, CUT_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_CUT_MESSAGE,PATH_FACTOR,SNC_PATH_MESSAGE>::reparametrize_snc_factors()
{
    const double primal_cost_before = this->lp_->EvaluatePrimal();
    read_in_mcf_costs(true);
    mcf_->solve();
   if(diagnostics())  std::cout << "mcf cost = " << mcf_->objective() << "\n";
    write_back_mcf_costs();
    const double primal_cost_after = this->lp_->EvaluatePrimal();
    if(diagnostics()) std::cout << "primal cost before = " << primal_cost_before << ", primal cost after = " << primal_cost_after << "\n";
    assert(std::abs(primal_cost_before - primal_cost_after) <= 1e-6);
 //   sncDebug();
}


}
