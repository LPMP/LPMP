#pragma once

#include "LP.h"
#include "solver.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "ldp_triangle_factor.hxx"
#include "MCF-SSP/mcf_ssp.hxx"
#include <unordered_map>
#include<andres/graph/components.hxx>
#include<map>

#include <memory>

namespace LPMP {

template<class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
class lifted_disjoint_paths_constructor
{
public:
    using FMC = FACTOR_MESSAGE_CONNECTION;
    template<typename SOLVER>
    lifted_disjoint_paths_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

    //void construct(const lifted_disjoint_paths_instance& i);
    void construct(const lifted_disjoint_paths::LdpInstance& instance);

    void ComputePrimal();

    size_t Tighten(const std::size_t nr_constraints_to_add);
    void pre_iterate() { reparametrize_snc_factors(); }

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
    void adjustTriangleLabels(size_t firstIndex=0);

    bool checkFeasibilityInSnc();
    bool checkFeasibilityLiftedInSnc();
    bool checkFeasibilityBaseInSnc();

    void sncDebug(size_t vertex,bool isOut);

    LP<FMC> *lp_;
    using mcf_solver_type = MCF::SSP<long, double>;
    std::unique_ptr<mcf_solver_type> mcf_; // minimum cost flow factor for base edges
    std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>> single_node_cut_factors_;
    std::vector<TRIANGLE_FACTOR_CONT*> triangle_factors_;
    std::vector<SINGLE_NODE_CUT_LIFTED_MESSAGE*> snc_lifted_messages_;
    std::vector<SNC_TRIANGLE_MESSAGE*> snc_triangle_messages_;
    std::vector<std::vector<std::unordered_set<size_t>>> usedTriangles;

};

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE, class SNC_TRIANGLE_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE, SNC_TRIANGLE_MESSAGE>::base_graph_node(const std::size_t mcf_node) const
{
    assert(mcf_node < mcf_->no_nodes());
    if(mcf_node == mcf_source_node())
        return base_graph_source_node();
    if (mcf_node == mcf_terminal_node())
        return base_graph_terminal_node();
    return mcf_node / 2;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::checkFeasibilityInSnc(){
    bool isFeasible=checkFeasibilityBaseInSnc();
    if(isFeasible){
        isFeasible=checkFeasibilityLiftedInSnc();
    }
    return isFeasible;
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::checkFeasibilityLiftedInSnc(){
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

    }

    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
        size_t vertex=i;
        auto* snc=single_node_cut_factors_[vertex][1]->get_factor();
        if(snc->isNodeActive()&&(snc->getPrimalBaseVertexID()==base_graph_terminal_node())){

            std::vector<bool> isOnPath(nr_nodes(),0);
            std::list<size_t> path;

            while(vertex!=base_graph_source_node()&&isFeasible){
                auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
                const std::vector<size_t>& liftedIDs= sncFactorOut->getLiftedIDs();
                const std::unordered_set<size_t>& activeInFactor=sncFactorOut->getPrimalLiftedIndices();
                for (int i = 0; i < liftedIDs.size(); ++i) {
                    if(isOnPath.at(liftedIDs.at(i))){
                        if(activeInFactor.count(i)==0){
                            isFeasible=false;
                            break;
                        }
                    }
                }
                if(!isFeasible) break;
                for(size_t activeV:activeInFactor){
                    size_t activeVertexID=liftedIDs.at(activeV);
                    if(!isOnPath.at(activeVertexID)){
                        isFeasible=false;
                        break;
                    }
                }
                if(!isFeasible) break;
                isOnPath[vertex]=1;
                path.push_front(vertex);
                auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
                vertex=sncFactorIn->getPrimalBaseVertexID();
            }
            if(isFeasible){
                for(size_t activeVertex:path){
                    auto* sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();

                    const std::vector<size_t>& liftedIDs= sncFactorIn->getLiftedIDs();
                    const std::unordered_set<size_t>& activeInFactor=sncFactorIn->getPrimalLiftedIndices();
                    for (int i = 0; i < liftedIDs.size(); ++i) {
                        if(isOnPath.at(liftedIDs.at(i))){
                            if(activeInFactor.count(i)==0){
                                isFeasible=false;
                                break;
                            }
                        }
                    }
                    if(!isFeasible) break;
                    for(size_t activeV:activeInFactor){
                        size_t activeVertexID=liftedIDs.at(activeV);
                        if(!isOnPath.at(activeVertexID)){
                            isFeasible=false;
                            break;
                        }
                    }
                    if(!isFeasible) break;
                }
            }
        }
    }

    return isFeasible;

}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::checkFeasibilityBaseInSnc(){

    //Check flow conservation
    bool isFeasible=true;

    std::unordered_map<size_t,size_t> predecessors; //Used later for lifted edges consistency
    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {

        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();

        if(sncFactorIn->isNodeActive()!=sncFactorOut->isNodeActive()){
            isFeasible=false;

        }
        else if(sncFactorIn->isNodeActive()&&sncFactorOut->isNodeActive()){
            size_t inputNode=sncFactorIn->getPrimalBaseVertexID();
            if(inputNode!=base_graph_source_node()){
                auto* sncFactorInputNodeOut=single_node_cut_factors_[inputNode][1]->get_factor();
                if(sncFactorInputNodeOut->getPrimalBaseVertexID()!=i){
                    isFeasible=false;

                }
            }

            size_t outputNode=sncFactorOut->getPrimalBaseVertexID();
            if(outputNode!=base_graph_terminal_node()){
                auto* sncFactorOuputNodeIn=single_node_cut_factors_[outputNode][0]->get_factor();
                if(sncFactorOuputNodeIn->getPrimalBaseVertexID()!=i){
                    isFeasible=false;

                }
            }
        }
    }


    return isFeasible;
}



template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::adjustLiftedLabels(){
//Assumes primal feasible solution w.r.t. base edges and node labels
    for (int i = 0; i < nr_nodes(); ++i) {
        size_t vertex=i;
        auto* snc=single_node_cut_factors_[vertex][1]->get_factor();
        if(snc->isNodeActive()&&(snc->getPrimalBaseVertexID()==base_graph_terminal_node())){

            std::vector<bool> isOnPath(nr_nodes(),0);
            isOnPath[vertex]=1;
            std::list<size_t> path;
           // path.push_front(vertex);

            while(vertex!=base_graph_source_node()){
                auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
                std::unordered_set<size_t> activeEndpointIndices;
                const std::vector<size_t>& liftedIDs= sncFactorOut->getLiftedIDs();
                for (int i = 0; i < liftedIDs.size(); ++i) {
                    if(isOnPath[liftedIDs.at(i)]) activeEndpointIndices.insert(i);
                }

                sncFactorOut->setPrimalLifted(activeEndpointIndices);
                isOnPath[vertex]=1;
                path.push_front(vertex);
                auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
                vertex=sncFactorIn->getPrimalBaseVertexID();
            }
            for(size_t activeVertex:path){
                auto* sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();
                std::unordered_set<size_t> activeEndpointIndices;
                const std::vector<size_t>& liftedIDs= sncFactorIn->getLiftedIDs();
                for (int i = 0; i < liftedIDs.size(); ++i) {
                    if(isOnPath[liftedIDs.at(i)]) activeEndpointIndices.insert(i);
                }
                sncFactorIn->setPrimalLifted(activeEndpointIndices);
            }
        }
        else if(!(snc->isNodeActive())){
            auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
            std::unordered_set<size_t> verticesOfActiveEdges;
            sncFactorOut->setPrimalLifted(verticesOfActiveEdges);
            auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
            sncFactorIn->setPrimalLifted(verticesOfActiveEdges);

        }
    }
}



template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::adjustTriangleLabels(size_t firstIndex){
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




template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::construct(const lifted_disjoint_paths::LdpInstance &instance)
{

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
        auto* outgoing_snc = lp_->template add_factor<SINGLE_NODE_CUT_FACTOR>(instance, i, true);
        outgoing_snc->get_factor()->initBaseCosts(0.5);
        outgoing_snc->get_factor()->initLiftedCosts(0.5);
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

    //lp_->add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>(left_snc, right_snc, left_node, right_node);

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
//            if(debug()){
//                std::cout<<"message left "<<vertex1<<", right"<<vertex2<<", ij "<<i<<","<<j<<std::endl;
//            }
        }
        //
        //        	for (int j = 0; j < instance.getGraph().numberOfEdgesFromVertex(vertex1); ++j) {
        //        		size_t vertex2=instance.getGraph().vertexFromVertex(vertex1,j);
        //        		auto right_snc=single_node_cut_factors_[vertex2][0];
        //        		size_t i=right_snc->get_factor()->getLiftedOrderToID(vertex1);
        //        		if(debug()) std::cout<<"message "<<vertex1<<", "<<vertex2<<": "<<i<<", "<<j<<std::endl;
        //        		auto * message=lp_->template add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE>(left_snc, right_snc, i, j);
        //        		snc_lifted_messages_.push_back(message);
        //
        //
        //			}
    }



    //        for (int i = 0; i < single_node_cut_factors_.size(); ++i) {
    //        	auto left_snc=single_node_cut_factors_[i][0];
    //        	for(const auto pair:left_snc->get_factor()->getLiftedCosts()){
    //        		size_t j=pair.first;
    //        		auto right_snc=single_node_cut_factors_[j][1];
    //
    //        		auto * message=lp_->template add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE>(left_snc, right_snc, i, j);
    //        		snc_lifted_messages_.push_back(message);
    //
    //
    //        	}
    //        }
    //        for (int i = 0; i < single_node_cut_factors_.size(); ++i) {
    //        	if(i==this->base_graph_source_node()||i==this->base_graph_terminal_node()) continue;
    //        	auto left_snc=single_node_cut_factors_[i][0];
    //        	auto right_snc=single_node_cut_factors_[i][1];
    //        	auto * message=lp_->template add_message<SNC_TRIANGLE_MESSAGE>(left_snc, right_snc, i);
    //        	SNC_TRIANGLE_MESSAGEs_.push_back(message);
    //        }

    //sncDebug(20,1);
    if(debug()) std::cout<<"messages added"<<std::endl;
    usedTriangles=std::vector<std::vector<std::unordered_set<size_t>>>(nr_nodes());
    for (int i = 0; i < nr_nodes(); ++i) {
        size_t nrIncomingEdges=instance.getGraph().numberOfEdgesToVertex(i);
        nrIncomingEdges+=instance.getGraphLifted().numberOfEdgesToVertex(i);
        usedTriangles[i]=std::vector<std::unordered_set<size_t>>(nrIncomingEdges);
    }
    //Tighten(200);
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::sncDebug(size_t vertex,bool isOut)
{
    std::cout<<"print node "<<vertex<<std::endl;
    auto* sncFactor=single_node_cut_factors_[vertex][isOut]->get_factor();
    const std::vector<double>& baseCosts=sncFactor->getBaseCosts();
    const std::vector<double>& liftedCosts=sncFactor->getLiftedCosts();

    std::cout<<"base costs"<<std::endl;

    for (int i = 0; i < baseCosts.size(); ++i) {
        std::cout<<i<<": ";
        std::cout<<(sncFactor->baseIDs[i])<<": "<<baseCosts[i];
        std::cout<<std::endl;
    }

    std::cout<<std::endl;
    std::cout<<"lifted costs"<<std::endl;
    for (int i = 0; i < liftedCosts.size(); ++i) {
        std::cout<<i<<": "<<(sncFactor->getLiftedIDs()[i])<<": ";
        std::cout<<liftedCosts[i];
        std::cout<<std::endl;
    }



    std::cout<<"getting all base min marginals"<<std::endl;

    sncFactor->getAllBaseMinMarginals();
    std::cout<<"isolated min marginals "<<std::endl;
    for (int i = 0; i < baseCosts.size(); ++i) {
        double value=sncFactor->getOneBaseEdgeMinMarginal(i);
        std::cout<<i<<": "<<"("<<(sncFactor->baseIDs[i])<<")"<<value<<std::endl;
    }



    std::cout<<"computing lifted min marginal "<<std::endl;
    double value=sncFactor->oneLiftedMinMarginal(1);
    //std::cout<<"lifted min marginal "<<value<<std::endl;
    //	sncFactor->getAllLiftedMinMarginals();

}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::ComputePrimal()
{
    read_in_mcf_costs();
    mcf_->solve();

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
                      std::cout<<pSNC->getBaseIDs()[vertex_index]<<", "<<graph_node<<std::endl;
                   // primalValue+=pSNC->EvaluatePrimal();
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
                      std::cout<<graph_node<<", "<<pSNC->getBaseIDs()[vertex_index]<<std::endl;
                   // primalValue+=pSNC->EvaluatePrimal();
                }
                ++vertex_index;
            }
        }
    }
    bool isFeasible=this->checkFeasibilityBaseInSnc();
    std::cout<<"checked feasibility: "<<isFeasible<<std::endl;
    assert(isFeasible);
    adjustLiftedLabels();
    isFeasible=this->checkFeasibilityLiftedInSnc();
    adjustTriangleLabels();
    assert(isFeasible);
    double primalValue=0;
    for (int i = 0; i < nr_nodes(); ++i) {

        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        primalValue+=sncFactorIn->EvaluatePrimal();
        primalValue+=sncFactorOut->EvaluatePrimal();


    }
    for (int i = 0; i < triangle_factors_.size(); ++i) {
        auto * trFactor=triangle_factors_[i]->get_factor();
        primalValue+=trFactor->EvaluatePrimal();

    }

    std::cout<<"primal value: "<<primalValue<<std::endl;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::Tighten(const std::size_t nr_constraints_to_add)
{
    //TODO: Remember triangles that have already been added!

    std::cout<<"TIGHTEN "<<nr_constraints_to_add<<std::endl;

    const lifted_disjoint_paths::LdpInstance &instance=single_node_cut_factors_[0][0]->get_factor()->getLdpInstance();
    const andres::graph::Digraph<>& baseGraph=instance.getGraph();
    const andres::graph::Digraph<>& liftedGraph=instance.getGraphLifted();

    std::vector<double> baseEdgeLabelsIn(baseGraph.numberOfEdges(),0);
    std::vector<double> liftedEdgeLabelsIn(liftedGraph.numberOfEdges(),0);
    std::vector<double> baseEdgeLabelsOut(baseGraph.numberOfEdges(),0);
    std::vector<double> liftedEdgeLabelsOut(liftedGraph.numberOfEdges(),0);

    std::vector<bool> baseEdgeUsed(baseGraph.numberOfEdges(),0);
    std::vector<bool> liftedEdgeUsed(liftedGraph.numberOfEdges(),0);

    //std::multimap<double,ldp_triangle_factor> candidateFactors;
    std::multimap<double,std::array<size_t,3>> candidates;


    for (size_t i = 0; i < nr_nodes(); ++i) {
        // std::cout<<"node "<<i<<std::endl;

        //TODO just half of the change for base!
        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        std::vector<double> minMarginalsIn=sncFactorIn->getAllBaseMinMarginals();
        std::vector<double> localBaseCostsIn=sncFactorIn->getBaseCosts();
        assert(baseGraph.numberOfEdgesToVertex(i)==minMarginalsIn.size());
        for (size_t j = 0; j < minMarginalsIn.size(); ++j) {
            minMarginalsIn[j]*=0.5;
            size_t edge=baseGraph.edgeToVertex(i,j);
            baseEdgeLabelsIn.at(edge)+=minMarginalsIn.at(j);
            localBaseCostsIn.at(j)-=minMarginalsIn.at(j);
        }

        //  std::cout<<"base min marginals in "<<i<<std::endl;



        std::vector<double> minMarginalsLiftedIn=sncFactorIn->getAllLiftedMinMarginals(&localBaseCostsIn);
        assert(liftedGraph.numberOfEdgesToVertex(i)==minMarginalsLiftedIn.size());
        for (size_t j = 0; j < minMarginalsLiftedIn.size(); ++j) {
            size_t edge=liftedGraph.edgeToVertex(i,j);
            liftedEdgeLabelsIn.at(edge)+=minMarginalsLiftedIn.at(j);
        }

        // std::cout<<"lifted min marginals in "<<i<<std::endl;

        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        std::vector<double> localBaseCostsOut=sncFactorOut->getBaseCosts();
        std::vector<double> minMarginalsOut=sncFactorOut->getAllBaseMinMarginals();
        assert(baseGraph.numberOfEdgesFromVertex(i)==minMarginalsOut.size());
        for (size_t j = 0; j < minMarginalsOut.size(); ++j) {
            minMarginalsOut[j]*=0.5;
            size_t edge=baseGraph.edgeFromVertex(i,j);
            baseEdgeLabelsOut.at(edge)+=minMarginalsOut.at(j);
            localBaseCostsOut.at(j)-=minMarginalsOut.at(j);
        }

        //  std::cout<<"base min marginals out "<<i<<std::endl;


        std::vector<double> minMarginalsLiftedOut=sncFactorOut->getAllLiftedMinMarginals(&localBaseCostsOut);
        assert(liftedGraph.numberOfEdgesFromVertex(i)==minMarginalsLiftedOut.size());
        for (size_t j = 0; j < minMarginalsLiftedOut.size(); ++j) {
            size_t edge=liftedGraph.edgeFromVertex(i,j);
            liftedEdgeLabelsOut.at(edge)+=minMarginalsLiftedOut.at(j);
        }

        // std::cout<<"lifted min marginals out "<<i<<std::endl;

    }

    //std::cout<<"edge scores obtained"<<std::endl;

    for (size_t vertex = 0; vertex < nr_nodes(); ++vertex) {
        size_t numberOfBaseOut=baseGraph.numberOfEdgesFromVertex(vertex);
        size_t numberOfBaseIn=baseGraph.numberOfEdgesToVertex(vertex);

        //base edges out
        for (size_t beOut = 0; beOut < baseGraph.numberOfEdgesFromVertex(vertex); ++beOut) {
            size_t vertexOut=baseGraph.vertexFromVertex(vertex,beOut);
            bool isStrongOut=instance.isStrongBase(vertex,vertexOut);
            size_t edgeOut=baseGraph.edgeFromVertex(vertex,beOut);
            double valueOut=baseEdgeLabelsOut.at(edgeOut)+baseEdgeLabelsIn.at(edgeOut);
            if((!isStrongOut&&valueOut>-eps)||(isStrongOut&&valueOut<eps&&valueOut>-eps )) continue;

            //base edges in
            for (size_t beIn = 0; beIn < baseGraph.numberOfEdgesToVertex(vertex); ++beIn) {
                if(usedTriangles.at(vertex).at(beIn).count(beOut)>0) continue;
                size_t vertexIn=baseGraph.vertexToVertex(vertex,beIn);
                bool isStrongIn=instance.isStrongBase(vertexIn,vertex);
                auto findEdge=liftedGraph.findEdge(vertexIn,vertexOut);
                if(findEdge.first){
                    size_t edgeIn=baseGraph.edgeToVertex(vertex,beIn);
                    double valueIn=baseEdgeLabelsIn.at(edgeIn)+baseEdgeLabelsOut.at(edgeIn);
                    double valueConnecting=liftedEdgeLabelsIn.at(findEdge.second)+liftedEdgeLabelsOut.at(findEdge.second);
                    // if((valueIn<eps&&valueIn>-eps)||(valueConnecting<eps&&valueConnecting>-eps)) continue;
                    if((valueIn<-eps&&valueConnecting>eps&&valueOut<-eps)||(valueIn<-eps&&valueConnecting<-eps&&valueOut>eps)||((valueIn>eps&&isStrongIn)&&valueConnecting<-eps&&valueOut<-eps)){
                        double improvement=std::min({std::abs(valueOut),std::abs(valueIn),std::abs(valueConnecting)});
                        std::array<size_t,3> toInsert={vertex,beIn,beOut};
                        candidates.emplace(improvement,toInsert);



                    }
                }

            }

            //lifted edges in
            for (size_t leIn = 0; leIn < liftedGraph.numberOfEdgesToVertex(vertex); ++leIn) {
                if(usedTriangles.at(vertex).at(leIn+numberOfBaseIn).count(beOut)>0) continue;
                size_t vertexIn=liftedGraph.vertexToVertex(vertex,leIn);

                auto findEdge=liftedGraph.findEdge(vertexIn,vertexOut);
                if(findEdge.first){
                    size_t edgeIn=liftedGraph.edgeToVertex(vertex,leIn);
                    double valueIn=liftedEdgeLabelsIn.at(edgeIn)+liftedEdgeLabelsOut.at(edgeIn);
                    double valueConnecting=liftedEdgeLabelsIn.at(findEdge.second)+liftedEdgeLabelsOut.at(findEdge.second);
                    if((valueIn<-eps&&valueConnecting>eps&&valueOut<-eps)||(valueIn<-eps&&valueConnecting<-eps&&valueOut>eps)||(valueIn>eps&&valueConnecting<-eps&&valueOut<-eps)){
                        double improvement=std::min({std::abs(valueOut),std::abs(valueIn),std::abs(valueConnecting)});

                        std::array<size_t,3> toInsert={vertex,leIn+numberOfBaseIn,beOut};
                        candidates.emplace(improvement,toInsert);
                    }

                }

            }

        }

        //lifted edges out
        for (size_t leOut = 0; leOut < liftedGraph.numberOfEdgesFromVertex(vertex); ++leOut) {
            size_t vertexOut=liftedGraph.vertexFromVertex(vertex,leOut);
            size_t edgeOut=liftedGraph.edgeFromVertex(vertex,leOut);
            double valueOut=liftedEdgeLabelsIn.at(edgeOut)+liftedEdgeLabelsOut.at(edgeOut);
            if(valueOut<eps&&valueOut>-eps) continue;

            //base edges in
            for (size_t beIn = 0; beIn < baseGraph.numberOfEdgesToVertex(vertex); ++beIn) {
                if(usedTriangles.at(vertex).at(beIn).count(leOut+numberOfBaseOut)>0) continue;
                size_t vertexIn=baseGraph.vertexToVertex(vertex,beIn);
                bool isStrongIn=instance.isStrongBase(vertexIn,vertex);
                auto findEdge=liftedGraph.findEdge(vertexIn,vertexOut);
                if(findEdge.first){
                    size_t edgeIn=baseGraph.edgeToVertex(vertex,beIn);
                    double valueIn=baseEdgeLabelsIn.at(edgeIn)+baseEdgeLabelsOut.at(edgeIn);
                    if(valueIn>-eps) continue;
                    double valueConnecting=liftedEdgeLabelsIn.at(findEdge.second)+liftedEdgeLabelsOut.at(findEdge.second);
                    if((valueIn<-eps&&valueConnecting>eps&&valueOut<-eps)||(valueIn<-eps&&valueConnecting<-eps&&valueOut>eps)||((valueIn>eps&&isStrongIn)&&valueConnecting<-eps&&valueOut<-eps)){
                        double improvement=std::min({std::abs(valueOut),std::abs(valueIn),std::abs(valueConnecting)});
                        std::array<size_t,3> toInsert={vertex,beIn,leOut+numberOfBaseOut};
                        candidates.emplace(improvement,toInsert);

                    }
                }

            }
            //lifted edges in
            for (size_t leIn = 0; leIn < liftedGraph.numberOfEdgesToVertex(vertex); ++leIn) {
                if(usedTriangles.at(vertex).at(leIn+numberOfBaseIn).count(leOut+numberOfBaseOut)>0) continue;
                size_t vertexIn=liftedGraph.vertexToVertex(vertex,leIn);
                auto findEdge=liftedGraph.findEdge(vertexIn,vertexOut);
                if(findEdge.first){
                    size_t edgeIn=liftedGraph.edgeToVertex(vertex,leIn);
                    double valueIn=liftedEdgeLabelsIn.at(edgeIn)+liftedEdgeLabelsOut.at(edgeIn);
                    double valueConnecting=liftedEdgeLabelsIn.at(findEdge.second)+liftedEdgeLabelsOut.at(findEdge.second);
                    if((valueOut<-eps&&valueIn<-eps&&valueConnecting>eps)||(valueOut<-eps&&valueIn>eps&&valueConnecting<-eps)||(valueOut>eps&&valueIn<-eps&&valueConnecting<-eps)){
                        double improvement=std::min({std::abs(valueOut),std::abs(valueIn),std::abs(valueConnecting)});
                        std::array<size_t,3> toInsert={vertex,leIn+numberOfBaseIn,leOut+numberOfBaseOut};
                        candidates.emplace(improvement,toInsert);

                    }

                }

            }


        }
    }

    std::cout<<"triangle factors candidates"<<std::endl;

    size_t counter=0;
    //for(auto iter=candidateFactors.rbegin();iter!=candidateFactors.rend()&&counter<nr_constraints_to_add;iter++){
    //    const ldp_triangle_factor& trFact=iter->second;
    //    //auto * triangleFactor =lp_->template add_factor<TRIANGLE_FACTOR_CONT>(trFact);
    //   // std::vector<double> costs={trFact.getEdgeCosts()[0],trFact.getEdgeCosts()[1],trFact.getEdgeCosts()[2]};
    //    std::array<double,3> costs={0,0,0};
    //    auto* triangleFactor = lp_->template add_factor<TRIANGLE_FACTOR_CONT>(trFact.getV1(),trFact.getV2(),trFact.getV3(),costs,trFact.isV1V2Base(),trFact.isV2V3Base());
    //    triangle_factors_.push_back(triangleFactor);
    //    counter++;
    //}

    double expectedImprovement=0;
    std::cout<<"candidate size "<<candidates.size()<<std::endl;


    size_t triangleFactorsOrigSize=triangle_factors_.size();
    for(auto iter=candidates.rbegin();iter!=candidates.rend()&&counter<nr_constraints_to_add;iter++){
       // std::cout<<"new candidate "<<std::endl;

        std::array<size_t,3>& triple =iter->second;
        std::array<double,3> costs={0,0,0};
        size_t v2=triple[0];
        assert(v2<nr_nodes());
        size_t v1;
        bool v1v2base=true;

        std::pair<bool,size_t> feV1V2;
        std::pair<bool,size_t> feV2V3;
        //  std::pair<bool,size_t> feV1V3;


        if(triple[1]>=baseGraph.numberOfEdgesToVertex(v2)){
            size_t liftedOrder=triple[1]-baseGraph.numberOfEdgesToVertex(v2);
            assert(liftedOrder<liftedGraph.numberOfEdgesToVertex(v2));
            v1=liftedGraph.vertexToVertex(v2,liftedOrder);
            v1v2base=false;
            feV1V2=liftedGraph.findEdge(v1,v2);
            assert(feV1V2.first);
            if(liftedEdgeUsed.at(feV1V2.second)){
               // std::cout<<"cont"<<std::endl;
                continue;
            }
            liftedEdgeUsed.at(feV1V2.second)=true;

        }
        else{
            v1=baseGraph.vertexToVertex(v2,triple[1]);
            feV1V2=baseGraph.findEdge(v1,v2);
            assert(feV1V2.first);
            if(baseEdgeUsed.at(feV1V2.second)){
              //  std::cout<<"cont"<<std::endl;
                continue;
            }
            baseEdgeUsed.at(feV1V2.second)=true;
        }

        assert(v2<nr_nodes());
        size_t v3;
        bool v2v3base=true;
        if(triple[2]>=baseGraph.numberOfEdgesFromVertex(v2)){
            size_t liftedOrder=triple[2]-baseGraph.numberOfEdgesFromVertex(v2);
            assert(liftedOrder<liftedGraph.numberOfEdgesFromVertex(v2));
            v3=liftedGraph.vertexFromVertex(v2,liftedOrder);
            v2v3base=false;
            feV2V3=liftedGraph.findEdge(v2,v3);
            assert(feV2V3.first);
            if(liftedEdgeUsed.at(feV2V3.second)){
              //  std::cout<<"cont"<<std::endl;
                continue;
            }
            liftedEdgeUsed.at(feV2V3.second)=true;
        }
        else{
            v3=baseGraph.vertexFromVertex(v2,triple[2]);
            feV2V3=baseGraph.findEdge(v2,v3);
            assert(feV2V3.first);
            if(baseEdgeUsed.at(feV2V3.second)){
               // std::cout<<"cont"<<std::endl;
                continue;
            }
            baseEdgeUsed.at(feV2V3.second)=true;

        }
       // std::cout<<"adding"<<std::endl;
        auto feV1V3=liftedGraph.findEdge(v1,v3);
        assert(feV1V3.first);
        if(liftedEdgeUsed.at(feV1V3.second)) continue;
        liftedEdgeUsed.at(feV1V3.second)=true;

        auto * sncFactorOutV1=single_node_cut_factors_[v1][1]->get_factor();
        auto * sncFactorInV2=single_node_cut_factors_[v2][0]->get_factor();
        auto * sncFactorOutV2=single_node_cut_factors_[v2][1]->get_factor();
        auto * sncFactorInV3=single_node_cut_factors_[v3][0]->get_factor();

        if(!v1v2base){

            size_t orderV2InV1=sncFactorOutV1->getLiftedIDToOrder(v2);
            size_t orderV1InV2=sncFactorInV2->getLiftedIDToOrder(v1);

            sncFactorOutV1->updateCostSimple(-liftedEdgeLabelsOut.at(feV1V2.second),orderV2InV1,true);
            sncFactorInV2->updateCostSimple(-liftedEdgeLabelsIn.at(feV1V2.second),orderV1InV2,true);
            costs[0]=liftedEdgeLabelsIn.at(feV1V2.second)+liftedEdgeLabelsOut.at(feV1V2.second);


        }
        else{

            size_t orderV2InV1=sncFactorOutV1->getBaseIDToOrder(v2);
            size_t orderV1InV2=sncFactorInV2->getBaseIDToOrder(v1);

            sncFactorOutV1->updateCostSimple(-baseEdgeLabelsOut.at(feV1V2.second),orderV2InV1,false);
            sncFactorInV2->updateCostSimple(-baseEdgeLabelsIn.at(feV1V2.second),orderV1InV2,false);
            costs[0]=baseEdgeLabelsIn.at(feV1V2.second)+baseEdgeLabelsOut.at(feV1V2.second);
        }


        if(!v2v3base){

            size_t orderV3InV2=sncFactorOutV2->getLiftedIDToOrder(v3);
            size_t orderV2InV3=sncFactorInV3->getLiftedIDToOrder(v2);


            sncFactorOutV2->updateCostSimple(-liftedEdgeLabelsOut.at(feV2V3.second),orderV3InV2,true);
            sncFactorInV3->updateCostSimple(-liftedEdgeLabelsIn.at(feV2V3.second),orderV2InV3,true);
            costs[1]=liftedEdgeLabelsIn.at(feV2V3.second)+liftedEdgeLabelsOut.at(feV2V3.second);
        }
        else{
            size_t orderV3InV2=sncFactorOutV2->getBaseIDToOrder(v3);
            size_t orderV2InV3=sncFactorInV3->getBaseIDToOrder(v2);


            sncFactorOutV2->updateCostSimple(-baseEdgeLabelsOut.at(feV2V3.second),orderV3InV2,false);
            sncFactorInV3->updateCostSimple(-baseEdgeLabelsIn.at(feV2V3.second),orderV2InV3,false);
            costs[1]=baseEdgeLabelsIn.at(feV2V3.second)+baseEdgeLabelsOut.at(feV2V3.second);
        }



        size_t orderV3InV1=sncFactorOutV1->getLiftedIDToOrder(v3);
        size_t orderV1InV3=sncFactorInV3->getLiftedIDToOrder(v1);


        sncFactorOutV1->updateCostSimple(-liftedEdgeLabelsOut.at(feV1V3.second),orderV3InV1,true);
        sncFactorInV3->updateCostSimple(-liftedEdgeLabelsIn.at(feV1V3.second),orderV1InV3,true);
        costs[2]=liftedEdgeLabelsIn.at(feV1V3.second)+liftedEdgeLabelsOut.at(feV1V3.second);

        assert(v3<nr_nodes());
        auto* triangleFactor = lp_->template add_factor<TRIANGLE_FACTOR_CONT>(v1,v2,v3,costs,v1v2base,v2v3base,instance);
        //    double lb=triangleFactor->LowerBound();
        //    double simpleLB=0;
        //    for (int i = 0; i < 3; ++i) {
        //        if(costs[i]<0){
        //            simpleLB+=costs[i];
        //        }
        //        else {
        //            simpleLB-=costs[i];
        //        }
        //    }



        expectedImprovement+=iter->first;
        triangle_factors_.push_back(triangleFactor);
        usedTriangles.at(triple[0]).at(triple[1]).insert(triple[2]);
        counter++;
    }

    std::cout<<"used triangles size "<<usedTriangles.size()<<std::endl;

    //TODO update cost in SNC factors that were used for creating the triangle factors!

    for (int i = triangleFactorsOrigSize; i < triangle_factors_.size(); ++i) {
        auto * triangleFactor=triangle_factors_[i];


        size_t v1=triangleFactor->get_factor()->getV1();
        size_t v2=triangleFactor->get_factor()->getV2();
        size_t v3=triangleFactor->get_factor()->getV3();

        assert(v1<nr_nodes()&&v2<nr_nodes()&&v3<nr_nodes());

        auto * sncOutV1=single_node_cut_factors_[v1][1];
        auto * sncOutV2=single_node_cut_factors_[v2][1];
        auto * sncInV2=single_node_cut_factors_[v2][0];
        auto * sncInV3=single_node_cut_factors_[v3][0];

        bool isv1v2Base=triangleFactor->get_factor()->isV1V2Base();
        bool isv2v3Base=triangleFactor->get_factor()->isV2V3Base();

        std::vector<size_t> trInd={0,2};
        size_t v2IndexInV1;
        if(isv1v2Base){
            v2IndexInV1=sncOutV1->get_factor()->getBaseIDToOrder(v2);
        }
        else{
            v2IndexInV1=sncOutV1->get_factor()->getLiftedIDToOrder(v2);
        }
        size_t v3IndexInV1=sncOutV1->get_factor()->getLiftedIDToOrder(v3);
        std::vector<size_t> vertices={v2IndexInV1,v3IndexInV1};
        std::vector<bool> isLifted={!isv1v2Base,true};

        auto * message1=lp_->template add_message<SNC_TRIANGLE_MESSAGE>(triangleFactor, sncOutV1, trInd, vertices,isLifted);
        snc_triangle_messages_.push_back(message1);

        trInd={0};
        size_t v1IndexInV2;
        if(isv1v2Base){
            v1IndexInV2=sncInV2->get_factor()->getBaseIDToOrder(v1);
        }
        else{
            v1IndexInV2=sncInV2->get_factor()->getLiftedIDToOrder(v1);
        }
        vertices={v1IndexInV2};
        isLifted={!isv1v2Base};
        auto * message2=lp_->template add_message<SNC_TRIANGLE_MESSAGE>(triangleFactor, sncInV2, trInd, vertices,isLifted);
        snc_triangle_messages_.push_back(message2);

        trInd={1};
        size_t v3IndexInV2;
        if(isv2v3Base){
            v3IndexInV2=sncOutV2->get_factor()->getBaseIDToOrder(v3);
        }
        else{
            v3IndexInV2=sncOutV2->get_factor()->getLiftedIDToOrder(v3);
        }
        vertices={v3IndexInV2};
        isLifted={!isv2v3Base};
        auto * message3=lp_->template add_message<SNC_TRIANGLE_MESSAGE>(triangleFactor, sncOutV2, trInd, vertices,isLifted);
        snc_triangle_messages_.push_back(message3);

        trInd={1,2};
        size_t v2IndexInV3;
        if(isv2v3Base){
            v2IndexInV3=sncInV3->get_factor()->getBaseIDToOrder(v2);
        }
        else{
            v2IndexInV3=sncInV3->get_factor()->getLiftedIDToOrder(v2);
        }
        size_t v1IndexInV3=sncInV3->get_factor()->getLiftedIDToOrder(v1);
        vertices={v2IndexInV3,v1IndexInV3};
        isLifted={!isv2v3Base,true};
        auto * message4=lp_->template add_message<SNC_TRIANGLE_MESSAGE>(triangleFactor, sncInV3, trInd, vertices,isLifted);
        snc_triangle_messages_.push_back(message4);

        //    auto * message=lp_->template add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE>(left_snc, right_snc, i, j);
        //    snc_lifted_messages_.push_back(message);



    }
    std::cout<<"tighten finished"<<std::endl;




    //    andres::graph::ComponentsBySearch<andres::graph::Digraph<>> components;
    //    struct ForLabels{
    //        ForLabels(const std::vector<double>& costs):myCosts(costs){}
    //        bool vertex(size_t i) const {return true;}
    //        bool edge(size_t i) const {
    //            return myCosts.at(i)<0;
    //        }
    //        const std::vector<double>& myCosts;
    //    };
    //    std::vector<double> costs(nr_nodes(),0);
    //    //TODO initialize costs with lifted min marginals (evt. base min marginals)
    //    andres::graph::Digraph<> graph; //TODO initialize this: Actually, it is needed lifted with some base edges. Maybe copy of lifted, add some edges?
    //    components.build(graph,costs);
    //    //for all edges with positive costs whose vertices belong to the same component: add triangle (priority queue)
    std::cout<<"tighten expected improvement "<<expectedImprovement<<", added constraints "<<counter<<std::endl;
    adjustTriangleLabels(triangleFactorsOrigSize);
    return counter;

}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::mcf_node_to_graph_node(std::size_t i) const
{
    assert(i < mcf_->no_nodes());
    if(i == mcf_source_node())
        return base_graph_source_node();
    if(i == mcf_terminal_node())
        return base_graph_terminal_node();
    else
        return i/2;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE, class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE, SNC_TRIANGLE_MESSAGE>::read_in_mcf_costs(const bool change_marginals)
{
    mcf_->reset_costs();

    const lifted_disjoint_paths::LdpInstance &instance=single_node_cut_factors_[0][0]->get_factor()->getLdpInstance();
    const andres::graph::Digraph<>& baseGraph=instance.getGraph();
    std::vector<std::unordered_map<size_t,double>> costsFromTriangles(nr_nodes());


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


    for(std::size_t i=0; i<nr_nodes(); ++i)
    {
        {
           // std::cout<<"incomming primal "<<std::endl;
            auto *incoming_snc = single_node_cut_factors_[i][0]->get_factor();
            const auto incoming_min_marg = incoming_snc->getAllBaseMinMarginals();
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
                //std::cout << "incoming min marginal edge " << i << "," << j << "\n";
                while (mcf_node_to_graph_node(mcf_->head(e)) != j)
                {
                    //std::cout << "mcf edge " << e << ": " << mcf_->tail(e) << "," << mcf_->head(e) << "; base edge: " << mcf_node_to_graph_node(mcf_->tail(e)) << "," << mcf_node_to_graph_node(mcf_->head(e)) << "\n";
                    assert(mcf_->tail(e) == incoming_mcf_node(i));
                    ++e;
                }
                const std::size_t start_node = mcf_->tail(e);
                assert(mcf_->tail(e) == incoming_mcf_node(i));
                const double m = incoming_min_marg[l];
               // std::cout<<i<<","<<j<<": "<<m<<std::endl;
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
                    incoming_snc->updateCostSimple(-incoming_min_marg[l], l, false);
                }
            }
        }

        {
           // std::cout<<"outgoing primal "<<std::endl;
            auto *outgoing_snc = single_node_cut_factors_[i][1]->get_factor();
            const auto outgoing_min_marg = outgoing_snc->getAllBaseMinMarginals();
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
               //  std::cout<<i<<","<<j<<": "<<m<<std::endl;
                assert(mcf_->lower_bound(e) == 0 && mcf_->upper_bound(e) == 1);
                mcf_->update_cost(e, m);
                auto iter=costsFromTriangles[i].find(j);
                if(iter!=costsFromTriangles[i].end()){
                    mcf_->update_cost(e,iter->second);
                }
            }
            if(change_marginals)
            {
                for (std::size_t l = 0; l < outgoing_min_marg.size(); ++l)
                {
                    outgoing_snc->updateCostSimple(-outgoing_min_marg[l], l, false);
                }
            }
        }
    }
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::write_back_mcf_costs()
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
                    //incoming_snc->updateCostSimple(- reduced_cost, base_graph_node(incoming_node), false);
                    //assert(vertex_index == incoming_snc->getBaseIDs().size()-1);
                    assert(incoming_snc->getBaseIDs()[vertex_index] == base_graph_node(incoming_node));
                    incoming_snc->updateCostSimple(- reduced_cost, vertex_index , false);
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
                    incoming_snc->updateCostSimple(-0.5 * reduced_cost, vertex_index, false);
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
                    //outgoing_snc->updateCostSimple(reduced_cost, base_graph_node(outgoing_node), false);
                    assert(outgoing_snc->getBaseIDs()[vertex_index] == base_graph_node(outgoing_node));
                    outgoing_snc->updateCostSimple(reduced_cost, vertex_index , false);
                    ++vertex_index;
                }
                else if (base_graph_node(outgoing_node) != i) // regular edge
                {
                    assert(flow == 0 || flow == 1);
                    if(flow == 0)
                        assert(reduced_cost >= -1e-8);
                    if (flow == 1)
                        assert(reduced_cost <= 1e-8);
                    //outgoing_snc->updateCostSimple(0.5 * reduced_cost, base_graph_node(outgoing_node), false);
                    assert(outgoing_snc->getBaseIDs()[vertex_index] == base_graph_node(outgoing_node));
                    outgoing_snc->updateCostSimple(0.5 * reduced_cost, vertex_index, false);
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

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR,class TRIANGLE_FACTOR_CONT, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_TRIANGLE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, TRIANGLE_FACTOR_CONT, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_TRIANGLE_MESSAGE>::reparametrize_snc_factors()
{
    const double primal_cost_before = this->lp_->EvaluatePrimal();
    read_in_mcf_costs(true);
    mcf_->solve();
    std::cout << "mcf cost = " << mcf_->objective() << "\n";
    write_back_mcf_costs();
    const double primal_cost_after = this->lp_->EvaluatePrimal();
    std::cout << "primal cost before = " << primal_cost_before << ", primal cost after = " << primal_cost_after << "\n";
    assert(std::abs(primal_cost_before - primal_cost_after) <= 1e-6);
}


}
