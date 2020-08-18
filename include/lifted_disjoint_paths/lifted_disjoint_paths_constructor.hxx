#pragma once

#include "LP.h"
#include "solver.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "MCF-SSP/mcf_ssp.hxx"
#include <unordered_map>

#include <memory>

namespace LPMP {

template<class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
class lifted_disjoint_paths_constructor
{
public:
    using FMC = FACTOR_MESSAGE_CONNECTION;
    template<typename SOLVER>
    lifted_disjoint_paths_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

    //void construct(const lifted_disjoint_paths_instance& i);
    void construct(const lifted_disjoint_paths::LdpInstance& instance);

    void ComputePrimal();

    void Tighten(const std::size_t nr_constraints_to_add);
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

    bool checkFeasibilityInSnc();
    bool checkFeasibilityLiftedInSnc();
    bool checkFeasibilityBaseInSnc();

    void sncDebug(size_t vertex,bool isOut);

    LP<FMC> *lp_;
    using mcf_solver_type = MCF::SSP<long, double>;
    std::unique_ptr<mcf_solver_type> mcf_; // minimum cost flow factor for base edges
    std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>> single_node_cut_factors_;
    std::vector<SINGLE_NODE_CUT_LIFTED_MESSAGE*> snc_lifted_messages_;
    std::vector<SNC_NODE_MESSAGE*> snc_node_messages_;

};

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE, class SNC_NODE_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE, SNC_NODE_MESSAGE>::base_graph_node(const std::size_t mcf_node) const
{
    assert(mcf_node < mcf_->no_nodes());
    if(mcf_node == mcf_source_node())
        return base_graph_source_node();
    if (mcf_node == mcf_terminal_node())
        return base_graph_terminal_node();
    return mcf_node / 2;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::checkFeasibilityInSnc(){
    bool isFeasible=checkFeasibilityBaseInSnc();
    if(isFeasible){
        isFeasible=checkFeasibilityLiftedInSnc();
    }
    return isFeasible;
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::checkFeasibilityLiftedInSnc(){
  //TODO check this
//    bool isFeasible=true;
//    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
//        size_t vertex=i;
//        const auto* snc=single_node_cut_factors_[vertex][1]->get_factor();
//        size_t primalBaseIndex=snc->getPrimalBase();
//        if(snc->baseIDs[primalBaseIndex]==base_graph_terminal_node()){

//            std::vector<bool> isOnPath(nr_nodes(),0);
//            isOnPath[vertex]=1;
//            std::list<size_t> path;
//            //path.push_front(vertex);

//            while(vertex!=base_graph_source_node()){
//                const auto* sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
//                const std::unordered_set<size_t>& indicesOfActiveEndpoints=sncFactorOut->getPrimalLiftedIndices();
//                std::vector<bool> activeVector(sncFactorOut->getLiftedCosts().size());
//                for(size_t v:indicesOfActiveEndpoints){
//                    activeVector[v]=1;
//                }
//                for (int j = 0; j < activeVector.size(); ++j) {
//                    size_t vertex2=sncFactorOut->getLiftedIDs()[j];
//                    if(activeVector[j]!=isOnPath[vertex2]){
//                        isFeasible=false;
//                        break;
//                    }

//                }

//                if(!isFeasible) break;


//                isOnPath[vertex]=1;
//                path.push_front(vertex);
//                auto* sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
//                size_t ind=sncFactorIn->getPrimalBase();
//                vertex=sncFactorIn.baseIDs[ind];

//            }
//            if(isFeasible){
//                for(size_t activeVertex:path){
//                    auto* sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();


//                    const std::unordered_set<size_t>& indicesOfActiveEndpoints=sncFactorIn.getPrimalLiftedIndices();
//                    std::vector<bool> activeVector(sncFactorIn->getLiftedCosts().size());
//                    for(size_t v:indicesOfActiveEndpoints){
//                        activeVector[v]=1;
//                    }
//                    for (int j = 0; j < activeVector.size(); ++j) {
//                        size_t vertex2=sncFactorIn->getLiftedIDs()[j];
//                        if(activeVector[j]!=isOnPath[vertex2]){
//                            isFeasible=false;
//                            break;
//                        }

//                    }

//                }
//            }
//        }
//    }
//    return isFeasible;

}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::checkFeasibilityBaseInSnc(){

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



template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::adjustLiftedLabels(){
 //TODO fix this
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
    }
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::construct(const lifted_disjoint_paths::LdpInstance &instance)
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

    //lp_->add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>(left_snc, right_snc, left_node, right_node);

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
    //        	auto * message=lp_->template add_message<SNC_NODE_MESSAGE>(left_snc, right_snc, i);
    //        	snc_node_messages_.push_back(message);
    //        }

    //sncDebug(20,1);
    if(debug()) std::cout<<"messages added"<<std::endl;
}


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::sncDebug(size_t vertex,bool isOut)
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


template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::ComputePrimal()
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
                   // primalValue+=pSNC->EvaluatePrimal();
                }
                ++vertex_index;
            }
        }
    }
    const bool isFeasible=this->checkFeasibilityBaseInSnc();
    std::cout<<"checked feasibility: "<<isFeasible<<std::endl;
    assert(isFeasible);
    adjustLiftedLabels();
    double primalValue=0;
    for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {

        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        primalValue+=sncFactorIn->EvaluatePrimal();
        primalValue+=sncFactorOut->EvaluatePrimal();


    }

    std::cout<<"primal value: "<<primalValue<<std::endl;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::Tighten(const std::size_t nr_constraints_to_add)
{
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::mcf_node_to_graph_node(std::size_t i) const
{
    assert(i < mcf_->no_nodes());
    if(i == mcf_source_node())
        return base_graph_source_node();
    if(i == mcf_terminal_node())
        return base_graph_terminal_node();
    else
        return i/2;
}

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE, class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE, SNC_NODE_MESSAGE>::read_in_mcf_costs(const bool change_marginals)
{
    mcf_->reset_costs();

    for(std::size_t i=0; i<nr_nodes(); ++i)
    {
        {
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
                assert(mcf_->lower_bound(e) == 0 && mcf_->upper_bound(e) == 1);
                mcf_->update_cost(e, m);
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

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::write_back_mcf_costs()
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

template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::reparametrize_snc_factors()
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
