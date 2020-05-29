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

        private:
            std::size_t mcf_node_to_graph_node(std::size_t i) const;
            void prepare_mcf_costs();
            void write_back_mcf_costs();
            void reparametrize_snc_factors();

            std::size_t nr_nodes() const { assert(single_node_cut_factors_.size() == (mcf_->no_nodes() - 2) / 2); return single_node_cut_factors_.size(); }
            std::size_t incoming_mcf_node(const std::size_t i) const { assert(i < nr_nodes()); return i*2; }
            std::size_t outgoing_mcf_node(const std::size_t i) const { assert(i < nr_nodes()); return i*2+1; }
            std::size_t mcf_source_node() const { return mcf_->no_nodes()-2; }
            std::size_t mcf_terminal_node() const { return mcf_->no_nodes()-1; }
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
    	bool isFeasible=true;
    	for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
    		size_t vertex=i;
    		const auto snc=single_node_cut_factors_[vertex][1]->get_factor();
    		size_t primalBaseIndex=snc->getPrimalBase();
    		if(snc->baseIDs[primalBaseIndex]==base_graph_terminal_node()){

    			std::vector<bool> isOnPath(nr_nodes(),0);
    			isOnPath[vertex]=1;
    			std::list<size_t> path;
    			path.push_front(vertex);

    			while(vertex!=base_graph_source_node()){
    				const auto& sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
    				const std::unordered_set<size_t>& indicesOfActiveEndpoints=sncFactorOut.getPrimalLiftedIndices();
    				std::vector<bool> activeVector(sncFactorOut.getLiftedCosts().size());
    				for(size_t v:indicesOfActiveEndpoints){
    					activeVector[v]=1;
    				}
    				for (int j = 0; j < activeVector.size(); ++j) {
    					size_t vertex2=sncFactorOut->liftedIDs[j];
    					if(activeVector[j]!=isOnPath[vertex2]){
    						isFeasible=false;
    						break;
    					}

    				}

    				if(!isFeasible) break;


    				isOnPath[vertex]=1;
    				path.push_front(vertex);
    				auto& sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
    				size_t ind=sncFactorIn.getPrimalBase();
    				vertex=sncFactorIn.baseIDs[ind];

    			}
    			if(isFeasible){
    				for(size_t activeVertex:path){
    					auto& sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();


    					const std::unordered_set<size_t>& indicesOfActiveEndpoints=sncFactorIn.getPrimalLiftedIndices();
    					std::vector<bool> activeVector(sncFactorIn.getLiftedCosts().size());
    					for(size_t v:indicesOfActiveEndpoints){
    						activeVector[v]=1;
    					}
    					for (int j = 0; j < activeVector.size(); ++j) {
    						size_t vertex2=sncFactorIn->liftedIDs[j];
    						if(activeVector[j]!=isOnPath[vertex2]){
    							isFeasible=false;
    							break;
    						}

    					}

    				}
    			}
    		}
    	}
    	return isFeasible;

    }


    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
    bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::checkFeasibilityBaseInSnc(){

    	//Check flow conservation
    	bool isFeasible=true;
    	std::unordered_map<size_t,size_t> predecessors; //Used later for lifted edges consistency
    	for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
    		auto *incoming_snc = single_node_cut_factors_[i][0]->get_factor();
    		const auto incoming_min_marg = incoming_snc->getAllBaseMinMarginals();

    		const auto sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
    		const auto sncFactorOut=single_node_cut_factors_[i][1]->get_factor();

    		if(sncFactorIn->isNodeActive()!=sncFactorOut->isNodeActive()){
    			isFeasible=false;

    		}
    		else if(sncFactorIn->isNodeActive()&&sncFactorOut->isNodeActive()){
    			size_t inputNode=sncFactorIn->getPrimalBaseVertexID();
    			if(inputNode!=base_graph_source_node()){
    				auto sncFactorInputNodeOut=single_node_cut_factors_[inputNode][1]->get_factor();
    				if(sncFactorInputNodeOut->getPrimalBaseVertexID()!=i){
    					isFeasible=false;

    				}
    			}

    			size_t outputNode=sncFactorOut->getPrimalBaseVertexID();
    			if(outputNode!=base_graph_terminal_node()){
    				auto sncFactorOuputNodeIn=single_node_cut_factors_[outputNode][0]->get_factor();
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

    	for (int i = 0; i < nr_nodes(); ++i) {
    		size_t vertex=i;
    		auto& snc=single_node_cut_factors_[vertex][1]->get_factor();
    		if(snc.getPrimalBaseVertexID()==base_graph_terminal_node()){

    			std::vector<bool> isOnPath(nr_nodes(),0);
    			isOnPath[vertex]=1;
    			std::list<size_t> path;
    			path.push_front(vertex);

    			while(vertex!=base_graph_source_node()){
    				auto& sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
    				std::unordered_set<size_t> activeEndpointIndices;
    				const std::vector<size_t>& liftedIDs= sncFactorOut->liftedIDs;
    				for (int i = 0; i < liftedIDs.size(); ++i) {
    					if(isOnPath[liftedIDs[i]]) activeEndpointIndices.insert(i);
					}

    				sncFactorOut.setPrimalLifted(activeEndpointIndices);
    				isOnPath[vertex]=1;
    				path.push_front(vertex);
    				auto& sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
    				vertex=sncFactorIn.getPrimalBaseVertexID();
    			}
    			for(size_t activeVertex:path){
    				auto& sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();
    				std::unordered_set<size_t> activeEndpointIndices;
    				const std::vector<size_t>& liftedIDs= sncFactorIn->liftedIDs;
    				for (int i = 0; i < liftedIDs.size(); ++i) {
    					if(isOnPath[liftedIDs[i]]) activeEndpointIndices.insert(i);
    				}
    				sncFactorIn.setPrimalLifted(activeEndpointIndices);
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
        std::cout << "source = " << base_graph_source << ", terminal = " << base_graph_terminal << "\n";


        const std::size_t nr_mcf_nodes = 2*nr_base_graph_nodes + 2; // source/terminal vertex + 2*ordinary vertices to ensure unit capacity vertices
        const std::size_t nr_mcf_edges = 3*nr_base_graph_nodes + instance.getGraph().numberOfEdges() + 1; // appearance/disappearance/uniqueness edge + connection edges + source/terminal edge
        mcf_ = std::make_unique<mcf_solver_type>(nr_mcf_nodes, nr_mcf_edges);

        const std::size_t mcf_source_node = nr_mcf_nodes - 2;
        const std::size_t mcf_terminal_node = nr_mcf_nodes - 1;
        auto incoming_edge = [](const std::size_t i) { return 2*i; };
        auto outgoing_edge = [](const std::size_t i) { return 2*i+1; };

        mcf_->add_edge(mcf_source_node, mcf_terminal_node, 0, nr_base_graph_nodes, 0.0);

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
               // std::cout << "base graph edge " << i << "," << j << "\n";
                if (j != base_graph_source && j != base_graph_terminal)
                {
                    assert(i<j);
                    const double c = instance.getEdgeScore(e);
                    mcf_->add_edge(outgoing_edge(i), incoming_edge(j), 0, 1, c);
                }
                assert(i < j);
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
           	auto left_snc=single_node_cut_factors_[vertex1][1];
           	auto left_snc_factor=left_snc->get_factor();
           //	std::cout<<"vertex 1 "<<vertex1<<std::endl;
           	for (int j = 0; j < left_snc_factor->liftedIDs.size(); ++j) {
				size_t vertex2=left_snc_factor->liftedIDs[j];
				//std::cout<<"vertex 2 "<<vertex2<<std::endl;
				assert(instance.getGraphLifted().findEdge(vertex1,vertex2).first);
				if(vertex2==this->base_graph_source_node()||vertex2==this->base_graph_source_node()) continue;
				auto right_snc=single_node_cut_factors_[vertex2][0];
				size_t i=right_snc->get_factor()->getLiftedOrderToID(vertex1);
	       		auto * message=lp_->template add_message<SINGLE_NODE_CUT_LIFTED_MESSAGE>(left_snc, right_snc, i, j);
	        		snc_lifted_messages_.push_back(message);
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
    	auto sncFactor=single_node_cut_factors_[vertex][isOut]->get_factor();
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
    		std::cout<<i<<": "<<(sncFactor->liftedIDs[i])<<": ";
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
        prepare_mcf_costs();
        mcf_->solve();
        double primalValue=0;
        for (int i = 0; i < mcf_->no_nodes(); ++i) {
			std::size_t first_out=mcf_->first_outgoing_arc(i);
			std::size_t graphNode=this->mcf_node_to_graph_node(i);
			bool foundOutput=false;
			for (int j = 0; j < mcf_->no_outgoing_arcs(i); ++j) {
				std::size_t edgeId=j+first_out;
				auto flow=mcf_->flow(edgeId);
				if(flow>0.5){
					foundOutput=true;

					std::size_t node=mcf_->head(edgeId);
					std::size_t graphNode2=mcf_node_to_graph_node(node);
					std::cout<<"active edge "<<graphNode<<", "<<graphNode2<<std::endl;
					if(graphNode!=graphNode2){
						if(i!=this->mcf_source_node()){
							auto pSNC=single_node_cut_factors_[graphNode][1]->get_factor();
							pSNC->setBaseEdgeActive(graphNode2);
							primalValue+=pSNC->EvaluatePrimal();
						}
						if(node!=this->mcf_terminal_node()){
							auto pSNC=single_node_cut_factors_[graphNode2][0]->get_factor();
							pSNC->setBaseEdgeActive(graphNode);
							primalValue+=pSNC->EvaluatePrimal();
						}
					}
				}

			}
			if(!foundOutput&&i!=mcf_terminal_node()&&i!=mcf_source_node()){
				std::cout<<"inactive node "<<graphNode<<std::endl;
				single_node_cut_factors_[graphNode][1]->get_factor()->setNoBaseEdgeActive();
				single_node_cut_factors_[graphNode][0]->get_factor()->setNoBaseEdgeActive();
			}
		}
        bool isFeasible=this->checkFeasibilityBaseInSnc();
        std::cout<<"checked feasibility: "<<isFeasible<<std::endl;
        assert(isFeasible);
        std::cout<<"primal value: "<<primalValue<<std::endl;



        // propagate base solution to single node cut factors
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
            

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::prepare_mcf_costs()
    {
        mcf_->reset_costs();

        for(std::size_t i=0; i<nr_nodes(); ++i)
        {
            {
                auto *incoming_snc = single_node_cut_factors_[i][0];
                const auto incoming_min_marg = incoming_snc->get_factor()->getAllBaseMinMarginals();
                std::cout << "incoming_min_marg size = " << incoming_min_marg.size() << "\n";
                std::vector<std::tuple<std::size_t, double>> incoming_min_margs_sorted;
                for (const auto [node, m] : incoming_min_marg)
                    incoming_min_margs_sorted.push_back({node, m});
                std::sort(incoming_min_margs_sorted.begin(), incoming_min_margs_sorted.end(), [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(b); });
                std::size_t e = mcf_->first_outgoing_arc(incoming_mcf_node(i));
                assert(incoming_min_margs_sorted.size() <= mcf_->no_outgoing_arcs(incoming_mcf_node(i)));
                for (std::size_t l = 0; l < incoming_min_margs_sorted.size(); ++l)
                {
                    const std::size_t j = std::get<0>(incoming_min_margs_sorted[l]);
                    assert(j != base_graph_terminal_node());
                    std::cout << "incoming min marginal edge " << i << "," << j << "\n";
                    while (mcf_node_to_graph_node(mcf_->head(e)) != j)
                    {
                        std::cout << "mcf edge " << e << ": " << mcf_->tail(e) << "," << mcf_->head(e) << "; base edge: " << mcf_node_to_graph_node(mcf_->tail(e)) << "," << mcf_node_to_graph_node(mcf_->head(e)) << "\n";
                        assert(mcf_->tail(e) == incoming_mcf_node(i));
                        ++e;
                    }
                    const std::size_t start_node = mcf_->tail(e);
                    assert(mcf_->tail(e) == incoming_mcf_node(i));
                    const double m = std::get<1>(incoming_min_margs_sorted[l]);
                    if (j != base_graph_source_node()){
                    	mcf_->update_cost(e, -m);
                    }
                    else{
                    	mcf_->update_cost(e, m);
                    }
                        //std::cout << "updates mcf edge " << e << ", reverse edge = " << mcf_->sister(e) << "\n";
                }
            }

            {
                auto *outgoing_snc = single_node_cut_factors_[i][1];
                const auto outgoing_min_marg = outgoing_snc->get_factor()->getAllBaseMinMarginals();
                std::cout << "outgoing_min_marg size = " << outgoing_min_marg.size() << "\n";
                std::vector<std::tuple<std::size_t, double>> outgoing_min_margs_sorted;
                for (const auto [node, m] : outgoing_min_marg)
                    outgoing_min_margs_sorted.push_back({node, m});
                std::sort(outgoing_min_margs_sorted.begin(), outgoing_min_margs_sorted.end(), [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(b); });
                std::size_t e = mcf_->first_outgoing_arc(outgoing_mcf_node(i));
                assert(outgoing_min_margs_sorted.size() <= mcf_->no_outgoing_arcs(outgoing_mcf_node(i)));
                for (std::size_t l = 0; l < outgoing_min_margs_sorted.size(); ++l)
                {
                    const std::size_t j = std::get<0>(outgoing_min_margs_sorted[l]);
                    assert(j != base_graph_source_node());
                    std::cout << "outgoing min marginal edge " << i << "," << std::get<0>(outgoing_min_margs_sorted[l]) << "\n";
                    while (mcf_node_to_graph_node(mcf_->head(e)) != std::get<0>(outgoing_min_margs_sorted[l]))
                    {
                        std::cout << "mcf edge: " << mcf_->tail(e) << "," << mcf_->head(e) << "; base edge: " << mcf_node_to_graph_node(mcf_->tail(e)) << "," << mcf_node_to_graph_node(mcf_->head(e)) << "\n";
                        assert(mcf_->tail(e) == outgoing_mcf_node(i));
                        ++e;
                    }
                    const std::size_t start_node = mcf_->tail(e);
                    assert(mcf_->tail(e) == outgoing_mcf_node(i));
                   // std::cout << "outgoing updates mcf edge " << e << ", reverse edge = " << mcf_->sister(e) << "\n";
                    mcf_->update_cost(e, std::get<1>(outgoing_min_margs_sorted[l]));
                }
            }
        }
    }

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::write_back_mcf_costs()
    {
        for(std::size_t i=0; i<nr_nodes(); ++i)
        {
            {
                auto *incoming_snc = single_node_cut_factors_[i][0]->get_factor();
            }
        } 
    }

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE,class SNC_NODE_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE,SNC_NODE_MESSAGE>::reparametrize_snc_factors()
    {
        prepare_mcf_costs();
        mcf_->solve();
        write_back_mcf_costs();
    } 

}
