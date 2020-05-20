#pragma once

#include "LP.h"
#include "solver.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "MCF-SSP/mcf_ssp.hxx"
#include <unordered_map>

#include <memory>

namespace LPMP {

    template<class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
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

    };

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
        bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::checkFeasibilityInSnc(){
    	bool isFeasible=checkFeasibilityBaseInSnc();
    	if(isFeasible){
    		isFeasible=checkFeasibilityLiftedInSnc();
    	}
    	return isFeasible;
    }


    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
        bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::checkFeasibilityLiftedInSnc(){
    	bool isFeasible=true;
    	for (int i = 0; i < nr_nodes(); ++i) {
    		size_t vertex=i;
    		const auto& snc=single_node_cut_factors_[vertex][1]->get_factor();
    		if(snc.getPrimalBase()==base_graph_terminal_node()){

    			std::vector<bool> isOnPath(nr_nodes(),0);
    			isOnPath[vertex]=1;
    			std::list<size_t> path;
    			path.push_front(vertex);

    			while(vertex!=base_graph_source_node()){
    				const auto& sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
    				const std::unordered_set<size_t>& activeEndpoints=sncFactorOut.getPrimalLifted();
    				const std::unordered_map<size_t,double>& liftedCosts= sncFactorOut.getLiftedCosts();
    				for(auto pair:liftedCosts){
    					size_t vertex2=pair.first;
    					if(isOnPath[vertex2]&&activeEndpoints.count(vertex2)==0){
    						isFeasible=false;
    						break;
    					}
    				}
    				for(size_t vertex2:activeEndpoints){
    					if(!isOnPath[vertex2]){
    						isFeasible=false;
    						break;
    					}
    				}
    				isOnPath[vertex]=1;
    				path.push_front(vertex);
    				auto& sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
    				vertex=sncFactorIn.getPrimalBase();
    			}
    			for(size_t activeVertex:path){
    				auto& sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();
    				const std::unordered_set<size_t>& activeEndpoints=sncFactorIn.getPrimalLifted();
    				const std::unordered_map<size_t,double>& liftedCosts= sncFactorIn.getLiftedCosts();
    				for(auto pair:liftedCosts){
    					size_t vertex2=pair.first;
    					if(isOnPath[vertex2]&&activeEndpoints.count(vertex2)==0){
    						isFeasible=false;
    						break;
    					}
    				}
    				for(size_t vertex2:activeEndpoints){
    					if(!isOnPath[vertex2]){
    						isFeasible=false;
    						break;
    					}
    				}
    			}
    		}
    	}
    	return isFeasible;

    }


    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    bool lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::checkFeasibilityBaseInSnc(){

    	//Check flow conservation
    	bool isFeasible=true;
    	std::unordered_map<size_t,size_t> predecessors; //Used later for lifted edges consistency
    	for (int i = 0; i < nr_nodes()&&isFeasible; ++i) {
    		auto *incoming_snc = single_node_cut_factors_[i][0]->get_factor();
    		const auto incoming_min_marg = incoming_snc->get_factor()->getAllBaseMinMarginals();

    		const auto sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
    		const auto sncFactorOut=single_node_cut_factors_[i][1]->get_factor();

    		if(sncFactorIn.isNodeActive()!=sncFactorOut.isNodeActive()){
    			isFeasible=false;

    		}
    		else if(sncFactorIn.isNodeActive()&&sncFactorOut.isNodeActive()){
    			size_t inputNode=sncFactorIn.getPrimalBase();
    			if(inputNode!=base_graph_source_node()){
    				auto& sncFactorInputNodeOut=single_node_cut_factors_[inputNode][1]->get_factor();
    				if(sncFactorInputNodeOut.getPrimalBase()!=i){
    					isFeasible=false;

    				}
    			}

    			size_t outputNode=sncFactorOut.getPrimalBase();
    			if(outputNode!=base_graph_terminal_node()){
    				auto& sncFactorOuputNodeIn=single_node_cut_factors_[outputNode][0]->get_factor();
    				if(sncFactorOuputNodeIn.getPrimalBase()!=i){
    					isFeasible=false;

    				}
    			}
    		}
    	}


    	return isFeasible;
    }



    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::adjustLiftedLabels(){

    	for (int i = 0; i < nr_nodes(); ++i) {
    		size_t vertex=i;
    		auto& snc=single_node_cut_factors_[vertex][1]->get_factor();
    		if(snc.getPrimalBase()==base_graph_terminal_node()){

    			std::vector<bool> isOnPath(nr_nodes(),0);
    			isOnPath[vertex]=1;
    			std::list<size_t> path;
    			path.push_front(vertex);

    			while(vertex!=base_graph_source_node()){
    				auto& sncFactorOut=single_node_cut_factors_[vertex][1]->get_factor();
    				std::unordered_set<size_t> activeEndpoints;
    				const std::unordered_map<size_t,double>& liftedCosts= sncFactorOut.getLiftedCosts();
    				for(auto pair:liftedCosts){
    					size_t vertex2=pair.first;
    					if(isOnPath[vertex2]){
    						activeEndpoints.insert(vertex2);
    					}
       				}
    				sncFactorOut.setPrimalLifted(activeEndpoints);
    				isOnPath[vertex]=1;
    				path.push_front(vertex);
    				auto& sncFactorIn=single_node_cut_factors_[vertex][0]->get_factor();
    				vertex=sncFactorIn.getPrimalBase();
    			}
    			for(size_t activeVertex:path){
    				auto& sncFactorIn=single_node_cut_factors_[activeVertex][0]->get_factor();
    				std::unordered_set<size_t> activeEndpoints;
    				const std::unordered_map<size_t,double>& liftedCosts= sncFactorIn.getLiftedCosts();
    				for(auto pair:liftedCosts){
    					size_t vertex2=pair.first;
    					if(isOnPath[vertex2]){
    						activeEndpoints.insert(vertex2);
    					}
    				}
    				sncFactorIn.setPrimalLifted(activeEndpoints);
    			}
    		}
    	}
    }

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::construct(const lifted_disjoint_paths::LdpInstance &instance)
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
                std::cout << "base graph edge " << i << "," << j << "\n";
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


        size_t nodeId=200;
      //  sncDebug(nodeId,1);
        //checkFeasibilityBaseInSnc();

    }


    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::sncDebug(size_t vertex,bool isOut)
    {
    	std::cout<<"print node "<<vertex<<std::endl;
    	auto sncFactor=single_node_cut_factors_[vertex][isOut]->get_factor();
    	const std::unordered_map<size_t, double>& baseCosts=sncFactor->getBaseCosts();
    	const std::unordered_map<size_t, double>& liftedCosts=sncFactor->getLiftedCosts();

    	std::cout<<"base costs"<<std::endl;
    	for(const auto pair:baseCosts){
    		std::cout<<pair.first<<": ";
    		std::cout<<pair.second;
    		std::cout<<std::endl;
    	}

    	std::cout<<"lifted costs"<<std::endl;
    	for(const auto pair:liftedCosts){
    		std::cout<<pair.first<<": ";
    		std::cout<<pair.second;
    		std::cout<<std::endl;
    	}


    	sncFactor->getAllBaseMinMarginals();
    	std::cout<<"isolated min marginals "<<std::endl;
    	for(const auto pair:baseCosts){
    		if(pair.first!=sncFactor->nodeID){

    			double value=sncFactor->getOneBaseEdgeMinMarginal(pair.first);
    			std::cout<<pair.first<<": "<<value<<std::endl;
    		}
    	}

    //	std::cout<<"computing lifted min marginal "<<std::endl;
    	//double value=sncFactor->oneLiftedMinMarginal(171);
    	//std::cout<<"lifted min marginal "<<value<<std::endl;
    	sncFactor->getAllLiftedMinMarginals();

    }


    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::ComputePrimal()
    {
        prepare_mcf_costs();
        mcf_->solve();

        // propagate base solution to single node cut factors
    }

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::Tighten(const std::size_t nr_constraints_to_add)
    {
    }

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::mcf_node_to_graph_node(std::size_t i) const
    {
        assert(i < mcf_->no_nodes());
        if(i == mcf_source_node())
            return base_graph_source_node();
        if(i == mcf_terminal_node())
            return base_graph_terminal_node();
        else
            return i/2; 
    }
            

    template <class FACTOR_MESSAGE_CONNECTION, class SINGLE_NODE_CUT_FACTOR, class SINGLE_NODE_CUT_LIFTED_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, SINGLE_NODE_CUT_FACTOR, SINGLE_NODE_CUT_LIFTED_MESSAGE>::prepare_mcf_costs()
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
                    if (j != base_graph_source_node())
                        mcf_->update_cost(e, -m);
                    else
                        mcf_->update_cost(e, m);
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


}
