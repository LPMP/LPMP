#pragma once

#include "LP.h"
#include "solver.hxx"
#include "lifted_disjoint_paths/ldp_instance.hxx"
#include "MCF-SSP/mcf_ssp.hxx"

#include <memory>

namespace LPMP {

    template<class FACTOR_MESSAGE_CONNECTION, class MCF_FACTOR, class SINGLE_NODE_CUT_FACTOR, class MCF_SINGLE_NODE_CUT_MESSAGE>
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
            std::size_t base_graph_source_node() const { return nr_nodes() + 1; }
            std::size_t base_graph_terminal_node() const { return nr_nodes() + 2; }

            LP<FMC> *lp_;
            using mcf_solver_type = MCF::SSP<long, double>;
            std::unique_ptr<mcf_solver_type> mcf_; // minimum cost flow factor for base edges
            std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>> single_node_cut_factors_;
    };


    template <class FACTOR_MESSAGE_CONNECTION, class MCF_FACTOR, class SINGLE_NODE_CUT_FACTOR, class MCF_SINGLE_NODE_CUT_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, MCF_FACTOR, SINGLE_NODE_CUT_FACTOR, MCF_SINGLE_NODE_CUT_MESSAGE>::construct(const lifted_disjoint_paths::LdpInstance &instance)
    {
        // first construct minimum cost flow factor for base edges
        const std::size_t nr_base_graph_nodes = instance.getGraph().numberOfVertices() - 2;
        const std::size_t base_graph_source = instance.getSourceNode();
        const std::size_t base_graph_terminal = instance.getTerminalNode();

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
            auto* outgoing_snc = lp_->template add_factor<SINGLE_NODE_CUT_FACTOR>(instance, i, true);
            single_node_cut_factors_.push_back({incoming_snc, outgoing_snc});
        }
    }

    template <class FACTOR_MESSAGE_CONNECTION, class MCF_FACTOR, class SINGLE_NODE_CUT_FACTOR, class MCF_SINGLE_NODE_CUT_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, MCF_FACTOR, SINGLE_NODE_CUT_FACTOR, MCF_SINGLE_NODE_CUT_MESSAGE>::ComputePrimal()
    {
        prepare_mcf_costs();
        mcf_->solve();

        // propagate base solution to single node cut factors
    }

    template <class FACTOR_MESSAGE_CONNECTION, class MCF_FACTOR, class SINGLE_NODE_CUT_FACTOR, class MCF_SINGLE_NODE_CUT_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, MCF_FACTOR, SINGLE_NODE_CUT_FACTOR, MCF_SINGLE_NODE_CUT_MESSAGE>::Tighten(const std::size_t nr_constraints_to_add)
    {
    }

    template <class FACTOR_MESSAGE_CONNECTION, class MCF_FACTOR, class SINGLE_NODE_CUT_FACTOR, class MCF_SINGLE_NODE_CUT_MESSAGE>
    std::size_t lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, MCF_FACTOR, SINGLE_NODE_CUT_FACTOR, MCF_SINGLE_NODE_CUT_MESSAGE>::mcf_node_to_graph_node(std::size_t i) const
    {
        assert(i < mcf_->no_nodes());
        if(i == mcf_source_node())
            return base_graph_source_node();
        if(i == mcf_terminal_node())
            return base_graph_terminal_node();
        else
            return i/2; 
    }
            

    template <class FACTOR_MESSAGE_CONNECTION, class MCF_FACTOR, class SINGLE_NODE_CUT_FACTOR, class MCF_SINGLE_NODE_CUT_MESSAGE>
    void lifted_disjoint_paths_constructor<FACTOR_MESSAGE_CONNECTION, MCF_FACTOR, SINGLE_NODE_CUT_FACTOR, MCF_SINGLE_NODE_CUT_MESSAGE>::prepare_mcf_costs()
    {
       mcf_->reset_costs();

        for(std::size_t i=0; i<nr_nodes(); ++i)
        {
            auto *incoming_snc = single_node_cut_factors_[i][0];
            const auto incoming_min_marg = incoming_snc->get_factor()->adjustCostsAndSendMessages();
            std::cout << "incoming_min_marg size = " << incoming_min_marg.size() << "\n";
            //for (std::size_t l = 0; l < incoming_min_marg.size(); ++l, ++e)
            //{
            //    const std::size_t start_node = mcf_->tail(e);
            //    assert(mcf_->tail(e) == i);
            //    assert(mcf_->cost(e) == 0.0);
            //    assert(incoming_min_marg.count(l) > 0);
            //    mcf_->update_cost(e, incoming_min_marg.find(l)->second);
            //}
            std::vector<std::tuple<std::size_t,double>> incoming_min_margs_sorted;
            for(const auto [node,m] : incoming_min_marg)
                incoming_min_margs_sorted.push_back({node, m});
            std::sort(incoming_min_margs_sorted.begin(), incoming_min_margs_sorted.end(), [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(b); });
            std::size_t e = mcf_->first_outgoing_arc(incoming_mcf_node(i));
            assert(incoming_min_margs_sorted.size() <= mcf_->no_outgoing_arcs(i));
            for (std::size_t l = 0; l < incoming_min_margs_sorted.size(); ++l)
            {
                std::cout << "incoming min marginal edge " << i << "," << std::get<0>(incoming_min_margs_sorted[l]) << "\n";
                while (mcf_node_to_graph_node(mcf_->head(e)) != std::get<0>(incoming_min_margs_sorted[l]))
                {
                    assert(mcf_->tail(e) == i);
                    assert(mcf_->cost(e) == 0.0);
                    ++e;
                }
                const std::size_t start_node = mcf_->tail(e);
                assert(mcf_->tail(e) == i);
                assert(mcf_->cost(e) == 0.0);
            }

            auto *outgoing_snc = single_node_cut_factors_[i][1];
            const auto outgoing_min_marg = outgoing_snc->get_factor()->adjustCostsAndSendMessages();
            std::cout << "outgoing_min_marg size = " << outgoing_min_marg.size() << "\n";
        }
    }


}
