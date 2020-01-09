#pragma once

#include <vector>
#include <array>
#include <random>
#include <unordered_set>

#include "mrf/max_flow_instance.hxx"
#include "mrf/binary_MRF_instance.hxx"
#include "max_cut/max_cut_instance.hxx"

namespace LPMP {

std::vector<std::array<std::size_t,2>> generate_random_graph(const std::size_t no_nodes, const std::size_t no_edges, std::random_device& rd)
{
    std::mt19937 gen{rd()};

    std::vector<std::array<std::size_t,2>> edges;
    // generate edges randomly. Check if they exist already and if yes, skip
    std::uniform_int_distribution<std::size_t> ed(0,no_nodes-1);
    std::size_t no_generated_edges = 0;
    std::unordered_set<std::array<std::size_t,2>> edge_map;
    for(std::size_t c=0; c<no_edges; ++c) {
        auto i = ed(gen);
        auto j = ed(gen);
        if(i == j) { continue; }
        if(i>j) { std::swap(i,j); }
        if(edge_map.count({i,j})) { continue; }

        edge_map.insert({i,j});
        edges.push_back({i,j});
    }   
    return edges; 
}

max_cut_instance generate_random_max_cut_instance(const std::size_t no_nodes, const std::size_t no_edges, std::random_device& rd)
{
    const auto edges = generate_random_graph(no_nodes, no_edges, rd);

    max_cut_instance output;

    std::mt19937 gen{rd()};
    std::normal_distribution ud(0.0,5.0);

    for(auto e : edges)
        output.add_edge(e[0], e[1], ud(gen));

    return output;
}

binary_MRF_instance generate_random_binary_MRF_instance(const std::size_t no_nodes, const std::size_t no_edges, std::random_device& rd)
{
    const auto edges = generate_random_graph(no_nodes, no_edges, rd);

    binary_MRF_instance output;
    output.unaries.reserve(no_nodes);
    output.pairwise_potentials.reserve(edges.size()); 

    std::mt19937 gen{rd()};
    std::normal_distribution ud(0.0,5.0);

    for(std::size_t i=0; i<no_nodes; ++i) {
        output.unaries.push_back({ud(gen), ud(gen)});
    }   

    for(auto e : edges) {
        binary_MRF_instance::binary_pairwise_potential p(e[0], e[1], {{{ud(gen), ud(gen)}, {ud(gen), ud(gen)}}});
        output.pairwise_potentials.push_back(p);
    }

    return output;
}

binary_Potts_instance generate_random_binary_Potts_instance(const std::size_t no_nodes, const std::size_t no_edges, std::random_device& rd)
{
    const auto edges = generate_random_graph(no_nodes, no_edges, rd);

    binary_Potts_instance output;
    output.unaries.reserve(no_nodes);
    output.pairwise_potentials.reserve(edges.size()); 

    std::mt19937 gen{rd()};
    std::normal_distribution ud(0.0,5.0);

    for(std::size_t i=0; i<no_nodes; ++i) {
        output.unaries.push_back({ud(gen), ud(gen)});
    }   

    for(auto e : edges) {
        output.pairwise_potentials.push_back({e[0], e[1], ud(gen)});
    }

    return output;
}

max_flow_instance generate_random_qpbo_instance(const std::size_t no_nodes, const std::size_t no_edges, std::random_device& rd, const bool submodular = false)
{
    max_flow_instance output;
    output.source = 0;
    output.terminal = 1;
    output.no_nodes = 2 + 2*no_nodes;

    std::mt19937 gen{rd()};
    std::bernoulli_distribution bd(0.5);
    std::normal_distribution ud(0.0,5.0);

    // generate edges from source and to terminal
    for(std::size_t i=0; i<no_nodes; ++i) {
        const auto cost = ud(gen);
        if(cost > 0) { // let label 0 have cost
            output.add_arc(output.source, 2+i, 0.5*cost);
            output.add_arc(2+no_nodes+i, output.terminal, 0.5*cost);
        } else if(cost < 0) { // let label 1 have -cost
            output.add_arc(output.source, 2+no_nodes+i, -0.5*cost);
            output.add_arc(2+i, output.terminal, -0.5*cost); 
        }
    }   

    // generate random edges representing pairwise potentials
    const auto edges = generate_random_graph(no_nodes, no_edges, rd);
    for(const auto& e : edges) {
        const auto i = e[0];
        const auto j = e[1];
        assert(i<j);
        const auto cost_1 = std::abs(ud(gen));
        const auto cost_2 = std::abs(ud(gen));
        const bool potential_submodular = submodular ? 1 : bd(gen);
        if(potential_submodular) {
            output.add_arc(2+i, 2+j, cost_1);
            output.add_arc(2+j, 2+i, cost_2);
            output.add_arc(2+no_nodes+i, 2+no_nodes+j, cost_2);
            output.add_arc(2+no_nodes+j, 2+no_nodes+i, cost_1);
        } else {
            output.add_arc(2+i, 2+no_nodes+j, cost_1);
            output.add_arc(2+no_nodes+j, 2+i, cost_2);
            output.add_arc(2+no_nodes+i, 2+j, cost_2);
            output.add_arc(2+j, 2+no_nodes+i, cost_1); 
        }
    }

    return output;
}

max_flow_instance generate_random_graph_cut_instance(const std::size_t no_nodes, const std::size_t no_edges, std::random_device& rd)
{
    max_flow_instance output;
    output.source = 0;
    output.terminal = 1;

    std::mt19937 gen{rd()};
    std::bernoulli_distribution bd(0.5);
    std::normal_distribution ud(0.0,5.0);

    // generate edges from source and to terminal
    for(std::size_t i=0; i<no_nodes; ++i) {
        const auto cap = std::abs(ud(gen));
        if(bd(gen)) {
            output.add_arc(output.source, 2+i, cap);
        } else {
            output.add_arc(2+i, output.terminal, cap);
        }
    }   

    // generate random edges representing pairwise potentials
    const auto edges = generate_random_graph(no_nodes, no_edges, rd);
    for(const auto& e : edges) {
        const auto i = e[0];
        const auto j = e[1];
        output.add_arc(2+i, 2+j, std::abs(ud(gen)));
        output.add_arc(2+j, 2+i, std::abs(ud(gen)));
    } 

    return output;
}

} // namespace LPMP
