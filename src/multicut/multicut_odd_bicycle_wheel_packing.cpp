#include "multicut/multicut_odd_bicycle_wheel_packing.h"
#include "bipartite_graph_helper.hxx"
#include "two_dimensional_variable_array.hxx"
#include <array>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cassert>

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    bool nodes_valid(const std::array<std::size_t,2> axle_nodes, const std::array<std::size_t,2> wheel_nodes)
    {
        if(axle_nodes[0] > axle_nodes[1] || wheel_nodes[0] > wheel_nodes[1]) return false;
        if(axle_nodes[0] == wheel_nodes[0] || axle_nodes[0] == wheel_nodes[1]) return false;
        if(axle_nodes[1] == wheel_nodes[0] || axle_nodes[1] == wheel_nodes[1]) return false;
        return true; 
    }

    // axle edge index, wheel edge index
    std::array<std::size_t,2> get_axle_wheel_edge_indices(const std::array<std::size_t,2> axle_nodes, const std::array<std::size_t,2> wheel_nodes)
    {
        assert(nodes_valid(axle_nodes, wheel_nodes));
        // quadruplet edge order: 01,02,12,03,13,23
        std::array<std::size_t,4> nodes;
        nodes = {axle_nodes[0], axle_nodes[1], wheel_nodes[0], wheel_nodes[1]};
        if(std::is_sorted(nodes.begin(), nodes.end())) 
            return {0,5};
        nodes = {axle_nodes[0], wheel_nodes[0], axle_nodes[1], wheel_nodes[1]};
        if(std::is_sorted(nodes.begin(), nodes.end())) 
            return {1,4};
        nodes = {axle_nodes[0], wheel_nodes[0], wheel_nodes[1], axle_nodes[1]};
        if(std::is_sorted(nodes.begin(), nodes.end())) 
            return {3,2};
        nodes = {wheel_nodes[0], axle_nodes[0], axle_nodes[1], wheel_nodes[1]};
        if(std::is_sorted(nodes.begin(), nodes.end())) 
            return {2,3};
        nodes = {wheel_nodes[0], axle_nodes[0], wheel_nodes[1], axle_nodes[1]};
        if(std::is_sorted(nodes.begin(), nodes.end())) 
            return {4,1};
        nodes = {wheel_nodes[0], wheel_nodes[1], axle_nodes[0], axle_nodes[1]};
        assert(std::is_sorted(nodes.begin(), nodes.end()));
        return {5,0};
    }

    // compute difference between minimal costs such that exactly one edge incident to center node is cut against cost when when zero or two incident to it are cut
    double compute_quadruplet_th(const std::array<std::size_t,2> axle_nodes, const std::array<std::size_t,2> wheel_nodes, const multicut_quadruplet_factor& t)
    {
        assert(nodes_valid(axle_nodes, wheel_nodes));
        double participating_cases = std::numeric_limits<double>::infinity();
        double min_other_cases = 0.0;
        const auto [axle_edge_index, wheel_edge_index] = get_axle_wheel_edge_indices(axle_nodes, wheel_nodes);
        auto update_costs = [&](const std::bitset<6> labeling, const double cost) {
            if(labeling[axle_edge_index] && labeling[wheel_edge_index])
                participating_cases = std::min(cost, participating_cases);
            else
                min_other_cases = std::min(cost, min_other_cases); 
        };
        t.for_each_labeling(update_costs);
        return participating_cases - min_other_cases;
    }

    void reparametrize_quadruplet(multicut_quadruplet_factor& t, const std::array<std::size_t,2> axle_nodes, const std::array<std::size_t,2> wheel_nodes, const double weight)
    {
        assert(weight >= 0.0);
        const auto [axle_edge_index, wheel_edge_index] = get_axle_wheel_edge_indices(axle_nodes, wheel_nodes);
        auto update_costs = [&](const std::bitset<6> labeling, double& cost) {
            if(labeling[axle_edge_index] && labeling[wheel_edge_index])
                cost -= weight;
        };
        t.for_each_labeling(update_costs);
    }

    odd_bicycle_wheel_packing compute_multicut_odd_bicycle_wheel_packing_impl(const quadruplet_multicut_instance& input, const bool record_odd_bicycle_wheels)
    {
        odd_bicycle_wheel_packing obwp;

        // iterate over all edges of the instance which will be axles of odd bicycle wheels.
        struct quadruple_item {
            std::array<std::size_t,2> other_nodes;
            double weight;
            multicut_quadruplet_factor* f;
        };
        std::vector<std::size_t> no_quadruplets_per_axle;
        std::vector<std::array<std::size_t,2>> quadruplet_axles;
        std::unordered_map<std::array<std::size_t,2>, std::size_t> nodes_to_edge;
        for(const auto& q : input.quadruplets()) {
            auto add_edge = [&](const std::size_t i, const std::size_t j) {
                if(nodes_to_edge.count({i,j}) == 0) {
                    quadruplet_axles.push_back({i,j});
                    no_quadruplets_per_axle.push_back(0);
                    nodes_to_edge.insert( std::make_pair(std::array<std::size_t,2>{i,j}, nodes_to_edge.size()) );
                }
                const std::size_t edge_no = nodes_to_edge.find({i,j})->second;
                ++no_quadruplets_per_axle[edge_no];
            };

            add_edge(q[0], q[1]);
            add_edge(q[0], q[2]);
            add_edge(q[0], q[3]);
            add_edge(q[1], q[2]);
            add_edge(q[1], q[3]);
            add_edge(q[2], q[3]);
        }

        two_dim_variable_array<quadruple_item> quadruplets_per_axle(no_quadruplets_per_axle.begin(), no_quadruplets_per_axle.end());

        std::vector<multicut_quadruplet_factor> quadruplets;
        quadruplets.reserve(input.quadruplets().size());
        std::fill(no_quadruplets_per_axle.begin(), no_quadruplets_per_axle.end(), 0);
        for(const auto& q : input.quadruplets()) {
            quadruplets.push_back(q.cost);

            auto add_quadruplet_factor_pointer = [&](const std::array<std::size_t,2> axle_nodes, const std::array<std::size_t,2> other_nodes, multicut_quadruplet_factor& f) {
                assert(nodes_valid(axle_nodes, other_nodes));
                assert(nodes_to_edge.count(axle_nodes) > 0);
                const std::size_t e = nodes_to_edge[axle_nodes];
                quadruplets_per_axle[e][no_quadruplets_per_axle[e]].other_nodes = other_nodes;
                quadruplets_per_axle[e][no_quadruplets_per_axle[e]++].f = &f;
            };
            add_quadruplet_factor_pointer({q[0], q[1]}, {q[2], q[3]}, quadruplets.back());
            add_quadruplet_factor_pointer({q[0], q[2]}, {q[1], q[3]}, quadruplets.back());
            add_quadruplet_factor_pointer({q[0], q[3]}, {q[1], q[2]}, quadruplets.back());
            add_quadruplet_factor_pointer({q[1], q[2]}, {q[0], q[3]}, quadruplets.back());
            add_quadruplet_factor_pointer({q[1], q[3]}, {q[0], q[2]}, quadruplets.back());
            add_quadruplet_factor_pointer({q[2], q[3]}, {q[0], q[1]}, quadruplets.back());
        }

        double lower_bound = input.lower_bound();

        compressed_bipartite_graph_helper<quadruple_item> bfs_helper(input.no_nodes());
        assert(quadruplets_per_axle.size() == quadruplet_axles.size());
        for(std::size_t axle_no=0; axle_no<quadruplets_per_axle.size(); ++axle_no) {
            const std::array<std::size_t,2> axle_nodes = quadruplet_axles[axle_no];

            bfs_helper.construct_compressed_bipartite_graph(quadruplets_per_axle[axle_no].begin(), quadruplets_per_axle[axle_no].end(), 
                    [](const quadruple_item& q) {
                    return q.other_nodes;
                    },
                    [&](const quadruple_item& q) {
                    const double cost = compute_quadruplet_th(axle_nodes, q.other_nodes, *q.f);
                    return quadruple_item{q.other_nodes, cost, q.f}; 
                    } );

            for(std::size_t cycle_length=2; cycle_length<=bfs_helper.get_graph().no_nodes(); ++cycle_length) {
                for(std::size_t ci=0; ci<bfs_helper.no_compressed_nodes(); ++ci) {

                    double cycle_cap = std::numeric_limits<double>::infinity();
                    auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const quadruple_item& e) { cycle_cap = std::min(cycle_cap, e.weight); };

                    auto mask_small_edges = [cycle_length](const std::size_t i, const std::size_t j, const quadruple_item& e, const std::size_t distance) { 
                        if(e.weight <= tolerance) return false;
                        if(distance >= cycle_length) return false;
                        return true;
                    };

                    auto cycle = bfs_helper.get_bfs().find_path(ci, ci+bfs_helper.no_compressed_nodes(), mask_small_edges, cycle_capacity);
                    if(cycle.size() == 0) 
                        continue;
                    assert(cycle_cap >= -1e8);

                    std::transform(cycle.begin(), cycle.end(), cycle.begin(), [&](const std::size_t val) { return val % bfs_helper.no_compressed_nodes(); });
                    cycle = find_subcycle(cycle);

                    lower_bound += cycle_cap;

                    // substract weight from quadruplets and update underlying quadruplets
                    for(std::size_t c=1; c<cycle.size(); ++c) {
                        const std::size_t ci = cycle[c-1];
                        const std::size_t cj = cycle[c];
                        bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes()).weight -= cycle_cap;
                        bfs_helper.get_graph().edge(cj+bfs_helper.no_compressed_nodes(), ci).weight -= cycle_cap;
                        bfs_helper.get_graph().edge(cj, ci+bfs_helper.no_compressed_nodes()).weight -= cycle_cap;
                        bfs_helper.get_graph().edge(ci+bfs_helper.no_compressed_nodes(), cj).weight -= cycle_cap;
                        
                        assert(bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes()).weight >= 0.0);
                        assert(bfs_helper.get_graph().edge(cj+bfs_helper.no_compressed_nodes(), ci).weight >= 0.0);
                        assert(bfs_helper.get_graph().edge(cj, ci+bfs_helper.no_compressed_nodes()).weight >= 0.0);
                        assert(bfs_helper.get_graph().edge(ci+bfs_helper.no_compressed_nodes(), cj).weight >= 0.0);

                        const std::size_t w1 = bfs_helper.compressed_to_original_node(ci);
                        const std::size_t w2 = bfs_helper.compressed_to_original_node(cj);
                        auto& e = bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes());
                        reparametrize_quadruplet( *e.f, axle_nodes, {std::min(w1,w2), std::max(w1,w2)}, cycle_cap);
                    }
                    cycle.resize(cycle.size()-1); 
                    bfs_helper.compressed_path_to_original(cycle);
                    if(record_odd_bicycle_wheels)
                        obwp.add_odd_bicycle_wheel(axle_nodes, cycle.begin(), cycle.end(), cycle_cap);
                }
            }
        }
        return obwp;
    }

void multicut_odd_bicycle_wheel_packing(const quadruplet_multicut_instance& input)
{
    compute_multicut_odd_bicycle_wheel_packing_impl(input, false);
}
odd_bicycle_wheel_packing compute_multicut_odd_bicycle_wheel_packing(const quadruplet_multicut_instance& input)
{
    return compute_multicut_odd_bicycle_wheel_packing_impl(input, true);
}

}
