#include "max_cut/max_cut_odd_bicycle_wheel_packing.h"
#include "bipartite_graph_helper.hxx"
#include "two_dimensional_variable_array.hxx"
#include "max_cut/max_cut_factors_messages.h"
#include <array>
#include <vector>
#include <unordered_map>
#include <variant>
#include <limits>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <variant>

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    double compute_wheel_cost(max_cut_quadruplet_factor& q, const std::size_t axle_edge_index)
    {
        assert(axle_edge_index < 6);
        double min_participating_cases = std::numeric_limits<double>::infinity();
        double min_other_cases = 0.0;
        auto update_costs = [&](const std::bitset<6> labeling, const double cost) {
            const std::size_t no_cut_edges = labeling.count();
            assert(no_cut_edges <= 4);
            if(no_cut_edges >= 3 && labeling[axle_edge_index])
                min_participating_cases = std::min(cost, min_participating_cases);
            else
                min_other_cases = std::min(cost, min_other_cases); 
        };
        q.for_each_labeling(update_costs);
        assert(min_participating_cases < std::numeric_limits<double>::max());
        return -(min_participating_cases - min_other_cases);  // TODO: sign?
    }

    // 2 - (x(0,1) + x(0,2) + x(1,2)) = 0
    double compute_wheel_cost(const max_cut_triplet_factor& t)
    {
        double min_participating_cases = std::numeric_limits<double>::infinity();
        double min_other_cases = 0.0;
        auto update_costs = [&](const std::bitset<3> labeling, const double cost) {
            assert(labeling.count() == std::size_t(labeling[0]) + std::size_t(labeling[1]) + std::size_t(labeling[2]));
            const std::size_t no_cut_edges = labeling.count();
            assert(no_cut_edges <= 2);
            if(no_cut_edges == 2)
                min_participating_cases = std::min(cost, min_participating_cases);
            else
                min_other_cases = std::min(cost, min_other_cases); 
        };
        t.for_each_labeling(update_costs);
        return +(min_other_cases - min_participating_cases); // TODO: sign

    }
    double compute_wheel_cost(const max_cut_triplet_factor& t1, const max_cut_triplet_factor& t2)
    {
        return std::max(0.0, compute_wheel_cost(t1)) + std::max(0.0, compute_wheel_cost(t2));
        //return std::min(0.0, compute_wheel_cost(t1)) + std::max(0.0, compute_wheel_cost(t2));
        //return std::min(compute_wheel_cost(t1), compute_wheel_cost(t2));
        //return std::max(compute_wheel_cost(t1), compute_wheel_cost(t2)); // TODO: what is the best way to compute costs from two triplets?
    }

    double compute_wheel_cost(const std::array<std::size_t,2>& axle_nodes, const double axle_edge_cost, const std::array<std::size_t,3>& nodes1, const max_cut_triplet_factor& t1, const std::array<std::size_t,3>& nodes2, const max_cut_triplet_factor& t2)
    {
        // edges 01,02,12,03,13,23
        max_cut_quadruplet_factor q;

        std::array<std::size_t,4> all_nodes;
        all_nodes[0] = nodes1[0];
        all_nodes[1] = nodes1[1];
        all_nodes[2] = nodes1[2];
        if(std::find(nodes1.begin(), nodes1.end(), nodes2[0]) == nodes1.end())
            all_nodes[3] = nodes2[0];
        if(std::find(nodes1.begin(), nodes1.end(), nodes2[1]) == nodes1.end())
            all_nodes[3] = nodes2[1];
        if(std::find(nodes1.begin(), nodes1.end(), nodes2[2]) == nodes1.end())
            all_nodes[3] = nodes2[2];
        std::sort(all_nodes.begin(), all_nodes.end());

        using axle_msg_type = std::variant<max_cut_edge_quadruplet_message_0,max_cut_edge_quadruplet_message_1,max_cut_edge_quadruplet_message_2,max_cut_edge_quadruplet_message_3,max_cut_edge_quadruplet_message_4,max_cut_edge_quadruplet_message_5>;
        auto [axle_msg, axle_edge_index] = [&]() -> std::pair<axle_msg_type,std::size_t> {
            if(axle_nodes == std::array<std::size_t,2>{all_nodes[0], all_nodes[1]})
                return {max_cut_edge_quadruplet_message_0{}, 0};
            if(axle_nodes == std::array<std::size_t,2>{all_nodes[0], all_nodes[2]})
                return {max_cut_edge_quadruplet_message_1{}, 1};
            if(axle_nodes == std::array<std::size_t,2>{all_nodes[1], all_nodes[2]})
                return {max_cut_edge_quadruplet_message_2{}, 2};
            if(axle_nodes == std::array<std::size_t,2>{all_nodes[0], all_nodes[3]})
                return {max_cut_edge_quadruplet_message_3{}, 3};
            if(axle_nodes == std::array<std::size_t,2>{all_nodes[1], all_nodes[3]})
                return {max_cut_edge_quadruplet_message_4{}, 4};
            assert(axle_nodes[0] == all_nodes[2] && axle_nodes[1] == all_nodes[3]);
            return {max_cut_edge_quadruplet_message_5{}, 5};
        }();

        std::visit([&](auto m) {
                array<double,1> msg_val{0.0};
                max_cut_edge_factor edge_factor;
                edge_factor[0] = axle_edge_cost;
                m.send_message_to_right(edge_factor, msg_val, 1.0);
                m.RepamRight(q, -msg_val);
                }, axle_msg);

        using msg_type = std::variant<
            max_cut_triplet_quadruplet_message_012,
            max_cut_triplet_quadruplet_message_013,
            max_cut_triplet_quadruplet_message_023,
            max_cut_triplet_quadruplet_message_123>;
        auto get_triplet_quadruplet_msg = [&](const std::array<std::size_t,3>& triplet_nodes, const std::array<std::size_t,4>& quadruplet_nodes) -> msg_type
        {
            if(triplet_nodes[0] == quadruplet_nodes[0] && triplet_nodes[1] == quadruplet_nodes[1] && triplet_nodes[2] == quadruplet_nodes[2])
                return max_cut_triplet_quadruplet_message_012{};
            if(triplet_nodes[0] == quadruplet_nodes[0] && triplet_nodes[1] == quadruplet_nodes[1] && triplet_nodes[2] == quadruplet_nodes[3])
                return max_cut_triplet_quadruplet_message_013{};
            if(triplet_nodes[0] == quadruplet_nodes[0] && triplet_nodes[1] == quadruplet_nodes[2] && triplet_nodes[2] == quadruplet_nodes[3])
                return max_cut_triplet_quadruplet_message_023{};
            assert(triplet_nodes[0] == quadruplet_nodes[1] && triplet_nodes[1] == quadruplet_nodes[2] && triplet_nodes[2] == quadruplet_nodes[3]);
            return max_cut_triplet_quadruplet_message_123{}; 
        };

        auto m1 = get_triplet_quadruplet_msg(nodes1, all_nodes);
        auto m2 = get_triplet_quadruplet_msg(nodes2, all_nodes);

        std::visit( [&](auto msg) 
        { 
                auto msg_val = typename decltype(msg)::msg_val_type{};
                msg.send_message_to_right(t1, msg_val, 1.0);
                msg.RepamRight(q, -msg_val);
        }, m1);
        std::visit( [&](auto msg) 
        { 
                auto msg_val = typename decltype(msg)::msg_val_type{};
                msg.send_message_to_right(t1, msg_val, 1.0);
                msg.RepamRight(q, -msg_val);
                }, m2);

        return compute_wheel_cost(q, axle_edge_index);
    }

    void reparametrize_triplet(max_cut_triplet_factor& t, const double weight)
    {
        assert(weight >= 0.0);
        const double repam_weight = 0.0;//std::max(0.0, std::min(weight, compute_wheel_cost(t)));
        assert(repam_weight >= 0.0);
        //assert(weight <= compute_wheel_cost(t));

        auto update_costs = [&](const std::bitset<3> labeling, double& cost) {
            const std::size_t no_cut_edges = labeling.count();
            assert(no_cut_edges == std::size_t(labeling[0]) + std::size_t(labeling[1]) + std::size_t(labeling[2]));
            assert(no_cut_edges <= 2);
            if(no_cut_edges == 2)
                cost += repam_weight; // TODO: cases?
        };

        t.for_each_labeling(update_costs);
    }

    void reparametrize_triplet(max_cut_triplet_factor& t0, max_cut_triplet_factor& t1, const double weight)
    {
        throw std::runtime_error("kwas");
        assert(weight >= 0.0);
        const double w0 = compute_wheel_cost(t0);
        const double w1 = compute_wheel_cost(t1);
        assert(weight <= w0 + w1);

        double repam_weight;
        auto update_costs = [&](const std::bitset<3> labeling, double& cost) {
            const std::size_t no_cut_edges = std::size_t(labeling[0]) + std::size_t(labeling[1]) + std::size_t(labeling[2]);
            assert(no_cut_edges <= 2);
            if(no_cut_edges == 2) // TODO: cases?
                cost += repam_weight;
        };

        repam_weight =  weight * w0/(w0+w1);
        t0.for_each_labeling(update_costs);
        assert(compute_wheel_cost(t0) >= 0.0);
        repam_weight =  weight * w1/(w0+w1);
        t1.for_each_labeling(update_costs);
        assert(compute_wheel_cost(t1) >= 0.0);
    }

    odd_bicycle_wheel_packing compute_max_cut_odd_bicycle_wheel_packing_impl(const triplet_max_cut_instance& input, const bool record_odd_bicycle_wheels)
    {
        odd_bicycle_wheel_packing obwp;

        // iterate over all edges of the instance which will be axles of odd bicycle wheels.
        struct triplet_item {
            std::array<std::size_t,2> other_nodes;
            double weight;
            max_cut_triplet_factor* f;
        };
        auto triplet_item_sort = [](const auto& t1, const auto& t2) {
            return std::lexicographical_compare(t1.other_nodes.begin(), t1.other_nodes.end(), t2.other_nodes.begin(), t2.other_nodes.end());
        };

        std::unordered_map<std::array<std::size_t,2>, std::size_t> edge_nodes_to_edge_index;
        std::vector<std::size_t> no_triplets_per_node(input.no_nodes(), 0);

        std::vector<std::size_t> no_triplets_per_edge(input.no_edges(), 0);
        // TODO: below edge_index is used but we then knwo that edge_index already exists. Optimize!
        auto edge_index = [&](const std::size_t i, const std::size_t j) -> std::size_t {
            if(edge_nodes_to_edge_index.count({i,j}) == 0)
                edge_nodes_to_edge_index.insert(std::make_pair(std::array<std::size_t,2>{i,j}, edge_nodes_to_edge_index.size()));
            const std::size_t edge_index = edge_nodes_to_edge_index.find({i,j})->second; 
            return edge_index;
        };

        for(const auto& q : input.triplets()) {
            no_triplets_per_node[q[0]]++;
            no_triplets_per_node[q[1]]++;
            no_triplets_per_node[q[2]]++;

            auto add_triplet_to_edge = [&](const std::size_t i, const std::size_t j) {
                no_triplets_per_edge[edge_index(i,j)]++;
            };

            add_triplet_to_edge(q[0], q[1]);
            add_triplet_to_edge(q[0], q[2]);
            add_triplet_to_edge(q[1], q[2]);
        }

        two_dim_variable_array<triplet_item> triplets_per_node(no_triplets_per_node.begin(), no_triplets_per_node.end());
        std::fill(no_triplets_per_node.begin(), no_triplets_per_node.end(), 0);
        std::vector<max_cut_triplet_factor> triplets;
        triplets.reserve(input.triplets().size());

        using msg_variant = std::variant<max_cut_edge_triplet_message_0, max_cut_edge_triplet_message_1, max_cut_edge_triplet_message_2>;
        struct triplet_edge_item {
            msg_variant msg;
            max_cut_triplet_factor* f;
        };
        two_dim_variable_array<triplet_edge_item> triplets_per_edge(no_triplets_per_edge.begin(), no_triplets_per_edge.end());
        std::fill(no_triplets_per_edge.begin(), no_triplets_per_edge.end(), 0);

        for(const auto& q : input.triplets()) {
            triplets.push_back(q.cost);
            triplets_per_node[q[0]][no_triplets_per_node[q[0]]++] = {q[1], q[2], 0.0, &triplets.back()};
            triplets_per_node[q[1]][no_triplets_per_node[q[1]]++] = {q[0], q[2], 0.0, &triplets.back()};
            triplets_per_node[q[2]][no_triplets_per_node[q[2]]++] = {q[0], q[1], 0.0, &triplets.back()};

            auto add_triplet_to_edge = [&](const std::array<std::size_t,2> edge, const std::size_t other_node, max_cut_triplet_factor* t) {
                assert(t != nullptr);
                const std::size_t edge_idx = edge_index(edge[0], edge[1]);
                const auto msg = [&]() -> msg_variant {
                    assert(other_node != edge[0] && other_node != edge[1]);
                    assert(edge[0] < edge[1]);
                    if(other_node < edge[0])
                        return max_cut_edge_triplet_message_0{};
                    if(other_node < edge[1])
                        return max_cut_edge_triplet_message_1{};
                    return max_cut_edge_triplet_message_2{};
                }();
                triplets_per_edge[edge_idx][no_triplets_per_edge[edge_idx]++] = {msg, t};
            };

            add_triplet_to_edge({q[0], q[1]}, q[2], &triplets.back());
            add_triplet_to_edge({q[0], q[2]}, q[1], &triplets.back());
            add_triplet_to_edge({q[1], q[2]}, q[0], &triplets.back());
        }

        for(std::size_t e=0; e<triplets_per_edge.size(); ++e) {
            assert(no_triplets_per_edge[e] == triplets_per_edge[e].size());
            for(std::size_t idx=0; idx<triplets_per_edge[e].size(); ++idx) {
                assert(triplets_per_edge(e,idx).f != nullptr);
            }
        }
        
        // sort for fast intersection later
        for(std::size_t i=0; i<triplets_per_node.size(); ++i)
            std::sort(triplets_per_node[i].begin(), triplets_per_node[i].end(), triplet_item_sort);

        double lower_bound = input.lower_bound();
        std::cout << "Initial lower bound for odd bicycle wheel packing = " << lower_bound << "\n";

        struct bfs_item {
            std::array<max_cut_triplet_factor*,2> f;
            double weight;
        };
        compressed_bipartite_graph_helper<bfs_item> bfs_helper(input.no_nodes());
        for(const auto& e : input.edges()) {
            // check whether edge would be on. For this, compute over all triplets that share this edge
            // TODO: collect all edge activations and sort
            const std::size_t edge_idx = edge_index(e[0], e[1]); 
            auto marg_func = [](const double val, triplet_edge_item& t) { 
                const double triplet_val = std::visit( [&](auto msg) {
                array<double,1> msg_val{0.0};   
                msg.send_message_to_left(*t.f, msg_val, 1.0);
                return msg_val[0];
                }, t.msg);
                return val + triplet_val;
            };
            //std::cout << e.cost[0] << "," << std::accumulate(triplets_per_edge[edge_idx].begin(), triplets_per_edge[edge_idx].end(), 0.0, marg_func) << "\n";
            double edge_cost = e.cost[0] + std::accumulate(triplets_per_edge[edge_idx].begin(), triplets_per_edge[edge_idx].end(), 0.0, marg_func);

            if(edge_cost >= -tolerance)
                continue;
            //std::cout  << e.cost[0] << "," << std::accumulate(triplets_per_edge[edge_idx].begin(), triplets_per_edge[edge_idx].end(), 0.0, marg_func) << "\n";

            const std::array<std::size_t,2> axle_nodes {e[0], e[1]};

            struct triplet_intersection_type {
                std::array<std::size_t,2> wheel_nodes;
                std::array<std::size_t,2> triplet_item_index;
            };
            std::vector<triplet_intersection_type> triplet_intersection; // TODO: make more efficient by preallocating
            set_intersection_merge(
                    triplets_per_node[axle_nodes[0]].begin(), triplets_per_node[axle_nodes[0]].end(),
                    triplets_per_node[axle_nodes[1]].begin(), triplets_per_node[axle_nodes[1]].end(),
                    std::back_inserter(triplet_intersection), 
                    triplet_item_sort,
                    [&](triplet_item& t1, triplet_item& t2) -> triplet_intersection_type {
                        assert(t1.other_nodes[0] == t1.other_nodes[0]);
                        assert(t1.other_nodes[1] == t1.other_nodes[1]);
                        const std::size_t t1_index = std::distance(&triplets_per_node[axle_nodes[0]][0], &t1);
                        const std::size_t t2_index = std::distance(&triplets_per_node[axle_nodes[1]][0], &t2);
                        assert(t1_index < triplets_per_node[axle_nodes[0]].size());
                        assert(t2_index < triplets_per_node[axle_nodes[1]].size());
                        return {t1.other_nodes, t1_index, t2_index};
                    });

            // filter triplet_intersections by checking whether they have cost >= tolerance
            auto triplet_intersection_end = std::remove_if(triplet_intersection.begin(), triplet_intersection.end(), 
                    [&](const triplet_intersection_type& t) {
                        //const double cost = compute_wheel_cost(triplets[t.triplet_item_index[0]], triplets[t.triplet_item_index[1]]);
                    std::array<std::size_t,3> triplet_nodes1 = {t.wheel_nodes[0], t.wheel_nodes[1], axle_nodes[0]};
                    std::sort(triplet_nodes1.begin(), triplet_nodes1.end());
                    std::array<std::size_t,3> triplet_nodes2 = {t.wheel_nodes[0], t.wheel_nodes[1], axle_nodes[1]};
                    std::sort(triplet_nodes2.begin(), triplet_nodes2.end());
                    const double cost = compute_wheel_cost(axle_nodes, edge_cost, triplet_nodes1, triplets[t.triplet_item_index[0]], triplet_nodes2, triplets[t.triplet_item_index[1]]);
                        return cost < tolerance;
                        }
                    );
            triplet_intersection.resize(std::distance(triplet_intersection.begin(), triplet_intersection_end));

            bfs_helper.construct_compressed_bipartite_graph(triplet_intersection.begin(), triplet_intersection.end(), 
                    [](const triplet_intersection_type& t) {
                        return t.wheel_nodes;
                    },
                    [&](const triplet_intersection_type& t) {
                    std::array<std::size_t,3> triplet_nodes1 = {t.wheel_nodes[0], t.wheel_nodes[1], axle_nodes[0]};
                    std::sort(triplet_nodes1.begin(), triplet_nodes1.end());
                    std::array<std::size_t,3> triplet_nodes2 = {t.wheel_nodes[0], t.wheel_nodes[1], axle_nodes[1]};
                    std::sort(triplet_nodes2.begin(), triplet_nodes2.end());
                    const double cost = compute_wheel_cost(axle_nodes, edge_cost, triplet_nodes1, triplets[t.triplet_item_index[0]], triplet_nodes2, triplets[t.triplet_item_index[1]]);
                        //const double cost = compute_wheel_cost(triplets[t.triplet_item_index[0]], triplets[t.triplet_item_index[1]]);
                        return bfs_item{{&triplets[t.triplet_item_index[0]], &triplets[t.triplet_item_index[1]]}, cost}; 
                    } );

            const std::array<std::size_t,10> cycle_lengths{2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
            for(const std::size_t cycle_length : cycle_lengths) {
                for(std::size_t ci=0; ci<bfs_helper.no_compressed_nodes(); ++ci) {
                    if(!bfs_helper.get_union_find().connected(ci, ci+bfs_helper.no_compressed_nodes()))
                        continue;

                    while(edge_cost < -tolerance) {
                        double cycle_cap = -edge_cost;
                        assert(cycle_cap > 0.0);
                        auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const bfs_item& e) { cycle_cap = std::min(cycle_cap, e.weight); };

                        auto mask_small_edges = [cycle_length](const std::size_t i, const std::size_t j, const bfs_item& e, const std::size_t distance) { 
                            if(e.weight <= tolerance) return false;
                            if(distance >= cycle_length) return false;
                            return true;
                        };

                        auto cycle = bfs_helper.get_bfs().find_path(ci, ci+bfs_helper.no_compressed_nodes(), mask_small_edges, cycle_capacity);
                        if(cycle.size() == 0) 
                            break;
                        assert(cycle_cap >= tolerance);

                        std::transform(cycle.begin(), cycle.end(), cycle.begin(), [&](const std::size_t val) { return val % bfs_helper.no_compressed_nodes(); });
                        cycle = find_subcycle(cycle);

                        lower_bound += cycle_cap;
                        edge_cost += cycle_cap;
                        assert(edge_cost <= tolerance);

                        // substract weight from triplets and update underlying triplets
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
                            // TODO: or 1.0?
                            //reparametrize_triplet(*e.f[0], *e.f[1], cycle_cap);
                            //reparametrize_triplet(*e.f[0], cycle_cap);
                            //reparametrize_triplet(*e.f[1], cycle_cap);
                        }
                        cycle.resize(cycle.size()-1); 
                        bfs_helper.compressed_path_to_original(cycle);
                        if(record_odd_bicycle_wheels)
                            obwp.add_odd_bicycle_wheel(axle_nodes, cycle.begin(), cycle.end(), cycle_cap);
                    }
                }
            }
        }

        std::cout << "odd bicycle_wheel_packing found " << obwp.no_odd_bicycle_wheels() << " wheels\n";
        std::cout << "lower bound after after packing odd bicycle wheels = " << lower_bound << "\n";

        return obwp;
    }

void max_cut_odd_bicycle_wheel_packing(const triplet_max_cut_instance& input)
{
    compute_max_cut_odd_bicycle_wheel_packing_impl(input, false);
}
odd_bicycle_wheel_packing compute_max_cut_odd_bicycle_wheel_packing(const triplet_max_cut_instance& input)
{
    return compute_max_cut_odd_bicycle_wheel_packing_impl(input, true);
}

}

