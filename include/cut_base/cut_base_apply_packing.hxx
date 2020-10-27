#pragma once

#include <array>
#include <cassert>
#include <algorithm>
#include <variant>
#include <optional>
#include <iostream>
#include <unordered_map>
#include "cut_base_packing.h"
#include "cut_base_instance.hxx"

namespace LPMP {

    namespace {
        using edge_type = std::array<std::size_t,2>;
        using triplet_type = std::array<std::size_t,3>; 
        using quadruplet_type = std::array<std::size_t,4>; 
        using quintuplet_type = std::array<std::size_t,5>; 
    }

    template<typename ITERATOR, typename GET_EDGE_FUNC, typename GET_TRIPLET_FUNC, typename GET_MSG_FUNC>
        void triangulate_cycle(ITERATOR cycle_begin, ITERATOR cycle_end, GET_EDGE_FUNC get_edge_func, GET_TRIPLET_FUNC get_triplet_func, GET_MSG_FUNC get_msg_func, const double cycle_weight, const bool pack_cycle = true)
        {
            assert(std::distance(cycle_begin, cycle_end) >= 3);
            assert(cycle_weight > 0.0); // TODO: current version works without cycle weight. lower bound increase should be at least as big as cycle weight, though.
            // TODO: possibly get normal form for cycle by rotating
            const std::size_t first_node = *cycle_begin;

            enum class wheel_edge_type {wheel, incoming_chord, outgoing_chord};
            using msg_variant = decltype(get_msg_func(edge_type{},triplet_type{}));
            using msg_type = std::pair<msg_variant, wheel_edge_type>;

            auto get_msg_type_func = [&](const auto cycle_iterator, const edge_type e_nodes, const triplet_type t_nodes) -> msg_type {
                auto msg = get_msg_func(e_nodes, t_nodes);
                if(cycle_begin+1 != cycle_iterator && e_nodes[0] == first_node && e_nodes[1] == *cycle_iterator) {
                    return {msg, wheel_edge_type::incoming_chord};
                } else if(cycle_iterator+1 != cycle_end && e_nodes[0] == first_node && e_nodes[1] == *(cycle_iterator+1)) {
                    return {msg, wheel_edge_type::outgoing_chord}; 
                } else {
                    return {msg, wheel_edge_type::wheel}; 
                }
            };

            // TODO: make msg types in lambda const auto&
            auto send_message_to_triplet = [&](auto& edge, auto& triplet, auto& msg) {
                typename std::remove_reference_t<decltype(msg)>::msg_val_type msg_val;
                std::fill(msg_val.begin(), msg_val.end(), 0.0);
                msg.send_message_to_right(edge, msg_val, 1.0);
                msg.RepamLeft(edge, msg_val);
                std::transform(msg_val.begin(), msg_val.end(), msg_val.begin(), [](const auto& x) { return -x; });
                msg.RepamRight(triplet, msg_val);
            };

            auto send_message_to_edge = [&](auto& edge, auto& triplet, auto& msg) {
                typename std::remove_reference_t<decltype(msg)>::msg_val_type msg_val;
                std::fill(msg_val.begin(), msg_val.end(), 0.0);
                msg.send_message_to_left(triplet, msg_val, 1.0);
                msg.RepamRight(triplet, msg_val);
                std::transform(msg_val.begin(), msg_val.end(), msg_val.begin(), [](const auto& x) { return -x; });
                msg.RepamLeft(edge, msg_val);
            };

            auto send_msg_forward = [&](const wheel_edge_type et, auto& edge, auto& triplet, auto& msg) {
                if(et == wheel_edge_type::wheel || et == wheel_edge_type:: outgoing_chord) {
                    std::visit([&](auto& msg) {
                            send_message_to_triplet(edge, triplet, msg);
                    }, msg);
                } else {
                    assert(et == wheel_edge_type:: incoming_chord);
                    std::visit([&](auto& msg) {
                            send_message_to_edge(edge, triplet, msg);
                    }, msg);
                } 
            }; 

            // triangulate cycle: first put all reparametrization into triplet factors and make edges zero
            for(auto it=cycle_begin+1; it+1!=cycle_end; ++it) {
                if(first_node != *it && first_node != *(it+1)) {
                    triplet_type nodes({first_node, *it, *(it+1)});
                    std::sort(nodes.begin(), nodes.end());

                    auto& t = get_triplet_func(nodes);
                    using edge_type = decltype(get_edge_func(std::array<std::size_t,2>{}));
                    edge_type edge_01 = get_edge_func({nodes[0], nodes[1]});
                    edge_type edge_02 = get_edge_func({nodes[0], nodes[2]});
                    edge_type edge_12 = get_edge_func({nodes[1], nodes[2]});
                    auto [edge_01_msg, edge_01_type] = get_msg_type_func(it, {nodes[0], nodes[1]}, nodes);
                    auto [edge_02_msg, edge_02_type] = get_msg_type_func(it, {nodes[0], nodes[2]}, nodes);
                    auto [edge_12_msg, edge_12_type] = get_msg_type_func(it, {nodes[1], nodes[2]}, nodes);
                    if(pack_cycle) {
                        send_msg_forward(edge_01_type, edge_01, t, edge_01_msg);
                        send_msg_forward(edge_02_type, edge_02, t, edge_02_msg);
                        send_msg_forward(edge_12_type, edge_12, t, edge_12_msg);
                    }
                } else {
                    assert(false);
                }
            }

            auto send_msg_backward = [&](const wheel_edge_type et, auto& edge, auto& triplet, auto& msg) {
                if(et == wheel_edge_type::wheel || et == wheel_edge_type:: outgoing_chord) {
                    std::visit([&](auto& msg) {
                            send_message_to_edge(edge, triplet, msg);
                    }, msg);
                } else {
                    assert(et == wheel_edge_type:: incoming_chord);
                    std::visit([&](auto& msg) {
                            send_message_to_triplet(edge, triplet, msg);
                    }, msg);
                } 
            }; 

            // go over triplets backwards now and put as much weight back on the edges
            if(pack_cycle) {
                for(auto it=cycle_end-2; it!=cycle_begin; --it) {
                    triplet_type nodes({first_node, *it, *(it+1)});
                    std::sort(nodes.begin(), nodes.end());

                    auto& t = get_triplet_func(nodes);
                    using edge_type = decltype(get_edge_func(std::array<std::size_t,2>{}));
                    edge_type edge_01 = get_edge_func({nodes[0], nodes[1]});
                    edge_type edge_02 = get_edge_func({nodes[0], nodes[2]});
                    edge_type edge_12 = get_edge_func({nodes[1], nodes[2]});
                    auto [edge_01_msg, edge_01_type] = get_msg_type_func(it, {nodes[0], nodes[1]}, nodes);
                    auto [edge_02_msg, edge_02_type] = get_msg_type_func(it, {nodes[0], nodes[2]}, nodes);
                    auto [edge_12_msg, edge_12_type] = get_msg_type_func(it, {nodes[1], nodes[2]}, nodes);
                    send_msg_backward(edge_01_type, edge_01, t, edge_01_msg);
                    send_msg_backward(edge_02_type, edge_02, t, edge_02_msg);
                    send_msg_backward(edge_12_type, edge_12, t, edge_12_msg);

                }
            }
        }

    // template methods for packing multicut/max-cut instances/Lagrange decompositions
    template<typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2, typename TRIPLET_INSTANCE>
        TRIPLET_INSTANCE pack_cut_base_instance(const cut_base_instance& input, const cycle_packing& cp)
        {
            std::cout << "pack " << cp.no_cycles() << " cycles\n";
            using edge_factor = typename TRIPLET_INSTANCE::edge_factor;
            using triplet_factor = typename TRIPLET_INSTANCE::triplet_factor;
            std::unordered_map<triplet_type, triplet_factor> triplets;
            std::unordered_map<edge_type, edge_factor> edges;

            for(const auto& e : input.edges())
                edges.insert(std::make_pair(edge_type{e[0],e[1]}, e.cost));

            for(std::size_t c=0; c<cp.no_cycles(); ++c) {
                auto [cycle_begin, cycle_end] = cp.get_cycle(c);
                const double cycle_weight = cp.get_cycle_weight(c);

                auto get_edge_func = [&](const edge_type nodes) -> edge_factor& {
                    assert(nodes[0] < nodes[1] && nodes[1] < input.edges().size());
                    if(edges.count(nodes) == 0)
                        edges.insert(std::make_pair(nodes, edge_factor{}));
                    return edges.find(nodes)->second;
                };

                auto get_triplet_func = [&](const triplet_type nodes) -> triplet_factor& {
                    assert(nodes[0] < nodes[1] && nodes[1] < nodes[2] && nodes[2] < input.edges().size());
                    if(triplets.count(nodes) == 0)
                        triplets.insert(std::make_pair(nodes, triplet_factor{}));
                    return triplets.find(nodes)->second; 
                };

                using msg_type = std::variant<UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;
                auto get_msg_func = [&](const edge_type e, const triplet_type t) -> msg_type {
                    if(t[0] == e[0] && t[1] == e[1]) {
                        return UNARY_TRIPLET_MESSAGE_0();
                    } else if(t[0] == e[0] && t[2] == e[1]) {
                        return UNARY_TRIPLET_MESSAGE_1();
                    } else {
                        assert(t[1] == e[0] && t[2] == e[1]);
                        return UNARY_TRIPLET_MESSAGE_2();
                    } 
                };

                triangulate_cycle(cycle_begin, cycle_end, get_edge_func, get_triplet_func, get_msg_func, cycle_weight);
            }

            // add all edges and triplets
            TRIPLET_INSTANCE output;
            for(const auto& e : edges)
                output.add_edge(e.first[0], e.first[1], e.second);
            for(const auto& t : triplets)
                output.add_triplet({t.first[0], t.first[1], t.first[2]}, t.second);

            assert(output.normalized());

            std::cout << "triplet multicut instance has " << output.edges().size() << " edges\n";
            std::cout << "triplet multicut instance has " << output.triplets().size() << " triplets\n";
            std::cout << "triplet multicut instance lower bound = " << output.lower_bound() << "\n";
            return output;
        }

    template<typename TRIPLET_INSTANCE, typename QUADRUPLET_INSTANCE, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        QUADRUPLET_INSTANCE pack_triplet_cut_base_instance(const TRIPLET_INSTANCE& input, const odd_wheel_packing& owp)
        {
            std::cout << "pack " << owp.no_odd_wheels() << " odd wheels\n";
            using triplet_factor = typename TRIPLET_INSTANCE::triplet_factor;
            using quadruplet_factor = typename QUADRUPLET_INSTANCE::quadruplet_factor;
            std::unordered_map<triplet_type, triplet_factor> triplets;
            std::unordered_map<quadruplet_type, quadruplet_factor> quadruplets;

            for(const auto& t : input.triplets())
                triplets.insert(std::make_pair(triplet_type{t[0], t[1], t[2]}, t.cost));

            for(std::size_t c=0; c<owp.no_odd_wheels(); ++c) {
                auto [cycle_begin, cycle_end] = owp.get_cycle(c);
                const std::size_t center_node = owp.get_center_node(c);
                const double odd_wheel_weight = owp.get_odd_wheel_weight(c);

                auto get_triplet_func = [&](const edge_type edge_nodes) -> triplet_factor& {
                    triplet_type nodes{edge_nodes[0], edge_nodes[1], center_node};
                    std::sort(nodes.begin(), nodes.end());
                    if(triplets.count(nodes) == 0)
                        triplets.insert(std::make_pair(nodes, triplet_factor{}));
                    return triplets.find(nodes)->second; 
                };

                auto get_quadruplet_func = [&](const triplet_type triplet_nodes) -> quadruplet_factor& {
                    quadruplet_type nodes{triplet_nodes[0], triplet_nodes[1], triplet_nodes[2], center_node};
                    std::sort(nodes.begin(), nodes.end());
                    if(quadruplets.count(nodes) == 0)
                        quadruplets.insert(std::make_pair(nodes, quadruplet_factor{}));
                    return quadruplets.find(nodes)->second;
                };

                using msg_type = std::variant<TRIPLET_QUADRUPLET_MESSAGE_012,TRIPLET_QUADRUPLET_MESSAGE_013,TRIPLET_QUADRUPLET_MESSAGE_023,TRIPLET_QUADRUPLET_MESSAGE_123>;
                auto get_msg_func = [&](const edge_type e, const triplet_type t) -> msg_type {
                    triplet_type ec{e[0], e[1], center_node};
                    std::sort(ec.begin(), ec.end());

                    quadruplet_type tc{t[0], t[1], t[2], center_node};
                    std::sort(tc.begin(), tc.end());

                    if(ec[0] == tc[0] && ec[1] == tc[1] && ec[2] == tc[2])
                        return TRIPLET_QUADRUPLET_MESSAGE_012();
                    if(ec[0] == tc[0] && ec[1] == tc[1] && ec[2] == tc[3])
                        return TRIPLET_QUADRUPLET_MESSAGE_013();
                    if(ec[0] == tc[0] && ec[1] == tc[2] && ec[2] == tc[3])
                        return TRIPLET_QUADRUPLET_MESSAGE_023();
                    assert(ec[0] == tc[1] && ec[1] == tc[2] && ec[2] == tc[3]);
                    return TRIPLET_QUADRUPLET_MESSAGE_123();
                };

                triangulate_cycle(cycle_begin, cycle_end, get_triplet_func, get_quadruplet_func, get_msg_func, odd_wheel_weight);
            }

            // add all edges, triplets and quadruplets
            QUADRUPLET_INSTANCE output;
            for(const auto& e : input.edges())
                output.add_edge(e[0], e[1], e.cost);
            for(const auto& t : triplets)
                output.add_triplet({t.first[0], t.first[1], t.first[2]}, t.second);
            for(const auto& q : quadruplets)
                output.add_quadruplet({q.first[0], q.first[1], q.first[2], q.first[3]}, q.second);

            assert(output.normalized());

            std::cout << "quadruplet multicut instance has " << output.edges().size() << " edges\n";

            return output;
        }

    template<typename QUADRUPLET_INSTANCE, typename QUINTUPLET_INSTANCE, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234>
        QUINTUPLET_INSTANCE pack_quadruplet_cut_base_instance(const QUADRUPLET_INSTANCE& input, const odd_bicycle_wheel_packing& obwp)
        {
            assert(false);

            QUINTUPLET_INSTANCE output;
            return output;
        }

} // namespace LPMP
