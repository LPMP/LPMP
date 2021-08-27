#include <vector>
#include <queue>
#include <numeric>
#include <cassert>
#include "asymmetric_multiway_cut/asymmetric_multiway_cut_gaec.h"
#include "dynamic_graph.hxx"
#include "union_find.hxx"
#include <iostream>

namespace LPMP {

    asymmetric_multiway_cut_labeling asymmetric_multiway_cut_gaec(const asymmetric_multiway_cut_instance& instance)
    {
        const size_t nr_nodes = instance.nr_nodes();
        const size_t nr_edges = instance.nr_edges();
        const size_t nr_labels = instance.nr_labels();
        constexpr static size_t no_label = std::numeric_limits<size_t>::max();

        std::vector<double> label_costs(nr_nodes*nr_labels);
        for(size_t i=0; i<nr_nodes; ++i)
            for(size_t l=0; l<nr_labels; ++l)
                label_costs[i*nr_labels + l] = instance.node_costs(i,l);

        union_find partition(nr_nodes);

        auto label_cost = [&](const size_t i, const size_t l) -> double& {
            assert(i < nr_nodes && l < nr_labels);
            const size_t i_c = partition.find(i);
            return label_costs[i_c*nr_labels + l]; 
        };

        // first, label each node with its most probable class.
        std::vector<size_t> node_labels(nr_nodes, no_label);
        for(size_t i=0; i<nr_nodes; ++i)
        {
            double min_label_cost = std::numeric_limits<double>::infinity();
            size_t min_label = 0;
            for(size_t l=0; l<nr_labels; ++l)
            {
                if(instance.node_costs(i,l) <= min_label_cost)
                {
                    min_label_cost = instance.node_costs(i,l);
                    min_label = l;
                } 
            } 
            node_labels[i] = min_label;
        }
        
        // second, merge individual nodes and switch labels for optimal cost decrease
        struct edge_type {
            double cost;
            size_t stamp;
        };
        dynamic_graph<edge_type> g(instance.edge_costs.edges().begin(), instance.edge_costs.edges().end(), [](const auto& e) -> edge_type { return {e.cost, 0}; });

        struct edge_type_q : public std::array<size_t,2> {
            double cost;
            size_t label;
            size_t stamp;
        };

        auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost > e2.cost; }; 
        std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)> Q(pq_cmp);

        std::vector<std::pair<std::array<size_t,2>, edge_type>> insert_candidates; // vector stores elements to be added later. if we first remove a node and then add edges, we will reuse the space of the deleted edges. This gives a slight, but real performance improvement.

        auto compute_edge_cost = [&](const double mc_edge_cost, const size_t i, const size_t j) -> std::tuple<double,size_t> {
            // cost of join nodes, assign partition of joint minimum cost
            const size_t i_c = partition.find(i);
            const size_t j_c = partition.find(j);
            double min_label_cost = std::numeric_limits<double>::infinity();
            size_t min_label = 0;
            for(size_t l=0; l<nr_labels; ++l)
            {
                if(label_costs[i_c*nr_labels + l] + label_costs[j_c*nr_labels + l] <= min_label_cost)
                {
                    min_label_cost = label_costs[i_c*nr_labels + l] + label_costs[j_c*nr_labels + l];
                    min_label = l;
                } 
            }
            const double sep_cost = mc_edge_cost + label_costs[i_c*nr_labels + node_labels[i_c]] + label_costs[j_c*nr_labels + node_labels[j_c]];
            const double join_cost = min_label_cost;
            return {join_cost - sep_cost, min_label};
        };

        double min_join_cost = std::numeric_limits<double>::infinity();
        double max_join_cost = -std::numeric_limits<double>::infinity();
        for(const auto& e : instance.edge_costs.edges())
        {
            const auto [join_cost, join_label] = compute_edge_cost(e.cost, e[0], e[1]);
            min_join_cost = std::min(min_join_cost, join_cost);
            max_join_cost = std::max(max_join_cost, join_cost);
            assert(partition.find(e[0]) == e[0]);
            assert(partition.find(e[1]) == e[1]);
            if(join_cost <= 0.0)
                Q.push(edge_type_q{e[0], e[1], join_cost, join_label, 0});
        }

        while(!Q.empty()) {
            const edge_type_q e_q = Q.top();
            Q.pop();
            const size_t i = e_q[0];
            const size_t j = e_q[1];

            if(!g.edge_present(i,j))
                continue;
            const auto& e = g.edge(i,j);
            if(e_q.stamp < e.stamp)
                continue;
            if(e_q.cost >= 0.0)
                break;

            const size_t c_i = partition.find(i);
            const size_t c_j = partition.find(j);
            partition.merge(i,j);
            const size_t c_ij = partition.find(i);
            for(size_t l=0; l<nr_labels; ++l)
                label_costs[c_ij*nr_labels + l] = label_costs[c_i*nr_labels + l] + label_costs[c_j*nr_labels + l]; 
            node_labels[c_ij] = e_q.label;

            const auto [stable_node, merge_node] = [&]() -> std::array<size_t,2> {
                if(g.no_edges(i) < g.no_edges(j))
                    return {j,i};
                else
                    return {i,j};
            }();

            for(size_t edge_index=g.first_outgoing_edge_index(merge_node); edge_index!=decltype(g)::no_next_edge; edge_index=g.next_outgoing_edge_index(edge_index)) {
                const size_t head = g.head(edge_index);
                if(head == stable_node)
                    continue;
                auto& p = g.edge(merge_node,head);
                if(g.edge_present(stable_node, head)) {
                    // update costs and new minimum label
                    auto& pp = g.edge(stable_node, head);
                    pp.cost += p.cost;
                    pp.stamp++;

                    const auto [join_cost, join_label] = compute_edge_cost(pp.cost, stable_node, head);
                    if(join_cost <= 0.0)
                    {
                        Q.push(edge_type_q{stable_node, head, join_cost, join_label, pp.stamp});
                    }
                } else {
                    const auto [join_cost, join_label] = compute_edge_cost(p.cost, stable_node, head);
                    if(join_cost <= 0.0)
                        Q.push(edge_type_q{stable_node, head, join_cost, join_label, 0});
                    insert_candidates.push_back({{stable_node, head}, {p.cost, 0}});
                } 
            }
            g.remove_node(merge_node);
            for(const auto& e : insert_candidates)
                g.insert_edge(e.first[0], e.first[1], e.second);
            insert_candidates.clear();
        }

        // Merge non-partitionable classes:
        for(size_t e=0; e<nr_edges; ++e)
        {
            const size_t i = instance.edge_costs.edges()[e][0];
            const size_t j = instance.edge_costs.edges()[e][1];
            const size_t c_i = partition.find(i);
            const size_t c_j = partition.find(j);
            const size_t l_i = node_labels[c_i];
            const size_t l_j = node_labels[c_j];

            if(l_i == l_j && !instance.node_costs.partitionable(l_i) && !partition.connected(i,j))
                partition.merge(i, j);
        }

        // construct labeling
        asymmetric_multiway_cut_labeling labeling;

        for(size_t e=0; e<nr_edges; ++e)
        {
            const size_t i = instance.edge_costs.edges()[e][0];
            const size_t j = instance.edge_costs.edges()[e][1];
            if(partition.connected(i,j))
                labeling.edge_labels.push_back(0);
            else
                labeling.edge_labels.push_back(1); 
        }

        for(size_t i=0; i<nr_nodes; ++i)
        {
            const size_t c = partition.find(i);
            assert(node_labels[c] < nr_labels);
            labeling.node_labels.push_back(node_labels[c]);
            labeling.node_connected_components_ids.push_back(c);
        }

        assert(instance.feasible(labeling));
        return labeling;
    }
}
