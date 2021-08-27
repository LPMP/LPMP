#include <vector>
#include <queue>
#include <cassert>
#include <functional>
#include "multicut/multicut_balanced_edge_contraction.h"
#include "union_find.hxx"
#include "dynamic_graph.hxx"

namespace LPMP {
    multicut_edge_labeling multicut_balanced_edge_contraction(const multicut_instance& instance)
    {
        multicut_edge_labeling sol;
        std::vector<int> node_labels;
        std::tie(sol, node_labels) = multicut_balanced_edge_contraction_impl(instance);
        return sol;
    }

    std::tuple<multicut_edge_labeling, std::vector<int>> multicut_balanced_edge_contraction_impl(const multicut_instance& instance)
    {
        struct edge_type {
            double cost;
            std::size_t stamp;
        };

        dynamic_graph<edge_type> g(instance.edges().begin(), instance.edges().end(), [](const auto& e) -> edge_type { return {e.cost, 0}; });
        union_find partition(instance.no_nodes());

        struct edge_type_q : public std::array<std::size_t,2> {
            double cost;
            double balanced_cost;
            std::size_t stamp;
        };

        auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.balanced_cost < e2.balanced_cost; };
        std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)> Q(pq_cmp);

        auto compute_balanced_edge_cost = [&](const double mc_edge_cost, const size_t i, const size_t j) -> double {
            const size_t i_size = partition.no_elements(partition.find(i));
            const size_t j_size = partition.no_elements(partition.find(j));
            return mc_edge_cost / double(i_size + j_size);
        };

        std::vector<std::pair<std::array<std::size_t,2>, edge_type>> insert_candidates; // vector stores elements to be added later. if we first remove a node and then add edges, we will reuse the space of the deleted edges. This gives a slight, but real performance improvement.

        for(const auto& e : instance.edges())
            if(e.cost >= 0.0)
                Q.push(edge_type_q{e[0], e[1], e.cost, compute_balanced_edge_cost(e.cost, e[0], e[1]), 0});

        while(!Q.empty()) {
            const edge_type_q e_q = Q.top();
            Q.pop();
            const std::size_t i = e_q[0];
            const std::size_t j = e_q[1];

            if(!g.edge_present(i,j))
                continue;
            const auto& e = g.edge(i,j);
            if(e_q.stamp < e.stamp)
                continue;
            if(e.cost <= 0.0)
                break;

            partition.merge(i,j);

            const auto [stable_node, merge_node] = [&]() -> std::array<std::size_t,2> {
                if(g.no_edges(i) < g.no_edges(j))
                    return {j,i};
                else
                    return {i,j};
            }();

            for(std::size_t edge_index=g.first_outgoing_edge_index(merge_node); edge_index!=decltype(g)::no_next_edge; edge_index=g.next_outgoing_edge_index(edge_index)) {
                const std::size_t head  = g.head(edge_index);
                if(head == stable_node)
                    continue;
                auto& p = g.edge(merge_node,head);
                if(g.edge_present(stable_node, head)) {
                    auto& pp = g.edge(stable_node, head);
                    pp.cost += p.cost;
                    pp.stamp++;

                    if(pp.cost >= 0.0)
                        Q.push(edge_type_q{stable_node, head, pp.cost, compute_balanced_edge_cost(pp.cost, stable_node, head), pp.stamp});
                } else {
                    if(p.cost >= 0.0)
                        Q.push(edge_type_q{stable_node, head, p.cost, compute_balanced_edge_cost(p.cost, stable_node, head), 0});
                    insert_candidates.push_back({{stable_node, head}, {p.cost, 0}});
                } 
            }
            g.remove_node(merge_node);
            for(const auto& e : insert_candidates)
                g.insert_edge(e.first[0], e.first[1], e.second);
            insert_candidates.clear();
        }

        std::vector<int> node_connected_components_ids(instance.no_nodes());
        for(size_t i=0; i<instance.no_nodes(); ++i)
        {
            const size_t c = partition.find(i);
            node_connected_components_ids[i] = c;
        }
        return {multicut_edge_labeling(instance, partition), node_connected_components_ids}; 
    }

} // namespace LPMP 
