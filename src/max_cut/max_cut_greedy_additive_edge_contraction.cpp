#include <vector>
#include <queue>
#include <cassert>
#include <functional>
#include <algorithm>
#include "max_cut/max_cut_greedy_additive_edge_contraction.h"
#include "union_find.hxx"
#include "dynamic_graph.hxx"

namespace LPMP {

    max_cut_edge_labeling greedy_additive_edge_contraction(const max_cut_instance& instance)
    {
        std::cout << "graph connected: " << instance.graph_connected() << "\n";
        assert(instance.graph_connected());

        struct edge_type {
            double cost;
            std::size_t stamp;
        };

        dynamic_graph<edge_type> g(instance.edges().begin(), instance.edges().end(), [](const auto& e) -> edge_type { return {e.cost, 0}; });
        union_find partition(instance.no_nodes());

        struct edge_type_q : public std::array<std::size_t,2> {
            double cost;
            std::size_t stamp;
        };

        auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost < e2.cost; };
        std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)> Q(pq_cmp);

        // vector stores elements to be added later.
        // If we first remove a node and then add edges, we will reuse the space of the deleted edges.
        // This gives a slight, but consistent performance improvement.
        std::vector<std::pair<std::array<std::size_t,2>, edge_type>> insert_candidates; 

        for(const auto& e : instance.edges())
            if(e.cost >= 0.0)
                Q.push(edge_type_q{e[0], e[1], e.cost, 0});

        while(!Q.empty()) {
            const edge_type_q e_q = Q.top();
            Q.pop();
            const std::size_t i = e_q[0];
            const std::size_t j = e_q[1];
            assert(i != j);

            if(!g.edge_present(i,j))
                continue;
            const auto& e = g.edge(i,j);
            if(e_q.stamp < e.stamp)
                continue;
            if(partition.count() == 2)
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

                    Q.push(edge_type_q{stable_node, head, pp.cost, pp.stamp});
                } else {
                    Q.push(edge_type_q{stable_node, head, p.cost, 0});
                    insert_candidates.push_back({{stable_node, head}, {p.cost, 0}});
                } 
            }
            g.remove_node(merge_node);
            for(const auto& e : insert_candidates)
                g.insert_edge(e.first[0], e.first[1], e.second);
            insert_candidates.clear();
        }

        max_cut_edge_labeling sol(instance, partition);

        // if solution is worse than trivial one, return trivial solution
        if(instance.evaluate(sol) > 0.0)
            std::transform(sol.begin(), sol.end(), sol.begin(), [](const auto x) { return 0; });

        return sol;
    }

} // namespace LPMP 

