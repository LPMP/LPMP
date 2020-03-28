#include "max_cut/max_cut_sahni_gonzalez.h"
#include <algorithm>
#include <functional>
#include <vector>
#include <queue>
#include "graph.hxx"

namespace LPMP {

    struct node_type_q  {
        double cut_val_x, cut_val_y;
        std::size_t node;
        std::size_t stamp;
    };

    template<typename PQ_CMP>
    max_cut_node_labeling max_cut_sahni_gonzalez_impl(const max_cut_instance& instance, PQ_CMP pq_cmp)
    {
        assert(instance.no_nodes() >= 2);

        // first find edge with minimal weight and put endpoints in distinct components
        const auto min_edge = *std::min_element(instance.edges().begin(), instance.edges().end(), [](const auto e1, const auto e2) { return e1.cost[0] < e2.cost[0]; });

        const std::size_t x = min_edge[0];
        const std::size_t y = min_edge[1];
        max_cut_node_labeling partition(instance.no_nodes(), 0);
        partition[x] = 0;
        partition[y] = 1;

        std::vector<std::array<double,2>> cut_values (instance.no_nodes(), {0.0, 0.0});
        for(const auto& e : instance.edges()) {
            if(e[0] == x)
                cut_values[e[1]][0] += e.cost[0];
            if(e[1] == x)
                cut_values[e[0]][0] += e.cost[0];
            if(e[0] == y)
                cut_values[e[1]][1] += e.cost[0];
            if(e[1] == y)
                cut_values[e[0]][1] += e.cost[0];
        }

        // TODO: possibly a static graph will be sufficient, as we only need to mask edges?
        struct edge_type { double cost; bool active = true; };
        graph<edge_type> g(instance.edges().begin(), instance.edges().end(), 
                //[](const auto& e) -> std::array<std::size_t,2> { return {e[0], e[1]}; },
                [](const auto& e) -> edge_type { return {e.cost}; }
                ); 

        std::priority_queue<node_type_q, std::vector<node_type_q>, decltype(pq_cmp)> Q(pq_cmp);
        std::vector<std::size_t> stamps(instance.no_nodes(), 0);

        auto get_cut_value = [&](const std::size_t i, const std::size_t y) {
            if(g.edge_present(i,y))
                return g.edge(i,y).cost;
            else
                return 0.0;
        };

        for(std::size_t i=0; i<instance.no_nodes(); ++i) {
            if(i != x && i != y) { // or only add incident edges to x or y?
                const double cut_val_x = cut_values[i][0];
                const double cut_val_y = cut_values[i][1];
                Q.push({cut_val_x, cut_val_y, i, 0});
            }
        }

        while(!Q.empty()) {
            const node_type_q i_q = Q.top();
            Q.pop();
            const std::size_t i = i_q.node;
            const std::size_t stamp = i_q.stamp;
            assert(stamp < std::numeric_limits<std::size_t>::max());
            if(stamp < stamps[i])
                continue;
            if(i_q.cut_val_x < i_q.cut_val_y)
                partition[i] = 1;
            else
                partition[i] = 0;
            const std::size_t partition_index = i_q.cut_val_x < i_q.cut_val_y ? 1 : 0;
            const std::size_t partition_node = partition_index == 0 ? x : y;

            // update cut scores
            //for(std::size_t edge_index=g.first_outgoing_edge_index(i); edge_index!=decltype(g)::no_next_edge; edge_index=g.next_outgoing_edge_index(edge_index)) {
            for(auto edge_it=g.begin(i); edge_it!=g.end(i); ++edge_it) {
                if(!edge_it->edge().active) continue;
                const std::size_t j = edge_it->head();
                if(j == x || j == y)
                    continue;
                auto& p = edge_it->edge();
                cut_values[j][partition_index] += p.cost;
                Q.push({cut_values[j][0], cut_values[j][1], j, ++stamps[j]});
                p.active = false;
                edge_it->sister().edge().active = false;
            }

            // remove node i
            //g.remove_node(i);
            stamps[i] = std::numeric_limits<std::size_t>::max(); // node i will not be considered anymore
        }

        return partition;
    }

    max_cut_node_labeling max_cut_sahni_gonzalez_1(const max_cut_instance& instance)
    {
        auto pq_cmp = [](const node_type_q& e1, const node_type_q& e2) { return std::min(e1.cut_val_x, e1.cut_val_y) > std::min(e2.cut_val_x, e2.cut_val_y); };
        return max_cut_sahni_gonzalez_impl(instance, pq_cmp); 
    }
    max_cut_node_labeling max_cut_sahni_gonzalez_2(const max_cut_instance& instance)
    {
        auto pq_cmp = [](const node_type_q& e1, const node_type_q& e2) { return std::max(e1.cut_val_x, e1.cut_val_y) > std::max(e2.cut_val_x, e2.cut_val_y); };
        return max_cut_sahni_gonzalez_impl(instance, pq_cmp); 
    }
    max_cut_node_labeling max_cut_sahni_gonzalez_3(const max_cut_instance& instance)
    {
        auto pq_cmp = [](const node_type_q& e1, const node_type_q& e2) { return std::abs(e1.cut_val_x - e1.cut_val_y) < std::abs(e2.cut_val_x - e2.cut_val_y); };
        return max_cut_sahni_gonzalez_impl(instance, pq_cmp); 
    }


} // namespace LPMP
