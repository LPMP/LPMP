#pragma once

#include "two_dimensional_variable_array.hxx"
#include "pseudo_peripheral_node.hxx"
#include "permutation.hxx"
#include <limits>
#include <vector>
#include <queue>

namespace LPMP {

    template<typename ADJACENCY_GRAPH>
    permutation bfs_ordering(const ADJACENCY_GRAPH& adj)
    {
        struct node
        {
            node(std::size_t index, std::size_t dist, double avg_adj_dist) : 
                index_(index), dist_(dist), avg_adj_dist_(avg_adj_dist)
            {}

            bool operator <(node const & other) const
            {
                if (dist_ != other.dist_)
                    return dist_ < other.dist_;
                else
                    return avg_adj_dist_ < other.avg_adj_dist_;
            }

            std::size_t index_;
            std::size_t dist_;
            double avg_adj_dist_;
        };

        permutation ordering(adj.size());

        std::queue<std::size_t> Q;
        std::vector<std::size_t> dist(adj.size(), std::numeric_limits<std::size_t>::max());
        std::vector<char> seen(adj.size(), 0);
        std::vector<char> is_terminal(adj.size(), 0);
        std::size_t pos = adj.size()-1;

        std::vector<size_t> pseudo_peripheral_nodes = find_pseudo_peripheral_nodes(adj);

        // loop over connected components
        for (size_t s : pseudo_peripheral_nodes)
        {
            Q.push(s);
            seen[s] = 1;
            dist[s] = 0;
            std::vector<std::size_t> terminals;

            // forward run of BFS from source to determine distances and terminals of search tree
            while (!Q.empty())
            {
                const std::size_t i = Q.front();
                Q.pop();

                bool terminal = true;
                for(const std::size_t j : adj[i]) {
                    assert(j < dist.size() && i < dist.size());
                    if (dist[j] > dist[i] || is_terminal[j])
                        terminal = false;
                    if (seen[j])
                        continue;
                    seen[j] = 1;
                    dist[j] = dist[i] + 1;
                    Q.push(j);
                }

                if(terminal)
                {
                    terminals.push_back(i);
                    is_terminal[i] = 1;
                }
            }

            std::fill(seen.begin(), seen.end(), 0);
            std::priority_queue<node> dist_queue;
            for(const std::size_t t : terminals) {
                double avg_adj_dist = 0;
                for (const auto n : adj[t])
                    avg_adj_dist += dist[n];
                avg_adj_dist /= (double) adj[t].size();
                dist_queue.emplace(t, dist[t], avg_adj_dist);
                seen[t] = 1;
            }

            // backward run of BFS from all terminals, ordering vertices in decreasing distance from source
            while(!dist_queue.empty()) {

                const auto next = dist_queue.top();
                dist_queue.pop();
                const auto i = next.index_;
                // ordering.push_back(i);
                ordering[pos--] = i;

                for(const std::size_t j : adj[i]) {
                    if (seen[j])
                        continue;
                    double avg_adj_dist = 0;
                    for (const auto n : adj[j])
                        avg_adj_dist += dist[n];
                    avg_adj_dist /= (double) adj[j].size();
                    dist_queue.emplace(j, dist[j], avg_adj_dist);
                    seen[j] = 1;
                }
            }
        }

        // std::reverse(ordering.begin(), ordering.end());

        assert(ordering.size() == adj.size());
        assert(is_permutation(ordering.begin(), ordering.end()));
        return ordering; 
    }
}
