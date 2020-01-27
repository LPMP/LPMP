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
        permutation ordering(adj.size());

        const std::size_t s = find_pseudo_peripheral_node(adj);

        std::queue<std::size_t> Q;
        std::vector<char> seen(adj.size(), 0);
        std::vector<std::size_t> dist(adj.size(), std::numeric_limits<std::size_t>::max());
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
                if (dist[j] >= dist[i])
                    terminal = false;
                if (seen[j])
                    continue;
                seen[j] = 1;
                dist[j] = dist[i] + 1;
                Q.push(j);
            }

            if(terminal)
                terminals.push_back(i);
        }

        struct node
        {
            node(std::size_t index, std::size_t dist) : index_(index), dist_(dist)
            {}

            bool operator <(node const & other) const { return dist_ < other.dist_; }

            std::size_t index_;
            std::size_t dist_;
        };

        std::fill(seen.begin(), seen.end(), 0);
        std::priority_queue<node> dist_queue;
        for(const std::size_t t : terminals) {
            dist_queue.emplace(t, dist[t]);
            seen[t] = 1;
        }

        // backward run of BFS from all terminals, ordering vertices in decreasing distance from source
        std::size_t pos = adj.size()-1;
        while(!dist_queue.empty()) {

            // get next batch of nodes
            std::vector<node> batch;
            while (!dist_queue.empty())
            {
                const auto n = dist_queue.top();
                dist_queue.pop();
                batch.push_back(n);

            }
            for (const auto n : batch)
            {
                const auto i = n.index_;
                // ordering.push_back(i);
                ordering[pos--] = i;

                for(const std::size_t j : adj[i]) {
                    if (seen[j])
                        continue;
                    dist_queue.emplace(j, dist[j]);
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
