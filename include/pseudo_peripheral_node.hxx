#pragma once

#include "union_find.hxx"
#include <queue>
#include <vector>
#include <tuple>
#include <algorithm>

namespace LPMP {

template <typename ADJACENCY_GRAPH>
std::tuple<std::size_t, std::size_t> farthest_node(const ADJACENCY_GRAPH &adjacency, const std::size_t x)
{
    std::vector<char> visited(adjacency.size(), 0);
    return farthest_node(adjacency, x, visited, 0);
}

template <typename ADJACENCY_GRAPH, typename VISITED_VECTOR>
std::tuple<std::size_t, std::size_t> farthest_node(const ADJACENCY_GRAPH &adjacency, const std::size_t x, VISITED_VECTOR &visited, const std::size_t timestamp)
{
    assert(visited.size() == adjacency.size());
    assert(*std::max_element(visited.begin(), visited.end()) <= timestamp);

    std::size_t d = 0;
    struct queue_elem
    {
        std::size_t v;
        std::size_t d;
    };
    std::queue<queue_elem> Q;
    Q.push({x, 0});
    visited[x] = timestamp + 1;
    std::size_t farthest_node = x;
    std::size_t max_distance = 0;
    while (!Q.empty())
    {
        const auto [i, d] = Q.front();
        Q.pop();
        assert(visited[i] == timestamp + 1);
        visited[i] = timestamp + 2;
        if (d > max_distance)
        {
            max_distance = d;
            farthest_node = x;
        }
        for (const auto j : adjacency[i])
        {
            if (visited[j] <= timestamp)
            {
                Q.push({j, d + 1});
                visited[j] = timestamp + 1;
            }
        }
    }

    return {farthest_node, max_distance};
}

    template<typename ADJACENCY_GRAPH>
    std::size_t find_pseudo_peripheral_node(const ADJACENCY_GRAPH& adjacency)
    {
        std::size_t min_degree = adjacency[0].size();
        std::size_t x = 0;
        for(std::size_t i=0; i<adjacency.size(); ++i) {
            if(adjacency[i].size() < min_degree) {
                min_degree = adjacency[i].size();
                x = i;
            }
        }

        assert(x < adjacency.size());
        auto [y, d_y] = farthest_node(adjacency, x);
        auto [z, d_z] = farthest_node(adjacency, y);
        while(d_z > d_y) {
            std::swap(y,z);
            std::swap(d_z, d_y);
            std::tie(z, d_z) = farthest_node(adjacency,y);
        }
        return y; 
    }

    // return pseudo peripheral node of each connected component
    template<typename ADJACENCY_GRAPH>
    std::vector<std::size_t> find_pseudo_peripheral_nodes(const ADJACENCY_GRAPH& adjacency)
    {
        // compute connected components of graph and for each connected component determine node of minimum degree
        union_find uf(adjacency.size());
        for(std::size_t i=0; i<adjacency.size(); ++i)
            for (const std::size_t j : adjacency[i])
                uf.merge(i, j);

        struct min_degree_elem
        {
            std::size_t degree = std::numeric_limits<std::size_t>::max();
            std::size_t node = std::numeric_limits<std::size_t>::max();
        };
        std::vector<min_degree_elem> min_degree(adjacency.size());

        for(std::size_t i=0; i<adjacency.size(); ++i)
        {
            const std::size_t cc_id = uf.find(i);
            const std::size_t d = adjacency[i].size();
            if (adjacency[i].size() < min_degree[cc_id].degree)
            {
                min_degree[cc_id].degree = adjacency[i].size();
                min_degree[cc_id].node = i;
            }
        }

        std::vector<std::size_t> pseudo_peripheral_nodes;
        std::vector<std::size_t> visited(adjacency.size(), 0);
        std::size_t iter = 0;

        for (std::size_t i = 0; i < adjacency.size(); ++i)
        { 
            if(visited[i] != 0)
                continue;

            const std::size_t x = min_degree[uf.find(i)].node;

            assert(x < adjacency.size());

            auto [y, d_y] = farthest_node(adjacency, x, visited, 2*iter++);
            auto [z, d_z] = farthest_node(adjacency, y, visited, 2*iter++);
            while (d_z > d_y)
            {
                std::swap(y, z);
                std::swap(d_z, d_y);
                std::tie(z, d_z) = farthest_node(adjacency, y, visited, 2*iter++);
            }
            pseudo_peripheral_nodes.push_back(y);
        }

        return pseudo_peripheral_nodes;
    }

    template<typename ADJACENCY_GRAPH, typename NODE_ITERATOR>
    std::size_t find_pseudo_peripheral_node(const ADJACENCY_GRAPH& adjacency, NODE_ITERATOR node_begin, NODE_ITERATOR node_end)
    {
        //assert(std::distance(node_begin, node_end) > 0);
        std::size_t min_degree = adjacency[*node_begin].size();
        std::size_t x = *node_begin;
        for(auto node_it=node_begin; node_it!=node_end; ++node_it) {
            if(adjacency[*node_it].size() < min_degree) {
                min_degree = adjacency[*node_it].size();
                x = *node_it;
            }
        }

        assert(x < adjacency.size());
        auto [y, d_y] = farthest_node(adjacency, x);
        auto [z, d_z] = farthest_node(adjacency, y);
        while(d_z > d_y) {
            std::swap(y,z);
            std::swap(d_z, d_y);
            std::tie(z, d_z) = farthest_node(adjacency,y);
        }
        return y; 
    }

    /*
    struct iterator {
        std::size_t i;
        std::size_t operator*() const { return i; }
        void operator++() { ++i; }
        bool operator!=(const iterator o) { return o.i != this->i; }
        std::size_t operator-(const iterator o) const { return o.i - this->i; }
    }; 

    template<typename ADJACENCY_GRAPH>
    std::size_t find_pseudo_peripheral_node(const ADJACENCY_GRAPH& adjacency)
    {
        return find_pseudo_peripheral_node(adjacency, iterator({0}), iterator({adjacency.size()}));
    }
    */

} 
