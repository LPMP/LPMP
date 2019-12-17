#pragma once
#include "two_dimensional_variable_array.hxx"
#include <vector>
#include <queue>
#include <tuple>
#include <algorithm>

namespace LPMP {

    template<typename T>
    std::tuple<std::size_t, std::size_t> farthest_node(const two_dim_variable_array<T>& adjacency, const std::size_t x)
    {
        std::size_t d=0;
        struct queue_elem { std::size_t v; std::size_t d; };
        std::queue<queue_elem> Q;
        Q.push({x,0});
        std::vector<char> visited(adjacency.size(), 0);
        visited[x] = 1;
        std::size_t farthest_node = x;
        std::size_t max_distance = 0;
        while(!Q.empty()) {
            const auto [i,d] = Q.front();
            Q.pop();
            assert(visited[i] == 1);
            visited[i] = 2;
            if(d > max_distance) {
                max_distance = d;
                farthest_node = x;
            }
            for(const auto j : adjacency[i]) {
                if(visited[j] == 0) {
                    Q.push({j, d+1});
                    visited[j] = 1;
                } 
            }

        }

        return {farthest_node, max_distance};
    }

    template<typename T>
    std::size_t find_pseudo_peripheral_node(const two_dim_variable_array<T>& adjacency)
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

    template<typename T>
    std::vector<std::size_t> Cuthill_McKee(const two_dim_variable_array<T>& adjacency)
    {
        std::queue<std::size_t> Q;
        std::vector<std::size_t> result;
        result.reserve(adjacency.size());
        std::vector<std::size_t> remaining_degree;
        remaining_degree.reserve(adjacency.size());
        std::vector<char> visited(adjacency.size(), 0);

        for(std::size_t i=0; i<adjacency.size(); ++i) {
            remaining_degree.push_back(adjacency[i].size());
        }

        const std::size_t i = find_pseudo_peripheral_node(adjacency);
        result.push_back(i);
        Q.push(i);
        visited[i] = 1;

        while(!Q.empty()) {
            const std::size_t i = Q.front();
            Q.pop();
            assert(visited[i] == 1);
            visited[i] = 2;
            std::vector<std::size_t> a;
            for(const std::size_t x : adjacency[i]) {
                if(visited[x] == 0 && x != i) {
                    a.push_back(x); 
                }
                if(visited[x] != 2) {
                    assert(remaining_degree[i] > 0);
                    remaining_degree[i]--;
                    assert(remaining_degree[x] > 0);
                    remaining_degree[x]--;
                }
            }
            std::sort(a.begin(), a.end(), [&](const std::size_t x, const std::size_t y) { return remaining_degree[x] < remaining_degree[y]; });
            for(const auto x : a) {
                Q.push(x);
                result.push_back(x);
                visited[x] = 1;
            }
            assert(remaining_degree[i] == 0);
        } 

        if(result.size() != adjacency.size())
            throw std::runtime_error("Graph not connected.");

        return result;
    }

}
