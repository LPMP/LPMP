#pragma once

#include "two_dimensional_variable_array.hxx"
#include "pseudo_peripheral_node.hxx"
#include "permutation.hxx"
#include <vector>
#include <queue>
#include <tuple>
#include <algorithm>

namespace LPMP {

    template<typename T>
    permutation Cuthill_McKee(const two_dim_variable_array<T>& adjacency)
    {
        std::queue<std::size_t> Q;
        permutation result;
        result.reserve(adjacency.size());
        std::vector<std::size_t> remaining_degree;
        remaining_degree.reserve(adjacency.size());
        std::vector<char> visited(adjacency.size(), 0);

        for(std::size_t i=0; i<adjacency.size(); ++i) {
            remaining_degree.push_back(adjacency[i].size());
        }

        const auto pseudo_peripheral_nodes = find_pseudo_peripheral_nodes(adjacency);

        for (const std::size_t i : pseudo_peripheral_nodes)
        {
            result.push_back(i);
            Q.push(i);
            visited[i] = 1;

            while (!Q.empty())
            {
                const std::size_t i = Q.front();
                Q.pop();
                assert(visited[i] == 1);
                visited[i] = 2;
                std::vector<std::size_t> a;
                for (const std::size_t x : adjacency[i])
                {
                    if (visited[x] == 0 && x != i)
                    { // TODO: second check not needed!
                        a.push_back(x);
                    }
                    if (visited[x] != 2)
                    {
                        assert(remaining_degree[i] > 0);
                        remaining_degree[i]--;
                        assert(remaining_degree[x] > 0);
                        remaining_degree[x]--;
                    }
                }
                std::sort(a.begin(), a.end(), [&](const std::size_t x, const std::size_t y) { return remaining_degree[x] < remaining_degree[y]; });
                for (const auto x : a)
                {
                    Q.push(x);
                    result.push_back(x);
                    visited[x] = 1;
                }
                assert(remaining_degree[i] == 0);
            }
        }

        if(result.size() != adjacency.size())
            throw std::runtime_error("Graph not connected.");

        assert(is_permutation(result.begin(), result.end()));
        return result;
    }

}
