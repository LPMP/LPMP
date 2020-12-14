#pragma once

// Compute topological sorting of a DAG

#include <vector>
#include <array>
#include <stack>
#include <cassert>
#include "help_functions.hxx"
#include "two_dimensional_variable_array.hxx"

namespace LPMP {

namespace Topological_Sort {

// TODO: detect, when graph is not DAG
class Graph
{
   const std::size_t V;    // number of vertices
   std::vector<std::array<std::size_t,2>> edges;

public:
   inline Graph(const std::size_t V);
   inline void addEdge(std::size_t v, std::size_t w);
   inline std::vector<std::size_t> topologicalSort() const;
   inline bool sorting_valid(const std::vector<std::size_t>& ordering) const;
};
 
inline Graph::Graph(const std::size_t _V) : V(_V) {}

// add directed edge from v to w: this means v must come before w in sorting.
inline void Graph::addEdge(std::size_t v, std::size_t w)
{
   assert(v != w && v < V && w < V);
   edges.push_back({v,w});
}

inline bool Graph::sorting_valid(const std::vector<std::size_t>& ordering) const
{
   std::vector<std::size_t> inverse_ordering(ordering.size());
   for(std::size_t i=0; i<ordering.size(); ++i) {
      inverse_ordering[ordering[i]] = i; 
   }

   // check validity of sorting
   for(const auto e : edges) {
      const auto i = e[0];
      const auto j = e[1];
      assert(inverse_ordering[i] != inverse_ordering[j]);
      if(inverse_ordering[i] > inverse_ordering[j]) {
         return false; 
      }
   }
   return true; 
}

inline std::vector<std::size_t> topological_sort_adjacency_list(const two_dim_variable_array<size_t>& adj)
{
  constexpr unsigned char not_marked = 0;
  constexpr unsigned char temp_marked = 1;
  constexpr unsigned char perm_marked = 2;

  std::vector<unsigned char> mark(adj.size(),not_marked); 
  std::stack<std::pair<std::size_t,decltype(adj[0].begin())> > dfs;
  std::vector<std::size_t> post_order;
  post_order.reserve(adj.size());

  for(std::size_t i=0; i<adj.size(); ++i){
     if(mark[i] == not_marked) {
        dfs.push(std::make_pair(i, adj[i].begin()));
        while(!dfs.empty()){
           const std::size_t node = dfs.top().first;
           auto it = dfs.top().second;
           //if(mark[node] == perm_marked) {
           //   assert(false);
           //   throw std::runtime_error("graph not a dag");
           //}
           if(it == adj[node].begin()) { // first visit
              assert(mark[node] == not_marked);
              mark[node] = temp_marked;
           }
           while(it != adj[node].end() && mark[*it] != not_marked) {
              if(mark[*it] == temp_marked) {
                 throw std::runtime_error("graph not a dag");
              }
              ++it;
           }
           if(it != adj[node].end()) {
              std::size_t next = *it;
              ++it;
              dfs.top().second = it;
              dfs.push(std::make_pair(next,adj[next].begin()));
           } else {
              dfs.pop();
              mark[node] = perm_marked;
              post_order.push_back(node);
           }
        }
     }
  }

  assert(post_order.size() == adj.size());
  assert(LPMP::HasUniqueValues(post_order));
  //assert(sorting_valid(post_order));

  return post_order;

}

inline std::vector<std::size_t> Graph::topologicalSort() const
{
  std::vector<std::size_t> adjacency_list_size(V);
  for(const auto e : edges) {
     adjacency_list_size[e[1]]++;
  }

  two_dim_variable_array<std::size_t> adjacency_list(adjacency_list_size.begin(), adjacency_list_size.end());
  std::fill(adjacency_list_size.begin(), adjacency_list_size.end(), 0);

  // add edges in backward direction, since post_order will record nodes in backward topological sorting order.
  for(const auto e : edges) {
     adjacency_list(e[1],adjacency_list_size[e[1]]) = e[0];
     adjacency_list_size[e[1]]++;
  }

  return topological_sort_adjacency_list(adjacency_list);

}

} // namespace Topological_Sort

} // namespace LPMP
