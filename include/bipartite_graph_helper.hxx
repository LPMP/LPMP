#pragma once

#include <vector>
#include "graph.hxx"
#include "union_find.hxx"

namespace LPMP {

    // Given a list of edges, compress underying nodes to contiguous range and build bipartite graph out of them. 
    // Used in tightening for multicut (odd wheel, odd bicycle wheel)
    template<typename EDGE_INFORMATION>
        class compressed_bipartite_graph_helper {
            public:
                constexpr static std::size_t node_not_present = std::numeric_limits<std::size_t>::max();

                compressed_bipartite_graph_helper(const std::size_t max_no_nodes);

                // edge_func returns endpoints of edge.
                // edge_info_func populates the edge information
                template<typename EDGE_ITERATOR, typename EDGE_FUNC, typename EDGE_INFORMATION_FUNC>
                    void construct_compressed_bipartite_graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_FUNC edge_func, EDGE_INFORMATION_FUNC edge_info_func);
                std::size_t no_compressed_nodes() const { return compressed_to_orig_node.size(); }

                auto& get_bfs() { return bfs_; }
                auto& get_graph() { return g_; }
                union_find& get_union_find() { return uf_; }
                std::size_t compressed_to_original_node(const std::size_t c) const { assert(c < compressed_to_orig_node.size()); return compressed_to_orig_node[c]; }
                void compressed_path_to_original(std::vector<std::size_t>& p) const;

            private:
                std::vector<std::size_t> orig_to_compressed_nodes; // compresses node indices
                std::vector<std::size_t> compressed_to_orig_node; // compressed nodes to original
                struct compressed_edge : public std::array<std::size_t,2> { 
                    EDGE_INFORMATION e; 
                };
                std::vector<compressed_edge> compressed_bipartite_edges;

                graph<EDGE_INFORMATION> g_;
                bfs_data<decltype(g_)> bfs_; 
                union_find uf_;
        };

    template<typename EDGE_INFORMATION>
        compressed_bipartite_graph_helper<EDGE_INFORMATION>::compressed_bipartite_graph_helper(const std::size_t max_no_nodes)
        : orig_to_compressed_nodes(max_no_nodes, compressed_bipartite_graph_helper<EDGE_INFORMATION>::node_not_present),
        bfs_(g_)
        {}

    template<typename EDGE_INFORMATION>
        template<typename EDGE_ITERATOR, typename EDGE_FUNC, typename EDGE_INFORMATION_FUNC>
        void compressed_bipartite_graph_helper<EDGE_INFORMATION>::construct_compressed_bipartite_graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_FUNC edge_func, EDGE_INFORMATION_FUNC edge_info_func)
        {
            // clear originial nodes to compressed nodes
            for(const std::size_t orig_node : compressed_to_orig_node)
                orig_to_compressed_nodes[orig_node] = node_not_present;
            compressed_to_orig_node.clear();
            for(const auto v : orig_to_compressed_nodes) {
                assert(v == node_not_present);
            }

            auto add_compressed_node = [&](const std::size_t orig_node) {
                if(orig_to_compressed_nodes[orig_node] == node_not_present) {
                    orig_to_compressed_nodes[orig_node] = compressed_to_orig_node.size();
                    compressed_to_orig_node.push_back(orig_node);
                }
            };
            for(auto it=edge_begin; it!=edge_end; ++it) {
                add_compressed_node(edge_func(*it)[0]);
                add_compressed_node(edge_func(*it)[1]);
            }

            uf_.init(2*compressed_to_orig_node.size());
            const std::size_t no_compressed_nodes = compressed_to_orig_node.size();
            compressed_bipartite_edges.clear();
            for(auto it=edge_begin; it!=edge_end; ++it) {
                const auto edge_info = edge_info_func(*it);
                const std::size_t ci = orig_to_compressed_nodes[edge_func(*it)[0]];
                const std::size_t cj = orig_to_compressed_nodes[edge_func(*it)[1]];
                assert(ci != cj && std::max(ci,cj) < no_compressed_nodes);
                compressed_bipartite_edges.push_back(compressed_edge{ci, cj + no_compressed_nodes, edge_info});
                compressed_bipartite_edges.push_back(compressed_edge{ci + no_compressed_nodes, cj, edge_info});
                uf_.merge(ci, cj + no_compressed_nodes);
                uf_.merge(cj, ci + no_compressed_nodes);
            }

            auto copy_edge_info = [](const compressed_edge& e) { return e.e; };
            g_.construct(compressed_bipartite_edges.begin(), compressed_bipartite_edges.end(), copy_edge_info);
            //std::cout << "bipartite search graph has " << compressed_bipartite_edges.size() << " edges\n";
            bfs_.update_graph();
        }

    template<typename EDGE_INFORMATION>
        void compressed_bipartite_graph_helper<EDGE_INFORMATION>::compressed_path_to_original(std::vector<std::size_t>& p) const
        {
            for(auto& v : p)
                v = compressed_to_original_node(v % no_compressed_nodes());
        }

}
