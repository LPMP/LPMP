#pragma once

#include <vector>
#include <array>
#include <numeric>
#include <limits>
#include <algorithm>
#include <tsl/robin_map.h>
#include "hash_helper.hxx"

namespace LPMP {

    template<typename EDGE_INFORMATION>
    class dynamic_graph {
        //constexpr static std::size_t free_edge_store_factor = 2;  // for each node, we store the edges consecutively. For later edge additions, we leave some room for later input
        public:

        constexpr static std::size_t sister_pointer_inactive = std::numeric_limits<std::size_t>::max();
        constexpr static std::size_t no_next_edge = std::numeric_limits<std::size_t>::max();
        struct edge_type  {
            std::size_t head;
            std::size_t sister = sister_pointer_inactive;
            std::size_t prev_outgoing_edge = no_next_edge;
            std::size_t next_outgoing_edge = no_next_edge;
        };

		template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
            void construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func);

        dynamic_graph(const std::size_t no_nodes);

        constexpr static auto return_edge_op = [](const auto& edge) { return EDGE_INFORMATION{}; };

        template<typename EDGE_ITERATOR>
            dynamic_graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end)
            : dynamic_graph(edge_begin, edge_end, return_edge_op)
            {}

		template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
            dynamic_graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func);

		bool edge_present(const std::size_t i, const std::size_t j) const;
        EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j);
        const EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j) const;
        std::size_t tail(const std::size_t edge_index) const;
        std::size_t head(const std::size_t edge_index) const;
        std::array<std::size_t,2> endpoints(const std::size_t edge_index) const;

		std::size_t no_nodes() const { return nodes_.size(); }
		std::size_t no_edges(const std::size_t i) const { assert(i < no_nodes()); return nodes_[i].no_edges; }

        std::size_t first_outgoing_edge_index(const std::size_t i) const { return nodes_[i].first_outgoing_edge; }
        std::size_t next_outgoing_edge_index(const std::size_t e) const { return edges_[e].next_outgoing_edge; }

        void insert_edge(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& e);
        void remove_edge(const std::size_t e);
        void remove_node(const std::size_t i);

        private:
        std::vector<edge_type> edges_;
        struct node_item {
            std::size_t first_outgoing_edge = no_next_edge;
            std::size_t last_outgoing_edge = no_next_edge;
            std::size_t no_edges = 0;
        };
        std::vector<node_item> nodes_;
        std::size_t free_edge_list_ = no_next_edge;

        std::array<std::size_t,2> normal_edge(const std::size_t i, const std::size_t j) const { return {std::min(i,j), std::max(i,j)}; }
        tsl::robin_map<std::array<std::size_t,2>, EDGE_INFORMATION> edge_map_;

        std::size_t get_edge_space();

        void assert_edge_valid(const std::size_t edge_index, const std::size_t tail_node = std::numeric_limits<std::size_t>::max()) const;
    };

    template<typename EDGE_INFORMATION>
        template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
        void dynamic_graph<EDGE_INFORMATION>::construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func)
        {
            nodes_.clear();
            edges_.clear();
            edge_map_.clear();
            free_edge_list_ = no_next_edge;

            std::vector<std::size_t> adjacency_list_count;
			for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
				adjacency_list_count.resize(std::max({i+1,j+1,adjacency_list_count.size()}));
				adjacency_list_count[i]++;
				adjacency_list_count[j]++; 
			}

            nodes_.resize(adjacency_list_count.size());
            const std::size_t no_edges = std::accumulate(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
            edges_.reserve(no_edges);
            edge_map_.reserve(no_edges);

            for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
                insert_edge(i, j, edge_information_func(*edge_it));
            }
            /*
            std::vector<std::size_t> edge_offsets{0};
            edge_offsets.reserve(adjacency_list_count.size()+1);
            std::partial_sum(adjacency_list_count.begin(), adjacency_list_count.end(), std::back_inserter(edge_offsets));
            for(std::size_t i=0; i<adjacency_list_count.size(); ++i) {
                nodes_[i].first_outgoing_edge = edge_offsets[i];
                nodes_[i].last_outgoing_edge = edge_offsets[i];
                nodes_[i].no_edges = adjacency_list_count[i];
            }

            std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);

            for(std::size_t i=0; i<adjacency_list_count.size(); ++i) {
                edges_[edge_offsets[i]].prev_outgoing_edge = no_next_edge;
                for(std::size_t j=0; j<adjacency_list_count[i]-1; ++j) {
                    edges_[edge_offsets[i] + j].next_outgoing_edge = edge_offsets[i] + j + 1;
                    edges_[edge_offsets[i] + j + 1].prev_outgoing_edge = edge_offsets[i] + j;
                }
                edges_[edge_offsets[i] + adjacency_list_count[i] - 1].next_outgoing_edge = no_next_edge;
            }

            for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];

                edges_[ nodes_[i].last_outgoing_edge ] = edge_type{j, nodes_[j].last_outgoing_edge};
                edges_[ nodes_[j].last_outgoing_edge ] = edge_type{i, nodes_[i].last_outgoing_edge};
                //if(nodes_[i].first_outgoing_edge != nodes_[i].last_outgoing_edge) {
                 //   edges_[nodes_[i].last_outgoing_edge-1].next_outgoing_edge = nodes_[i].last_outgoing_edge;
                 //   edges_[nodes_[i].last_outgoing_edge].prev_outgoing_edge = nodes_[i].last_outgoing_edge-1;
                //}
                //if(nodes_[j].first_outgoing_edge != nodes_[j].last_outgoing_edge) {
                //    edges_[nodes_[j].last_outgoing_edge-1].next_outgoing_edge = nodes_[j].last_outgoing_edge;
                //    edges_[nodes_[j].last_outgoing_edge].prev_outgoing_edge = nodes_[j].last_outgoing_edge-1;
                //}
                nodes_[i].last_outgoing_edge++;
                nodes_[j].last_outgoing_edge++;

                assert(!edge_present(i,j));
                edge_map_.insert(std::make_pair(normal_edge(i,j), edge_information_func(*edge_it)));
            }
            for(std::size_t i=0; i<adjacency_list_count.size(); ++i)
                nodes_[i].last_outgoing_edge--;
                */

            for(std::size_t e=0; e<edges_.size(); ++e)
                assert_edge_valid(e, tail(e));
            for(std::size_t i=0; i<no_nodes(); ++i)
                for(std::size_t e=nodes_[i].first_outgoing_edge; e!=no_next_edge; e=next_outgoing_edge_index(e)) 
                    assert_edge_valid(e,i);
        }

    template<typename EDGE_INFORMATION>
        dynamic_graph<EDGE_INFORMATION>::dynamic_graph(const std::size_t no_nodes)
        {
            nodes_.clear();
            edges_.clear();
            edge_map_.clear();
            free_edge_list_ = no_next_edge;
            nodes_.resize(no_nodes);
        }

    template<typename EDGE_INFORMATION>
        template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
        dynamic_graph<EDGE_INFORMATION>::dynamic_graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func)
        {
            construct(edge_begin, edge_end, edge_information_func);
        }

    template<typename EDGE_INFORMATION>
        EDGE_INFORMATION& dynamic_graph<EDGE_INFORMATION>::edge(const std::size_t i, const std::size_t j)
        {
            assert(edge_present(i,j));
            //return edge_map_.find(normal_edge(i,j))->second;
            return edge_map_.find(normal_edge(i,j)).value();
        }

    template<typename EDGE_INFORMATION>
        const EDGE_INFORMATION& dynamic_graph<EDGE_INFORMATION>::edge(const std::size_t i, const std::size_t j) const
        {
            assert(edge_present(i,j));
            return edge_map_.find(normal_edge(i,j))->second;
        }

    template<typename EDGE_INFORMATION>
        std::size_t dynamic_graph<EDGE_INFORMATION>::tail(const std::size_t edge_index) const
        {
            assert(edge_index < edges_.size());
            return head(edges_[edge_index].sister); 
        }

    template<typename EDGE_INFORMATION>
        std::size_t dynamic_graph<EDGE_INFORMATION>::head(const std::size_t edge_index) const
        {
            assert(edge_index < edges_.size());
            return edges_[edge_index].head; 
        }

    template<typename EDGE_INFORMATION>
        std::array<std::size_t,2> dynamic_graph<EDGE_INFORMATION>::endpoints(const std::size_t edge_index) const
        {
            assert(edge_index < edges_.size());
            return {tail(edge_index), head(edge_index)};
        }

    template<typename EDGE_INFORMATION>
        bool dynamic_graph<EDGE_INFORMATION>::edge_present(const std::size_t i, const std::size_t j) const
        {
            assert(std::max(i,j) < no_nodes());
            assert(i != j);
            return edge_map_.count(normal_edge(i,j)) > 0;
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph<EDGE_INFORMATION>::insert_edge(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& edge_info)
        {
            assert(!edge_present(i,j));
            const std::size_t e = get_edge_space();
            const std::size_t e_sister = get_edge_space();
            //std::cout << "insert edge " << i << "," << j << " with indices " << e << "," << e_sister << "\n";

            edges_[e].head = j;
            edges_[e_sister].head = i;

            edges_[e].sister = e_sister;
            edges_[e_sister].sister = e;

            auto insert_list = [&](const std::size_t tail_node, const std::size_t edge_index) {
                assert(edge_index < edges_.size());
                if(nodes_[tail_node].last_outgoing_edge != no_next_edge) {
                    assert(edges_[nodes_[tail_node].last_outgoing_edge].next_outgoing_edge  == no_next_edge);
                    edges_[nodes_[tail_node].last_outgoing_edge].next_outgoing_edge = edge_index;
                }
                edges_[edge_index].prev_outgoing_edge = nodes_[tail_node].last_outgoing_edge;
                nodes_[tail_node].last_outgoing_edge = edge_index;
                if(nodes_[tail_node].first_outgoing_edge == no_next_edge) {
                    nodes_[tail_node].first_outgoing_edge = edge_index;
                } else {
                    assert(edges_[nodes_[tail_node].first_outgoing_edge].prev_outgoing_edge == no_next_edge); 
                }
                edges_[edge_index].next_outgoing_edge = no_next_edge;
            };

            insert_list(i, e);
            insert_list(j, e_sister);

            edge_map_.insert(std::make_pair(normal_edge(i,j), edge_info));

            nodes_[i].no_edges++;
            nodes_[j].no_edges++;

            assert_edge_valid(e, i);
            assert_edge_valid(e_sister, j);
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph<EDGE_INFORMATION>::remove_edge(const std::size_t e)
        {
            assert_edge_valid(e);
            const auto [i,j] = endpoints(e);
            //std::cout << "remove edge with endpoints " << i << "," << j << "\n";
            assert(edge_present(i,j));
            edge_map_.erase(normal_edge(i,j));
            const std::size_t e_sister = edges_[e].sister;

            auto set_prev_next = [&](const std::size_t tail_node, const std::size_t edge_index) {
                assert(tail(edge_index) == tail_node);
                if(edges_[edge_index].prev_outgoing_edge != no_next_edge)
                    edges_[edges_[edge_index].prev_outgoing_edge].next_outgoing_edge = edges_[edge_index].next_outgoing_edge;
                if(edges_[edge_index].next_outgoing_edge != no_next_edge)
                    edges_[edges_[edge_index].next_outgoing_edge].prev_outgoing_edge = edges_[edge_index].prev_outgoing_edge;
                if(nodes_[tail_node].first_outgoing_edge == edge_index) {
                    assert(edges_[edge_index].prev_outgoing_edge == no_next_edge);
                    nodes_[tail_node].first_outgoing_edge = edges_[edge_index].next_outgoing_edge;
                }
                if(nodes_[tail_node].last_outgoing_edge == edge_index) {
                    assert(edges_[edge_index].next_outgoing_edge == no_next_edge);
                    nodes_[tail_node].last_outgoing_edge = edges_[edge_index].prev_outgoing_edge;
                }
            };

            set_prev_next(i, e);
            set_prev_next(j, e_sister);

            edges_[e] = {std::numeric_limits<std::size_t>::max()};
            edges_[e_sister] = {std::numeric_limits<std::size_t>::max()};

            edges_[e_sister].next_outgoing_edge = free_edge_list_;
            edges_[e].next_outgoing_edge = e_sister;
            free_edge_list_ = e;

            assert(nodes_[i].no_edges > 0);
            nodes_[i].no_edges--;
            assert(nodes_[j].no_edges > 0);
            nodes_[j].no_edges--;
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph<EDGE_INFORMATION>::remove_node(const std::size_t i)
        {
            assert(i < no_nodes());
            while(no_edges(i) > 0) {
                const std::size_t e = first_outgoing_edge_index(i);
                remove_edge(e);
                //std::cout << "remove edge " << e << "," << edges_[e].sister << "\n";
            }
        }

    template<typename EDGE_INFORMATION>
        std::size_t dynamic_graph<EDGE_INFORMATION>::get_edge_space()
        {
            if(free_edge_list_ != no_next_edge) {
                const std::size_t e = free_edge_list_;
                edges_[e] = {std::numeric_limits<std::size_t>::max()};
                free_edge_list_ = edges_[free_edge_list_].next_outgoing_edge; 
                return e;
            } else {
                edges_.push_back({});
                return edges_.size()-1;
            } 
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph<EDGE_INFORMATION>::assert_edge_valid(const std::size_t edge_index, const std::size_t tail_node) const
        {
            assert(edge_index < edges_.size());
            assert(edges_[edge_index].sister != edge_index);
            assert(edges_[edge_index].sister < edges_.size());
            assert(edges_[edges_[edge_index].sister].sister == edge_index);
            if(tail_node != std::numeric_limits<std::size_t>::max()) {
                assert(tail(edge_index) == tail_node);
                assert(edges_[edge_index].head != tail_node);
                assert(tail_node < nodes_.size());
                assert(edges_[edges_[edge_index].sister].head == tail_node);
                if(nodes_[tail_node].first_outgoing_edge == edge_index) {
                    assert(edges_[edge_index].prev_outgoing_edge == no_next_edge);
                }
                if(nodes_[tail_node].last_outgoing_edge == edge_index) {
                    assert(edges_[edge_index].next_outgoing_edge == no_next_edge);
                }
            }
            assert(edge_present(tail(edge_index), head(edge_index)));
        }
}
