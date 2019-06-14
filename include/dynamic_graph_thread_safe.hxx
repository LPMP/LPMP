#pragma once

#include <vector>
#include <array>
#include <numeric>
#include <limits>
#include <algorithm>
#include <tsl/robin_map.h>
#include "hash_helper.hxx"

// This class provides a similar interface as dynamic_graph_thread_safe, but edges can be added/deleted concurrently if they do not share endpoints
namespace LPMP {

    template<typename EDGE_INFORMATION>
    class dynamic_graph_thread_safe {
        //constexpr static std::size_t free_edge_store_factor = 2;  // for each node, we store the edges consecutively. For later edge additions, we leave some room for later input
        public:

		template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
            void construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func);

        dynamic_graph_thread_safe(const std::size_t no_nodes);

        constexpr static auto return_edge_op = [](const auto& edge) { return EDGE_INFORMATION{}; };

        template<typename EDGE_ITERATOR>
            dynamic_graph_thread_safe(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end)
            : dynamic_graph_thread_safe(edge_begin, edge_end, return_edge_op)
            {}

		template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
            dynamic_graph_thread_safe(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func);

		bool edge_present(const std::size_t i, const std::size_t j) const;
        EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j);
        const EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j) const;

		std::size_t no_nodes() const { return edge_maps_.size(); }
		std::size_t no_edges(const std::size_t i) const { assert(i < no_nodes()); return edge_maps_[i].size(); }

        void insert_edge(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& e);
        void remove_edge(const std::size_t i, const std::size_t j);
        void remove_node(const std::size_t i);

        const auto edges(const std::size_t i) const { assert(i < no_nodes()); return edge_maps_[i]; }

        private:
        std::array<std::size_t,2> normal_edge(const std::size_t i, const std::size_t j) const { return {std::min(i,j), std::max(i,j)}; }
        std::vector< tsl::robin_map<std::size_t, EDGE_INFORMATION> > edge_maps_;

        void assert_edge_valid(const std::size_t edge_index, const std::size_t tail_node = std::numeric_limits<std::size_t>::max()) const;
    };

    template<typename EDGE_INFORMATION>
        template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
        void dynamic_graph_thread_safe<EDGE_INFORMATION>::construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func)
        {
            edge_maps_.clear();

            std::vector<std::size_t> adjacency_list_count;
			for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
				adjacency_list_count.resize(std::max({i+1,j+1,adjacency_list_count.size()}));
				adjacency_list_count[i]++;
				adjacency_list_count[j]++; 
			}

            edge_maps_.resize(adjacency_list_count.size());
            for(std::size_t i=0; i<adjacency_list_count.size(); ++i)
                edge_maps_[i].reserve(adjacency_list_count[i]);

            for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
				const auto i = (*edge_it)[0];
				const auto j = (*edge_it)[1];
                insert_edge(i, j, edge_information_func(*edge_it));
            }
        }

    template<typename EDGE_INFORMATION>
        dynamic_graph_thread_safe<EDGE_INFORMATION>::dynamic_graph_thread_safe(const std::size_t no_nodes)
        {
            edge_maps_.resize(no_nodes);
        }

    template<typename EDGE_INFORMATION>
        template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
        dynamic_graph_thread_safe<EDGE_INFORMATION>::dynamic_graph_thread_safe(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func)
        {
            construct(edge_begin, edge_end, edge_information_func);
        }

    template<typename EDGE_INFORMATION>
        EDGE_INFORMATION& dynamic_graph_thread_safe<EDGE_INFORMATION>::edge(const std::size_t i, const std::size_t j)
        {
            assert(edge_present(i,j));
            //return edge_maps_.find(normal_edge(i,j))->second;
            return edge_maps_[i].find(j).value();
        }

    template<typename EDGE_INFORMATION>
        const EDGE_INFORMATION& dynamic_graph_thread_safe<EDGE_INFORMATION>::edge(const std::size_t i, const std::size_t j) const
        {
            assert(edge_present(i,j));
            return edge_maps_[i].find(j)->second;
        }

    template<typename EDGE_INFORMATION>
        bool dynamic_graph_thread_safe<EDGE_INFORMATION>::edge_present(const std::size_t i, const std::size_t j) const
        {
            assert(std::max(i,j) < no_nodes());
            assert(i != j);
            return edge_maps_[i].count(j) > 0;
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph_thread_safe<EDGE_INFORMATION>::insert_edge(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& edge_info)
        {
            assert(!edge_present(i,j));
            edge_maps_[i].insert(std::make_pair(j, edge_info));
            edge_maps_[j].insert(std::make_pair(i, edge_info));
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph_thread_safe<EDGE_INFORMATION>::remove_edge(const std::size_t i, const std::size_t j)
        {
            assert(edge_present(i,j));
            edge_maps_[i].erase(j);
            edge_maps_[j].erase(i);
        }

    template<typename EDGE_INFORMATION>
        void dynamic_graph_thread_safe<EDGE_INFORMATION>::remove_node(const std::size_t i)
        {
            assert(i < no_nodes());
            for(const auto& e : edge_maps_[i])
                edge_maps_[e.first].erase(i);
            edge_maps_[i].clear();
        }
}

