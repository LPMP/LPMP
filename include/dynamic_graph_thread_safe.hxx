#pragma once

#include <vector>
#include <array>
#include <numeric>
#include <limits>
#include <algorithm>
#include <tsl/robin_map.h>
#include <taskflow/taskflow.hpp>
#include "hash_helper.hxx"

// This class provides a similar interface as dynamic_graph, but edges can be added/deleted concurrently if they do not share endpoints
namespace LPMP {

    template<typename EDGE_INFORMATION>
    class dynamic_graph_thread_safe {
        public:

            template<typename DEGREE_ITERATOR>
                void pre_allocate(DEGREE_ITERATOR degree_begin, DEGREE_ITERATOR degree_end);

            template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
                void construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA edge_information_func);

            dynamic_graph_thread_safe() {}
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

            void insert_arc(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& e);
            void insert_edge(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& e);
            void remove_edge(const std::size_t i, const std::size_t j);
            void remove_node(const std::size_t i);

            const auto edges(const std::size_t i) const { assert(i < no_nodes()); return edge_maps_[i]; }

        private:
            std::array<std::size_t,2> normal_edge(const std::size_t i, const std::size_t j) const { return {std::min(i,j), std::max(i,j)}; }
            std::vector< tsl::robin_map<std::size_t, EDGE_INFORMATION> > edge_maps_;
    };

    template<typename EDGE_INFORMATION>
        template<typename DEGREE_ITERATOR>
        void dynamic_graph_thread_safe<EDGE_INFORMATION>::pre_allocate(DEGREE_ITERATOR degree_begin, DEGREE_ITERATOR degree_end)
        {
            assert(std::distance(degree_begin, degree_end) == no_nodes());
            auto degree_iter = degree_begin;
            for(std::size_t i=0; i<no_nodes(); ++i, ++degree_iter) {
                const std::size_t degree = *degree_iter;
                edge_maps_[i].reserve(degree);
            }
        }

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

            // parallel implementation slower than sequential one
            /*
               tf::Executor executor;
               tf::Taskflow taskflow;

               auto [E_begin, E_end] = taskflow.parallel_for(std::size_t(0), adjacency_list_count.size(), std::size_t(1), [&](const std::size_t i) {
               edge_maps_[i].reserve(adjacency_list_count[i]);
               });

               std::vector<std::atomic<char>> edge_occupied(edge_maps_.size());
               std::fill(edge_occupied.begin(), edge_occupied.end(), false);

               auto mark = [&](const std::size_t i) {
               char mark_val = false;
               auto& mark_var = edge_occupied[i];
               while(!mark_var.compare_exchange_strong(mark_val, 1, std::memory_order_relaxed));
               };
               auto unmark = [&](const std::size_t i) {
               assert(edge_occupied[i] == true);
               edge_occupied[i].store(false);
               };

               auto [I_begin, I_end] = taskflow.parallel_for(edge_begin, edge_end, [&](auto edge) {
               const auto i = edge[0];
               const auto j = edge[1];
               mark(i);
               insert_arc(i, j, edge_information_func(edge));
               unmark(i);
               mark(j);
               insert_arc(j, i, edge_information_func(edge));
               unmark(j);
               });

               I_begin.gather(E_end);

               executor.run(taskflow);
               executor.wait_for_all();
             */
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
        void dynamic_graph_thread_safe<EDGE_INFORMATION>::insert_arc(const std::size_t i, const std::size_t j, const EDGE_INFORMATION& edge_info)
        {
            edge_maps_[i].insert(std::make_pair(j, edge_info));
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
