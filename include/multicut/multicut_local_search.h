#pragma once

#include <vector>
#include <array>
#include <tuple>
#include "multicut_instance.h"
#include "graph.hxx"

namespace LPMP {

    class multicut_local_search {
        public:
            multicut_local_search(const multicut_instance& instance, const multicut_node_labeling& labeling);
            void construct(const multicut_instance& instance, const multicut_node_labeling& labeling);

            double compute_cut_value(const std::size_t i, const std::size_t k) const;
            void move(const std::size_t i, const std::size_t new_label);
            double move_cost(const std::size_t i, const std::size_t new_label) const;
            double perform_1_swaps();
            double compute_2_swap_cost(const std::size_t i, const std::size_t j, const double edge_cost) const;
            double perform_2_swaps();
            std::pair<double, std::array<std::size_t,3>> compute_triangle_swap_cost(std::array<std::size_t,3> nodes, std::array<double,3> edge_costs) const;
            double perform_3_swaps();
            void perform_swaps();
            double perform_joins();

            multicut_node_labeling get_labeling() const;
            std::size_t empty_label() const;

        private:
            graph<double> g_;
            std::vector<double> total_edge_sum_;
            std::vector<double> same_component_edge_sum_;
            multicut_node_labeling labeling_;
            const multicut_instance& instance_;
            double lower_bound_;
            std::size_t empty_label_ = 0;
    }; 

}
