#pragma once

#include "multicut_instance.h"
#include "cut_base/cut_base_packing.h"
#include "atomic_helper.h"
#include <vector>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"

namespace LPMP {

	struct triangle_item {
        std::array<std::size_t, 3> nodes;
        std::array<CopyableAtomic<double>, 3> weights;
        std::array<std::size_t, 3> edge_indices;
        triangle_item(std::array<std::size_t, 3> n): nodes(n){ 
        	weights = {0.0, 0.0, 0.0};
        	edge_indices = {};
        };
        triangle_item(std::array<std::size_t, 3> n, std::array<double, 3> w, std::array<std::size_t, 3> i): nodes(n), edge_indices(i){ 
        	weights[0].store(w[0]);
        	weights[1].store(w[1]);
        	weights[2].store(w[2]);
        };
        void set_weights(double w1, double w2, double w3){
        	weights[0].store(w1);
        	weights[1].store(w2);
        	weights[2].store(w3);
        }
        void set_edge_indices(std::array<std::size_t, 3> indices){
        	edge_indices[0] = indices[0]; edge_indices[1] = indices[1]; edge_indices[2] = indices[2]; 
        }
    };

    struct edge_item {
        std::array<std::size_t, 2> nodes;
        CopyableAtomic<double> cost;
        std::vector<std::array<std::size_t, 2>> triangle_indices;
        edge_item(std::array<std::size_t, 2> n): nodes(n){ 
        	cost = 0.0;
        	triangle_indices = {};
        };
        edge_item(std::array<std::size_t, 2> n, double c, std::vector<std::array<std::size_t, 2>> i): nodes(n), triangle_indices(i){ 
        	cost.store(c);
        };
    };

    struct edge_t : public std::array<std::size_t,2> { double cost; };

	// implementation of the ICP algorithm from Lange et al's ICML18 algorithm.
	void multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads);
	cycle_packing compute_multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads);
	cycle_packing cycle_packing_triangulation_parallel(const multicut_instance& input, const int nr_threads, 
	        std::vector<triangle_item>& triangle_to_edge, std::vector<edge_item>& edge_to_triangle, std::vector<edge_t>& other_edges);

} // end namespace LPMP
