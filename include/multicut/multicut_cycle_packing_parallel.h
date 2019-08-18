#pragma once

#include "multicut_instance.h"
#include "cut_base/cut_base_packing.h"
#include "copyable_atomic.h"
#include <vector>

namespace LPMP {
	using atomic_edge_container = std::map<std::array<std::size_t,2>, CopyableAtomic<double>>;
	using multicut_triangle_factor = std::map<std::array<std::size_t,3>, std::array<CopyableAtomic<double>,3>>;
	using edge_to_triangle_map = std::unordered_map<std::array<std::size_t,2>, std::set<std::size_t>>;

	// implementation of the ICP algorithm from Lange et al's ICML18 algorithm.
	void multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads);
	cycle_packing compute_multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads);
	cycle_packing cycle_packing_triangulation_parallel(const multicut_instance& input, const int nr_threads, 
	        multicut_triangle_factor& T, edge_to_triangle_map& M, atomic_edge_container& all_edges);

} // end namespace LPMP
