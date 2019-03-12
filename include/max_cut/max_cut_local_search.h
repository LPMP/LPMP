#include "max_cut/max_cut_instance.hxx"
#include "graph.hxx"
#include <vector> 
#include <array>
#include <tuple>

namespace LPMP {

    class max_cut_local_search {
        public:
            max_cut_local_search(const max_cut_instance& instance, const max_cut_node_labeling& l);
            // compute cost change in objective that results from swapping labels of i
            double swap_1_cost(const std::size_t i) const;
            double swap_2_cost(const std::size_t i, const std::size_t j, const double edge_cost) const;
            std::tuple<double, std::array<std::size_t,3>> best_3_swap(const std::array<std::size_t,3>& nodes, const std::array<double,3>& edge_costs) const;

            void swap(const std::size_t i);

            double perform_1_swaps();
            double perform_2_swaps();
            double perform_3_swaps();
            double perform_swaps();

            max_cut_node_labeling get_labeling() const;

        private:
            graph<double> g;
            max_cut_node_labeling label;
            std::vector<double> cut_values;
            const max_cut_instance& instance_;
            double lower_bound;
    };
}
