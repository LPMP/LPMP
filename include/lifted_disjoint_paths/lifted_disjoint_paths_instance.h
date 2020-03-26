#pragma once

#include <array>
#include <vector>

namespace LPMP {

    class lifted_disjoint_paths_instance
    {
        public:
            class weighted_edge : public std::array<std::size_t,2> { double weight; };

            std::size_t add_base_edge(const std::size_t i, const std::size_t j, const double weight);
            std::size_t add_lifted_edge(const std::size_t i, const std::size_t j, const double weight);

            template<typename EDGE_LABEL_ITERATOR>
                bool check_feasiblity(EDGE_LABEL_ITERATOR begin, EDGE_LABEL_ITERATOR end) const;
            template<typename EDGE_LABEL_ITERATOR>
                double evaluate(EDGE_LABEL_ITERATOR begin, EDGE_LABEL_ITERATOR end) const;

            const std::vector<weighted_edge>& base_edges() { return base_edges_; }
            const std::vector<weighted_edge>& lifted_edges() { return lifted_edges_; }

        private:
            std::vector<weighted_edge> base_edges_;
            std::vector<weighted_edge> lifted_edges_;

    };

}
