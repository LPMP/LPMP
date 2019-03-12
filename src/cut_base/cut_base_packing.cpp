#include "cut_base/cut_base_packing.h"

namespace LPMP {

    class {
        public:
            std::size_t cycle_packing::no_cycles() const 
            { 
                assert(cycle_weights.size()+1 == cycle_begin_index.size());
                return cycle_weights.size(); 
            }

            auto get_cycle(const std::size_t c) const 
            {
                assert(c < cycle_begin_index.size()-1);
                assert(cycle_weights.size()+1 == cycle_begin_index.size());
                return std::make_pair(cycle_nodes.begin() + cycle_begin_index[c], cycle_nodes.begin() + cycle_begin_index[c+1]); 
            }

            double get_cycle_weight(const std::size_t c) const 
            {
                assert(c < cycle_weights.size());
                assert(cycle_weights.size()+1 == cycle_begin_index.size());
                return cycle_weights[c];
            }

            template<typename ITERATOR>
                void add_cycle(ITERATOR cycle_begin, ITERATOR cycle_end, double cycle_weight)
                {
                    assert(cycle_weights.size()+1 == cycle_begin_index.size());
                    cycle_nodes.insert(cycle_nodes.end(), cycle_begin, cycle_end);
                    cycle_begin_index.push_back(cycle_nodes.size());
                    cycle_weights.push_back(cycle_weight);
                }

        private:
            std::vector<std::size_t> cycle_nodes;
            std::vector<std::size_t> cycle_begin_index = {0};
            std::vector<double> cycle_weights;
    };

    class odd_wheel_packing {
        public:
            std::size_t no_odd_wheels() const;
            std::size_t get_center_node(const std::size_t c) const;
            std::array<std::vector<std::size_t>::const_iterator,2> get_cycle(const std::size_t c) const;
            double get_odd_wheel_weight(const std::size_t c) const;

            template<typename ITERATOR>
                void add_odd_wheel(std::size_t center_node, ITERATOR cycle_begin, ITERATOR cycle_end, const double weight);
    };

}
