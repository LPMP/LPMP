#pragma once

#include <vector>
#include <array>
#include <cassert>
#include <unordered_map>
#include <unordered_set>

namespace LPMP {

    // utility functinos for detecting and extracting subcycles found by cycle search
    inline bool has_subcycles(const std::vector<std::size_t>& cycle)
    {
        std::unordered_set<std::size_t> s;
        for(std::size_t c=0; c<cycle.size()-1; ++c) {
            if(s.count(cycle[c]) > 0)
                return true;
            s.insert(cycle[c]);
        }
        return false; 
    }

    // TODO: find more than one subcycle
    inline std::vector<std::size_t> find_subcycle(std::vector<std::size_t>& cycle)
    {
        assert(cycle.size() >= 4 && cycle[0] == cycle.back());
        assert(cycle.size() >= 4);
        assert(cycle[0] == cycle.back());

        std::unordered_map<std::size_t, std::size_t> first_occurence;
        first_occurence.reserve(cycle.size()-1);

        for(std::size_t c=0; c<cycle.size(); ++c) {
            if(first_occurence.count(cycle[c]) > 0) {
                const std::size_t p = first_occurence[cycle[c]];
                std::vector<std::size_t> subcycle;
                subcycle.reserve(c+1-p);
                std::copy(cycle.begin() + p, cycle.begin() + c+1, std::back_inserter(subcycle));
                assert(subcycle.size() >= 4 && subcycle[0] == subcycle.back());
                return subcycle;
            } else {
                first_occurence.insert({cycle[c], c});
            }
        }
        return cycle;
    }

    class cycle_packing {
        public:
            std::size_t no_cycles() const;
            std::array<std::vector<std::size_t>::const_iterator,2> get_cycle(const std::size_t c) const;
            double get_cycle_weight(const std::size_t c) const;

            template<typename ITERATOR>
                void add_cycle(ITERATOR cycle_begin, ITERATOR cycle_end, double cycle_weight);

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

        private:
            std::vector<std::size_t> cycle_nodes;
            struct odd_wheel_item {
                double weight;
                std::size_t center_node;
                std::size_t cycle_begin_index;
            };
            std::vector<odd_wheel_item> odd_wheel_items = {{std::numeric_limits<double>::infinity(), std::numeric_limits<std::size_t>::max(), 0}};
    };

    class odd_bicycle_wheel_packing {
        public:
            std::size_t no_odd_bicycle_wheels() const;
            std::array<std::size_t,2> get_axle(const std::size_t c) const;
            std::array<std::vector<std::size_t>::const_iterator,2> get_cycle(const std::size_t c) const;
            double get_odd_bicycle_wheel_weight(const std::size_t c) const;

            template<typename ITERATOR>
                void add_odd_bicycle_wheel(std::array<std::size_t,2> center_nodes, ITERATOR cycle_begin, ITERATOR cycle_end, const double weight);

        private:
            std::vector<std::size_t> cycle_nodes;
            struct odd_bicycle_wheel_item {
                double weight;
                std::array<std::size_t,2> center_nodes;
                std::size_t cycle_begin_index;
            }; 
            std::vector<odd_bicycle_wheel_item> odd_bicycle_wheel_items = {{std::numeric_limits<double>::infinity(), std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max(), 0}};
    };

    // implementation
    
    // cycle packing 
    inline std::size_t cycle_packing::no_cycles() const 
    { 
        assert(cycle_weights.size()+1 == cycle_begin_index.size());
        return cycle_weights.size(); 
    }

    inline std::array<std::vector<std::size_t>::const_iterator,2> cycle_packing::get_cycle(const std::size_t c) const 
    {
        assert(c < cycle_begin_index.size()-1);
        assert(cycle_weights.size()+1 == cycle_begin_index.size());
        return {cycle_nodes.begin() + cycle_begin_index[c], cycle_nodes.begin() + cycle_begin_index[c+1]};
    }

    inline double cycle_packing::get_cycle_weight(const std::size_t c) const 
    {
        assert(c < cycle_weights.size());
        assert(cycle_weights.size()+1 == cycle_begin_index.size());
        return cycle_weights[c];
    }

    template<typename ITERATOR>
        void cycle_packing::add_cycle(ITERATOR cycle_begin, ITERATOR cycle_end, double cycle_weight)
        {
            assert(std::distance(cycle_begin, cycle_end) >= 3);
            assert(cycle_weights.size()+1 == cycle_begin_index.size());
            cycle_nodes.insert(cycle_nodes.end(), cycle_begin, cycle_end);
            cycle_begin_index.push_back(cycle_nodes.size());
            cycle_weights.push_back(cycle_weight);
        }


    // odd wheel packing 
    inline std::size_t odd_wheel_packing::no_odd_wheels() const
    {
        return odd_wheel_items.size()-1;
    }

    inline std::size_t odd_wheel_packing::get_center_node(const std::size_t c) const
    {
        assert(c < no_odd_wheels());
        return odd_wheel_items[c].center_node; 
    }

    inline std::array<std::vector<std::size_t>::const_iterator,2> odd_wheel_packing::get_cycle(const std::size_t c) const 
    {
        assert(c < no_odd_wheels());
        return {cycle_nodes.begin() + odd_wheel_items[c].cycle_begin_index, cycle_nodes.begin() + odd_wheel_items[c+1].cycle_begin_index};
    } 

    inline double odd_wheel_packing::get_odd_wheel_weight(const std::size_t c) const 
    {
        assert(c < no_odd_wheels());
        return odd_wheel_items[c].weight;
    }

    template<typename ITERATOR>
        void odd_wheel_packing::add_odd_wheel(std::size_t center_node, ITERATOR cycle_begin, ITERATOR cycle_end, const double weight)
        {
            assert(std::distance(cycle_begin, cycle_end) >= 3);
            odd_wheel_items.back().center_node = center_node;
            odd_wheel_items.back().weight = weight;
            assert(odd_wheel_items.back().cycle_begin_index == cycle_nodes.size());
            cycle_nodes.insert(cycle_nodes.end(), cycle_begin, cycle_end);
            odd_wheel_items.push_back({std::numeric_limits<double>::infinity(), std::numeric_limits<std::size_t>::max(), cycle_nodes.size()});
        }

    // odd bicycle wheel packing 
    inline std::size_t odd_bicycle_wheel_packing::no_odd_bicycle_wheels() const
    {
        return odd_bicycle_wheel_items.size()-1;
    }

    inline std::array<std::size_t,2> odd_bicycle_wheel_packing::get_axle(const std::size_t c) const
    {
        assert(c < no_odd_bicycle_wheels());
        return odd_bicycle_wheel_items[c].center_nodes; 
    }

    inline std::array<std::vector<std::size_t>::const_iterator,2> odd_bicycle_wheel_packing::get_cycle(const std::size_t c) const 
    {
        assert(c < no_odd_bicycle_wheels());
        return {cycle_nodes.begin() + odd_bicycle_wheel_items[c].cycle_begin_index, cycle_nodes.begin() + odd_bicycle_wheel_items[c+1].cycle_begin_index};
    } 

    inline double odd_bicycle_wheel_packing::get_odd_bicycle_wheel_weight(const std::size_t c) const 
    {
        assert(c < no_odd_bicycle_wheels());
        return odd_bicycle_wheel_items[c].weight;
    }

    template<typename ITERATOR>
        void odd_bicycle_wheel_packing::add_odd_bicycle_wheel(std::array<std::size_t,2> axle, ITERATOR cycle_begin, ITERATOR cycle_end, const double weight)
        {
            assert(std::distance(cycle_begin, cycle_end) >= 3);
            odd_bicycle_wheel_items.back().center_nodes = axle;
            odd_bicycle_wheel_items.back().weight = weight;
            assert(odd_bicycle_wheel_items.back().cycle_begin_index == cycle_nodes.size());
            cycle_nodes.insert(cycle_nodes.end(), cycle_begin, cycle_end);
            odd_bicycle_wheel_items.push_back({std::numeric_limits<double>::infinity(), std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max(), cycle_nodes.size()});
        }
    
} // namespace LPMP
