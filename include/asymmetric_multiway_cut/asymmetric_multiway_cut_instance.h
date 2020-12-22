#pragma once

#include "multicut/multicut_instance.h"

namespace LPMP {

    class asymmetric_multiway_cut_node_costs
    {
        private:
            size_t nr_labels_ = 0;
            std::vector<double> costs;

        public:
            size_t nr_labels() const { return nr_labels_; }
            size_t nr_nodes() const { return costs.size()/nr_labels_; }
            double operator()(const size_t var, const size_t label) const; 
            template<typename ITERATOR>
                void push_back(ITERATOR cost_begin, ITERATOR cost_end);
    };

    inline double asymmetric_multiway_cut_node_costs::operator()(const size_t var, const size_t label) const
    {
        assert(var < nr_nodes());
        assert(label < nr_labels());
        return costs[var*nr_labels() + label];

    }

    class asymmetric_multiway_cut_labeling {
        public:
            std::vector<size_t> node_labels;
            std::vector<size_t> node_connected_components_ids;
            multicut_edge_labeling edge_labels; 
    };

    class asymmetric_multiway_cut_instance {
        public:
            asymmetric_multiway_cut_node_costs node_costs;
            multicut_instance edge_costs;
            size_t nr_nodes() const { return node_costs.nr_nodes(); }
            size_t nr_labels() const { return node_costs.nr_labels(); }
            size_t nr_edges() const { return edge_costs.no_edges(); }
            double evaluate(const asymmetric_multiway_cut_labeling& labeling) const;
            bool feasible(const asymmetric_multiway_cut_labeling& labeling) const;
    };

    ////////////////////
    // implementation //
    ////////////////////

    template<typename ITERATOR>
        void asymmetric_multiway_cut_node_costs::push_back(ITERATOR cost_begin, ITERATOR cost_end)
        {
            if(nr_labels_ == 0)
                nr_labels_ = std::distance(cost_begin, cost_end);

            assert(std::distance(cost_begin, cost_end) > 0);
            assert(std::distance(cost_begin, cost_end) == nr_labels_);

            for(auto it=cost_begin; it!=cost_end; ++it)
                costs.push_back(*it); 
        }
}
