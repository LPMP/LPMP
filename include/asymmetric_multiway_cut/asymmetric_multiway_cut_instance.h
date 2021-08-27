#pragma once

#include "multicut/multicut_instance.h"

namespace LPMP {

    class asymmetric_multiway_cut_node_costs
    {
        private:
            size_t nr_labels_ = 0;
            std::vector<double> costs;
            std::vector<bool> partitionable_class;

        public:
            size_t nr_labels() const { return nr_labels_; }
            size_t nr_nodes() const { return costs.size()/nr_labels_; }
            bool partitionable(const size_t i) const;
            double operator()(const size_t var, const size_t label) const; 
            template<typename ITERATOR>
                void push_back(ITERATOR cost_begin, ITERATOR cost_end);
            template<typename ITERATOR>
                void partitionable(ITERATOR partitionable_begin, ITERATOR partitionable_end);
    };

    inline double asymmetric_multiway_cut_node_costs::operator()(const size_t var, const size_t label) const
    {
        assert(var < nr_nodes());
        assert(label < nr_labels());
        return costs[var*nr_labels() + label];

    }

    inline bool asymmetric_multiway_cut_node_costs::partitionable(const size_t i) const
    {
        assert(i < nr_labels());
        if(partitionable_class.size() == 0)
            return true;
        assert(partitionable_class.size() == nr_labels());
        return partitionable_class[i]; 
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
            template<typename STREAM>
                void write(STREAM& s) const;
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

    template<typename ITERATOR>
        void asymmetric_multiway_cut_node_costs::partitionable(ITERATOR partitionable_begin, ITERATOR partitionable_end)
        {
            assert(nr_labels_ == 0 || std::distance(partitionable_begin, partitionable_end) == nr_labels_);
            nr_labels_ = std::distance(partitionable_begin, partitionable_end);
            partitionable_class = std::vector<bool>(partitionable_begin, partitionable_end); 
        }

    template<typename STREAM>
        void asymmetric_multiway_cut_instance::write(STREAM& s) const
        {           
            s << "ASYMMETRIC MULTIWAY CUT\n";
            s << "MULTICUT\n";
            for(const auto& e : edge_costs.edges())
                s << e[0] << " " << e[1] << " " << e.cost << "\n";

            s << "NODE COSTS\n";
            for(size_t i=0; i<nr_nodes(); ++i)
            {
                for(size_t l=0; l<nr_labels(); ++l)
                {
                    if(l > 0)
                        s << " ";
                    s << node_costs(i,l);
                }
                s << "\n";
            }
        }
}
