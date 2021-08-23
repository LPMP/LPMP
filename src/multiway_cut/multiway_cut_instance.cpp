#include "multiway_cut/multiway_cut_instance.h"

namespace LPMP {

    double multiway_cut_instance::evaluate(const multiway_cut_labeling& labeling) const
    {
        double mc_cost = 0.0;
        for(const auto e : edge_costs.edges())
            if(labeling[e[0]] != labeling[e[1]])
                mc_cost += e.cost; 

        double node_cost = 0.0;
        for(size_t i=0; i<nr_nodes(); ++i)
            node_cost += node_costs(i, labeling[i]);

        return mc_cost + node_cost; 
    }

    bool multiway_cut_instance::feasible(const multiway_cut_labeling& labeling) const
    {
        if(labeling.size() != nr_nodes())
            return false;

        // see whether nodes have admissible labels 
        for(size_t i=0; i<nr_nodes(); ++i)
            if(labeling[i] >= nr_labels())
                return false;
        
        return true;
    }

}

