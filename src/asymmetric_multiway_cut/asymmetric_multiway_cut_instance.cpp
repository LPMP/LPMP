#include "asymmetric_multiway_cut/asymmetric_multiway_cut_instance.h"

namespace LPMP {

    double asymmetric_multiway_cut_instance::evaluate(const asymmetric_multiway_cut_labeling& labeling) const
    {
        const double mc_cost = edge_costs.evaluate(labeling.edge_labels);

        double node_cost = 0.0;
        for(size_t i=0; i<nr_nodes(); ++i)
            node_cost += node_costs(i, labeling.node_labels[i]);

        return mc_cost + node_cost; 
    }

    bool asymmetric_multiway_cut_instance::feasible(const asymmetric_multiway_cut_labeling& labeling) const
    {
        if(labeling.edge_labels.size() != nr_edges())
            return false;
        if(labeling.node_labels.size() != nr_nodes())
            return false;

        if(edge_costs.feasible(labeling.edge_labels) == false)
            return false;

        // see whether nodes have admissible labels 
        for(size_t i=0; i<nr_nodes(); ++i)
            if(labeling.node_labels[i] >= nr_labels())
                return false;
        
        // see whether node labels are compatible with multicut
        for(size_t e=0; e<edge_costs.no_edges(); ++e)
        {
            const size_t i = edge_costs.edges()[e][0];
            const size_t j = edge_costs.edges()[e][1];

            if(!labeling.edge_labels[e] && labeling.node_labels[i] != labeling.node_labels[j])
                return false;
        }

        return true;
    }

}
