#pragma once

#include "cut_base/cut_base_instance.hxx"
#include "correlation_clustering_instance.h"
#include "multicut_factors.h"
#include "union_find.hxx"

namespace LPMP {

    class correlation_clustering_instance; // forward declaration
    class correlation_clustering_edge_labeling; // forward declaration
    class correlation_clustering_node_labeling; // forward declaration

    class multicut_edge_labeling; // forward declaration

    struct multicut_instance : public cut_base_instance {

        template<typename STREAM>
            void write_problem(STREAM& s) const;

        bool feasible(const multicut_edge_labeling& labeling) const;

        correlation_clustering_instance transform_to_correlation_clustering() const;
    };

    class multicut_node_labeling : public cut_base_node_labeling {
        public:
            using cut_base_node_labeling::cut_base_node_labeling;
            multicut_edge_labeling transform_to_edge_labeling(const multicut_instance& instance) const;
    };

    class multicut_edge_labeling : public cut_base_edge_labeling {
        public:
            using cut_base_edge_labeling::cut_base_edge_labeling;
            bool check_primal_consistency(const multicut_instance& input) const;
            multicut_node_labeling transform_to_node_labeling(const multicut_instance& instance) const;
            correlation_clustering_edge_labeling transform_to_correlation_clustering() const;
    };

    struct triplet_multicut_instance : public triplet_cut_base_instance<multicut_triplet_factor> {};

    struct quadruplet_multicut_instance : public quadruplet_cut_base_instance<multicut_quadruplet_factor, triplet_multicut_instance> {};

    // implementation

    template<typename STREAM>
        void multicut_instance::write_problem(STREAM& s) const
        {
            s << "MULTICUT\n";
            cut_base_instance::write_problem(s);
        }

    // TODO: put in separate cpp file
    ////////////////////
    // implementation //
    ////////////////////

    inline bool multicut_instance::feasible(const multicut_edge_labeling& labeling) const
    {
        if(labeling.size() != no_edges())
            return false;

        union_find uf(no_nodes());

        for(size_t e=0; e<no_edges(); ++e)
        {
            const size_t i = this->edges_[e][0];
            const size_t j = this->edges_[e][1];
            if(!labeling[e])
                uf.merge(i,j);
        }

        for(size_t e=0; e<no_edges(); ++e)
        {
            const size_t i = this->edges_[e][0];
            const size_t j = this->edges_[e][1];
            if(labeling[e] && uf.connected(i,j))
                return false;
            if(!labeling[e] && !uf.connected(i,j))
                return false;
        } 

        return true;
    }

} // namespace LPMP
