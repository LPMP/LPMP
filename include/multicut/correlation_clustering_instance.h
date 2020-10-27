#pragma once

#include "multicut_instance.h"
#include <algorithm>
#include <cassert>
#include "cut_base/cut_base_instance.hxx"
#include "union_find.hxx"

namespace LPMP {

    class multicut_instance; // forward declaration
    class multicut_edge_labeling; // forward declaration
    class multicut_node_labeling; // forward declaration
    class correlation_clustering_edge_labeling; // forward declaration

    class correlation_clustering_instance : public cut_base_instance {
        public:
            template<typename STREAM>
                void write_problem(STREAM& s) const;

            multicut_instance transform_to_multicut() const;
    }; 

    class correlation_clustering_node_labeling : public cut_base_node_labeling {
        public:
            using cut_base_node_labeling::cut_base_node_labeling;
            correlation_clustering_edge_labeling transform_to_edge_labeling(const correlation_clustering_instance& instance) const;
    };

    class correlation_clustering_edge_labeling : public cut_base_edge_labeling {
        public:
            using cut_base_edge_labeling::cut_base_edge_labeling;
            bool check_primal_consistency(const correlation_clustering_instance& input) const;
            correlation_clustering_node_labeling transform_to_node_labeling(const correlation_clustering_instance& instance) const;
            multicut_edge_labeling transform_to_multicut() const;
    };

    // implementation

    template<typename STREAM>
        void correlation_clustering_instance::write_problem(STREAM& s) const
        {
            s << "CORRELATION_CLUSTERING\n";
            this->write_problem(s);
        }

} // namespace LPMP
