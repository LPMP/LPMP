#pragma once

#include "graph_matching/matching_problem_input.h"
#include "correlation_clustering_instance.h"
#include <vector>
#include <memory>

namespace LPMP {

    class multigraph_matching_correlation_clustering_transform {
        public:
            multigraph_matching_correlation_clustering_transform();
            multigraph_matching_correlation_clustering_transform(std::shared_ptr<multigraph_matching_input> _mgm);

            void set_multigraph_matching_input(std::shared_ptr<multigraph_matching_input> _mgm);

            const correlation_clustering_instance& get_correlatino_clustering_instance() const { return *cc; }
            const multigraph_matching_input& get_multigraph_matching_instance() const { return *mgm; }

            multigraph_matching_input::labeling transform(const correlation_clustering_edge_labeling& cc_l) const;
            correlation_clustering_edge_labeling transform(const multigraph_matching_input::labeling& mgm_l) const;

        private:
            // functions to map indices from mgm to cc problem and vice versa
            std::size_t no_graphs() const; 
            std::size_t total_no_nodes() const;
            std::size_t no_nodes(const std::size_t i) const;

            std::size_t multigraph_matching_to_correlation_clustering_node(const std::size_t graph_no, const std::size_t node_no) const;
            std::array<std::size_t,2> multigraph_matching_from_correlation_clustering_node(const std::size_t i) const;
            std::size_t compute_no_graphs();
            void compute_no_nodes();
            void compute_node_offsets();

            void transform_multigraph_matching_to_correlation_clustering();

            std::vector<std::size_t> no_nodes_;
            std::vector<std::size_t> graph_node_offsets_;

            std::shared_ptr<correlation_clustering_instance> cc;
            std::shared_ptr<multigraph_matching_input> mgm;
    };

}
