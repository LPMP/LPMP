#pragma once

#include "graph_matching/graph_matching_instance.h"
#include "union_find.hxx"
#include <array>
#include <algorithm>

namespace LPMP {

    class multigraph_matching_instance;
    class multigraph_matching_labeling;

    class graph_size {
        public:
            template<typename SIZE_ITERATOR>
                graph_size(SIZE_ITERATOR size_begin, SIZE_ITERATOR size_end)
                : no_nodes_(size_begin, size_end)
                {
                    compute_node_offsets(); 
                }
            graph_size(const multigraph_matching_instance& instance);
            graph_size(const multigraph_matching_labeling& labeling);
            graph_size() {}

            std::size_t no_graphs() const;
            std::size_t total_no_nodes() const;
            std::size_t no_nodes(const std::size_t i) const;

            std::vector<std::size_t>& no_nodes();
            const std::vector<std::size_t>& no_nodes() const;

            std::size_t node_no(const std::size_t graph_no, const std::size_t node_no) const;
            std::array<std::size_t,2> graph_node_no(const std::size_t i) const;

        private:
            std::size_t compute_no_graphs(const multigraph_matching_instance& instance);
            void compute_no_nodes(const multigraph_matching_instance& instance);
            void compute_no_nodes(const multigraph_matching_labeling& labeling);
            void compute_node_offsets();

            std::vector<std::size_t> no_nodes_;
            std::vector<std::size_t> graph_node_offsets_;
    };
    class node_mapping;

    struct multigraph_matching_labeling_entry {
        std::size_t left_graph_no, right_graph_no;
        LAP_labeling labeling;
    }; 

    class multigraph_matching_labeling : public std::vector<multigraph_matching_labeling_entry> {
        public:
            using std::vector<multigraph_matching_labeling_entry>::vector;

            template<typename MATRIX>
                multigraph_matching_labeling(const MATRIX& m, const graph_size& mgm_size)
                {
                    assert(m.cols() == m.rows());
                    //for(std::size_t i=0; i<m.cols(); ++i) {
                    //   assert(m(i,i) == 1);
                    //}
                    for(std::size_t p=0; p<mgm_size.no_graphs(); ++p) {
                        for(std::size_t q=p+1; q<mgm_size.no_graphs(); ++q) {
                            // extract submatrix
                            const std::size_t first_col = mgm_size.node_no(q,0);
                            const std::size_t no_cols = mgm_size.no_nodes(q);
                            const std::size_t first_row = mgm_size.node_no(p,0);
                            const std::size_t no_rows = mgm_size.no_nodes(p);

                            auto gm_matching = m.block(first_row, first_col, no_rows, no_cols);
                            //auto gm_matching = m.block(first_col, first_row, no_cols, no_rows);
                            LAP_labeling gm_labeling(gm_matching);
                            this->push_back({p,q, gm_labeling});
                        }
                    }
                }

            bool check_primal_consistency() const;

            template<typename STREAM>
                void write_primal_matching(STREAM& s) const
                {
                    for(const auto& gm : *this) {
                        s << "graph matching " << gm.left_graph_no << " -> " << gm.right_graph_no << "\n";
                        gm.labeling.WritePrimal(s);
                    }
                }

            // write our primal matching as table with three columns: graph number, node number, cluster identifier
            template<typename STREAM>
                void write_primal_clustering(STREAM& s) const
                {
                    graph_size nm(*this);
                    union_find uf(nm.total_no_nodes());

                    s << "# ${graph_no}, ${node_no}, ${cluster_id}\n";
                    for(const auto& gm : *this)
                        for(std::size_t i=0; i<gm.labeling.size(); ++i)
                            if(gm.labeling[i] != std::numeric_limits<std::size_t>::max())
                                uf.merge( nm.node_no(gm.left_graph_no, i), nm.node_no(gm.right_graph_no, gm.labeling[i]) );

                    for(std::size_t i=0; i<nm.no_nodes().size(); ++i) {
                        const auto [graph_no, node_no] = nm.graph_node_no(i);
                        s << graph_no << ", " << node_no << ", " << uf.find(i) << "\n";
                    }
                }
    };

    struct multigraph_matching_instance_entry {
        std::size_t left_graph_no, right_graph_no;
        graph_matching_instance gm_instance;

    };

    class multigraph_matching_labeling;

    class multigraph_matching_instance : public std::vector<multigraph_matching_instance_entry> {

        double evaluate(const multigraph_matching_labeling& l) const;
    };



} // namespace LPMP
