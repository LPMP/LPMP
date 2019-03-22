#pragma once

#include <vector>
#include <limits>
#include <cassert>

namespace LPMP {

    class LAP_labeling; 

    class LAP_instance {
        public:

            void add_assignment(const std::size_t left_node, const std::size_t right_node, const double cost);

            template<typename MCF>
                void initialize_mcf(MCF& mcf, const double scaling = 1.0) const;

            std::size_t no_mcf_nodes() const;
            std::size_t no_mcf_edges() const;
            std::size_t no_left_nodes() const;
            std::size_t no_right_nodes() const; 

            void normalize();
            double evaluate(const LAP_labeling& l) const;

            const auto& assignments() const { return assignments_; }

        private:
            struct assignment {
                std::size_t left_node;
                std::size_t right_node;
                double cost;
            };
            std::vector<assignment> assignments_;

            std::size_t no_left_nodes_ = 0;
            std::size_t no_right_nodes_ = 0; 
    };

    class LAP_labeling : public std::vector<std::size_t> {
        public:
            using std::vector<std::size_t>::vector;

            template<typename MATRIX>
                LAP_labeling(const MATRIX& m);

            template<typename STREAM>
                void WritePrimal(STREAM& s) const;

            bool check_primal_consistency() const;

            std::size_t no_left_nodes() const { return this->size(); }

            std::size_t highest_matched_node() const;
    };

    class graph_matching_instance : public LAP_instance {
        public:
            double evaluate(const LAP_labeling& l) const;
            const auto& quadratic_terms() const { return quadratic_terms_; }
        private:
            struct quadratic {
                std::size_t assignment_1, assignment_2;
                double cost;
            };
            std::vector<quadratic> quadratic_terms_; 
    };

    template<typename MCF>
        void LAP_instance::initialize_mcf(MCF& mcf, const double scaling) const
    {
        const std::size_t no_nodes = no_left_nodes() + no_right_nodes();
        const std::size_t no_mcf_nodes = no_left_nodes() + no_right_nodes() + 2;
        const std::size_t no_mcf_edges = assignments().size() + no_left_nodes() + no_right_nodes() + 1;

        mcf = MCF(no_mcf_nodes, no_mcf_edges);

        for(const auto a : assignments())
            mcf.add_edge(a.left_node, no_left_nodes() + a.right_node, 0, 1, scaling*a.cost);

        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            mcf.add_edge(i, no_nodes+1, 0, 1, 0.0); // for non-assignment
            mcf.add_node_excess(i, 1);
        }

        for(std::size_t i=0; i<no_right_nodes(); ++i) {
            mcf.add_edge(no_nodes, no_left_nodes() + i, 0, 1, 0.0); // for non-assignment
            mcf.add_node_excess(no_left_nodes() + i, -1);
        }

        mcf.add_edge(no_nodes, no_nodes + 1, 0, no_nodes, 0.0);
        mcf.add_node_excess(no_nodes, no_right_nodes());
        mcf.add_node_excess(no_nodes+1, -no_left_nodes());

        mcf.order();
        mcf.solve();
    }

    template<typename MATRIX>
        LAP_labeling::LAP_labeling(const MATRIX& m)
        {
            // check that each row and columns has at most one 1 entry
            for(std::size_t i=0; i<m.cols(); ++i) {
                for(std::size_t j=0; j<m.rows(); ++j) {
                    assert(m(i,j) == 0 || m(i,j) == 1);
                }
            }
            for(std::size_t i=0; i<m.cols(); ++i) {
                assert(m.col(i).sum() <= 1);
            }
            for(std::size_t j=0; j<m.rows(); ++j) {
                assert(m.row(j).sum() <= 1);
            }
            this->resize(m.cols(), std::numeric_limits<std::size_t>::max());
            for(std::size_t i=0; i<m.cols(); ++i) {
                for(std::size_t j=0; j<m.rows(); ++j) {
                    if(m(i,j) == 1) {
                        assert((*this)[i] == std::numeric_limits<std::size_t>::max());
                        (*this)[i] = j;
                    } 
                }
            }
            assert(this->size() == m.cols());
        }

    template<typename STREAM>
        void LAP_labeling::WritePrimal(STREAM& s) const
        {
            for(std::size_t i=0; i<this->size(); ++i) {
                if((*this)[i] != std::numeric_limits<std::size_t>::max())
                    s << i << " -> " << (*this)[i] << "\n";
                else
                    s << i << " not matched\n";
            } 
        }

} // namespace LPMP
