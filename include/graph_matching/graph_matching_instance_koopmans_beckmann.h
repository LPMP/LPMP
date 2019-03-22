#pragma once

#include "graph_matching_instance.h"
#include <eigen3/Eigen/Eigen>

namespace LPMP {

    class graph_matching_instance_koopmans_beckmann {
        public:
            double evaluate(const LAP_labeling& labeling) const;
            template<typename MCF>
                void initialize_mcf(MCF& mcf) const;

            Eigen::MatrixXd L;
            Eigen::MatrixXd A,B;
    };

    inline double graph_matching_instance_koopmans_beckmann::evaluate(const LAP_labeling& labeling) const
    {
        assert(labeling.size() == L.rows());
        double linear_cost = 0.0;
        for(std::size_t i=0; i<L.size(); ++i) {
            assert(labeling[i] < L.rows());
            linear_cost += L(i,labeling[i]);
        }

        double quadratic_cost = 0.0;
        for(std::size_t i=0; i<L.size(); ++i) {
            for(std::size_t j=0; j<L.size(); ++j) {
                quadratic_cost += A(i,labeling[i]) * B(j,labeling[j]);
            }
        }

        return linear_cost + quadratic_cost;
    } 

    template<typename MCF>
        void graph_matching_instance_koopmans_beckmann::initialize_mcf(MCF& mcf) const
        {
            const std::size_t no_nodes = L.cols() + L.rows();
            const std::size_t no_edges = L.size();

            mcf = MCF(no_nodes, no_edges);

            for(std::size_t i=0; i<L.cols(); ++i)
                for(std::size_t j=0; j<L.rows(); ++j)
                    mcf.add_edge(i, L.cols() + j, 0, 1, L(i,j));

            for(std::size_t i=0; i<L.cols(); ++i)
                mcf.add_node_excess(i, 1);

            for(std::size_t j=0; j<L.rows(); ++j)
                mcf.add_node_excess(L.cols() + j, -1);

            mcf.order();
            mcf.solve();
        }

} // namespace LPMP
