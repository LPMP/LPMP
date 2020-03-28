#include <vector>
#include <algorithm>
#include <array>
#include <limits>
#include <eigen3/Eigen/Eigen>
#include "MCF-SSP/mcf_ssp.hxx"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    class multigraph_matching_frank_wolfe_universe_koopman_beckmans {
        public:
            template<typename MATRIX>
                bool feasible(const MATRIX& m) const;
            bool feasible() const { return feasible(M); }

            template<typename MATRIX>
                double evaluate(const MATRIX& m) const;
            double evaluate() const { return evaluate(M); }

            bool perform_fw_step(const std::size_t iter);
            void round();
            void solve();

            std::size_t universe_size() const;
            multigraph_matching_input::labeling get_solution() const; 

        private:
            std::array<std::size_t,2> linear_indices(const std::size_t p, const std::size_t q, const std::size_t i, const std::size_t j) const; // graph 1, graph 2, node in graph 1, node in graph 2

            Eigen::MatrixXd L;
            Eigen::MatrixXd M;
            std::vector<Eigen::MatrixXd> A;
            std::vector<MCF::SSP<int, double>> mcfs;
            multigraph_matching_input::graph_size gs;
            const multigraph_matching_koopman_beckmans_instance& instance_;
    };

    multigraph_matching_frank_wolfe_universe_koopman_beckmans::multigraph_matching_frank_wolfe_universe_koopman_beckmans(const multigraph_matching_koopman_beckmans_instance& instance, const std::size_t universe_size)
        : gs(instance),
        instance_(instance)
    {
        M = Eigen::MatrixXd::Zero(gs.total_no_nodes(), std::min(universe_size,gs.total_no_nodes()));
        L = Eigen::MatrixXd::Zero(gs.total_no_nodes(), gs.total_no_nodes());

        for(const auto& gm : instance) {
            const std::size_t p = gm.left_graph_no;
            const std::size_t q = gm.right_graph_no;
            for(const auto a : gm.gm_input.assignments) {
                const auto [l1, l2] = linear_indices(p, q, a.left_node, a.right_node);
                assert(L(l1,l2) == 0.0);
                L(l1,l2) = a.cost;
            }
        }

        mcfs.reserve(instance.size());
        for(const auto& gm : instance) {
            mcfs.push_back(MCF::SSP<int,double>());
            gm.gm_input.initialize_mcf(mcfs.back());
        }
    }

    std::array<std::size_t,2> multigraph_matching_frank_wolfe_universe::linear_indices(const std::size_t p, const std::size_t q, const std::size_t l, const std::size_t r) const
    {
        assert(p < gs.no_graphs() && q < gs.no_graphs());
        assert(l < gs.no_nodes(p) && r < gs.no_nodes(q));
        return {gs.node_no(p, l), gs.node_no(q,r)}; 
    }

    std::size_t multigraph_matching_frank_wolfe_universe::universe_size() const
    {
        return M.cols();
    }

    template<typename MATRIX>
    bool multigraph_matching_frank_wolfe_universe::feasible(const MATRIX& m) const
    {
        assert(m.cols() == M.cols() && m.rows() == M.rows());

        if((m.array() < -tolerance).any())
            return false;

        const auto r_sum = m.rowwise().sum();
        if( (r_sum.array() > 1.0 + tolerance).any() )
            return false;

        for(std::size_t p=0; p<gs.no_graphs(); ++p) {
            const auto p_slice = get_assignment_slice(p);
            const auto c_sum = p_slice.colwise().sum();
            if( (c_sum.array() > 1.0 + tolerance).any() )
                return false; 
        }

        return true;
    }

    template<typename MATRIX>
    double multigraph_matching_frank_wolfe_universe_koopman_beckmans::evaluate(const MATRIX& m) const
    {
        assert(m.cols() == M.cols() && m.rows() == M.rows());
        assert(feasible(m));

        const double linear_cost = (L.array() * (m*m.transpose()).array()).sum();

        double quadratic_cost = 0.0;
        for(std::size_t i=0; i<gs.no_graphs(); ++i) {
            const auto m_i = get_assignment_slice(m, i);
            for(std::size_t j=0; j<gs.no_graphs(); ++j) { // TODO: currently we could sum from i+1
                if(i == j)
                    continue;
                const auto m_j = get_assignment_slice(m, j);
                const Eigen::MatrixXd m_ij = m_i * m_j.transpose();
                const Eigen::Map<const Eigen::VectorXd> m_ij_v(m_ij.data(), m_ij.size());
                quadratic_cost += 0.0;
            }
        }

        return linear_cost + quadratic_cost; 
    }

}
