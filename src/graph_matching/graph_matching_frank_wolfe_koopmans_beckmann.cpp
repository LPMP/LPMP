#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching_instance_koopmans_beckmann.h"
#include "polynomial_functions.h"
#include <eigen3/Eigen/Eigen>
#include <array>
#include "MCF-SSP/mcf_ssp.hxx"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    class graph_matching_frank_wolfe_koopmans_beckmann {
        public:
            graph_matching_frank_wolfe_koopmans_beckmann(const graph_matching_instance_koopmans_beckmann& gm, const LAP_labeling& l = {});

            std::size_t no_left_nodes() const;
            std::size_t no_right_nodes() const;

            template<typename MATRIX>
                bool feasible(const MATRIX& m) const;
            bool feasible() const { return feasible(M); }
            template<typename MATRIX>
                double evaluate(const MATRIX& m) const;
            double evaluate() const { return evaluate(M); }
            bool perform_fw_step(const std::size_t iter);
            Eigen::SparseMatrix<double> round();
            void solve();

        private:
            template<typename MATRIX>
                void update_costs(const MATRIX& m);
            Eigen::SparseMatrix<double> read_solution() const;

            std::array<std::size_t,2> quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2);

            // objective: trace(A*M*B*M) + <L, M>
            Eigen::MatrixXd A,B; // quadratic costs
            Eigen::MatrixXd L; // linear costs
            Eigen::MatrixXd M; // current matching
            MCF::SSP<int, double> mcf; 
    };

    graph_matching_frank_wolfe_koopmans_beckmann::graph_matching_frank_wolfe_koopmans_beckmann(const graph_matching_instance_koopmans_beckmann& gm, const LAP_labeling& labeling)
        :
            L(gm.L),
            A(gm.A),
            B(gm.B)
    {
        M = Eigen::MatrixXd(A.cols(), A.rows());
        M.setOnes();
        M *= 1.0/M.cols();

        assert(feasible());
        if(labeling.size() > 0) {
            assert(std::abs(gm.evaluate(labeling) - evaluate()) <= tolerance);
        }

        gm.initialize_mcf(mcf);

        std::cout << "initialized problem\n";
    }

    std::size_t graph_matching_frank_wolfe_koopmans_beckmann::no_left_nodes() const
    {
        return L.rows();
    }
    std::size_t graph_matching_frank_wolfe_koopmans_beckmann::no_right_nodes() const
    {
        return L.cols();
    }

    std::array<std::size_t,2> graph_matching_frank_wolfe_koopmans_beckmann::quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2)
    {
        assert(l1 < no_left_nodes() && l2 < no_left_nodes());
        assert(r1 < no_right_nodes() && r2 < no_right_nodes());
        return {l1*no_right_nodes() + r1, l2*no_right_nodes() + r2};
    }

    template<typename MATRIX>
        bool graph_matching_frank_wolfe_koopmans_beckmann::feasible(const MATRIX& m) const
    {
        if((m.array() < -tolerance).any())
            return false;

        const auto r_sum = m.rowwise().sum();
        if( (r_sum.array() > 1.0 + tolerance).any() )
            return false;
        if( (r_sum.array() < 1.0 - tolerance).any() )
            return false;

        const auto c_sum = m.colwise().sum();
        if( (c_sum.array() > 1.0 + tolerance).any() )
            return false;
        if( (c_sum.array() < 1.0 - tolerance).any() )
            return false;

        return true;
    }

    template<typename MATRIX>
        double graph_matching_frank_wolfe_koopmans_beckmann::evaluate(const MATRIX& m) const
    {
        assert(feasible());

        const double linear_cost = (L.array()*m.array()).sum();

        const double quadratic_cost = (A*M*B*M).trace();
        
        return linear_cost + quadratic_cost;
    }

    template<typename MATRIX>
        void graph_matching_frank_wolfe_koopmans_beckmann::update_costs(const MATRIX& cost_m) 
        {
            mcf.reset_costs();
            //assert(cost_m.rows() == no_left_nodes() && cost_m.cols() == no_right_nodes());
            for(std::size_t i=0; i<no_left_nodes(); ++i) {
                auto e_first = mcf.first_outgoing_arc(i);
                for(std::size_t e_c=0; e_c<mcf.no_outgoing_arcs(i); ++e_c) {
                    const std::size_t e = e_first + e_c;
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    mcf.update_cost(e, cost_m(i,j));
                }
            }
        }

    // possibly sparse matrix does not accelerate
    Eigen::SparseMatrix<double> graph_matching_frank_wolfe_koopmans_beckmann::read_solution() const
    { 
        Eigen::SparseMatrix<double> sol(no_left_nodes(),no_right_nodes());
        sol.reserve(Eigen::VectorXi::Constant(sol.cols(), 1));
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            auto e_first = mcf.first_outgoing_arc(i);
            for(std::size_t e_c=0; e_c<mcf.no_outgoing_arcs(i); ++e_c) {
                const std::size_t e = e_first + e_c;
                assert(e < mcf.no_edges());
                assert(mcf.flow(e) == 0 || mcf.flow(e) == 1);
                if(mcf.flow(e) == 1) {
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    sol.insert(i,j) = 1.0;
                }
            }
        }
        assert(feasible(Eigen::MatrixXd(sol)));
        return sol; 
    }

    bool graph_matching_frank_wolfe_koopmans_beckmann::perform_fw_step(const std::size_t iter)
    {
        assert(feasible());

        // compute gradient
        const Eigen::MatrixXd g = L + (B*M*A + A*M*B).transpose();

        // min_x <x,g> s.t. x in matching polytope
        update_costs(g);
        mcf.solve();
        const auto s = read_solution();
        const Eigen::MatrixXd d(s-M); // update direction
        const double gap = -(d.array()*g.array()).sum();
        if(iter%20 == 0)
            std::cout << "gap = " << gap << "\n";
        if(gap < 1e-4) {
            std::cout << "gap smaller than threshold, terminate!\n";
            return false;
        }

        // optimal step size
        double linear_term = (L.array()*d.array()).sum() + (A*M*B*d).trace() + (A*d*B*M).trace();
        const double quadratic_term = (A*d*B*d).trace();
        const double gamma = std::min(quadratic_function(quadratic_term, linear_term, 0.0).minimum(), 1.0);

        // diminishing step size
        //const double gamma = 2.0/(iter+3.0);

        M += gamma * d;

        if(iter%20 == 0)
            std::cout << "iteration = " << iter << ", optimal step size = " << gamma << ", quadratic term = " << quadratic_term << ", linear term = " << linear_term << "\n";

        assert(gamma >= 0.0);
        return true;
    }

    Eigen::SparseMatrix<double> graph_matching_frank_wolfe_koopmans_beckmann::round()
    {
        update_costs(-M);
        mcf.solve();
        return read_solution();
    }

    void graph_matching_frank_wolfe_koopmans_beckmann::solve()
    {
        for(std::size_t iter=0; iter<400; ++iter) {
            if(!perform_fw_step(iter))
                break;
            if(iter%20 == 0)
                std::cout << "objective = " << evaluate() << "\n";
        }

        M = round();
        std::cout << "Rounded solution cost = " << evaluate() << "\n";
    }
}
