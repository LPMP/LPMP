#include "graph_matching/matching_problem_input.h"
#include <eigen3/Eigen/Eigen>
#include <array>
#include "MCF-SSP/mcf_ssp.hxx"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    class graph_matching_frank_wolfe {
        public:
            graph_matching_frank_wolfe(const graph_matching_input& gm, const graph_matching_input::labeling& l = {});

            std::size_t no_left_nodes() const;
            std::size_t no_right_nodes() const;

            bool feasible() const;
            double evaluate() const;
            bool perform_fw_step(const std::size_t iter);
            void round();
            void solve();
        private:
            template<typename MATRIX>
                void update_costs(const MATRIX& m);
            Eigen::MatrixXd read_solution() const;

            std::array<std::size_t,2> quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2);

            // objective: M.vec().transpose() * Q * M.vec() + <L, M>
            Eigen::MatrixXd Q; // quadratic costs
            Eigen::MatrixXd L; // linear costs
            Eigen::MatrixXd M; // current matching
            MCF::SSP<int, double> mcf; 
    };

    graph_matching_frank_wolfe::graph_matching_frank_wolfe(const graph_matching_input& gm, const graph_matching_input::labeling& labeling)
        :
            Q(Eigen::MatrixXd::Zero(gm.no_left_nodes*gm.no_right_nodes, gm.no_left_nodes*gm.no_right_nodes)),
            L(Eigen::MatrixXd::Zero(gm.no_left_nodes, gm.no_right_nodes)),
            M(Eigen::MatrixXd::Zero(gm.no_left_nodes, gm.no_right_nodes))
    {
        for(const auto a : gm.assignments)
            L(a.left_node, a.right_node) = a.cost; 

        for(const auto& q : gm.quadratic_terms) {
            const std::size_t l1 = gm.assignments[q.assignment_1].left_node;
            const std::size_t l2 = gm.assignments[q.assignment_2].left_node;
            const std::size_t r1 = gm.assignments[q.assignment_1].right_node;
            const std::size_t r2 = gm.assignments[q.assignment_2].right_node;

            const auto q_idx = quadratic_indices(l1,l2,r1,r2);
            assert(Q(q_idx[0], q_idx[1]) == 0.0);
            Q(q_idx[0], q_idx[1]) = q.cost;
        } 

        assert(labeling.size() <= M.rows());
        for(std::size_t l=0; l<labeling.size(); ++l) {
            if(labeling[l] != std::numeric_limits<std::size_t>::max()) {
                assert(labeling[l] < no_right_nodes());
                M(l, labeling[l]) = 1.0;
            }
        }

        assert(feasible());
        if(labeling.size() > 0) {
            assert(std::abs(gm.evaluate(labeling) - evaluate()) <= tolerance);
        } else {
            assert(std::abs(evaluate()) <= tolerance);
        }

        gm.initialize_mcf(mcf);
    }

    std::size_t graph_matching_frank_wolfe::no_left_nodes() const
    {
        return L.rows();
    }
    std::size_t graph_matching_frank_wolfe::no_right_nodes() const
    {
        return L.cols();
    }
    std::array<std::size_t,2> graph_matching_frank_wolfe::quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2)
    {
        assert(l1 < no_left_nodes() && l2 < no_left_nodes());
        assert(r1 < no_right_nodes() && r2 < no_right_nodes());
        return {l1*no_right_nodes() + r1, l2*no_right_nodes() + r2};
    }

    bool graph_matching_frank_wolfe::feasible() const
    {
        if((M.array() < -tolerance).any())
            return false;

        const auto r_sum = M.rowwise().sum();
        if( (r_sum.array() > 1.0 + tolerance).any() )
            return false;

        const auto c_sum = M.colwise().sum();
        if( (c_sum.array() > 1.0 + tolerance).any() )
            return false;

        return true;
    }
    double graph_matching_frank_wolfe::evaluate() const
    {
        assert(feasible());

        const double linear_cost = (L.array()*M.array()).sum();

        const Eigen::Map<const Eigen::VectorXd> M_v(M.data(), M.size());
        const double quadratic_cost = (M_v.transpose() * Q * M_v)(0,0);
        
        return linear_cost + quadratic_cost;
    }

    template<typename MATRIX>
        void graph_matching_frank_wolfe::update_costs(const MATRIX& cost_m) 
        {
            mcf.reset_costs();
            assert(cost_m.rows() == no_left_nodes() && cost_m.cols() == no_right_nodes());
            for(std::size_t i=0; i<no_left_nodes(); ++i) {
                auto e_first = mcf.first_outgoing_arc(i);
                for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                    const std::size_t e = e_first + e_c;
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    mcf.update_cost(e, cost_m(i,j));
                }
            }
        };

    // possibly use sparse matrix here
    Eigen::MatrixXd graph_matching_frank_wolfe::read_solution() const
    { 
        Eigen::MatrixXd sol = Eigen::MatrixXd::Zero(no_left_nodes(), no_right_nodes());
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            auto e_first = mcf.first_outgoing_arc(i);
            for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                const std::size_t e = e_first + e_c;
                assert(mcf.flow(e) == 0 || mcf.flow(e) == 1);
                if(mcf.flow(e) == 1) {
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    sol(i,j) = 1.0;
                }
            }
        }
        return sol; 
    };
    bool graph_matching_frank_wolfe::perform_fw_step(const std::size_t iter)
    {
        assert(feasible());

        // compute gradient
        const Eigen::Map<const Eigen::VectorXd> M_v(M.data(), M.size());
        const Eigen::Map<const Eigen::VectorXd> L_v(L.data(), L.size());
        const Eigen::VectorXd g_v = L_v + Q*M_v + Q.transpose()*M_v;
        const Eigen::Map<const Eigen::MatrixXd> g(g_v.data(), M.rows(), M.cols());

        // min_x <x,g> s.t. x in matching polytope
        update_costs(g);
        mcf.solve();
        const auto s = read_solution();
        const Eigen::MatrixXd d = s-M; // update direction
        const Eigen::Map<const Eigen::VectorXd> d_v(d.data(), d.size()); 
        const double gap = (d.array() * (-g.array())).sum();
        std::cout << "gap = " << gap << "\n";
        if(gap < 1e-4) {
            std::cout << "gap smaller than threshold, terminate!\n";
            return false;
        }

        // optimal step size
        const double linear_term = (d.array()*L.array()).sum() + (d_v.transpose()*Q*M_v)(0,0) + (M_v.transpose()*Q*d_v)(0,0);
        const double quadratic_term = ( d_v.transpose()*Q*d_v )(0,0);
        //const double gamma = std::min(-linear_term / (2.0*quadratic_term), 1.0);
        //assert(gamma >= 0.0);
        //const double gamma = std::min(1.0, gap/10000); // TODO: estimate correctly
        const double gamma = 2.0/(iter+3.0);
        std::cout << "iteration = " << iter << ", optimal step size = " << gamma << ", quadratic term = " << quadratic_term << ", linear term = " << linear_term << "\n";
        M += gamma * d;

        return true;
    }

    void graph_matching_frank_wolfe::round()
    {
        update_costs(-M);
        mcf.solve();
        M = read_solution();
    }

    void graph_matching_frank_wolfe::solve()
    {
        for(std::size_t iter=0; iter<200; ++iter) {
            if(!perform_fw_step(iter))
                break;
            std::cout << "objective = " << evaluate() << "\n";
        }

        round();
        std::cout << "Rounded solution cost = " << evaluate() << "\n";
    }

}
