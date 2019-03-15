#include "graph_matching/matching_problem_input.h"
#include <eigen3/Eigen/Eigen>
#include <array>
#include "MCF-SSP/mcf_ssp.hxx"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    class graph_matching_frank_wolfe {
        public:
            graph_matching_frank_wolfe(const graph_matching_input& gm, const graph_matching_input::labeling& l = {});

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
            std::size_t linear_coeff(const std::size_t i, const std::size_t j) const { assert(i<no_left_nodes()); assert(j<no_right_nodes()); return i*no_right_nodes() + j; }

            template<typename MATRIX>
                void update_costs(const MATRIX& m);
            Eigen::SparseVector<double> read_solution() const;

            std::array<std::size_t,2> quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2);

            // objective: M.vec().transpose() * Q * M.vec() + <L, M>
            const std::size_t no_left_nodes_, no_right_nodes_;
            QUADRATIC_COST_TYPE Q; // quadratic costs
            ASSIGNMENT_VECTOR_TYPE L; // linear costs
            ASSIGNMENT_VECTOR_TYPE M; // current matching
            MCF::SSP<int, double> mcf; 
    };

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::graph_matching_frank_wolfe(const graph_matching_input& gm, const graph_matching_input::labeling& labeling)
        :
            no_left_nodes_(gm.no_left_nodes),
            no_right_nodes_(gm.no_right_nodes),
            Q(gm.no_left_nodes*gm.no_right_nodes, gm.no_left_nodes*gm.no_right_nodes),
            L(gm.no_left_nodes*gm.no_right_nodes),
            M(gm.no_left_nodes*gm.no_right_nodes)
    {
        Q.setZero();
        L.setZero();
        M.setZero();

        // TODO: reserve size for Q and L is they are sparse
        L.reserve(gm.assignments.size());
        Q.reserve(gm.quadratic_terms.size());

        for(const auto a : gm.assignments) {
            if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>) {
                assert(L(linear_coeff(a.left_node, a.right_node)) == 0.0); 
                L(linear_coeff(a.left_node, a.right_node)) = a.cost; 
            } else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>) {
                assert(L.coeff(linear_coeff(a.left_node, a.right_node)) == 0.0); 
                L.insert(linear_coeff(a.left_node, a.right_node)) = a.cost; 
            } else
                static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
        }

        for(const auto& q : gm.quadratic_terms) {
            const std::size_t l1 = gm.assignments[q.assignment_1].left_node;
            const std::size_t l2 = gm.assignments[q.assignment_2].left_node;
            const std::size_t r1 = gm.assignments[q.assignment_1].right_node;
            const std::size_t r2 = gm.assignments[q.assignment_2].right_node;

            const auto q_idx = quadratic_indices(l1,l2,r1,r2);
            if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::MatrixXd>) {
                assert(Q(q_idx[0], q_idx[1]) == 0.0);
                Q(q_idx[0], q_idx[1]) = q.cost;
            } else if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>) {
                assert(Q.coeff(q_idx[0], q_idx[1]) == 0.0);
                Q.insert(q_idx[0], q_idx[1]) = q.cost;
            } else
                static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
        } 

        if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::MatrixXd>)
            Q = (0.5*(Q + Q.transpose())).eval();
        else if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>)
            Q = (0.5*(Q + Eigen::SparseMatrix<double>(Q.transpose()))).eval(); 
        else
            static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");

        assert(labeling.size() <= M.rows());
        for(std::size_t l=0; l<labeling.size(); ++l) {
            if(labeling[l] != std::numeric_limits<std::size_t>::max()) {
                assert(labeling[l] < no_right_nodes());
                if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>)
                    M(linear_coeff(l, labeling[l])) = 1.0;
                else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>)
                    M.insert(linear_coeff(l, labeling[l])) = 1.0;
                else
                    static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
            }
        }

        assert(feasible());
        if(labeling.size() > 0) {
            assert(std::abs(gm.evaluate(labeling) - evaluate()) <= tolerance);
        } else {
            assert(std::abs(evaluate()) <= tolerance);
        }

        gm.initialize_mcf(mcf);

        std::cout << "initialized problem\n";

        if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseMatrix<double>>)
            L.makeCompressed();

        if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>)
            Q.makeCompressed();
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::size_t graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::no_left_nodes() const
    {
        return no_left_nodes_;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::size_t graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::no_right_nodes() const
    {
        return no_right_nodes_;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::array<std::size_t,2> graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2)
    {
        assert(l1 < no_left_nodes() && l2 < no_left_nodes());
        assert(r1 < no_right_nodes() && r2 < no_right_nodes());
        return {l1*no_right_nodes() + r1, l2*no_right_nodes() + r2};
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        bool graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::feasible(const MATRIX& m) const
    {
        return false;
        /*
        if((m.array() < -tolerance).any())
            return false;

        const auto r_sum = m.rowwise().sum();
        if( (r_sum.array() > 1.0 + tolerance).any() )
            return false;

        const auto c_sum = m.colwise().sum();
        if( (c_sum.array() > 1.0 + tolerance).any() )
            return false;

        return true;
        */
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        double graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::evaluate(const MATRIX& m) const
    {
        assert(feasible());

        const double linear_cost = L.dot(m); // (L.array()*m.array()).sum();

        const double quadratic_cost = m.cwiseProduct(Q*m).sum();
        
        return linear_cost + quadratic_cost;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        void graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::update_costs(const MATRIX& cost_m) 
        {
            mcf.reset_costs();
            //assert(cost_m.rows() == no_left_nodes() && cost_m.cols() == no_right_nodes());
            for(std::size_t i=0; i<no_left_nodes(); ++i) {
                auto e_first = mcf.first_outgoing_arc(i);
                for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                    const std::size_t e = e_first + e_c;
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>)
                        mcf.update_cost(e, cost_m(linear_coeff(i,j)));
                    else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>)
                        mcf.update_cost(e, cost_m.coeff(linear_coeff(i,j)));
                    else
                        static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
                }
            }
        }

    // possibly sparse matrix does not accelerate
    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    Eigen::SparseVector<double> graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::read_solution() const
    { 
        Eigen::SparseVector<double> sol(no_left_nodes()*no_right_nodes());
        sol.reserve(Eigen::VectorXi::Constant(sol.cols(), 1));
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            auto e_first = mcf.first_outgoing_arc(i);
            for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                const std::size_t e = e_first + e_c;
                assert(mcf.flow(e) == 0 || mcf.flow(e) == 1);
                if(mcf.flow(e) == 1) {
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    sol.insert(linear_coeff(i,j)) = 1.0;
                }
            }
        }
        return sol; 
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    bool graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::perform_fw_step(const std::size_t iter)
    {
        assert(feasible());

        // compute gradient
        const ASSIGNMENT_VECTOR_TYPE g(L + 2*Q*M);
        //const Eigen::Map<const Eigen::MatrixXd> g(g_v.data(), M.rows(), M.cols());

        // min_x <x,g> s.t. x in matching polytope
        update_costs(g);
        mcf.solve();
        const auto s = read_solution();
        const ASSIGNMENT_VECTOR_TYPE d(s-M); // update direction
        const double gap = -d.cwiseProduct(g).sum();
        if(iter%20 == 0)
            std::cout << "gap = " << gap << "\n";
        if(gap < 1e-4) {
            std::cout << "gap smaller than threshold, terminate!\n";
            return false;
        }

        // optimal step size
        double linear_term = d.cwiseProduct(L).sum() + 2.0*d.cwiseProduct(Q*M).sum();
        const double quadratic_term = d.cwiseProduct(Q*d).sum();
        const double gamma = [&]() {
            if(quadratic_term > 0.0)
                return std::min(-linear_term / (2.0*quadratic_term), 1.0);
            else
                return 1.0;
        }();

        // diminishing step size
        //const double gamma = 2.0/(iter+3.0);

        M += gamma * d;

        if(iter%20 == 0)
            std::cout << "iteration = " << iter << ", optimal step size = " << gamma << ", quadratic term = " << quadratic_term << ", linear term = " << linear_term << "\n";

        return true;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    Eigen::SparseMatrix<double> graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::round()
    {
        update_costs(ASSIGNMENT_VECTOR_TYPE(-M));
        mcf.solve();
        return read_solution();
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    void graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::solve()
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

    using graph_matching_frank_wolfe_dense = graph_matching_frank_wolfe<Eigen::VectorXd, Eigen::MatrixXd>;
    using graph_matching_frank_wolfe_sparse = graph_matching_frank_wolfe<Eigen::SparseVector<double>, Eigen::SparseMatrix<double>>;

}
