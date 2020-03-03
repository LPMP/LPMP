#pragma once

#include <variant>
#include <memory>
#include "graph_matching/matching_problem_input.h"

namespace LPMP {

    struct graph_matching_frank_wolfe_options {
        std::size_t max_iter = 1000;
        double gap = 1e-4;
    };

    class graph_matching_frank_wolfe_dense;
    class graph_matching_frank_wolfe_sparse;

    class graph_matching_frank_wolfe {
        public:
            graph_matching_frank_wolfe(const graph_matching_input& gm, const graph_matching_input::labeling& l = {}, const graph_matching_frank_wolfe_options o = {});
            graph_matching_input::labeling solve(); 
        private:
            std::variant<graph_matching_frank_wolfe_dense*, graph_matching_frank_wolfe_sparse*> solver;
            //std::variant<std::unique_ptr<graph_matching_frank_wolfe_dense>, std::unique_ptr<graph_matching_frank_wolfe_sparse>> solver;
    };
}


/*
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
            ASSIGNMENT_VECTOR_TYPE construct_matching_matrix(const graph_matching_input::labeling& l) const;

            std::size_t no_left_nodes() const;
            std::size_t no_right_nodes() const;

            template<typename MATRIX>
                bool feasible(const MATRIX& m) const;
            bool feasible() const { return feasible(M); }
            template<typename MATRIX>
                double evaluate(const MATRIX& m) const;
            double evaluate() const { return evaluate(M); }
            bool perform_fw_step(const std::size_t iter);
            Eigen::SparseVector<double> round();
            graph_matching_input::labeling get_solution();
            void solve();

        private:
            std::size_t linear_coeff(const std::size_t i, const std::size_t j) const;
            std::array<std::size_t,2> matrix_coeffs(const std::size_t x) const;

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
            double constant_ = 0.0;
    };

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::graph_matching_frank_wolfe(const graph_matching_input& gm, const graph_matching_input::labeling& labeling)
        :
            no_left_nodes_(gm.no_left_nodes),
            no_right_nodes_(gm.no_right_nodes),
            Q((gm.no_left_nodes+1)*(gm.no_right_nodes+1), (gm.no_left_nodes+1)*(gm.no_right_nodes+1)),
            L((gm.no_left_nodes+1)*(gm.no_right_nodes+1)),
            M((gm.no_left_nodes+1)*(gm.no_right_nodes+1)),
            constant_(gm.get_constant())
    {
        Q.setZero();
        L.setZero();
        M.setZero();

        // linear terms
        for(const auto a : gm.assignments) {
            if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>) {
                assert(L(linear_coeff(a.left_node, a.right_node)) == 0.0); 
                L(linear_coeff(a.left_node, a.right_node)) = a.cost; 
            } else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>) {
                // TODO: insert with triplets as for quadratic
                assert(L.coeff(linear_coeff(a.left_node, a.right_node)) == 0.0); 
                L.insert(linear_coeff(a.left_node, a.right_node)) = a.cost; 
            } else
                static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
        }

        // quadratic terms
        if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::MatrixXd>) {
            for(const auto& q : gm.quadratic_terms) {
                const std::size_t l1 = gm.assignments[q.assignment_1].left_node;
                const std::size_t l2 = gm.assignments[q.assignment_2].left_node;
                const std::size_t r1 = gm.assignments[q.assignment_1].right_node;
                const std::size_t r2 = gm.assignments[q.assignment_2].right_node;

                const auto q_idx = quadratic_indices(l1,l2,r1,r2);
                //if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::MatrixXd>) {
                assert(Q(q_idx[0], q_idx[1]) == 0.0);
                Q(q_idx[0], q_idx[1]) = q.cost;
                //} else if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>) {
                //    assert(Q.coeff(q_idx[0], q_idx[1]) == 0.0);
                //    Q.insert(q_idx[0], q_idx[1]) = q.cost;
                //} else
                //    static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
            } 
        } else if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>) {
            std::vector< Eigen::Triplet<double> > quadratic_terms_list;
            quadratic_terms_list.reserve(gm.quadratic_terms.size());
            for(const auto& q : gm.quadratic_terms) {
                const std::size_t l1 = gm.assignments[q.assignment_1].left_node;
                const std::size_t l2 = gm.assignments[q.assignment_2].left_node;
                const std::size_t r1 = gm.assignments[q.assignment_1].right_node;
                const std::size_t r2 = gm.assignments[q.assignment_2].right_node;
                const auto q_idx = quadratic_indices(l1,l2,r1,r2);
                quadratic_terms_list.push_back({q_idx[0], q_idx[1], q.cost});
            }
            Q.setFromTriplets(quadratic_terms_list.begin(), quadratic_terms_list.end());
        } else {
            static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
        } 


        if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::MatrixXd>)
            Q = (0.5*(Q + Q.transpose())).eval();
        else if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>)
            Q = (0.5*(Q + Eigen::SparseMatrix<double>(Q.transpose()))).eval(); 
        else
            static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");

        M = construct_matching_matrix(labeling);

        assert(feasible());
        if(labeling.size() > 0) {
            assert(std::abs(gm.evaluate(labeling) - evaluate(M)) <= tolerance);
        }

        gm.initialize_mcf(mcf);

        std::cout << "initialized problem\n";

        //if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>)
        //    L.makeCompressed();

        if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>)
            Q.makeCompressed();
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
        ASSIGNMENT_VECTOR_TYPE 
        graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::construct_matching_matrix(const graph_matching_input::labeling& labeling) const
        {
            ASSIGNMENT_VECTOR_TYPE M((no_left_nodes()+1)*(no_right_nodes()+1));
            if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>) {
                for(std::size_t i=0; i<no_left_nodes(); ++i)
                    M(linear_coeff(i, graph_matching_input::no_assignment)) = 1.0;
                for(std::size_t i=0; i<no_right_nodes(); ++i)
                    M(linear_coeff(graph_matching_input::no_assignment, i)) = 1.0;
                for(std::size_t l=0; l<labeling.size(); ++l) {
                    if(labeling[l] != graph_matching_input::no_assignment) {
                        assert(labeling[l] < no_right_nodes());
                        M(linear_coeff(l, labeling[l])) = 1.0;
                        M(linear_coeff(l, graph_matching_input::no_assignment)) = 0.0;
                        M(linear_coeff(graph_matching_input::no_assignment, labeling[l])) = 0.0;
                    }
                }
            } else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>) {
                M.reserve(no_left_nodes() + no_right_nodes());
                for(std::size_t i=0; i<no_left_nodes(); ++i)
                    M.insert(linear_coeff(i, graph_matching_input::no_assignment)) = 1.0;
                for(std::size_t i=0; i<no_right_nodes(); ++i)
                    M.insert(linear_coeff(graph_matching_input::no_assignment, i)) = 1.0;
                for(std::size_t l=0; l<labeling.size(); ++l) {
                    if(labeling[l] != graph_matching_input::no_assignment) {
                        assert(labeling[l] < no_right_nodes());
                        M.insert(linear_coeff(l, labeling[l])) = 1.0;
                        M.coeffRef(linear_coeff(l, graph_matching_input::no_assignment)) = 0.0;
                        M.coeffRef(linear_coeff(graph_matching_input::no_assignment, labeling[l])) = 0.0;
                    }
                }
                //std::cout << "M sum = " << M.sum() << "\n";
            } else {
                static_assert("Only Eigen::VectorXd and Eigen::SparseVector<double> allowed as template types"); 
            }

            assert(feasible(M));
            
            return M;
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
        assert(l1 < no_left_nodes() || l1 == graph_matching_input::no_assignment);
        assert(l2 < no_left_nodes() || l2 == graph_matching_input::no_assignment);
        assert(r1 < no_right_nodes() || r1 == graph_matching_input::no_assignment);
        assert(r2 < no_right_nodes() || r2 == graph_matching_input::no_assignment);
        const std::size_t L1 = std::min(l1, no_left_nodes());
        const std::size_t L2 = std::min(l2, no_left_nodes());
        const std::size_t R1 = std::min(r1, no_right_nodes());
        const std::size_t R2 = std::min(r2, no_right_nodes());
        return {L1*(no_right_nodes()+1) + R1, L2*(no_right_nodes()+1) + R2};
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        bool graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::feasible(const MATRIX& m) const
    {
        if constexpr(std::is_same_v<MATRIX,Eigen::VectorXd>) {
            if((m.array() < -tolerance).any())
                return false;

            const auto r_sum = m.rowwise().sum();
            for(std::size_t i=0; i<r_sum.size()-1; ++i)
                if(std::abs(r_sum[i] - 1.0) > tolerance)
                    return false;

            const auto c_sum = m.colwise().sum();
            for(std::size_t i=0; i<c_sum.size()-1; ++i)
                if(std::abs(c_sum[i] - 1.0) > tolerance)
                    return false;

            return true;
        } else if constexpr(std::is_same_v<MATRIX,Eigen::SparseVector<double>>) {
            {
                for(std::size_t i=0; i<no_left_nodes(); ++i) {
                    double marg = 0.0;
                    for(std::size_t j=0; j<no_right_nodes(); ++j) {
                        if(m.coeff(linear_coeff(i,j)) < 0.0)
                            return false;
                        marg += m.coeff(linear_coeff(i,j)); 
                    }
                    marg += m.coeff(linear_coeff(i,graph_matching_input::no_assignment)); 
                    if(std::abs(marg - 1.0) > tolerance)
                        return false;
                }

                for(std::size_t j=0; j<no_right_nodes(); ++j) {
                    double marg = 0.0;
                    for(std::size_t i=0; i<no_left_nodes(); ++i) {
                        marg += m.coeff(linear_coeff(i,j)); 
                    }
                    marg += m.coeff(linear_coeff(graph_matching_input::no_assignment, j)); 
                    if(std::abs(marg - 1.0) > tolerance)
                        return false;
                }
            }

            return true;
        } else {
            std::runtime_error("Matrix type not supported.");
            return false;
        }
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        double graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::evaluate(const MATRIX& m) const
    {
        assert(feasible());

        const double linear_cost = L.dot(m); // (L.array()*m.array()).sum();

        const double quadratic_cost = m.cwiseProduct(Q*m).sum();
        
        return linear_cost + quadratic_cost + constant_;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        void graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::update_costs(const MATRIX& cost_m) 
        {
            mcf.reset_costs();
            //assert(cost_m.rows() == no_left_nodes() && cost_m.cols() == no_right_nodes());
            // set assignment costs
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
            // set non-assignment costs
            for(std::size_t i=0; i<no_left_nodes(); ++i) {
                const std::size_t e = mcf.first_outgoing_arc(i) + mcf.no_outgoing_arcs(i) - 1;
                assert(mcf.head(e) == no_left_nodes() + no_right_nodes() + 1);
                if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>)
                    mcf.update_cost(e, cost_m(linear_coeff(i, graph_matching_input::no_assignment)));
                else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>)
                    mcf.update_cost(e, cost_m.coeff(linear_coeff(i, graph_matching_input::no_assignment)));
                else
                    static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
            }

            for(std::size_t i=0; i<no_right_nodes(); ++i) {
                const std::size_t e = mcf.first_outgoing_arc(no_left_nodes() + i) + mcf.no_outgoing_arcs(no_left_nodes() + i) - 1;
                assert(mcf.head(e) == no_left_nodes() + no_right_nodes());
                // TODO: possibly use negative cost here?
                if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>)
                    mcf.update_cost(e, -cost_m(linear_coeff(graph_matching_input::no_assignment, i)));
                else if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>)
                    mcf.update_cost(e, -cost_m.coeff(linear_coeff(graph_matching_input::no_assignment, i)));
                else
                    static_assert("Only Eigen::MatrixXd and Eigen::SparseMatrix<double> allowed as template types");
            }
        }

    // possibly sparse matrix does not accelerate
    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    Eigen::SparseVector<double> graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::read_solution() const
    { 
        Eigen::SparseVector<double> sol((no_left_nodes()+1)*(no_right_nodes()+1));
        sol.reserve(no_left_nodes() + no_right_nodes());
        // assignment
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            const auto e_first = mcf.first_outgoing_arc(i);
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
        // non-assignment
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            const std::size_t e = mcf.first_outgoing_arc(i) + mcf.no_outgoing_arcs(i) - 1;
            if(mcf.flow(e) == 1) {
                sol.insert(linear_coeff(i, graph_matching_input::no_assignment)) = 1.0;
            }
        }

        for(std::size_t i=0; i<no_right_nodes(); ++i) {
            const std::size_t e = mcf.first_outgoing_arc(no_left_nodes() + i) + mcf.no_outgoing_arcs(no_left_nodes() + i) - 1;
            assert(mcf.flow(e) == 0 || mcf.flow(e) == -1);
            if(mcf.flow(e) == -1) {
                sol.insert(linear_coeff(graph_matching_input::no_assignment, i)) = 1.0;
            }
        }

        assert(feasible(sol));
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
        assert(feasible(M));

        if(iter%20 == 0)
            std::cout << "iteration = " << iter << ", optimal step size = " << gamma << ", quadratic term = " << quadratic_term << ", linear term = " << linear_term << "\n";

        return true;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    Eigen::SparseVector<double> graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::round()
    {
        update_costs(ASSIGNMENT_VECTOR_TYPE(-M));
        mcf.solve();
        return read_solution();
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
        graph_matching_input::labeling graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::get_solution()
        {
            const Eigen::SparseVector<double> sol_matrix = round();
            assert(sol_matrix.size() == (no_left_nodes()+1)*(no_right_nodes()+1));
            graph_matching_input::labeling sol(no_left_nodes(), graph_matching_input::no_assignment);

            for(Eigen::SparseVector<double>::InnerIterator iter(sol_matrix); iter; ++iter){
                assert(iter.value() == 1.0);
                const auto [i,j] = matrix_coeffs(iter.index());
                if(i < no_left_nodes()) {
                    assert(sol[i] == graph_matching_input::no_assignment);
                    sol[i] = j;
                } else {
                    assert(i == graph_matching_input::no_assignment);
                }
            }
            return sol;
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

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::size_t graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::linear_coeff(const std::size_t i, const std::size_t j) const 
    { 
        assert(i<no_left_nodes() || i == graph_matching_input::no_assignment);
        assert(j<no_right_nodes() || j == graph_matching_input::no_assignment);
        const std::size_t I = std::min(i, no_left_nodes());
        const std::size_t J = std::min(j, no_right_nodes());
        return I*(no_right_nodes()+1) + J; 
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::array<std::size_t,2> graph_matching_frank_wolfe<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::matrix_coeffs(const std::size_t x) const
    {
        assert(x < (no_left_nodes()+1)*(no_right_nodes()+1));
        const std::size_t i = x/(no_right_nodes()+1);
        const std::size_t j = x%(no_right_nodes()+1);
        const std::size_t ii = i < no_left_nodes() ? i : graph_matching_input::no_assignment;
        const std::size_t jj = j < no_right_nodes() ? j : graph_matching_input::no_assignment;
        assert(x == linear_coeff(ii,jj));
        return {ii, jj};
    }

    using graph_matching_frank_wolfe_dense = graph_matching_frank_wolfe<Eigen::VectorXd, Eigen::MatrixXd>;
    using graph_matching_frank_wolfe_sparse = graph_matching_frank_wolfe<Eigen::SparseVector<double>, Eigen::SparseMatrix<double>>;

}
*/
