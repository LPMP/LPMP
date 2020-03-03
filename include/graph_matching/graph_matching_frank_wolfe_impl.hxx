#pragma once

#include "graph_matching/matching_problem_input.h"
#include <Eigen/Eigen>
#include <array>
#include "MCF-SSP/mcf_ssp.hxx"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    class graph_matching_frank_wolfe_impl {
        public:
            graph_matching_frank_wolfe_impl(const graph_matching_input& gm, const graph_matching_input::labeling& l,  const graph_matching_frank_wolfe_options o);
            void set_options(const graph_matching_frank_wolfe_options o) { options = o; }
            ASSIGNMENT_VECTOR_TYPE construct_matching_matrix(const graph_matching_input::labeling& l) const;

            std::size_t no_left_nodes() const;
            std::size_t no_right_nodes() const;

            template<typename MATRIX>
                bool feasible(const MATRIX& m) const;
            bool feasible() const { return feasible(M); }
            template <typename MATRIX>
            double evaluate_linear(const MATRIX &m) const;
            template <typename MATRIX>
            double evaluate_quadratic(const MATRIX &m) const;
            template <typename MATRIX>
            double evaluate(const MATRIX &m) const;
            double evaluate() const { return evaluate(M); }
            bool perform_fw_step(const std::size_t iter);
            //Eigen::SparseVector<double> round();
            void round(Eigen::SparseVector<double>& sol);
            graph_matching_input::labeling get_solution();
            void solve();

        private:
            std::size_t linear_coeff(const std::size_t i, const std::size_t j) const;
            std::array<std::size_t,2> matrix_coeffs(const std::size_t x) const;

            template<typename MATRIX>
                void update_costs(const MATRIX& m);
            //Eigen::SparseVector<double> read_solution() const;
            void read_solution(Eigen::SparseVector<double>& sol) const;

            std::array<std::size_t,2> quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2) const;

            // objective: M.vec().transpose() * Q * M.vec() + <L, M>
            const std::size_t no_left_nodes_, no_right_nodes_;
            QUADRATIC_COST_TYPE Q; // quadratic costs
            ASSIGNMENT_VECTOR_TYPE L; // linear costs
            ASSIGNMENT_VECTOR_TYPE M; // current matching
            MCF::SSP<int, double> mcf;
            double constant_ = 0.0;
            graph_matching_frank_wolfe_options options;
    };

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::graph_matching_frank_wolfe_impl(const graph_matching_input& gm, const graph_matching_input::labeling& labeling, const graph_matching_frank_wolfe_options o)
        :
            no_left_nodes_(gm.no_left_nodes),
            no_right_nodes_(gm.no_right_nodes),
            //Q((gm.no_left_nodes+1)*(gm.no_right_nodes+1), (gm.no_left_nodes+1)*(gm.no_right_nodes+1)),
            //L((gm.no_left_nodes+1)*(gm.no_right_nodes+1)),
            //M((gm.no_left_nodes+1)*(gm.no_right_nodes+1)),
            constant_(gm.get_constant()),
            options(o)
    {
        Q = QUADRATIC_COST_TYPE((gm.no_left_nodes+1)*(gm.no_right_nodes+1), (gm.no_left_nodes+1)*(gm.no_right_nodes+1));
        L = ASSIGNMENT_VECTOR_TYPE((gm.no_left_nodes + 1) * (gm.no_right_nodes + 1));
        M = ASSIGNMENT_VECTOR_TYPE((gm.no_left_nodes + 1) * (gm.no_right_nodes + 1));
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
                //assert(Q(q_idx[0], q_idx[1]) == 0.0);
                Q(q_idx[0], q_idx[1]) += q.cost;
            } 
            assert(Q((no_left_nodes()+1)*(no_right_nodes()+1)-1, (no_left_nodes()+1)*(no_right_nodes()+1)-1) == 0.0);
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
            assert(std::abs(evaluate_linear(M) - gm.evaluate_linear(labeling)) <= 1e-8);
            assert(std::abs(evaluate_quadratic(M) - gm.evaluate_quadratic(labeling)) <= 1e-8);
            assert(std::abs(gm.evaluate(labeling) - evaluate(M)) <= tolerance);
        }

        gm.initialize_mcf(mcf); 

        //if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::SparseVector<double>>)
        //    L.makeCompressed();

        if constexpr(std::is_same_v<QUADRATIC_COST_TYPE, Eigen::SparseMatrix<double>>)
            Q.makeCompressed();
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
        ASSIGNMENT_VECTOR_TYPE 
        graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::construct_matching_matrix(const graph_matching_input::labeling& labeling) const
        {
            ASSIGNMENT_VECTOR_TYPE M((no_left_nodes()+1)*(no_right_nodes()+1));
            if constexpr(std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd>) {
                M = Eigen::VectorXd::Constant(M.size(), 0.0);
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
    std::size_t graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::no_left_nodes() const
    {
        return no_left_nodes_;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::size_t graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::no_right_nodes() const
    {
        return no_right_nodes_;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::array<std::size_t,2> graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::quadratic_indices(const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2) const
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
        bool graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::feasible(const MATRIX& m) const
    {
        return true;
        // TODO: currently not working
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
        double graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::evaluate_linear(const MATRIX& m) const
    {
        return constant_ + L.dot(m);
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        double graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::evaluate_quadratic(const MATRIX& m) const
    {
        return (m.transpose() * Q * m).sum();
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        double graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::evaluate(const MATRIX& m) const
    {
        assert(feasible());
        return evaluate_linear(m) + evaluate_quadratic(m);

        const double linear_cost = L.dot(m); // (L.array()*m.array()).sum();

        const double quadratic_cost = (m.transpose() * Q * m).sum();//m.cwiseProduct(Q*m).sum();
        if constexpr (std::is_same_v<ASSIGNMENT_VECTOR_TYPE, Eigen::VectorXd> && std::is_same_v<QUADRATIC_COST_TYPE, Eigen::MatrixXd>)
        {
            double q_cost = 0.0;
            auto index_trafo = [](const std::size_t idx, const std::size_t no_nodes) {
                if(idx < no_nodes) return idx;
                else return graph_matching_input::no_assignment;
            };
            for (std::size_t i1 = 0; i1 <= no_left_nodes(); ++i1)
                for (std::size_t i2 = 0; i2 <= no_left_nodes(); ++i2)
                    for (std::size_t j1 = 0; j1 <= no_right_nodes(); ++j1)
                        for (std::size_t j2 = 0; j2 <= no_right_nodes(); ++j2) {
                            const std::size_t I1 = index_trafo(i1,no_left_nodes());
                            const std::size_t I2 = index_trafo(i2,no_left_nodes());
                            const std::size_t J1 = index_trafo(j1,no_right_nodes());
                            const std::size_t J2 = index_trafo(j2,no_right_nodes());
                            assert(std::abs(M(linear_coeff(I1,J1))) < 1e-8 || std::abs(M(linear_coeff(I1,J1)) - 1.0) < 1e-8);
                            if (M(linear_coeff(I1,J1)) >= 0.5  && M(linear_coeff(I2,J2)) >= 0.5)
                            {
                                const auto q_idx = quadratic_indices(I1, I2, J1, J2);
                                q_cost += Q(q_idx[0], q_idx[1]);
                            }
                        }

            //assert(std::abs(quadratic_cost - q_cost) <= 1e-8);
            return linear_cost + q_cost + constant_;
        }
        return linear_cost + quadratic_cost + constant_;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    template<typename MATRIX>
        void graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::update_costs(const MATRIX& cost_m) 
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
    void graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::read_solution(Eigen::SparseVector<double>& sol) const
    { 
        Eigen::SparseVector<double> tmp((no_left_nodes()+1)*(no_right_nodes()+1));
        tmp.reserve(no_left_nodes() + no_right_nodes());
        // assignment
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            const auto e_first = mcf.first_outgoing_arc(i);
            for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                const std::size_t e = e_first + e_c;
                assert(mcf.flow(e) == 0 || mcf.flow(e) == 1);
                if(mcf.flow(e) == 1) {
                    assert(mcf.head(e) >= no_left_nodes());
                    const std::size_t j = mcf.head(e) - no_left_nodes();
                    tmp.insert(linear_coeff(i,j)) = 1.0;
                }
            }
        }
        // non-assignment
        for(std::size_t i=0; i<no_left_nodes(); ++i) {
            const std::size_t e = mcf.first_outgoing_arc(i) + mcf.no_outgoing_arcs(i) - 1;
            if(mcf.flow(e) == 1) {
                tmp.insert(linear_coeff(i, graph_matching_input::no_assignment)) = 1.0;
            }
        }

        for(std::size_t i=0; i<no_right_nodes(); ++i) {
            const std::size_t e = mcf.first_outgoing_arc(no_left_nodes() + i) + mcf.no_outgoing_arcs(no_left_nodes() + i) - 1;
            assert(mcf.flow(e) == 0 || mcf.flow(e) == -1);
            if(mcf.flow(e) == -1) {
                tmp.insert(linear_coeff(graph_matching_input::no_assignment, i)) = 1.0;
            }
        }

        assert(feasible(tmp));
        sol = tmp;
        //return sol; 
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    bool graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::perform_fw_step(const std::size_t iter)
    {
        assert(feasible());

        // compute gradient
        const ASSIGNMENT_VECTOR_TYPE g(L + 2*Q*M);
        //const Eigen::Map<const Eigen::MatrixXd> g(g_v.data(), M.rows(), M.cols());

        // min_x <x,g> s.t. x in matching polytope
        update_costs(g);
        mcf.solve();
        Eigen::SparseVector<double> s;
        read_solution(s);
        //const auto s = read_solution();
        const ASSIGNMENT_VECTOR_TYPE d(s-M); // update direction
        const double gap = -d.cwiseProduct(g).sum();
        //if(iter%20 == 0)
            //std::cout << "gap = " << gap << "\n";
        if(gap < options.gap) {
            //std::cout << "gap smaller than threshold, terminate!\n";
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

        M += (gamma * d).eval(); // TODO: needed?
        assert(feasible(M));

        //if(iter%20 == 0)
        //    std::cout << "iteration = " << iter << ", optimal step size = " << gamma << ", quadratic term = " << quadratic_term << ", linear term = " << linear_term << "\n";

        return true;
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    void graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::round(Eigen::SparseVector<double>& sol)
    {
        update_costs(ASSIGNMENT_VECTOR_TYPE(-M));
        mcf.solve();
        read_solution(sol);
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
        graph_matching_input::labeling graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::get_solution()
        {
            Eigen::SparseVector<double> sol_matrix;
            round(sol_matrix);
            //const Eigen::SparseVector<double> sol_matrix = round();
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
    void graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::solve()
    {
        for(std::size_t iter=0; iter<options.max_iter; ++iter) {
            if(!perform_fw_step(iter))
                break;
            //if(iter%20 == 0)
            //    std::cout << "objective = " << evaluate() << "\n";
        }

        Eigen::SparseVector<double> M_temp;
        round(M_temp);
        M = ASSIGNMENT_VECTOR_TYPE(M_temp);
        //M = round();
        //std::cout << "Rounded solution cost = " << evaluate() << "\n";
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::size_t graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::linear_coeff(const std::size_t i, const std::size_t j) const 
    { 
        assert(i<no_left_nodes() || i == graph_matching_input::no_assignment);
        assert(j<no_right_nodes() || j == graph_matching_input::no_assignment);
        const std::size_t I = std::min(i, no_left_nodes());
        const std::size_t J = std::min(j, no_right_nodes());
        return I*(no_right_nodes()+1) + J; 
    }

    template<typename ASSIGNMENT_VECTOR_TYPE, typename QUADRATIC_COST_TYPE>
    std::array<std::size_t,2> graph_matching_frank_wolfe_impl<ASSIGNMENT_VECTOR_TYPE, QUADRATIC_COST_TYPE>::matrix_coeffs(const std::size_t x) const
    {
        assert(x < (no_left_nodes()+1)*(no_right_nodes()+1));
        const std::size_t i = x/(no_right_nodes()+1);
        const std::size_t j = x%(no_right_nodes()+1);
        const std::size_t ii = i < no_left_nodes() ? i : graph_matching_input::no_assignment;
        const std::size_t jj = j < no_right_nodes() ? j : graph_matching_input::no_assignment;
        assert(x == linear_coeff(ii,jj));
        return {ii, jj};
    }

    struct graph_matching_frank_wolfe_dense : public graph_matching_frank_wolfe_impl<Eigen::VectorXd, Eigen::MatrixXd> 
    {
        using graph_matching_frank_wolfe_impl<Eigen::VectorXd, Eigen::MatrixXd>::graph_matching_frank_wolfe_impl;;
    };

    struct  graph_matching_frank_wolfe_sparse : public graph_matching_frank_wolfe_impl<Eigen::SparseVector<double>, Eigen::SparseMatrix<double>>
    {
        using graph_matching_frank_wolfe_impl<Eigen::SparseVector<double>, Eigen::SparseMatrix<double>>::graph_matching_frank_wolfe_impl;
    };

}
