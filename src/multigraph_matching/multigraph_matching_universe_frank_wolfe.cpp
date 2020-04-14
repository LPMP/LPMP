#include "graph_matching/matching_problem_input.h"
#include <vector>
#include <algorithm>
#include <array>
#include <limits>
#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>
#include "MCF-SSP/mcf_ssp.hxx"
#include "commutation_matrix.h"
#include "polynomial_functions.h"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    class multigraph_matching_frank_wolfe_universe {
        public:
            multigraph_matching_frank_wolfe_universe(const multigraph_matching_input& instance, const std::size_t universe_size = std::numeric_limits<std::size_t>::max());
            multigraph_matching_frank_wolfe_universe(const multigraph_matching_input &instance, const multigraph_matching_input::labeling &l);

                template <typename MATRIX>
                bool feasible(const MATRIX &m) const;
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
            template<typename MATRIX>
                void update_costs(const MATRIX& m);
            Eigen::SparseMatrix<double> read_solution() const;

            std::array<std::size_t,2> linear_indices(const std::size_t p, const std::size_t q, const std::size_t i, const std::size_t j) const; // graph 1, graph 2, node in graph 1, node in graph 2
            std::array<std::size_t,2> quadratic_indices(const std::size_t p, const std::size_t q, const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2) const; // graph 1, graph 2, nodes in graph 1, nodes in graph 2

            auto get_cost(const std::size_t p, const std::size_t q) { return Q.block(q_gs.node_no(p,0), q_gs.node_no(q,0), q_gs.no_nodes(p), q_gs.no_nodes(q)); }
            const auto get_cost(const std::size_t p, const std::size_t q) const { return Q.block(q_gs.node_no(p,0), q_gs.node_no(q,0), q_gs.no_nodes(p), q_gs.no_nodes(q)); }

            template<typename MATRIX>
                auto get_assignment_slice(MATRIX& m, const std::size_t p) const { return m.block(gs.node_no(p,0), 0, gs.no_nodes(p), universe_size()); }
            auto get_assignment_slice(const std::size_t p) { return get_assignment_slice(M, p); }
            const auto get_assignment_slice(const std::size_t p) const { return get_assignment_slice(M, p); }

            Eigen::MatrixXd Q;
            Eigen::MatrixXd L;
            Eigen::MatrixXd M;
            std::vector<MCF::SSP<int, double>> mcfs;
            multigraph_matching_input::graph_size gs;
            multigraph_matching_input::graph_size q_gs;

            const multigraph_matching_input& instance_;
    };

    multigraph_matching_frank_wolfe_universe::multigraph_matching_frank_wolfe_universe(const multigraph_matching_input& instance, const multigraph_matching_input::labeling& l)
        : gs(instance),
        instance_(instance)
    {
        // determine universe size from labeling and transform present labeling into universe one.
        throw std::runtime_error("not implemented yet.");
    }

    multigraph_matching_frank_wolfe_universe::multigraph_matching_frank_wolfe_universe(const multigraph_matching_input& instance, const std::size_t universe_size)
        : gs(instance),
        instance_(instance)
    {
        // TODO: can be better solved with transformed iterators
        std::vector<std::size_t> squared_node_no;
        std::transform(gs.no_nodes().begin(), gs.no_nodes().end(), std::back_inserter(squared_node_no), [](const std::size_t x) { return x*x; });
        q_gs = multigraph_matching_input::graph_size(squared_node_no.begin(), squared_node_no.end());

        M = Eigen::MatrixXd::Zero(gs.total_no_nodes(), std::min(universe_size,gs.total_no_nodes()));
        L = Eigen::MatrixXd::Zero(gs.total_no_nodes(), gs.total_no_nodes());
        Q = Eigen::MatrixXd::Zero(q_gs.total_no_nodes(), q_gs.total_no_nodes());

        for(const auto& gm : instance) {
            const std::size_t p = gm.left_graph_no;
            const std::size_t q = gm.right_graph_no;
            for(const auto a : gm.gm_input.assignments) {
                const auto [l1, l2] = linear_indices(p, q, a.left_node, a.right_node);
                assert(L(l1,l2) == 0.0);
                L(l1,l2) = a.cost;
            }
        }

        for(const auto& gm : instance) {
            const std::size_t left_graph = gm.left_graph_no;
            const std::size_t right_graph = gm.right_graph_no;
            for(const auto q : gm.gm_input.quadratic_terms) {
                const std::size_t l1 = gm.gm_input.assignments[q.assignment_1].left_node;
                const std::size_t l2 = gm.gm_input.assignments[q.assignment_2].left_node;
                const std::size_t r1 = gm.gm_input.assignments[q.assignment_1].right_node;
                const std::size_t r2 = gm.gm_input.assignments[q.assignment_2].right_node;

                const auto q_idx = quadratic_indices(left_graph, right_graph, l1, l2, r1, r2);
                assert(Q(q_idx[0], q_idx[1]) == 0.0);
                Q(q_idx[0], q_idx[1]) = q.cost;
            }
        }

        //Q = 0.5*(Q + Q.transpose()); // TODO: check if algorithm works on non-transposed the same

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

    std::array<std::size_t,2> multigraph_matching_frank_wolfe_universe::quadratic_indices(const std::size_t p, const std::size_t q, const std::size_t l1, const std::size_t l2, const std::size_t r1, const std::size_t r2) const
    {
        assert(p < gs.no_graphs() && q < gs.no_graphs());
        assert(l1 < gs.no_nodes(p) && l2 < gs.no_nodes(p));
        assert(r1 < gs.no_nodes(q) && r2 < gs.no_nodes(q));

        const std::array<std::size_t,2> offset = {l1*gs.no_nodes(q) + r1, l2*gs.no_nodes(q) + r2};
        return {q_gs.node_no(p,0) + offset[0], q_gs.node_no(q,0) + offset[1]}; 
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
    double multigraph_matching_frank_wolfe_universe::evaluate(const MATRIX& m) const
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
                const auto c = get_cost(j,i);
                quadratic_cost += m_ij_v.transpose() * c * m_ij_v; 
            }
        }

        return linear_cost + quadratic_cost; 
    }

    template<typename MATRIX>
        void multigraph_matching_frank_wolfe_universe::update_costs(const MATRIX& cost_m) 
        {
            assert(cost_m.rows() == M.rows() && cost_m.cols() == M.cols());
#pragma omp parallel for
            for(std::size_t i=0; i<mcfs.size(); ++i) {
                mcfs[i].reset_costs();
            }
            assert(cost_m.rows() == gs.total_no_nodes() && cost_m.cols() == universe_size());
            for(std::size_t p=0; p<gs.no_graphs(); ++p) {
                auto& mcf = mcfs[p];
                const auto block = cost_m.block(gs.node_no(p,0), 0, gs.no_nodes(p), universe_size());

                for(std::size_t i=0; i<gs.no_nodes(p); ++i) {
                    auto e_first = mcf.first_outgoing_arc(i);
                    for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                        const std::size_t e = e_first + e_c;
                        assert(mcf.head(e) >= gs.no_nodes(p));
                        const std::size_t j = mcf.head(e) - gs.no_nodes(p);
                        mcf.update_cost(e, block(i,j));
                    }
                }
            }
        }

    // possibly use sparse matrix here
    Eigen::SparseMatrix<double> multigraph_matching_frank_wolfe_universe::read_solution() const
    { 
        Eigen::SparseMatrix<double> sol(M.rows(), M.cols());
        sol.reserve(Eigen::VectorXi::Constant(sol.cols(), gs.no_graphs()));
        for(std::size_t p=0; p<gs.no_graphs(); ++p) {
            const auto& mcf = mcfs[p];

            for(std::size_t i=0; i<gs.no_nodes(p); ++i) {
                auto e_first = mcf.first_outgoing_arc(i);
                for(std::size_t e_c=0; e_c+1<mcf.no_outgoing_arcs(i); ++e_c) {
                    const std::size_t e = e_first + e_c;
                    assert(mcf.flow(e) == 0 || mcf.flow(e) == 1);
                    if(mcf.flow(e) == 1) {
                        assert(mcf.head(e) >= gs.no_nodes(p));
                        const std::size_t j = mcf.head(e) - gs.no_nodes(p);
                        sol.insert(gs.node_no(p,i),j) = 1.0;
                    }
                }
            }
        }
        return sol; 
    }

    bool multigraph_matching_frank_wolfe_universe::perform_fw_step(const std::size_t iter)
    {
        assert(feasible());

        // compute gradient
        // (i) linear term
        const auto l_g = (L + L.transpose())*M;
        
        // (ii) quadratic term
        Eigen::MatrixXd q_g(M.rows(), M.cols());
        q_g.setZero();

#pragma omp parallel for
        for(std::size_t i=0; i<gs.no_graphs(); ++i) {
            auto g_i = get_assignment_slice(q_g, i);
            const Eigen::MatrixXd M_i = get_assignment_slice(M, i);
            const Eigen::Map<const Eigen::VectorXd> M_i_v(M_i.data(), M_i.size());

            for(std::size_t j=0; j<gs.no_graphs(); ++j) {
                if(i == j) 
                    continue;
                const Eigen::MatrixXd M_j = get_assignment_slice(M, j);
                const auto Id = identity_matrix(gs.no_nodes(i));
                const auto c = get_cost(j,i);
                const Eigen::MatrixXd Q = Eigen::kroneckerProduct(M_j.transpose(), Id) * (c+c.transpose()) * Eigen::kroneckerProduct(M_j, Id);
                const Eigen::VectorXd q_diff_v = Q * M_i_v;
                const Eigen::Map<const Eigen::MatrixXd> q_diff(q_diff_v.data(), g_i.rows(), g_i.cols());
                g_i += q_diff;
            }

            for(std::size_t j=0; j<gs.no_graphs(); ++j) {
                if(i == j)
                    continue;
                const auto K = commutation_matrix(M_i.cols(), M_i.rows());
                const Eigen::MatrixXd M_j = get_assignment_slice(M, j);
                const auto Id = identity_matrix(gs.no_nodes(i));
                const auto c = get_cost(i,j);
                const Eigen::MatrixXd Q = K.transpose() * Eigen::kroneckerProduct(Id, M_j.transpose()) * (c + c.transpose()) * Eigen::kroneckerProduct(Id, M_j) * K; 
                const Eigen::VectorXd q_diff_v = Q * M_i_v;
                const Eigen::Map<const Eigen::MatrixXd> q_diff(q_diff_v.data(), g_i.rows(), g_i.cols());
                g_i += q_diff; 
            } 
        }

        const auto g = l_g + q_g; // TODO: coefficient before quadratic gradient part: 0.5?

        // min_x <x,g> s.t. x in matching polytope
        update_costs(g);

#pragma omp parallel for
        for(std::size_t i=0; i<mcfs.size(); ++i) {
            mcfs[i].solve();
        }

        const auto s = read_solution();
        Eigen::MatrixXd d = -M; // update direction
        d += s;
        const Eigen::Map<const Eigen::VectorXd> d_v(d.data(), d.size()); 
        const double gap = (d.array() * (-g.array())).sum();
        if(iter % 20 == 0)
            std::cout << "gap = " << gap << "\n";
        //if(gap < 1e-4) {
        //    std::cout << "gap smaller than threshold, terminate!\n";
        //    return false;
        //}


        // optimal step size
        // quartic function has 3 possible optima and at most two minima. Since we need to cut off at 1.0, we evaluate the minima in (0,1] and evaluate at 1 as well
        double linear_term = (L.array() * (M*d.transpose() + d*M.transpose()).array()).sum();
        double quadratic_term = (L.array() * (d*d.transpose()).array()).sum();
        double cubic_term = 0.0;
        double quartic_term = 0.0;
        for(std::size_t i=0; i<gs.no_graphs(); ++i) {
            const auto m_i = get_assignment_slice(M, i);
            const auto d_i = get_assignment_slice(d, i);
            for(std::size_t j=0; j<gs.no_graphs(); ++j) {
                if(i == j)
                    continue;
                const auto m_j = get_assignment_slice(M, j);
                const auto d_j = get_assignment_slice(d, j);
                const Eigen::MatrixXd m_i_m_j = m_i * m_j.transpose();
                const Eigen::Map<const Eigen::VectorXd> m_i_m_j_v(m_i_m_j.data(), m_i_m_j.size());
                const Eigen::MatrixXd m_i_d_j = m_i * d_j.transpose();
                const Eigen::Map<const Eigen::VectorXd> m_i_d_j_v(m_i_d_j.data(), m_i_d_j.size());
                const Eigen::MatrixXd d_i_m_j = d_i * m_j.transpose();
                const Eigen::Map<const Eigen::VectorXd> d_i_m_j_v(d_i_m_j.data(), d_i_m_j.size());
                const Eigen::MatrixXd d_i_d_j = d_i * d_j.transpose();
                const Eigen::Map<const Eigen::VectorXd> d_i_d_j_v(d_i_d_j.data(), d_i_d_j.size());
                const auto c = get_cost(j,i);
                
                quartic_term += d_i_d_j_v.transpose() * c * d_i_d_j_v;

                cubic_term += m_i_d_j_v.transpose() * c * d_i_d_j_v;
                cubic_term += d_i_m_j_v.transpose() * c * d_i_d_j_v;
                cubic_term += d_i_d_j_v.transpose() * c * m_i_d_j_v;
                cubic_term += d_i_d_j_v.transpose() * c * d_i_m_j_v;

                quadratic_term += m_i_m_j_v.transpose() * c * d_i_d_j_v;
                quadratic_term += m_i_d_j_v.transpose() * c * m_i_d_j_v;
                quadratic_term += m_i_d_j_v.transpose() * c * d_i_m_j_v; 
                quadratic_term += d_i_m_j_v.transpose() * c * m_i_d_j_v; 
                quadratic_term += d_i_m_j_v.transpose() * c * d_i_m_j_v; 
                quadratic_term += d_i_d_j_v.transpose() * c * m_i_m_j_v;

                linear_term += d_i_m_j_v.transpose() * c * m_i_m_j_v;
                linear_term += m_i_d_j_v.transpose() * c * m_i_m_j_v; 
                linear_term += m_i_m_j_v.transpose() * c * d_i_m_j_v; 
                linear_term += m_i_m_j_v.transpose() * c * m_i_d_j_v;
            }
        }

        quartic_function f(quartic_term, cubic_term, quadratic_term, linear_term, 0.0);
        auto minima = f.minima();
        double best_step_size = 1.0;
        double best_val = f.evaluate(1.0);
        for(const double x : minima) {
            if(x > 0.0 && x != std::numeric_limits<double>::infinity()) {
                if(f.evaluate(x) < best_val) {
                    best_val = f.evaluate(x);
                    best_step_size = x;
                }
            }
        }
        const double gamma = 0.1*best_step_size;

        //const double gamma = 2.0/(0.1*double(iter)+3.0);

        //if(iter % 20 == 0)
            std::cout << "iteration = " << iter << ", step size = " << gamma << "\n";
        M += gamma * d;

        return true;
    }

    void multigraph_matching_frank_wolfe_universe::round()
    {
        update_costs(-M);

#pragma omp parallel for
        for(std::size_t i=0; i<mcfs.size(); ++i) {
            mcfs[i].solve();
        }

        M = read_solution();
    }

    void multigraph_matching_frank_wolfe_universe::solve()
    {
        for(std::size_t iter=0; iter<200; ++iter) {
            if(!perform_fw_step(iter))
                break;
            if(iter % 20 == 0)
                std::cout << "objective = " << evaluate() << "\n";
        }

        round();

        for(std::size_t iter=0; iter<200; ++iter) {
            if(!perform_fw_step(iter))
                break;
            if(iter % 20 == 0)
                std::cout << "objective = " << evaluate() << "\n";
        }

        round();

        std::cout << "Rounded solution cost = " << evaluate() << "\n";
        const auto s = get_solution();
        std::cout << "Exported solution cost = " << instance_.evaluate(s) << "\n";
    }

    multigraph_matching_input::labeling multigraph_matching_frank_wolfe_universe::get_solution() const
    {
        multigraph_matching_input::labeling l;
        const Eigen::MatrixXd MM = M * M.transpose();
        // assume that solution is rounded
        for(const auto& gm : instance_) {
            multigraph_matching_input::multigraph_matching_labeling_entry gm_l;
            const std::size_t p = gm.left_graph_no;
            const std::size_t q = gm.right_graph_no;
            gm_l.left_graph_no = p;
            gm_l.right_graph_no = q;
            gm_l.labeling.resize(gs.no_nodes(p), std::numeric_limits<std::size_t>::max());
            //const auto M_pq = MM.block(gs.node_no(p,0), gs.node_no(1,0), gs.no_nodes(p), gs.no_nodes(q));
            const auto M_pq = get_assignment_slice(p) * get_assignment_slice(q).transpose();
            for(std::size_t i=0; i<gs.no_nodes(p); ++i) {
                for(std::size_t j=0; j<gs.no_nodes(q); ++j) {
                    assert(std::abs(M_pq(i,j)) <= tolerance || std::abs(M_pq(i,j) - 1.0) <= tolerance);
                    if(M_pq(i,j) > 0.5) { 
                        assert(gm_l.labeling[i] == std::numeric_limits<std::size_t>::max());
                        gm_l.labeling[i] = j;
                    }
                }
            }

            l.push_back(gm_l);
            assert(gm_l.labeling.check_primal_consistency());
        }
        return l;
    }

}
