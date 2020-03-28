#pragma once

#include "bdd_solver_interface.h"
#include "bdd_branch_node.h"
#include "bdd_storage.h"
#include "cuthill-mckee.h"
#include "bfs_ordering.hxx"
#include "minimum_degree_ordering.hxx"
#include "two_dimensional_variable_array.hxx"
#include "ILP_input.h"
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <limits>
#include <chrono> // for now
#include <iostream>
#include "tclap/CmdLine.h"

namespace LPMP {

    struct bdd_anisotropic_diffusion_options {
        bdd_anisotropic_diffusion_options(int argc, char** argv);
        bdd_anisotropic_diffusion_options() {}
        enum class order{ Cuthill_McKee, BFS, minimum_degree, given_order } order = order::minimum_degree;
        enum class direction { isotropic, anisotropic };
        std::vector<direction> directions = {direction::anisotropic};
        std::vector<double> residuals = {0.0}; // round robin for residual
    };

    class bdd_anisotropic_diffusion : public bdd_solver_interface {
        public:

            void set_options(const bdd_anisotropic_diffusion_options& o) { options = o; }

            void init(const ILP_input& instance);
            void init();
            template<typename ITERATOR>
                void set_costs(ITERATOR begin, ITERATOR end);

            template<typename BDD_VARIABLES_ITERATOR>
                void add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR var_begin, BDD_VARIABLES_ITERATOR var_end, Cudd& bdd_mgr); 

            void forward_run();
            double backward_run();

            void anisotropic_diffusion_iteration();
            void anisotropic_diffusion_forward();
            void anisotropic_diffusion_backward();
            void iteration() { anisotropic_diffusion_iteration(); }

            double lower_bound() { return backward_run(); }

            template<typename SOL_ITERATOR>
                bool check_feasibility(SOL_ITERATOR sol_begin, SOL_ITERATOR sol_end) const;
            template<typename SOL_ITERATOR>
                double evaluate(SOL_ITERATOR sol_begin, SOL_ITERATOR sol_end) const;

            std::size_t nr_bdds() const { return bdd_variable_delimiters.size(); }
            std::size_t nr_variables() const { return Lagrange_multipliers.size(); }
            std::size_t nr_variables(const std::size_t bdd_nr) const { assert(bdd_nr < nr_bdds()); return bdd_variable_delimiters[bdd_nr].size(); }

            template<typename STREAM>
                void export_dot(STREAM& s) const;

            template<typename STREAM>
                void export_state_dot(STREAM& s) const;

        private:
            void init_bdd_branch_nodes();

            std::array<std::size_t,2> bdd_branch_node_range(const std::size_t bdd_nr) const;
            std::array<std::size_t,2> bdd_branch_node_range(const std::size_t bdd_nr, const std::size_t variable_index) const;
            void forward_step(const std::size_t bdd_nr, const std::size_t variable_index);
            void backward_step(const std::size_t bdd_nr, const std::size_t variable_index);

            std::size_t bdd_branch_node_index(const bdd_branch_node_opt* bdd) const;
            std::size_t bdd_branch_node_index(const bdd_branch_node_opt& bdd) const { return bdd_branch_node_index(&bdd); }

            std::vector<std::size_t> bdd_branch_node_variables() const;

            two_dim_variable_array<std::size_t> bdd_adjacency() const;
            permutation given_bdd_order() const;
            permutation find_cuthill_mkkee_order() const;
            permutation find_minimum_degree_bdd_ordering() const;
            permutation find_bfs_bdd_ordering() const;

            bdd_anisotropic_diffusion_options::direction current_direction() const;
            double current_residual() const;

            bdd_storage bdd_storage_;

            std::vector<bdd_branch_node_opt> bdd_branch_nodes_;
            std::vector<std::size_t> bdd_delimiters; // TODO: still needed?

            struct Lagrange_multiplier { std::size_t bdd_nr; std::size_t bdd_branch_node_offset; double Lagrange_multiplier; };
            two_dim_variable_array<Lagrange_multiplier> Lagrange_multipliers; 

            constexpr static std::size_t no_Lagrange_multipler = std::numeric_limits<std::size_t>::max();
            struct bdd_variable_delimiter { std::size_t variable; std::size_t bdd_branch_node_offset; std::size_t Lagrange_multipliers_index; };
            two_dim_variable_array<bdd_variable_delimiter> bdd_variable_delimiters;

            bdd_anisotropic_diffusion_options options;
            std::size_t iteration_nr = 0;
    };

    bdd_anisotropic_diffusion_options::bdd_anisotropic_diffusion_options(int argc, char** argv)
        {
            TCLAP::CmdLine cmd("Command line parser for bdd anisotropic diffusion options", ' ', " ");
            TCLAP::MultiArg<std::string> directions_arg("d","direction","propagation direction",false,"{isotropic|anisotropic}");
            cmd.add(directions_arg);
            TCLAP::ValueArg<std::string> order_arg("o","order","order of bdd subproblems",false,"bfs","{bfs|minimum-degree|cuthill-mckee|given}");
            cmd.add(order_arg);
            TCLAP::MultiArg<double> residual_arg("r","residual","residual weight in message passing",false,"[0,1]");
            cmd.add(residual_arg); 

            cmd.parse(argc, argv);

            if(directions_arg.getValue().size() > 0) {
                for(const auto d : directions_arg.getValue()) {
                    if(d == "anisotropic")
                        directions.push_back(direction::anisotropic);
                    else if(d == "isotropic")
                        directions.push_back(direction::isotropic);
                    else
                        throw std::runtime_error("direction not recognized");
                }
            }

            if(order_arg.getValue() == "bfs")
                order = order::BFS;
            else if(order_arg.getValue() == "minimum-degree")
                order = order::minimum_degree;
            else if(order_arg.getValue() == "cuthill-mckee")
                order = order::Cuthill_McKee;
            else if(order_arg.getValue() == "given")
                order = order::given_order;
            else
                throw std::runtime_error("order argument not recognized");

            if(residual_arg.getValue().size() == 0)
                residuals = {0.0};
            else 
                residuals = residual_arg.getValue(); 
        }


    two_dim_variable_array<std::size_t> bdd_anisotropic_diffusion::bdd_adjacency() const
    {
        std::vector<std::vector<std::size_t>> bdd_of_variable(bdd_storage_.nr_variables());
        for(std::size_t bdd_nr=0; bdd_nr<bdd_storage_.nr_bdds(); ++bdd_nr) {
            const std::size_t first_bdd_node = bdd_storage_.bdd_delimiters()[bdd_nr];
            const std::size_t last_bdd_node = bdd_storage_.bdd_delimiters()[bdd_nr+1];
            for(std::size_t bdd_node=first_bdd_node; bdd_node<last_bdd_node; ++bdd_node) {
                const std::size_t v = bdd_storage_.bdd_nodes()[bdd_node].variable;
                bdd_of_variable[v].push_back(bdd_nr);
            }
        }

        std::unordered_set<std::array<std::size_t,2>> bdd_adjacency; 
        for(std::size_t v=0; v<bdd_of_variable.size(); ++v) {
            // clear duplicates
            auto& bdd_nrs = bdd_of_variable[v];
            bdd_nrs.erase(std::unique(bdd_nrs.begin(), bdd_nrs.end()), bdd_nrs.end());
            for(std::size_t bdd_idx_1=0; bdd_idx_1<bdd_nrs.size(); ++bdd_idx_1) {
                for(std::size_t bdd_idx_2=bdd_idx_1+1; bdd_idx_2<bdd_nrs.size(); ++bdd_idx_2) {
                    assert(bdd_nrs[bdd_idx_1] < bdd_nrs[bdd_idx_2]);
                    bdd_adjacency.insert({bdd_nrs[bdd_idx_1], bdd_nrs[bdd_idx_2]});
                }
            }
        }

        std::vector<std::size_t> adjacency_size(bdd_storage_.nr_bdds(),0);
        for(const auto [i,j] : bdd_adjacency) {
            adjacency_size[i]++;
            adjacency_size[j]++;
        }

        two_dim_variable_array<std::size_t> adjacency(adjacency_size.begin(), adjacency_size.end());
        std::fill(adjacency_size.begin(), adjacency_size.end(), 0);
        for(const auto e : bdd_adjacency) {
            const auto [i,j] = e;
            assert(i<j);
            assert(adjacency_size[i] < adjacency[i].size());
            assert(adjacency_size[j] < adjacency[j].size());
            adjacency(i, adjacency_size[i]++) = j;
            adjacency(j, adjacency_size[j]++) = i;
        }

        return adjacency;
    }

    std::vector<std::size_t> bdd_anisotropic_diffusion::bdd_branch_node_variables() const
    {
        std::vector<std::size_t> vars(bdd_branch_nodes_.size(), std::numeric_limits<std::size_t>::max());
        for(std::size_t v=0; v<nr_variables(); ++v) {
            for(std::size_t bdd_index=0; bdd_index=Lagrange_multipliers[v].size(); ++bdd_index) {
                const auto& first_bdd = bdd_branch_nodes_[Lagrange_multipliers(v,bdd_index).bdd_branch_node_offset];
                for(std::size_t bdd_node_index = Lagrange_multipliers(v,bdd_index).bdd_branch_node_offset; ; ++bdd_node_index) {
                    const auto& bdd = bdd_branch_nodes_[bdd_node_index];
                    if(first_bdd.variable_cost == bdd.variable_cost)
                        vars[bdd_node_index] = v;
                    else
                        break; 
                } 
            } 
        } 

        for(const std::size_t x : vars) {
            assert(x != std::numeric_limits<std::size_t>::max());
        }

        return vars;
    }


    permutation bdd_anisotropic_diffusion::given_bdd_order() const
    {
        permutation order;
        for(std::size_t i=0; i<bdd_storage_.nr_bdds(); ++i)
            order.push_back(i);
        return order;
    }

    permutation bdd_anisotropic_diffusion::find_cuthill_mkkee_order() const
    {
        const two_dim_variable_array<std::size_t> bdd_adj = bdd_adjacency();
        return Cuthill_McKee(bdd_adj);
    }

    permutation bdd_anisotropic_diffusion::find_minimum_degree_bdd_ordering() const
    {
        const two_dim_variable_array<std::size_t> bdd_adj = bdd_adjacency();
        return minimum_degree_ordering(bdd_adj);
    }

    permutation bdd_anisotropic_diffusion::find_bfs_bdd_ordering() const
    {
        const two_dim_variable_array<std::size_t> bdd_adj = bdd_adjacency();
        auto order = bfs_ordering(bdd_adj);
        return order;
    }

    void bdd_anisotropic_diffusion::init(const ILP_input& instance)
    {
        bdd_storage_ = bdd_storage(instance);
        init_bdd_branch_nodes();
        set_costs(instance.objective().begin(), instance.objective().end()); 
    }

    template<typename BDD_VARIABLES_ITERATOR>
        void bdd_anisotropic_diffusion::add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR var_begin, BDD_VARIABLES_ITERATOR var_end, Cudd& bdd_mgr)
        {
            bdd_storage_.add_bdd(bdd_mgr, bdd, var_begin, var_end);
        }

    void bdd_anisotropic_diffusion::init()
    {
        init_bdd_branch_nodes();
    }

    template<typename ITERATOR>
        void bdd_anisotropic_diffusion::set_costs(ITERATOR begin, ITERATOR end)
        {
            assert(std::distance(begin, end) <= nr_variables());
            for(std::size_t v=0; v<std::distance(begin,end); ++v) {
                const double cost = *(begin+v)/Lagrange_multipliers[v].size();
                for(std::size_t j=0; j<Lagrange_multipliers[v].size(); ++j) {
                    Lagrange_multipliers(v,j).Lagrange_multiplier = cost;
                }
            }
            for(std::size_t v=std::distance(begin, end); v<Lagrange_multipliers.size(); ++v) {
                for(std::size_t j=0; j<Lagrange_multipliers[v].size(); ++j) {
                    Lagrange_multipliers(v,j).Lagrange_multiplier = 0.0;
                }
            }
        }

    void bdd_anisotropic_diffusion::init_bdd_branch_nodes()
    {
        bdd_branch_nodes_.clear();
        bdd_branch_nodes_.resize(bdd_storage_.bdd_nodes().size());

        bdd_delimiters.clear();
        bdd_delimiters.reserve(bdd_storage_.nr_bdds()+1);
        bdd_delimiters.push_back(0);

        const auto bdd_order = [&]() {
            if(options.order == bdd_anisotropic_diffusion_options::order::given_order)
                return given_bdd_order();
            if(options.order == bdd_anisotropic_diffusion_options::order::Cuthill_McKee)
                return find_cuthill_mkkee_order();
            if(options.order == bdd_anisotropic_diffusion_options::order::minimum_degree)
                return find_minimum_degree_bdd_ordering();
            if(options.order == bdd_anisotropic_diffusion_options::order::BFS)
                return find_bfs_bdd_ordering();
            throw std::runtime_error("order not known");
        }();
        //for(const auto x : bdd_order)
        //    std::cout << x << ",";
        //std::cout << "\n";

        std::vector<std::size_t> bdd_variable_delimiter_size(bdd_storage_.nr_bdds(), 0);
        std::vector<std::vector<bdd_variable_delimiter>> bdd_variable_delimiters_tmp;
        struct Lagrange_multiplier_tmp { std::size_t bdd_nr; std::size_t bdd_branch_node_offset; };
        std::vector<std::vector<Lagrange_multiplier_tmp>> Lagrange_multipliers_tmp(bdd_storage_.nr_variables());
        std::vector<std::size_t> Lagrange_multipliers_size(bdd_storage_.nr_variables(), 0);
        std::vector<std::size_t> bdd_branch_node_variable(bdd_storage_.bdd_nodes().size(), std::numeric_limits<std::size_t>::max());

        // distribute bdd nodes according to (i) bdd nr in order and (ii) variable the bdd node corresponds to
        std::size_t cur_bdd_branch_nodes_offset = 0;
        for(std::size_t i=0; i<bdd_order.size(); ++i) {
            const std::size_t bdd_nr = bdd_order[i];
            const std::size_t first_bdd_node = bdd_storage_.bdd_delimiters()[bdd_nr];
            const std::size_t last_bdd_node = bdd_storage_.bdd_delimiters()[bdd_nr+1];
            bdd_delimiters.push_back(last_bdd_node);

            std::unordered_map<std::size_t,std::size_t> bdd_node_counter_per_var;
            struct bdd_node_counter { std::size_t variable; std::size_t nr_bdd_nodes; };
            std::vector<bdd_node_counter> bdd_variable_nr_sorted;

            for(std::size_t bdd_storage_node_idx = first_bdd_node; bdd_storage_node_idx<last_bdd_node; ++bdd_storage_node_idx) {
                const auto& bdd_storage_node = bdd_storage_.bdd_nodes()[bdd_storage_node_idx];
                bdd_node_counter_per_var[bdd_storage_node.variable]++;
            }
            for(const auto nr_nodes_of_var : bdd_node_counter_per_var)
                bdd_variable_nr_sorted.push_back({std::get<0>(nr_nodes_of_var), std::get<1>(nr_nodes_of_var)});
            bdd_node_counter_per_var.clear();
            std::sort(bdd_variable_nr_sorted.begin(), bdd_variable_nr_sorted.end(), [](const bdd_node_counter& a, const bdd_node_counter& b) { return a.variable < b.variable; });
            std::size_t offset = cur_bdd_branch_nodes_offset;
            //std::cout << "bdd offset = " << offset << "\n";
            bdd_variable_delimiters_tmp.push_back({});
            for(auto& c : bdd_variable_nr_sorted) {
                bdd_node_counter_per_var.insert({c.variable, offset});
                bdd_variable_delimiters_tmp.back().push_back({c.variable, offset, Lagrange_multipliers_tmp[c.variable].size()});
                Lagrange_multipliers_tmp[c.variable].push_back({i,bdd_variable_delimiters_tmp.back().size()-1});
                offset += c.nr_bdd_nodes;
            }
            cur_bdd_branch_nodes_offset += offset - cur_bdd_branch_nodes_offset;

            bdd_variable_delimiter_size[i] = bdd_node_counter_per_var.size();

            for(const auto [v, nr_bdd_nodes] : bdd_variable_nr_sorted)
                Lagrange_multipliers_size[v]++;

            std::unordered_map<std::size_t, bdd_branch_node_opt*> bdd_storage_index_to_branch_node_index;

            auto bdd_branch_node_index_from_bdd_storage_index = [&](const std::size_t idx) -> bdd_branch_node_opt* {
                if(idx == bdd_storage::bdd_node::terminal_0)
                    return bdd_branch_node_opt::terminal_0();
                else if(idx == bdd_storage::bdd_node::terminal_1)
                    return bdd_branch_node_opt::terminal_1();
                assert(bdd_storage_index_to_branch_node_index.count(idx) > 0);
                return bdd_storage_index_to_branch_node_index.find(idx)->second; 
            };

            auto new_branch_node_index = [&](std::size_t idx) -> std::size_t {
                assert(bdd_storage_index_to_branch_node_index.count(idx) == 0);
                const bdd_storage::bdd_node bdd_storage_node = bdd_storage_.bdd_nodes()[idx];
                const std::size_t variable = bdd_storage_node.variable;
                const std::size_t bdd_branch_node_index = bdd_node_counter_per_var.find(variable)->second;
                bdd_node_counter_per_var.find(variable)->second++;
                // assert(bdd_branch_nodes_[bdd_branch_node_index].is_initial_state());
                bdd_storage_index_to_branch_node_index.insert({idx, &bdd_branch_nodes_[bdd_branch_node_index]});
                return bdd_branch_node_index; 
            };

            for(std::size_t bdd_storage_node_idx = first_bdd_node; bdd_storage_node_idx<last_bdd_node; ++bdd_storage_node_idx) {
                const auto& bdd_storage_node = bdd_storage_.bdd_nodes()[bdd_storage_node_idx];

                const std::size_t branch_node_index = new_branch_node_index(bdd_storage_node_idx);
                auto& bdd_node = bdd_branch_nodes_[branch_node_index];
                assert(bdd_branch_node_variable[branch_node_index] == std::numeric_limits<std::size_t>::max());
                bdd_branch_node_variable[branch_node_index] = bdd_storage_node.variable;

                auto* bdd_low_outgoing = bdd_branch_node_index_from_bdd_storage_index(bdd_storage_node.low);
                bdd_node.low_outgoing = bdd_low_outgoing;
                if(!bdd_branch_node_opt::is_terminal(bdd_low_outgoing)) {
                    bdd_node.next_low_incoming = bdd_low_outgoing->first_low_incoming;
                    bdd_low_outgoing->first_low_incoming = &bdd_node;
                }

                auto* bdd_high_outgoing = bdd_branch_node_index_from_bdd_storage_index(bdd_storage_node.high);
                bdd_node.high_outgoing = bdd_high_outgoing;
                if(!bdd_branch_node_opt::is_terminal(bdd_high_outgoing)) {
                    bdd_node.next_high_incoming = bdd_high_outgoing->first_high_incoming;
                    bdd_high_outgoing->first_high_incoming = &bdd_node;
                }
            }
        }

        // initialize bdd_variable_delimiters and Lagrange_multipliers
        bdd_variable_delimiters.resize(bdd_variable_delimiter_size.begin(), bdd_variable_delimiter_size.end());
        std::fill(bdd_variable_delimiter_size.begin(), bdd_variable_delimiter_size.end(), 0);

        Lagrange_multipliers.resize(Lagrange_multipliers_size.begin(), Lagrange_multipliers_size.end());
        std::fill(Lagrange_multipliers_size.begin(), Lagrange_multipliers_size.end(), 0);

        for(std::size_t i=0; i<bdd_variable_delimiters_tmp.size(); ++i) {
            for(std::size_t j=0; j<bdd_variable_delimiters_tmp[i].size(); ++j) {
                bdd_variable_delimiters(i,j) = bdd_variable_delimiters_tmp[i][j];
            } 
        }
        for(std::size_t i=0; i<Lagrange_multipliers_tmp.size(); ++i) {
            for(std::size_t j=0; j<Lagrange_multipliers_tmp[i].size(); ++j) {
                Lagrange_multipliers(i,j).bdd_nr = Lagrange_multipliers_tmp[i][j].bdd_nr;
                Lagrange_multipliers(i,j).bdd_branch_node_offset = Lagrange_multipliers_tmp[i][j].bdd_branch_node_offset;
            }
        }
        /*
        for(std::size_t i=0; i<bdd_storage_.nr_bdds(); ++i) {
            assert(bdd_variable_delimiters[i].size() == bdd_variable_delimiters_tmp[i].size());
            const std::size_t bdd_nr = bdd_order[i];
            const std::size_t first_bdd_node = bdd_storage_.bdd_delimiters()[bdd_nr];
            const std::size_t last_bdd_node = bdd_storage_.bdd_delimiters()[bdd_nr+1];
            std::size_t prev_variable = std::numeric_limits<std::size_t>::max();
            std::size_t c = 0;
            for(std::size_t bdd_node_idx = first_bdd_node; bdd_node_idx<last_bdd_node; ++bdd_node_idx) {
                const auto& bdd_instruction = bdd_branch_nodes_[bdd_node_idx];
                const std::size_t variable = bdd_branch_node_variable[bdd_node_idx];
                assert(variable != std::numeric_limits<std::size_t>::max());
                if(variable != prev_variable) {
                    prev_variable = variable;
                    auto& current_bdd_var_delimiter = bdd_variable_delimiters(i,c);
                    auto& current_Lagrange_multiplier = Lagrange_multipliers(variable, Lagrange_multipliers_size[variable]);

                    // TODO: not correct. Store contiguous in bdd_branch_nodes_
                    current_bdd_var_delimiter.bdd_branch_node_offset = bdd_variable_delimiters_tmp[i][c];//bdd_node_idx;
                    current_bdd_var_delimiter.Lagrange_multipliers_index = Lagrange_multipliers_size[variable];
                    current_bdd_var_delimiter.variable = variable;

                    current_Lagrange_multiplier.bdd_nr = i;
                    current_Lagrange_multiplier.bdd_branch_node_offset = bdd_node_idx;

                    c++;
                    Lagrange_multipliers_size[variable]++;
                }
            }
            assert(c == bdd_variable_delimiters[i].size());
        }
        */

        // set bdd branch instruction costs
        for(std::size_t i=0; i<nr_bdds(); ++i) {
            for(std::size_t j=0; j<nr_variables(i); ++j) {
                const auto [first_bdd_node, last_bdd_node] = bdd_branch_node_range(i, j);
                std::size_t variable = bdd_branch_node_variable[first_bdd_node];
                const auto& current_bdd_var_delimiter = bdd_variable_delimiters(i,j);
                double* Lagrange_mult = &Lagrange_multipliers(variable, current_bdd_var_delimiter.Lagrange_multipliers_index).Lagrange_multiplier;
                for(std::size_t bdd_node_idx = first_bdd_node; bdd_node_idx<last_bdd_node; ++bdd_node_idx) {
                    bdd_branch_nodes_[bdd_node_idx].variable_cost = Lagrange_mult;
                }
            }
        }

        for(const auto& bdd : bdd_branch_nodes_) {
            check_bdd_branch_node(bdd);
        }

        // check whether variables in same bdd and same variable share same Lagrange multiplier
        for(std::size_t i=0; i<nr_bdds(); ++i) {
            for(std::size_t j=0; j<nr_variables(i); ++j) {
                const auto [first_bdd_node, last_bdd_node] = bdd_branch_node_range(i, j);
                for(std::size_t bdd_node_idx = first_bdd_node+1; bdd_node_idx<last_bdd_node; ++bdd_node_idx) {
                    const auto& prev_bdd_branch_node = bdd_branch_nodes_[bdd_node_idx-1];
                    const auto& cur_bdd_branch_node = bdd_branch_nodes_[bdd_node_idx];
                    assert(prev_bdd_branch_node.variable_cost == cur_bdd_branch_node.variable_cost);
                }
            }
        }
    }

    void bdd_anisotropic_diffusion::forward_run()
    {
        for(std::size_t i=0; i<bdd_branch_nodes_.size(); ++i) {
            auto& bdd_node = bdd_branch_nodes_[i];
            bdd_node.forward_step();
        }
    }

    double bdd_anisotropic_diffusion::backward_run()
    {
        for(std::ptrdiff_t i=bdd_branch_nodes_.size()-1; i>=0; --i) {
            auto& bdd_node = bdd_branch_nodes_[i];
            bdd_node.backward_step();
        } 

        double lb = 0.0;
        for(std::size_t bdd_nr=0; bdd_nr<nr_bdds(); ++bdd_nr) {
            const auto& first_bdd_instr = bdd_branch_nodes_[bdd_variable_delimiters(bdd_nr,0).bdd_branch_node_offset];
            assert(first_bdd_instr.m == first_bdd_instr.cost_from_terminal());
            lb += first_bdd_instr.m;
        }

        return lb;
    }

    std::array<std::size_t,2> bdd_anisotropic_diffusion::bdd_branch_node_range(const std::size_t bdd_nr) const
    {
        assert(bdd_nr < nr_bdds());
        const std::size_t first_bdd_node_index = bdd_variable_delimiters(bdd_nr, 0).bdd_branch_node_offset;
        const std::size_t last_bdd_node_index = [&]() {
            if(bdd_nr+1 == nr_bdds()) {
                return bdd_branch_nodes_.size();
            } else {
                return bdd_variable_delimiters(bdd_nr+1,0).bdd_branch_node_offset;
            } 
        }();

        return {first_bdd_node_index, last_bdd_node_index}; 
    }

    std::array<std::size_t,2> bdd_anisotropic_diffusion::bdd_branch_node_range(const std::size_t bdd_nr, const std::size_t variable_index) const
    {
        assert(bdd_nr < nr_bdds());
        assert(variable_index < nr_variables(bdd_nr));

        const std::size_t first_bdd_node_index = bdd_variable_delimiters(bdd_nr, variable_index).bdd_branch_node_offset;
        const std::size_t last_bdd_node_index = [&]() {
            if(bdd_nr+1 == nr_bdds() && variable_index +1 == nr_variables(bdd_nr)) {
                return bdd_branch_nodes_.size();
            } else if(bdd_nr+1 < nr_bdds() && variable_index +1 == nr_variables(bdd_nr)) {
                return bdd_variable_delimiters(bdd_nr+1,0).bdd_branch_node_offset;
            } else {
                return bdd_variable_delimiters(bdd_nr,variable_index+1).bdd_branch_node_offset;
            } 
        }();

        return {first_bdd_node_index, last_bdd_node_index}; 
    }

    bdd_anisotropic_diffusion_options::direction bdd_anisotropic_diffusion::current_direction() const
    {
        return options.directions[iteration_nr%options.directions.size()];
    }

    double bdd_anisotropic_diffusion::current_residual() const
    {
        return options.residuals[iteration_nr%options.residuals.size()]; 
    }

    void bdd_anisotropic_diffusion::forward_step(const std::size_t bdd_nr, const std::size_t variable_index)
    {
        const auto [first_bdd_node_index, last_bdd_node_index] = bdd_branch_node_range(bdd_nr, variable_index);

        // compute min-marginal
        std::array<double,2> marginal = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        for(std::size_t i=first_bdd_node_index; i<last_bdd_node_index; ++i) {
            auto& bdd_node = bdd_branch_nodes_[i];
            bdd_node.forward_step();
            const auto [m0, m1] = bdd_node.min_marginal();
            marginal[0] = std::min(marginal[0], m0);
            marginal[1] = std::min(marginal[1], m1); 
        }

        // distribute min-marginals
        const std::size_t Lagrange_multipliers_index = bdd_variable_delimiters(bdd_nr, variable_index).Lagrange_multipliers_index;
        const std::size_t variable = bdd_variable_delimiters(bdd_nr, variable_index).variable;
        assert(Lagrange_multipliers_index < Lagrange_multipliers[variable].size());
        const double subgradient = marginal[1] < marginal[0] ? 1.0 : 0.0;
        const double r = current_residual();
        if(current_direction() == bdd_anisotropic_diffusion_options::direction::anisotropic) {
            const double delta = (1.0-r)*(marginal[1] - marginal[0]) / (Lagrange_multipliers[variable].size() - 1 - Lagrange_multipliers_index);
            //const double delta = (marginal[1] - marginal[0]) / std::max(Lagrange_multipliers_index,Lagrange_multipliers[variable].size() - 1 - Lagrange_multipliers_index);
            //const double delta = -0.01*subgradient + 0.99*(marginal[1] - marginal[0]) / std::max(Lagrange_multipliers_index, Lagrange_multipliers[variable].size() - (Lagrange_multipliers_index - 1));
            //const double delta = (marginal[1] - marginal[0]) / std::max(Lagrange_multipliers_index, Lagrange_multipliers[variable].size() - (Lagrange_multipliers_index - 1));
            //std::cout << "forward marginal for bdd " << bdd_nr << ", variable index " << variable_index << " = " << marginal[1] << "," << marginal[0] << ", marginal diff = " << (marginal[1] - marginal[0]) << ", delta = " << delta << "\n";
            const double subgradient_stepsize = 0.0*std::max( 0.5*std::abs(delta), 0.1);
            for(std::size_t l=Lagrange_multipliers_index+1; l<Lagrange_multipliers[variable].size(); ++l) {
                //std::cout << "l = " << l << "\n";
                assert(std::isfinite(delta));
                Lagrange_multipliers(variable, l).Lagrange_multiplier += delta - subgradient_stepsize*subgradient;
                Lagrange_multipliers(variable, Lagrange_multipliers_index).Lagrange_multiplier -= delta - subgradient_stepsize*subgradient;;
            }
        } else if(current_direction() == bdd_anisotropic_diffusion_options::direction::isotropic) {
            const double delta = (1.0-r)*(marginal[1] - marginal[0]) / (Lagrange_multipliers[variable].size()-1);
            for(std::size_t l=0; l<Lagrange_multipliers[variable].size(); ++l) {
                if(l != Lagrange_multipliers_index) {
                    assert(std::isfinite(delta));
                    Lagrange_multipliers(variable, l).Lagrange_multiplier += delta;
                    Lagrange_multipliers(variable, Lagrange_multipliers_index).Lagrange_multiplier -= delta;
                }
            }
        }
    }

    void bdd_anisotropic_diffusion::backward_step(const std::size_t bdd_nr, const std::size_t variable_index)
    {
        const auto [first_bdd_node_index, last_bdd_node_index] = bdd_branch_node_range(bdd_nr, variable_index);

        // compute min-marginal
        std::array<double,2> marginal = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        for(std::size_t i=first_bdd_node_index; i<last_bdd_node_index; ++i) {
            auto& bdd_node = bdd_branch_nodes_[i];
            const auto [m0, m1] = bdd_node.min_marginal();
            marginal[0] = std::min(marginal[0], m0);
            marginal[1] = std::min(marginal[1], m1); 
        }

        // distribute min-marginals
        const std::size_t Lagrange_multipliers_index = bdd_variable_delimiters(bdd_nr, variable_index).Lagrange_multipliers_index;
        const std::size_t variable = bdd_variable_delimiters(bdd_nr, variable_index).variable;
        //std::cout << "Lagrange multiplier index = " << Lagrange_multipliers_index << ", #Lagrange multipliers = " << Lagrange_multipliers[variable].size() << "\n";
        const double r = current_residual();
        if(current_direction() == bdd_anisotropic_diffusion_options::direction::anisotropic) {
            const double subgradient = marginal[1] < marginal[0] ? 1.0 : 0.0;
            const double delta = (1.0-r)*(marginal[1] - marginal[0]) / Lagrange_multipliers_index;
            //const double delta = (marginal[1] - marginal[0]) / std::max(Lagrange_multipliers_index,Lagrange_multipliers[variable].size() - 1 - Lagrange_multipliers_index);
            //const double delta = -0.01*subgradient + 0.99*(marginal[1] - marginal[0]) / std::max(Lagrange_multipliers_index, Lagrange_multipliers[variable].size() - (Lagrange_multipliers_index - 1));
            //const double delta = (marginal[1] - marginal[0]) / std::max(Lagrange_multipliers_index, Lagrange_multipliers[variable].size() - (Lagrange_multipliers_index - 1));
            //assert(std::isfinite(delta));
            //std::cout << "backward marginal for bdd " << bdd_nr << ", variable index " << variable_index << " = " << marginal[1] << "," << marginal[0] << ", marginal diff = " << (marginal[1] - marginal[0]) << ", delta = " << delta << "\n";
            const double subgradient_stepsize = 0.0*std::max( 0.5*std::abs(delta), 0.1);
            for(std::size_t l=0; l<Lagrange_multipliers_index; ++l) {
                assert(std::isfinite(delta));
                Lagrange_multipliers(variable, l).Lagrange_multiplier += delta - subgradient_stepsize*subgradient; 
                Lagrange_multipliers(variable, Lagrange_multipliers_index).Lagrange_multiplier -= delta - subgradient_stepsize*subgradient;
            }
        } else if(current_direction() == bdd_anisotropic_diffusion_options::direction::isotropic) {
            const double delta = (1.0-r)*(marginal[1] - marginal[0]) / (Lagrange_multipliers[variable].size() - 1);
            for(std::size_t l=0; l<Lagrange_multipliers[variable].size(); ++l) {
                if(l != Lagrange_multipliers_index) {
                    assert(std::isfinite(delta));
                    Lagrange_multipliers(variable, l).Lagrange_multiplier += delta;
                    Lagrange_multipliers(variable, Lagrange_multipliers_index).Lagrange_multiplier -= delta;
                }
            } 
        }

        for(std::size_t i=first_bdd_node_index; i<last_bdd_node_index; ++i) {
            auto& bdd_node = bdd_branch_nodes_[i];
            bdd_node.backward_step();
        }
    }

    void bdd_anisotropic_diffusion::anisotropic_diffusion_iteration()
    {
        const auto begin_time = std::chrono::steady_clock::now();
        anisotropic_diffusion_forward();
        const auto after_forward = std::chrono::steady_clock::now(); 
        std::cout << "forward pass took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(after_forward - begin_time).count() << " milliseconds\n";
        anisotropic_diffusion_backward();
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "backward pass took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - after_forward).count() << " milliseconds\n";
        ++iteration_nr;
    }

    void bdd_anisotropic_diffusion::anisotropic_diffusion_forward()
    {
        for(std::size_t bdd_nr=0; bdd_nr<nr_bdds(); ++bdd_nr) {
            const auto [first_bdd_index,last_bdd_index] = bdd_branch_node_range(bdd_nr);
            for(std::ptrdiff_t i=last_bdd_index-1; i>=std::ptrdiff_t(first_bdd_index); --i) {
                auto& bdd_node = bdd_branch_nodes_[i];
                bdd_node.backward_step();
            } 
            for(std::size_t variable_index=0; variable_index<bdd_variable_delimiters[bdd_nr].size(); ++variable_index) {
                forward_step(bdd_nr, variable_index);
            }
        }
    }

    void bdd_anisotropic_diffusion::anisotropic_diffusion_backward()
    {
        for(std::ptrdiff_t bdd_nr=nr_bdds()-1; bdd_nr>=0; --bdd_nr) {
            const auto [first_bdd_index,last_bdd_index] = bdd_branch_node_range(bdd_nr);
            for(std::size_t i=first_bdd_index; i<last_bdd_index; ++i) {
                auto& bdd_node = bdd_branch_nodes_[i];
                bdd_node.forward_step();
            }
            for(std::ptrdiff_t variable_index=bdd_variable_delimiters[bdd_nr].size()-1; variable_index>=0; --variable_index) {
                backward_step(bdd_nr, variable_index);
            }
        }
    }

    std::size_t bdd_anisotropic_diffusion::bdd_branch_node_index(const bdd_branch_node_opt* bdd) const
    {
        assert(bdd >= &bdd_branch_nodes_[0]);
        const std::size_t i = bdd - &bdd_branch_nodes_[0];
        assert(i < bdd_branch_nodes_.size());
        return i; 
    }

    template<typename SOL_ITERATOR>
        bool bdd_anisotropic_diffusion::check_feasibility(SOL_ITERATOR sol_begin, SOL_ITERATOR sol_end) const
        {
            assert(std::distance(sol_begin, sol_end) == nr_variables());

            std::vector<char> bdd_branch_node_marks(bdd_branch_nodes_.size(), 0);

            for(std::size_t v=0; v<nr_variables(); ++v) {
                const char val = *(sol_begin+v);
                assert(val ==0 || val == 1);
                for(std::size_t bdd_index=0; bdd_index<Lagrange_multipliers[v].size(); ++bdd_index) {
                    const auto& first_bdd = bdd_branch_nodes_[Lagrange_multipliers(v,bdd_index).bdd_branch_node_offset];
                    if(first_bdd.is_first())
                        bdd_branch_node_marks[bdd_index] = 1;
                    for(std::size_t bdd_node_index = Lagrange_multipliers(v,bdd_index).bdd_branch_node_offset; ; ++bdd_node_index) {
                        const auto& bdd = bdd_branch_nodes_[bdd_node_index];
                        if(bdd_node_index < bdd_branch_nodes_.size() && first_bdd.variable_cost == bdd.variable_cost) {
                            if(bdd_branch_node_marks[bdd_node_index] == 1) {
                                const auto& bdd = bdd_branch_nodes_[bdd_node_index];
                                const auto* bdd_next_index = [&]() {
                                    if(val == false)
                                        return bdd.low_outgoing;
                                    else 
                                        return bdd.high_outgoing;
                                }();

                                if(bdd_next_index == bdd_branch_node_opt::terminal_0())
                                    return false;
                                if(bdd_next_index == bdd_branch_node_opt::terminal_1()) {
                                } else { 
                                    bdd_branch_node_marks[ bdd_branch_node_index(bdd_next_index) ] = 1;
                                }
                            }
                        } else {
                            break;
                        }
                    } 
                }
            }
            return true;
        }

template<typename SOL_ITERATOR>
double bdd_anisotropic_diffusion::evaluate(SOL_ITERATOR sol_begin, SOL_ITERATOR sol_end) const
{
    if(!check_feasibility(sol_begin, sol_end))
                return std::numeric_limits<double>::infinity();
            
            // evaluate cost
            double cost = 0.0;
            for(std::size_t v=0; v<nr_variables(); ++v) {
                const char val = *(sol_begin+v);
                assert(val == 0 || val == 1);
                if(val == 1)
                    for(std::size_t j=0; j<Lagrange_multipliers[v].size(); ++j) {
                        cost += Lagrange_multipliers(v,j).Lagrange_multiplier;
                    }
            }

            return cost;
        }

    template<typename STREAM>
        void bdd_anisotropic_diffusion::export_dot(STREAM& s) const
        {
            s << "digraph bdd_anisotropic_diffusion {\n";
            std::vector<char> bdd_node_visited(bdd_branch_nodes_.size(), false);
            std::size_t cur_bdd_index = 0;
            for(const auto& bdd : bdd_branch_nodes_) {
                const std::size_t i = bdd_branch_node_index(bdd);
                if(bdd_node_visited[i])
                    continue;
                bdd_node_visited[i] = true;
                auto get_node_string = [&](const bdd_branch_node_opt* bdd) -> std::string {
                    if(bdd == bdd_branch_node_opt::terminal_0())
                        return "false";
                    if(bdd == bdd_branch_node_opt::terminal_1())
                        return "true";
                    return std::to_string(bdd_branch_node_index(bdd));
                };

                s << i << " -> " << get_node_string(bdd.low_outgoing) << " [label=\"0\"];\n";
                s << i << " -> " << get_node_string(bdd.high_outgoing) << " [label=\"1\"];\n";
                const auto [L_1, L_2] = Lagrange_multipliers.indices(reinterpret_cast<Lagrange_multiplier*>(bdd.variable_cost - 2));

                s << i << " -> " << L_1 * 1000 + L_2 + 1000000 << " [label=\"" << L_1 << "," << L_2 << "\"];\n";
            }

            s << "}\n"; 
        }

template<typename STREAM>
void bdd_anisotropic_diffusion::export_state_dot(STREAM& s) const
{
    auto node_nr = [&](const std::size_t bdd_nr, const std::size_t variable_index) {
        return variable_index*nr_bdds() + bdd_nr;
    };
    s << "digraph bdd_anisotropic_diffusion_state {\n";
    for(std::size_t bdd_nr=0; bdd_nr<nr_bdds(); ++bdd_nr) {
        //s << "subgraph cluster_" << bdd_nr << " {\n";
        for(std::size_t variable_index=0; variable_index<nr_variables(bdd_nr); ++variable_index) {
            const std::size_t Lagrange_multipliers_index = bdd_variable_delimiters(bdd_nr, variable_index).Lagrange_multipliers_index;
            const std::size_t variable = bdd_variable_delimiters(bdd_nr, variable_index).variable;
            const auto [first_bdd_node_index, last_bdd_node_index] = bdd_branch_node_range(bdd_nr, variable_index);
            const bdd_branch_node_opt& bdd = bdd_branch_nodes_[first_bdd_node_index];
            const auto [m0, m1] = bdd.min_marginal_debug();
            s << node_nr(bdd_nr, variable) << " [label=\"m0=" << m0 << ", m1=" << m1 << ", L=" << *bdd.variable_cost << "\" pos=\"" << bdd_nr << "," << variable << "!\"];\n"; 
            if(Lagrange_multipliers_index > 0) {
                const std::size_t prev_bdd_nr = Lagrange_multipliers(variable,Lagrange_multipliers_index-1).bdd_nr;
                s << node_nr(prev_bdd_nr, variable) << " -> " << node_nr(bdd_nr, variable) << ";\n";
            }
        }
        //s << "graph[style=dotted];\n}\n";
    } 
    s << "}\n";

}

}
