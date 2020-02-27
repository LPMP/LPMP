#pragma once

#include "bdd_solver_interface.h"
#include "two_dimensional_variable_array.hxx"
#include "cuddObj.hh"
#include "ILP_input.h"
#include "convert_pb_to_bdd.h"
#include "bdd_variable.h"
#include "bdd_branch_node.h"
#include <cassert>
#include <vector>
#include <stack>
#include <unordered_map>
#include <tsl/robin_map.h>
#include <numeric>
#include <chrono> // for now
#include "bdd_storage.h"
#include "tclap/CmdLine.h"

namespace LPMP {

    ////////////////////////////////////////////////////
    // Base Class
    ////////////////////////////////////////////////////

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    class bdd_base {
        public:
            bdd_base() {}
            bdd_base(const bdd_base&) = delete; // no copy constructor because of pointers in bdd_branch_node

            template<typename BDD_VARIABLES_ITERATOR>
                void add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end, Cudd& bdd_mgr);

            void init(const ILP_input& input);
            void init();

            std::size_t nr_variables() const { return bdd_variables_.size(); }
            std::size_t nr_bdds() const { return bdd_variables_.size()-1; }
            std::size_t nr_bdds(const std::size_t var) const { assert(var<nr_variables()); return bdd_variables_[var].size(); }

        protected:
            std::size_t bdd_branch_node_index(const BDD_BRANCH_NODE* bdd) const;
            std::size_t bdd_branch_node_index(const BDD_BRANCH_NODE& bdd) const { return bdd_branch_node_index(&bdd); }
            std::size_t variable_index(const BDD_BRANCH_NODE& bdd) const;
            std::size_t variable_index(const BDD_VARIABLE& bdd_var) const;

            std::size_t first_variable_of_bdd(const BDD_VARIABLE& bdd_var) const;
            std::size_t last_variable_of_bdd(const BDD_VARIABLE& bdd_var) const;
            bool first_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const;
            bool last_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const;

            std::vector<BDD_BRANCH_NODE> bdd_branch_nodes_;
            two_dim_variable_array<BDD_VARIABLE> bdd_variables_;

        private:
            void init_branch_nodes();

            bdd_storage bdd_storage_;
    };



  /*
   bdd_min_marginal_averaging_options::bdd_min_marginal_averaging_options(int argc, char** argv)
        {
            TCLAP::CmdLine cmd("Command line parser for bdd min marginal averation options", ' ', " ");
            TCLAP::ValueArg<std::string> reweighting_arg("o","order","reweighting scheme",false,"","{classic|SRMP}");
            cmd.add(reweighting_arg);

            cmd.parse(argc, argv);

            if(reweighting_arg.getValue() == "classic")
                averaging_type = bdd_min_marginal_averaging_options::averaging_type::classic;
            else if(reweighting_arg.getValue() == "SRMP")
                averaging_type = bdd_min_marginal_averaging_options::averaging_type::SRMP;
            else
                throw std::runtime_error("direction not recognized");
        }
        */
    
    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename BDD_VARIABLES_ITERATOR>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end, Cudd& bdd_mgr)
    {
        bdd_storage_.add_bdd(bdd_mgr, bdd, bdd_vars_begin, bdd_vars_end); 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init(const ILP_input& input)
    {
        Cudd bdd_mgr;
        bdd_storage_ = bdd_storage(input, bdd_mgr);
        init_branch_nodes();
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init()
    {
        init_branch_nodes();
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init_branch_nodes()
    {
        // allocate datastructures holding bdd instructions
        
        // helper vectors used throughout initialization
        std::vector<std::size_t> nr_bdd_nodes_per_variable(bdd_storage_.nr_variables(), 0);
        std::vector<std::size_t> nr_bdds_per_variable(bdd_storage_.nr_variables(), 0);

        // TODO: std::unordered_set might be faster (or tsl::unordered_set)
        std::vector<char> variable_counted(bdd_storage_.nr_variables(), 0);
        std::vector<std::size_t> variables_covered;

        auto variable_covered = [&](const std::size_t v) {
            assert(v < bdd_storage_.nr_variables());
            assert(variable_counted[v] == 0 || variable_counted[v] == 1);
            return variable_counted[v] == 1;
        };
        auto cover_variable = [&](const std::size_t v) {
            assert(!variable_covered(v));
            variable_counted[v] = 1;
            variables_covered.push_back(v);
        };
        auto uncover_variables = [&]() {
            for(const std::size_t v : variables_covered) {
                assert(variable_covered(v));
                variable_counted[v] = 0;
            }
            variables_covered.clear(); 
        };

        for(const auto& bdd_node : bdd_storage_.bdd_nodes()) {
            ++nr_bdd_nodes_per_variable[bdd_node.variable];
            // if we have additional gap nodes, count them too
            const auto& next_high = bdd_storage_.bdd_nodes()[bdd_node.high];
            //std::cout << bdd_node.variable << " ";
        }
        //std::cout << "\n";
        assert(bdd_storage_.bdd_nodes().size() == std::accumulate(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end(), 0));

        //std::cout << "cover variables\n";
        for(std::size_t bdd_index=0; bdd_index<bdd_storage_.bdd_delimiters().size()-1; ++bdd_index) {
            for(std::size_t i=bdd_storage_.bdd_delimiters()[bdd_index]; i<bdd_storage_.bdd_delimiters()[bdd_index+1]; ++i) {
                const std::size_t bdd_variable = bdd_storage_.bdd_nodes()[i].variable;
                if(!variable_covered(bdd_variable)) {
                    cover_variable(bdd_variable);
                    ++nr_bdds_per_variable[bdd_variable];
                }
            }
            //std::cout << "bdd index = " << bdd_index << "\n";
            uncover_variables();
        }
        //std::cout << "uncover variables\n";
        
        std::vector<std::size_t> bdd_offset_per_variable;
        bdd_offset_per_variable.reserve(bdd_storage_.nr_variables());
        bdd_offset_per_variable.push_back(0);
        std::partial_sum(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end()-1, std::back_inserter(bdd_offset_per_variable));
        assert(bdd_offset_per_variable.size() == bdd_storage_.nr_variables());
        assert(bdd_offset_per_variable.back() + nr_bdd_nodes_per_variable.back() == bdd_storage_.bdd_nodes().size());

        bdd_branch_nodes_.resize(bdd_storage_.bdd_nodes().size());
        bdd_variables_.resize(nr_bdds_per_variable.begin(), nr_bdds_per_variable.end());

        for(std::size_t v=0; v<nr_bdds_per_variable.size(); ++v) {
            //std::cout << "v = " << v << ", nr bdd nodes per var = " << nr_bdd_nodes_per_variable[v] << ", offset = " << bdd_offset_per_variable[v] << "\n";
        }

        // fill branch instructions into datastructures and set pointers
        std::fill(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end(), 0); // counter for bdd instruction per variable
        std::fill(nr_bdds_per_variable.begin(), nr_bdds_per_variable.end(), 0); // counter for bdd per variable

        //std::unordered_map<std::size_t, BDD_BRANCH_NODE*> stored_bdd_node_index_to_bdd_address;
        tsl::robin_map<std::size_t, BDD_BRANCH_NODE*> stored_bdd_node_index_to_bdd_address;
        //stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, BDD_BRANCH_NODE*>(bdd_storage::bdd_node::terminal_0, bdd_branch_node_terminal_0));
        //stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, BDD_BRANCH_NODE*>(bdd_storage::bdd_node::terminal_1, bdd_branch_node_terminal_1));


        std::size_t c = 0; // check if everything is read contiguously
        for(std::size_t bdd_index=0; bdd_index<bdd_storage_.nr_bdds(); ++bdd_index) {
            uncover_variables(); 
            stored_bdd_node_index_to_bdd_address.clear();
            stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, BDD_BRANCH_NODE*>(bdd_storage::bdd_node::terminal_0, BDD_BRANCH_NODE::terminal_0()));
            stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, BDD_BRANCH_NODE*>(bdd_storage::bdd_node::terminal_1, BDD_BRANCH_NODE::terminal_1()));
            BDD_VARIABLE* next_bdd_var = nullptr;

            //std::cout << "bdd index = " << bdd_index << "\n";
            const std::size_t first_stored_bdd_node = bdd_storage_.bdd_delimiters()[bdd_index]; 
            const std::size_t last_stored_bdd_node = bdd_storage_.bdd_delimiters()[bdd_index+1];
            //std::cout << "bdd delimiter = " << bdd_storage_.bdd_delimiters()[bdd_index+1] << "\n";
            for(std::size_t stored_bdd_node_index=first_stored_bdd_node; stored_bdd_node_index<last_stored_bdd_node; ++stored_bdd_node_index, ++c) {
                assert(c == stored_bdd_node_index);
                //std::cout << "stored bdd node index = " << stored_bdd_node_index << "\n";
                const auto& stored_bdd = bdd_storage_.bdd_nodes()[stored_bdd_node_index];
                const std::size_t v = stored_bdd.variable;
                //std::cout << "bdd variable = " << v << ", bdd variable offset = " << bdd_offset_per_variable[v] << "\n";

                auto& bdd_var = bdd_variables_(v, nr_bdds_per_variable[v]);
                if(!variable_covered(v)) {
                    // assert(bdd_var.is_initial_state());
                    cover_variable(v);

                    bdd_var.first_node_index = bdd_offset_per_variable[v];
                    bdd_var.last_node_index = bdd_offset_per_variable[v];
                    //std::cout << "bdd level offset for var = " << v << ", bdd index = " << bdd_index << " = " << stored_bdd_node_index << "\n";
                    // TODO: remove, not needed anymore
                    if(next_bdd_var != nullptr) {
                        //assert(next_bdd_var > &bdd_var);
                        bdd_var.next = next_bdd_var;
                        next_bdd_var->prev = &bdd_var;
                    }

                    next_bdd_var = &bdd_var;
                } else {
                    // assert(!bdd_var.is_initial_state());
                }

                const std::size_t bdd_branch_nodes_index = bdd_var.last_node_index; 
                bdd_var.last_node_index++;

                BDD_BRANCH_NODE& bdd = bdd_branch_nodes_[bdd_branch_nodes_index];
                // assert(bdd.is_initial_state());
                //std::cout << "address = " << &bdd << "\n";

                stored_bdd_node_index_to_bdd_address.insert({stored_bdd_node_index, &bdd});

                assert(stored_bdd_node_index_to_bdd_address.count(stored_bdd.low) > 0);
                BDD_BRANCH_NODE& bdd_low = *(stored_bdd_node_index_to_bdd_address.find(stored_bdd.low)->second);
                bdd.low_outgoing = &bdd_low;
                assert(bdd.low_outgoing != nullptr);
                if(!BDD_BRANCH_NODE::is_terminal(&bdd_low)) {
                    bdd.next_low_incoming = bdd_low.first_low_incoming;
                    bdd_low.first_low_incoming = &bdd;
                }

                assert(stored_bdd_node_index_to_bdd_address.count(stored_bdd.high) > 0);
                BDD_BRANCH_NODE& bdd_high = *(stored_bdd_node_index_to_bdd_address.find(stored_bdd.high)->second);
                bdd.high_outgoing = &bdd_high;
                if(!BDD_BRANCH_NODE::is_terminal(&bdd_high)) {
                    bdd.next_high_incoming = bdd_high.first_high_incoming;
                    bdd_high.first_high_incoming = &bdd;
                }

                //if(!bdd.low_outgoing->is_terminal()) { assert(variable_index(bdd) < variable_index(*bdd.low_outgoing)); }
                //if(!bdd.high_outgoing->is_terminal()) { assert(variable_index(bdd) < variable_index(*bdd.high_outgoing)); }
                check_bdd_branch_node(bdd, v+1 == nr_variables(), v == 0);
                //check_bdd_branch_instruction_level(bdd_level, v+1 == nr_variables(), v == 0);
                ++nr_bdd_nodes_per_variable[v]; 
                ++bdd_offset_per_variable[v];
            }

            for(const std::size_t v : variables_covered) {
                ++nr_bdds_per_variable[v];
            }
        }

        for(std::size_t v=0; v<nr_bdds_per_variable.size(); ++v) {
            assert(nr_bdds_per_variable[v] == bdd_variables_[v].size());
        }

        for(const auto& bdd : bdd_branch_nodes_) {
            check_bdd_branch_node(bdd);
        }

        for(std::size_t v=0; v<nr_variables(); ++v) {
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(v); ++bdd_index) {
                // TODO: reactive check_bdd_branch_variable
                //check_bdd_branch_instruction_level(bdd_branch_instruction_levels(v,bdd_index));
            }
        }

        // TODO: clear bdd_storage
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    std::size_t bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::bdd_branch_node_index(const BDD_BRANCH_NODE* bdd) const
    {
        assert(bdd >= &bdd_branch_nodes_[0]);
        const std::size_t i = bdd - &bdd_branch_nodes_[0];
        assert(i < bdd_branch_nodes_.size());
        return i; 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    std::size_t bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::variable_index(const BDD_BRANCH_NODE& bdd) const
    {
        const std::size_t i = bdd_branch_node_index(&bdd);

        std::size_t lb = 0;
        std::size_t ub = nr_variables()-1;
        std::size_t v = nr_variables()/2;
        while (! (i >= bdd_variables_(v, 0).first_node_index && i < bdd_variables_[v].back().last_node_index))
        {
            if (i > bdd_variables_(v, 0).first_node_index)
                lb = v+1;
            else
                ub = v-1;
            v = (lb+ub)/2;
        }
        return v;
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    std::size_t bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::variable_index(const BDD_VARIABLE& bdd_var) const
    {
        assert(&bdd_var >= &bdd_variables_(0,0));
        assert(&bdd_var <= &bdd_variables_.back().back());
        std::size_t lb = 0;
        std::size_t ub = nr_variables()-1;
        std::size_t v = nr_variables()/2;
        while(! (&bdd_var >= &bdd_variables_(v,0) && &bdd_var <= &bdd_variables_[v].back()) ) {
            if( &bdd_var > &bdd_variables_(v,0) ) {
                lb = v+1;
            } else {
                ub = v-1;
            }
            v = (ub+lb)/2; 
        }
        assert(&bdd_var >= &bdd_variables_(v,0) && &bdd_var <= &bdd_variables_[v].back());
        return v; 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    std::size_t bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::first_variable_of_bdd(const BDD_VARIABLE& bdd_var) const
    {
        BDD_VARIABLE const* p = &bdd_var;
        while(p->prev != nullptr)
            p = p->prev;
        return variable_index(*p);
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    std::size_t bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::last_variable_of_bdd(const BDD_VARIABLE& bdd_var) const
    {
        BDD_VARIABLE const* p = &bdd_var;
        while(p->next != nullptr)
            p = p->next;
        return variable_index(*p); 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    bool bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::last_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const
    {
        const BDD_VARIABLE& bdd_var = bdd_variables_(var, bdd_index);
        return bdd_var.next == nullptr; 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    bool bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::first_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const
    {
        const BDD_VARIABLE& bdd_var = bdd_variables_(var, bdd_index);
        return bdd_var.prev == nullptr; 
    }


    ////////////////////////////////////////////////////
    // Min-Marginal Averaging Solver
    ////////////////////////////////////////////////////

    struct bdd_min_marginal_averaging_options {
        bdd_min_marginal_averaging_options(int argc, char** argv);
        bdd_min_marginal_averaging_options() {}
        enum class averaging_type {classic, SRMP} averaging_type = averaging_type::classic;
    };

    class bdd_min_marginal_averaging : public bdd_base<bdd_variable_mma, bdd_branch_node_opt>, public bdd_solver_interface {
        public:

            bdd_min_marginal_averaging() {}
            bdd_min_marginal_averaging(const bdd_min_marginal_averaging&) = delete; // no copy constructor because of pointers in bdd_branch_node

            void init(const ILP_input& input);
            void init();

            template<typename ITERATOR>
                void set_costs(ITERATOR begin, ITERATOR end);

            template<typename ITERATOR>
                bool check_feasibility(ITERATOR var_begin, ITERATOR var_end) const; 
            template<typename ITERATOR>
                double evaluate(ITERATOR var_begin, ITERATOR var_end) const; 

            void forward_step(const std::size_t var, const std::size_t bdd_index);
            void backward_step(const std::size_t var, const std::size_t bdd_index);
            void forward_step(const std::size_t var);
            void backward_step(const std::size_t var);

            void backward_run(); // also used to initialize
            void forward_run();

            double lower_bound() { return lower_bound_; }
            double compute_lower_bound();
            double lower_bound_backward(const std::size_t var, const std::size_t bdd_index);

            double compute_upper_bound(const std::vector<char> & primal_solution) const;
            std::vector<double> total_min_marginals();

            void min_marginal_averaging_iteration();
            void min_marginal_averaging_forward();
            void min_marginal_averaging_backward();
            void min_marginal_averaging_forward_SRMP();
            void min_marginal_averaging_backward_SRMP();
            void min_marginal_averaging_iteration_SRMP();

            void iteration();

            const bdd_branch_node_opt& get_bdd_branch_node(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const;

            template<typename STREAM>
                void export_dot(STREAM& s) const;

            void set_options(const bdd_min_marginal_averaging_options o) { options = o; }

        protected:
            void init_costs();

            std::array<double,2> min_marginal(const std::size_t var, const std::size_t bdd_index) const;
            template<typename ITERATOR>
                static std::array<double,2> average_marginals(ITERATOR marginals_begin, ITERATOR marginals_end);
            template<typename ITERATOR>
                std::pair<std::array<double,2>, bool> average_marginals_forward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const;
            template<typename ITERATOR>
                std::pair<std::array<double,2>, bool> average_marginals_backward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const;
            void set_marginal(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals);
            void set_marginal_forward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals, const bool default_avg);
            void set_marginal_backward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals, const bool default_avg);

            std::vector<double> costs_; 
            double lower_bound_ = -std::numeric_limits<double>::infinity();
            bdd_min_marginal_averaging_options options;
    };


   bdd_min_marginal_averaging_options::bdd_min_marginal_averaging_options(int argc, char** argv)
        {
            TCLAP::CmdLine cmd("Command line parser for bdd min marginal averation options", ' ', " ");
            TCLAP::ValueArg<std::string> reweighting_arg("o","order","reweighting scheme",false,"","{classic|SRMP}");
            cmd.add(reweighting_arg);

            cmd.parse(argc, argv);

            if(reweighting_arg.getValue() == "classic")
                averaging_type = bdd_min_marginal_averaging_options::averaging_type::classic;
            else if(reweighting_arg.getValue() == "SRMP")
                averaging_type = bdd_min_marginal_averaging_options::averaging_type::SRMP;
            else
                throw std::runtime_error("direction not recognized");
        }

    void bdd_min_marginal_averaging::init_costs()
    {
        for (size_t var = 0; var < nr_variables(); var++)
        for (size_t bdd_index = 0; bdd_index < nr_bdds(var); bdd_index++)
        {
            auto & bdd_var = bdd_variables_(var, bdd_index);
            for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                bdd_branch_nodes_[node_index].variable_cost = & bdd_var.cost;
        }
        costs_.resize(nr_variables(), std::numeric_limits<double>::infinity());
    }

    void bdd_min_marginal_averaging::init(const ILP_input& input)
    {
        bdd_base<bdd_variable_mma, bdd_branch_node_opt>::init(input);
        init_costs();
        set_costs(input.objective().begin(), input.objective().end());
    }

    void bdd_min_marginal_averaging::init()
    {
        bdd_base<bdd_variable_mma, bdd_branch_node_opt>::init();
        init_costs();
    }

    void bdd_min_marginal_averaging::forward_step(const std::size_t var, const std::size_t bdd_index)
    {
        assert(var < bdd_variables_.size());
        assert(bdd_index < bdd_variables_[var].size());

        auto& bdd_var = bdd_variables_(var,bdd_index);
        assert(var != 0 || bdd_var.prev == nullptr);

        // iterate over all bdd nodes and make forward step
        const std::size_t first_node_index = bdd_var.first_node_index;
        const std::size_t last_node_index = bdd_var.last_node_index;
        for(std::size_t i=first_node_index; i<last_node_index; ++i) {
            //std::cout << "forward step for var = " << var << ", bdd_index = " << bdd_index << ", bdd_node_index = " << i << "\n";
            bdd_branch_nodes_[i].forward_step();
        }
    }

    void bdd_min_marginal_averaging::forward_step(const std::size_t var)
    {
        assert(var < bdd_variables_.size());
        for(std::size_t bdd_index = 0; bdd_index<bdd_variables_[var].size(); ++bdd_index)
            forward_step(var, bdd_index);
    }

    void bdd_min_marginal_averaging::backward_step(const std::size_t var, const std::size_t bdd_index)
    {
        assert(var < bdd_variables_.size());
        assert(bdd_index < bdd_variables_[var].size());

        auto& bdd_var = bdd_variables_(var,bdd_index);
        assert(var+1 != nr_variables() || bdd_var.next == nullptr);

        // iterate over all bdd nodes and make forward step
        const std::ptrdiff_t first_node_index = bdd_var.first_node_index;
        const std::ptrdiff_t last_node_index = bdd_var.last_node_index;
        //std::cout << "no bdd nodes for var " << var << " = " << last_node_index - first_node_index << "\n";
        //std::cout << "bdd branch instruction offset = " << first_node_index << "\n";
        for(std::ptrdiff_t i=last_node_index-1; i>=first_node_index; --i) {
            check_bdd_branch_node(bdd_branch_nodes_[i], var+1 == nr_variables(), var == 0);
            bdd_branch_nodes_[i].backward_step();
        }
    }

    void bdd_min_marginal_averaging::backward_step(const std::size_t var)
    {
        assert(var < bdd_variables_.size());
        for(std::size_t bdd_index = 0; bdd_index<bdd_variables_[var].size(); ++bdd_index)
            backward_step(var, bdd_index);
    }

    std::array<double,2> bdd_min_marginal_averaging::min_marginal(const std::size_t var, const std::size_t bdd_index) const
    {
        assert(var < nr_variables());
        assert(bdd_index < nr_bdds(var));
        std::array<double,2> m = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        const auto& bdd_var = bdd_variables_(var,bdd_index);
        for(std::size_t bdd_node_index=bdd_var.first_node_index; bdd_node_index<bdd_var.last_node_index; ++bdd_node_index) {
            //std::cout << "min marginal for var = " << var << ", bdd_index = " << bdd_index << ", bdd_node_index = " << bdd_node_index << "\n";
            const auto& bdd = bdd_branch_nodes_[bdd_node_index];
            const auto [m0,m1] = bdd.min_marginal();
            m[0] = std::min(m[0], m0);
            m[1] = std::min(m[1], m1);
        }
        assert(std::isfinite(m[0]));
        assert(std::isfinite(m[1]));
        return m;
    }

    void bdd_min_marginal_averaging::set_marginal(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals)
    {
        assert(var < nr_variables());
        assert(bdd_index < nr_bdds(var));
        assert(min_marginals == min_marginal(var,bdd_index));
        assert(std::isfinite(marginals[0]) && std::isfinite(marginals[1]));
        assert(std::isfinite(min_marginals[0]) && std::isfinite(min_marginals[1]));

        auto& bdd_var = bdd_variables_(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        assert(std::isfinite(marginal_diff));
        assert(std::isfinite(marginal_diff_target));
        bdd_var.cost += -marginal_diff + marginal_diff_target; 
    }

    void bdd_min_marginal_averaging::set_marginal_forward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals, const bool default_avg)
    {
        auto& bdd_var = bdd_variables_(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        assert(std::isfinite(marginal_diff));
        assert(std::isfinite(marginal_diff_target));
        if (default_avg)
        {
            bdd_var.cost += -marginal_diff + marginal_diff_target;
        }
        else if(last_variable_of_bdd(var, bdd_index)) {
            bdd_var.cost -= marginal_diff;
        } else {
            assert(std::isfinite(marginal_diff_target));
            bdd_var.cost += -marginal_diff + marginal_diff_target; 
        }

    }
    void bdd_min_marginal_averaging::set_marginal_backward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals, const bool default_avg)
    {
        auto& bdd_var = bdd_variables_(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        assert(std::isfinite(marginal_diff));
        assert(std::isfinite(marginal_diff_target));
        if (default_avg)
        {
            bdd_var.cost += -marginal_diff + marginal_diff_target;
        }
        else if(first_variable_of_bdd(var, bdd_index)) {
            bdd_var.cost -= marginal_diff;
        } else {
            bdd_var.cost += -marginal_diff + marginal_diff_target; 
            assert(std::isfinite(marginal_diff_target));
            bdd_var.cost += -marginal_diff + marginal_diff_target; 
        }
    }

    void bdd_min_marginal_averaging::backward_run()
    {
        for(long int var=nr_variables()-1; var>=0; --var) {
            backward_step(var); 
        } 
    }

    double bdd_min_marginal_averaging::compute_lower_bound()
    {
        double lb = 0.0;
        for(long int var=nr_variables()-1; var>=0; --var)
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index)
                lb += lower_bound_backward(var,bdd_index);
        return lb;
    }

    double bdd_min_marginal_averaging::lower_bound_backward(const std::size_t var, const std::size_t bdd_index)
    {
        const auto& bdd_var = bdd_variables_(var,bdd_index);
        if(first_variable_of_bdd(var, bdd_index)) {
            assert(bdd_var.nr_bdd_nodes() == 1);
            const auto& bdd = get_bdd_branch_node(var, bdd_index, 0);
            return bdd.m;
        } else {
            return 0.0;
        }
    }

    void bdd_min_marginal_averaging::forward_run()
    {
        for(std::size_t var=0; var<nr_variables(); ++var) {
            forward_step(var); 
        }
    }

    double bdd_min_marginal_averaging::compute_upper_bound(const std::vector<char> & primal_solution) const
    {
        if (primal_solution.size() != nr_variables())
            return std::numeric_limits<double>::infinity();
        if (!check_feasibility(primal_solution.begin(), primal_solution.end()))
            return std::numeric_limits<double>::infinity();
        else
            return evaluate(primal_solution.begin(), primal_solution.end());
    }

    std::vector<double> bdd_min_marginal_averaging::total_min_marginals()
    {
        std::vector<double> total_min_marginals;
        std::vector<double> current_marginals;
        size_t nr_conflicts = 0;
        for(std::size_t var=0; var<nr_variables(); ++var)
        {
            auto sign = [](const double x) -> int {
                if (x > 0.0) return 1;
                if (x < 0.0) return -1;
                return 0;
            };
            current_marginals.clear();
            double total_min_marg = 0;
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index)
            {
                forward_step(var,bdd_index);
                std::array<double,2> min_marg = min_marginal(var,bdd_index);
                total_min_marg += (min_marg[1] - min_marg[0]);
                current_marginals.push_back(min_marg[1] - min_marg[0]);
            }

            total_min_marginals.push_back(total_min_marg);

            // for(std::size_t i=0; i+1<current_marginals.size(); ++i) 
            //     if(sign(current_marginals[i]) != sign(current_marginals[i+1])) {
            //         total_min_marginals.back() = 0.0;
            //         ++nr_conflicts;
            //         continue;
            //     }
        }
        // std::cout << "#conflicts in min-marginals = " << nr_conflicts << "\n";
        backward_run();
        return total_min_marginals;
    }

    void bdd_min_marginal_averaging::min_marginal_averaging_iteration()
    {
        const auto begin_time = std::chrono::steady_clock::now();
        min_marginal_averaging_forward();
        const auto after_forward = std::chrono::steady_clock::now();
        std::cout << "forward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(after_forward - begin_time).count() << " ms, " << std::flush;
        const auto before_backward = std::chrono::steady_clock::now();
        min_marginal_averaging_backward();
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "backward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - before_backward).count() << " ms, " << std::flush;
    }

    template<typename ITERATOR>
        std::array<double,2> bdd_min_marginal_averaging::average_marginals(ITERATOR marginals_begin, ITERATOR marginals_end)
        {
            std::array<double,2> average_marginal = {0.0, 0.0};
            for(auto m_iter=marginals_begin; m_iter!=marginals_end; ++m_iter) {
                average_marginal[0] += (*m_iter)[0];
                average_marginal[1] += (*m_iter)[1];
            }
            const double no_marginals = std::distance(marginals_begin, marginals_end);
            average_marginal[0] /= no_marginals;
            average_marginal[1] /= no_marginals;

            assert(std::isfinite(average_marginal[0]));
            assert(std::isfinite(average_marginal[1]));
            return average_marginal;
        }

    template<typename ITERATOR>
        std::pair<std::array<double,2>, bool> bdd_min_marginal_averaging::average_marginals_forward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const
        {
            assert(nr_bdds(var) == std::distance(marginals_begin, marginals_end));
            std::array<double,2> average_marginal = {0.0, 0.0};
            std::size_t nr_averaged_marginals = 0;
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    average_marginal[0] += (*(marginals_begin+bdd_index))[0];
                    average_marginal[1] += (*(marginals_begin+bdd_index))[1];
                if(!last_variable_of_bdd(var, bdd_index))
                    ++nr_averaged_marginals;
            }
            // if no BDD satisfies forward condition, resort to averaging over all BDDs
            bool default_avg = false;
            if (nr_averaged_marginals == 0)
            {
                nr_averaged_marginals = nr_bdds(var);
                default_avg = true;
            }

            average_marginal[0] /= double(nr_averaged_marginals);
            average_marginal[1] /= double(nr_averaged_marginals);

            return std::make_pair(average_marginal, default_avg);
        }

    template<typename ITERATOR>
        std::pair<std::array<double,2>, bool> bdd_min_marginal_averaging::average_marginals_backward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const
        {
            assert(nr_bdds(var) == std::distance(marginals_begin, marginals_end));
            std::array<double,2> average_marginal = {0.0, 0.0};
            std::size_t nr_averaged_marginals = 0;
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    average_marginal[0] += (*(marginals_begin+bdd_index))[0];
                    average_marginal[1] += (*(marginals_begin+bdd_index))[1];
                if(!first_variable_of_bdd(var, bdd_index))
                    ++nr_averaged_marginals;
            }
            // if no BDD satisfies forward condition, resort to averaging over all BDDs
            bool default_avg = false;
            if (nr_averaged_marginals == 0)
            {
                nr_averaged_marginals = nr_bdds(var);
                default_avg = true;
            }

            average_marginal[0] /= double(nr_averaged_marginals);
            average_marginal[1] /= double(nr_averaged_marginals);

            assert(std::isfinite(average_marginal[0]));
            assert(std::isfinite(average_marginal[1]));
            return std::make_pair(average_marginal, default_avg);
        }

    // min marginal averaging
    void bdd_min_marginal_averaging::min_marginal_averaging_forward()
    {
        std::vector<std::array<double,2>> min_marginals;
        for(std::size_t var=0; var<nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                forward_step(var,bdd_index);
                min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
            } 
        }
    }

    void bdd_min_marginal_averaging::min_marginal_averaging_forward_SRMP()
    {
        std::vector<std::array<double,2>> min_marginals;
        for(std::size_t var=0; var<nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                forward_step(var,bdd_index);
                min_marginals.push_back(min_marginal(var,bdd_index));
            }


            const auto average_marginal = average_marginals_forward_SRMP(min_marginals.begin(), min_marginals.end(), var);
            const std::array<double,2> avg_marg = average_marginal.first;
            const bool default_averaging = average_marginal.second;

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                set_marginal_forward_SRMP(var,bdd_index,avg_marg,min_marginals[bdd_index], default_averaging);
            } 
        }
    }


    void bdd_min_marginal_averaging::min_marginal_averaging_backward()
    {
        double lb = 0.0;
        std::vector<std::array<double,2>> min_marginals;

        for(long int var=nr_variables()-1; var>=0; --var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            if(min_marginals.size() > 1) {
                //std::cout << "min marginals in backward pass: ";
                //for(const auto m : min_marginals)
                //    std::cout << "(" << m[0] << "," << m[1] << "), ";
                //std::cout << "\n";
            }

            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
                backward_step(var, bdd_index);
                lb += lower_bound_backward(var,bdd_index);
            }
        }

        lower_bound_ = lb; 
    }

    void bdd_min_marginal_averaging::min_marginal_averaging_backward_SRMP()
    {
        double lb = 0.0;
        std::vector<std::array<double,2>> min_marginals;
        for(long int var=nr_variables()-1; var>=0; --var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            const auto average_marginal = average_marginals_backward_SRMP(min_marginals.begin(), min_marginals.end(), var);
            const std::array<double,2> avg_marg = average_marginal.first;
            const bool default_averaging = average_marginal.second;

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                set_marginal_backward_SRMP(var,bdd_index,avg_marg,min_marginals[bdd_index], default_averaging);
                backward_step(var, bdd_index);
                lb += lower_bound_backward(var,bdd_index);
            }
        }

        lower_bound_ = lb; 
    }

    void bdd_min_marginal_averaging::min_marginal_averaging_iteration_SRMP()
    {
        const auto begin_time = std::chrono::steady_clock::now();
        min_marginal_averaging_forward_SRMP();
        const auto after_forward = std::chrono::steady_clock::now();
        std::cout << "forward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(after_forward - begin_time).count() << " ms, " << std::flush;
        const auto before_backward = std::chrono::steady_clock::now();
        //min_marginal_averaging_backward_SRMP();
        min_marginal_averaging_backward();
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "backward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - before_backward).count() << " ms, " << std::flush;
    }


    const bdd_branch_node_opt& bdd_min_marginal_averaging::get_bdd_branch_node(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const
    {
        assert(var < nr_variables());
        assert(bdd_index < bdd_variables_[var].size());
        assert(bdd_node_index < bdd_variables_(var,bdd_index).nr_bdd_nodes());
        return bdd_branch_nodes_[ bdd_variables_(var,bdd_index).first_node_index + bdd_node_index ];
    }


    //const bdd_branch_instruction& bdd_min_marginal_averaging::get_bdd_branch_instruction(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const
   
    void bdd_min_marginal_averaging::iteration()
    {
        if(options.averaging_type == bdd_min_marginal_averaging_options::averaging_type::classic)
            min_marginal_averaging_iteration();
        else if(options.averaging_type == bdd_min_marginal_averaging_options::averaging_type::SRMP)
            min_marginal_averaging_iteration_SRMP();
        else
            assert(false);
    }

    template<typename ITERATOR>
        void bdd_min_marginal_averaging::set_costs(ITERATOR begin, ITERATOR end)
        {
            // TODO: remove costs_ array
            std::fill(costs_.begin(), costs_.end(), 0.0);
            assert(std::distance(begin,end) <= nr_variables());
            std::copy(begin, end, costs_.begin());
            //std::fill(costs_.begin() + std::distance(begin, end), costs_.end(), 0.0);

            // distribute costs to bdds uniformly
            for(std::size_t v=0; v<nr_variables(); ++v) {
                assert(bdd_variables_[v].size() > 0);
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(v); ++bdd_index) {
                    assert(nr_bdds(v) > 0);
                    const double cost = v < std::distance(begin,end) ? *(begin+v)/nr_bdds(v) : 0.0; 
                    bdd_variables_(v,bdd_index).cost = cost;
                    assert(!std::isnan(bdd_variables_(v,bdd_index).cost));
                }
            }
            for(const auto& bdd : bdd_branch_nodes_) {
                assert(!std::isnan(*bdd.variable_cost));
            }
            backward_run();
        }

    template<typename ITERATOR>
        bool bdd_min_marginal_averaging::check_feasibility(ITERATOR var_begin, ITERATOR var_end) const
        {
            assert(std::distance(var_begin, var_end) == nr_variables());

            std::vector<char> bdd_nbranch_node_marks(bdd_branch_nodes_.size(), 0);

            std::size_t var = 0;
            for(auto var_iter=var_begin; var_iter!=var_end; ++var_iter, ++var) {
                const bool val = *(var_begin+var);
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    const auto& bdd_var = bdd_variables_(var, bdd_index);
                    if(bdd_var.is_first_bdd_variable()) {
                        bdd_nbranch_node_marks[bdd_var.first_node_index] = 1;
                    }
                    for(std::size_t bdd_node_index=bdd_var.first_node_index; bdd_node_index<bdd_var.last_node_index; ++bdd_node_index) {
                        if(bdd_nbranch_node_marks[bdd_node_index] == 1) {
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
                                bdd_nbranch_node_marks[ bdd_branch_node_index(bdd_next_index) ] = 1;
                            }
                        }
                    }
                }
            }

            return true;
        }

    template<typename ITERATOR>
        double bdd_min_marginal_averaging::evaluate(ITERATOR var_begin, ITERATOR var_end) const
        {
            assert(std::distance(var_begin, var_end) == nr_variables());

            if(!check_feasibility(var_begin, var_end))
                return std::numeric_limits<double>::infinity();

            double cost = 0.0;
            std::size_t var = 0;
            for(auto var_iter=var_begin; var_iter!=var_end; ++var_iter, ++var)
                cost += *var_iter * costs_[var];
            return cost;
        }

    template<typename STREAM>
        void bdd_min_marginal_averaging::export_dot(STREAM& s) const
        {
            return bdd_storage_.export_dot(s);
            s << "digraph bdd_min_marginal_averaging {\n";

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
            }

            s << "}\n"; 
        }

    ////////////////////////////////////////////////////
    // Variable Fixing
    ////////////////////////////////////////////////////

    struct log_entry
    {
        log_entry(char * var_value)
        : var_value_(var_value) {}
        log_entry(bdd_branch_node_fix * source, bdd_branch_node_fix * target, bool high)
        : source_(source), target_(target), high_(high) {}

        void restore();

        char * var_value_ = nullptr;

        bdd_branch_node_fix * source_;
        bdd_branch_node_fix * target_;
        bool high_;
    };

    void log_entry::restore()
    {
        if (var_value_ != nullptr)
        {
            *var_value_ = 2;
            return;    
        }

        assert(target_ != bdd_branch_node_fix::terminal_0());
        if (high_)
        {
            assert(source_->high_outgoing == bdd_branch_node_fix::terminal_0());
            assert(source_->prev_high_incoming == nullptr);
            assert(source_->next_high_incoming == nullptr);
            source_->high_outgoing = target_;
            source_->bdd_var->nr_feasible_high_arcs++;
            if (bdd_branch_node_fix::is_terminal(target_))
                return;
            if (target_->first_high_incoming != nullptr)
                target_->first_high_incoming->prev_high_incoming = source_;
            source_->next_high_incoming = target_->first_high_incoming;
            target_->first_high_incoming = source_;
        }
        else
        {
            assert(source_->low_outgoing == bdd_branch_node_fix::terminal_0());
            assert(source_->prev_low_incoming == nullptr);
            assert(source_->next_low_incoming == nullptr);
            source_->low_outgoing = target_;
            source_->bdd_var->nr_feasible_low_arcs++;
            if (bdd_branch_node_fix::is_terminal(target_))
                return;
            if (target_->first_low_incoming != nullptr)
                target_->first_low_incoming->prev_low_incoming = source_;
            source_->next_low_incoming = target_->first_low_incoming;
            target_->first_low_incoming = source_;
        }
    }

    class bdd_fixing : public bdd_base<bdd_variable_fix, bdd_branch_node_fix> {
        public:
            bool fix_variables(const std::vector<size_t> & indices, const std::vector<char> & values);

            void init(const ILP_input& input);
            void init();

            const std::vector<char> & primal_solution() const { return primal_solution_; }

        private:
            void init_pointers();

            bool fix_variables_recursive(const std::vector<size_t> & vars, const std::vector<char> & vals, const size_t index);
            bool fix_variable(const size_t var, const char value);
            bool is_fixed(const size_t var) const;

            bool remove_all_incoming_arcs(bdd_branch_node_fix & bdd_node);
            void remove_all_outgoing_arcs(bdd_branch_node_fix & bdd_node);
            void remove_outgoing_low_arc(bdd_branch_node_fix & bdd_node);
            void remove_outgoing_high_arc(bdd_branch_node_fix & bdd_node);

            // void restore_arc(const log_entry & entry);
            void revert_changes(const size_t target_log_size);

            std::vector<char> primal_solution_;
            std::stack<log_entry, std::vector<log_entry>> log_;
    };

    void bdd_fixing::init_pointers()
    {
        for (size_t var = 0; var < nr_variables(); var++)
        {
            for (size_t bdd_index = 0; bdd_index < nr_bdds(var); bdd_index++)
            {
                auto & bdd_var = bdd_variables_(var, bdd_index);
                bdd_var.nr_feasible_low_arcs = 0;
                bdd_var.nr_feasible_high_arcs = 0;
                bdd_var.variable_index = var;
                for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                {
                    auto & bdd_node = bdd_branch_nodes_[node_index];
                    if (bdd_node.low_outgoing != bdd_branch_node_fix::terminal_0())
                        bdd_var.nr_feasible_low_arcs++;
                    if (bdd_node.high_outgoing != bdd_branch_node_fix::terminal_0())
                        bdd_var.nr_feasible_high_arcs++;

                    bdd_node.bdd_var = & bdd_var;

                    auto * low_incoming = bdd_node.first_low_incoming;
                    while (low_incoming != nullptr && low_incoming->next_low_incoming != nullptr)
                    {
                        low_incoming->next_low_incoming->prev_low_incoming = low_incoming;
                        low_incoming = low_incoming->next_low_incoming;
                    }
                    auto * high_incoming = bdd_node.first_high_incoming;
                    while (high_incoming != nullptr && high_incoming->next_high_incoming != nullptr)
                    {
                        high_incoming->next_high_incoming->prev_high_incoming = high_incoming;
                        high_incoming = high_incoming->next_high_incoming;
                    }
                }
            }
        }
        primal_solution_.resize(nr_variables(), 2); 
    }

    void bdd_fixing::init(const ILP_input& input)
    {
        bdd_base<bdd_variable_fix, bdd_branch_node_fix>::init(input);
        init_pointers();
    }

    void bdd_fixing::init()
    {
        bdd_base<bdd_variable_fix, bdd_branch_node_fix>::init();
        init_pointers();
    }

    bool bdd_fixing::fix_variable(const std::size_t var, const char value)
    {
        assert(0 <= value && value <= 1);
        assert(primal_solution_.size() == nr_variables());
        assert(var < primal_solution_.size());

        // check if variable is already fixed
        if (primal_solution_[var] == value)
            return true;
        else if (is_fixed(var))
            return false;

        // mark variable as fixed
        primal_solution_[var] = value;
        const log_entry entry(&primal_solution_[var]);
        log_.push(entry);
        std::vector<std::pair<size_t, char>> restrictions;

        for (size_t bdd_index = 0; bdd_index < nr_bdds(var); bdd_index++)
        {
            auto & bdd_var = bdd_variables_(var, bdd_index);
            for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
            {
                // auto & bdd_node = bdd_branch_nodes_[node_index];
                auto & bdd_node = bdd_branch_nodes_.at(node_index);

                // skip isolated branch nodes
                if (bdd_node.is_first() && bdd_node.is_dead_end())
                    continue;

                if (value == 1)
                    remove_outgoing_low_arc(bdd_node);
                if (value == 0)
                    remove_outgoing_high_arc(bdd_node);

                // restructure parents if node is dead-end
                if (bdd_node.is_dead_end())
                {
                    if (!remove_all_incoming_arcs(bdd_node))
                        return false;
                }
            }

            // check if other variables in BDD are now restricted
            auto * cur = bdd_var.prev;
            while (cur != nullptr)
            {
                if (is_fixed(cur->variable_index))
                {
                    cur = cur->prev;
                    continue;
                }
                if (cur->nr_feasible_low_arcs == 0)
                    restrictions.emplace_back(cur->variable_index, 1);
                if (cur->nr_feasible_high_arcs == 0)
                    restrictions.emplace_back(cur->variable_index, 0);
                cur = cur->prev;
            }
            cur = bdd_var.next;
            while (cur != nullptr)
            {
                if (is_fixed(cur->variable_index))
                {
                    cur = cur->next;
                    continue;
                }
                if (cur->nr_feasible_low_arcs == 0)
                    restrictions.emplace_back(cur->variable_index, 1);
                if (cur->nr_feasible_high_arcs == 0)
                    restrictions.emplace_back(cur->variable_index, 0);
                cur = cur->next;
            }
        }

        // fix implied restrictions
        for (auto & restriction : restrictions)
        {
            if (!fix_variable(restriction.first, restriction.second))
                return false;
        }

        return true;
    }

    bool bdd_fixing::remove_all_incoming_arcs(bdd_branch_node_fix & bdd_node)
    {
        if (bdd_node.is_first())
            return false;
        // low arcs
        {
            auto * cur = bdd_node.first_low_incoming;
            while (cur != nullptr)
            {
                // log change
                assert(cur->low_outgoing == &bdd_node);
                auto * temp = cur;
                const log_entry entry(cur, &bdd_node, false);
                log_.push(entry);
                // remove arc
                cur->low_outgoing = bdd_branch_node_fix::terminal_0();
                assert(cur->bdd_var != nullptr);
                cur->bdd_var->nr_feasible_low_arcs--;
                bdd_node.first_low_incoming = cur->next_low_incoming;
                cur = cur->next_low_incoming;
                // remove list pointers
                if (cur != nullptr)
                {
                    assert(cur->prev_low_incoming != nullptr);
                    cur->prev_low_incoming->next_low_incoming = nullptr;
                    cur->prev_low_incoming = nullptr;    
                }
                // recursive call if parent is dead-end
                if (temp->is_dead_end())
                {
                    if (!remove_all_incoming_arcs(*temp))
                        return false;
                }
            }
        }
        // high arcs
        {
            auto * cur = bdd_node.first_high_incoming;
            while (cur != nullptr)
            {
                assert(cur->high_outgoing == &bdd_node);
                auto * temp = cur;
                const log_entry entry(cur, &bdd_node, true);
                log_.push(entry);
                cur->high_outgoing = bdd_branch_node_fix::terminal_0();
                assert(cur->bdd_var != nullptr);
                cur->bdd_var->nr_feasible_high_arcs--;
                bdd_node.first_high_incoming = cur->next_high_incoming; 
                cur = cur->next_high_incoming; 
                if (cur != nullptr)
                {
                    assert(cur->prev_low_incoming != nullptr);
                    cur->prev_high_incoming->next_high_incoming = nullptr;
                    cur->prev_high_incoming = nullptr;    
                }
                if (temp->is_dead_end())
                {
                    if (!remove_all_incoming_arcs(*temp))
                        return false;
                }
            }      
        }
        return true;
    }

    void bdd_fixing::remove_all_outgoing_arcs(bdd_branch_node_fix & bdd_node)
    {
        remove_outgoing_low_arc(bdd_node);
        remove_outgoing_high_arc(bdd_node);
    }

    void bdd_fixing::remove_outgoing_low_arc(bdd_branch_node_fix & bdd_node)
    {
        if (!bdd_branch_node_fix::is_terminal(bdd_node.low_outgoing))
        {
            // change pointers
            if (bdd_node.prev_low_incoming == nullptr)
                bdd_node.low_outgoing->first_low_incoming = bdd_node.next_low_incoming;
            else
                bdd_node.prev_low_incoming->next_low_incoming = bdd_node.next_low_incoming;
            if (bdd_node.next_low_incoming != nullptr)
                bdd_node.next_low_incoming->prev_low_incoming = bdd_node.prev_low_incoming;
            bdd_node.prev_low_incoming = nullptr;
            bdd_node.next_low_incoming = nullptr;
            // recursive call if child node is unreachable
            if (bdd_node.low_outgoing->is_first())
                remove_all_outgoing_arcs(*bdd_node.low_outgoing);
        }
        if (bdd_node.low_outgoing != bdd_branch_node_fix::terminal_0())
        {
            // log change
            const log_entry entry(&bdd_node, bdd_node.low_outgoing, false);
            log_.push(entry);
            // remove arc
            bdd_node.low_outgoing = bdd_branch_node_fix::terminal_0();
            bdd_node.bdd_var->nr_feasible_low_arcs--;
        } 
    }

    void bdd_fixing::remove_outgoing_high_arc(bdd_branch_node_fix & bdd_node)
    {
        if (!bdd_branch_node_fix::is_terminal(bdd_node.high_outgoing))
        {
            if (bdd_node.prev_high_incoming == nullptr)
                bdd_node.high_outgoing->first_high_incoming = bdd_node.next_high_incoming;
            else
                bdd_node.prev_high_incoming->next_high_incoming = bdd_node.next_high_incoming;
            if (bdd_node.next_high_incoming != nullptr)
                bdd_node.next_high_incoming->prev_high_incoming = bdd_node.prev_high_incoming;
            bdd_node.prev_high_incoming = nullptr;
            bdd_node.next_high_incoming = nullptr;
            if (bdd_node.high_outgoing->is_first())
                remove_all_outgoing_arcs(*bdd_node.high_outgoing);
        }
        if (bdd_node.high_outgoing != bdd_branch_node_fix::terminal_0())
        {
            const log_entry entry(&bdd_node, bdd_node.high_outgoing, true);
            log_.push(entry);
            bdd_node.high_outgoing = bdd_branch_node_fix::terminal_0();
            bdd_node.bdd_var->nr_feasible_high_arcs--;
        }
    }

    bool bdd_fixing::fix_variables(const std::vector<size_t> & variables, const std::vector<char> & values)
    {
        assert(primal_solution_.size() == nr_variables());
        assert(variables.size() == values.size());

        std::fill(primal_solution_.begin(), primal_solution_.end(), 2);

        struct var_fix
        {
            var_fix(const size_t log_size, const size_t index, const char val)
            : log_size_(log_size), index_(index), val_(val) {}
            
            const size_t log_size_;
            const size_t index_;
            const char val_;
        };

        std::stack<var_fix, std::vector<var_fix>> variable_fixes;
        variable_fixes.emplace(log_.size(), 0, 1-values[0]);
        variable_fixes.emplace(log_.size(), 0, values[0]);

        size_t nfixes = 0;
        size_t max_fixes = nr_variables();
        std::cout << "\nExpanded " << nfixes << " out of " << max_fixes << " search tree nodes.." << std::flush;

        while (!variable_fixes.empty())
        {
            nfixes++;
            std::cout << "\rExpanded " << nfixes << " out of " << max_fixes << " search tree nodes.." << std::flush;
            if (nfixes > max_fixes)
                return false;

            auto fix = variable_fixes.top();
            variable_fixes.pop();
            size_t index = fix.index_;

            revert_changes(fix.log_size_);
            bool feasible = fix_variable(variables[index], fix.val_);

            if (!feasible)
                continue;

            while (is_fixed(variables[index]))
            {
                index++;
                if (index >= variables.size())
                    return true;
            }

            variable_fixes.emplace(log_.size(), index, 1-values[index]);
            variable_fixes.emplace(log_.size(), index, values[index]);
        }

        return false;
    }

    bool bdd_fixing::is_fixed(const size_t var) const
    {
        assert(primal_solution_.size() == nr_variables());
        assert(var < primal_solution_.size());
        return primal_solution_[var] < 2;
    }

    void bdd_fixing::revert_changes(const size_t target_log_size)
    {
        while (log_.size() > target_log_size)
        {
            log_.top().restore();
            log_.pop();
        }
    }



    ////////////////////////////////////////////////////
    // Solver with Final Rounding
    ////////////////////////////////////////////////////

    class bdd_opt : public bdd_solver_interface {
        public:

            bdd_opt() {}
            bdd_opt(const bdd_opt&) = delete; // no copy constructor because of pointers in bdd_branch_node

            void init(const ILP_input& input);
            double lower_bound() { return bdd_mma_.lower_bound(); }
            void iteration() { bdd_mma_.iteration(); }
            void mma_iteration() { bdd_mma_.min_marginal_averaging_iteration(); }
            void srmp_iteration() { bdd_mma_.min_marginal_averaging_iteration_SRMP(); }

            bool fix_variables();
            double compute_upper_bound() { return bdd_mma_.compute_upper_bound(bdd_fix_.primal_solution()); }
            std::vector<char> primal_solution() const { return bdd_fix_.primal_solution(); }

            void set_options(const bdd_min_marginal_averaging_options o) { bdd_mma_.set_options(o); }

        private:
            
            bdd_min_marginal_averaging bdd_mma_;
            bdd_fixing bdd_fix_;
    };

    void bdd_opt::init(const ILP_input& input)
    {
        bdd_mma_.init(input);
        bdd_mma_.compute_lower_bound();
        bdd_fix_.init(input);
    }


    bool bdd_opt::fix_variables()
    {
        mma_iteration();
        std::vector<double> total_min_marginals = bdd_mma_.total_min_marginals();
        std::vector<size_t> variables;
        for (size_t i = 0; i < bdd_mma_.nr_variables(); i++)
            variables.push_back(i);

        // TODO change to more sensible ordering
        auto order = [&](const size_t a, const size_t b)
        {
            if (total_min_marginals[a] >= 0.0 && total_min_marginals[b] >= 0.0)
                return (total_min_marginals[a]) > (total_min_marginals[b]);
            return (total_min_marginals[a]) < (total_min_marginals[b]);
        };
        std::sort(variables.begin(), variables.end(), order);

        std::vector<char> values;
        for (size_t i = 0; i < variables.size(); i++)
        {
            const char val = (total_min_marginals[variables[i]] < 0) ? 1 : 0;
            values.push_back(val);
        }

        return bdd_fix_.fix_variables(variables, values);
    }

}
