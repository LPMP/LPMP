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

            void forward_step(const std::size_t var, const std::size_t bdd_index);
            void backward_step(const std::size_t var, const std::size_t bdd_index);
            void forward_step(const std::size_t var);
            void backward_step(const std::size_t var);

            void backward_run(); // also used to initialize
            void forward_run();

            const BDD_VARIABLE &get_bdd_variable(const std::size_t var, const std::size_t bdd_index) const;
            const BDD_BRANCH_NODE &get_bdd_branch_node(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const;

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
        // last parameter for recording primal fixations
        bdd_storage_ = bdd_storage(input, bdd_mgr, false);
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
                //check_bdd_branch_node(bdd, v+1 == nr_variables(), v == 0); // TODO: cannot be enabled because not all pointers are set yet.
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

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::forward_step(const std::size_t var, const std::size_t bdd_index)
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

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::forward_step(const std::size_t var)
    {
        assert(var < bdd_variables_.size());
        for(std::size_t bdd_index = 0; bdd_index<bdd_variables_[var].size(); ++bdd_index)
            forward_step(var, bdd_index);
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::backward_step(const std::size_t var, const std::size_t bdd_index)
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

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::backward_step(const std::size_t var)
    {
        assert(var < bdd_variables_.size());
        for(std::size_t bdd_index = 0; bdd_index<bdd_variables_[var].size(); ++bdd_index)
            backward_step(var, bdd_index);
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::backward_run()
    {
        for(std::ptrdiff_t var=nr_variables()-1; var>=0; --var)
            backward_step(var); 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::forward_run()
    {
        for(std::size_t var=0; var<nr_variables(); ++var)
            forward_step(var); 
    }

    template <typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    const BDD_VARIABLE &bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::get_bdd_variable(const std::size_t var, const std::size_t bdd_index) const
    {
        assert(var < nr_variables());
        assert(bdd_index < bdd_variables_[var].size());
        return bdd_variables_(var, bdd_index);
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    const BDD_BRANCH_NODE& bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::get_bdd_branch_node(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const
    {
        assert(var < nr_variables());
        assert(bdd_index < bdd_variables_[var].size());
        assert(bdd_node_index < bdd_variables_(var,bdd_index).nr_bdd_nodes());
        return bdd_branch_nodes_[ bdd_variables_(var,bdd_index).first_node_index + bdd_node_index ];
    }

    ////////////////////////////////////////////////////
    // Min-Marginal Averaging Solver
    ////////////////////////////////////////////////////

    struct bdd_min_marginal_averaging_options {
        bdd_min_marginal_averaging_options(int argc, char** argv);
        bdd_min_marginal_averaging_options() {}
        enum class averaging_type {classic, SRMP} averaging_type = averaging_type::classic;
        enum class variable_order {input, bfs, cuthill, mindegree} variable_order = variable_order::input;
    };

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    class bdd_mma_base;

    typedef bdd_mma_base<bdd_variable_mma, bdd_branch_node_opt> bdd_min_marginal_averaging;

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    class bdd_mma_base : public bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>, public bdd_solver_interface
    {
    public:
        bdd_mma_base() {}
        bdd_mma_base(const bdd_mma_base &) = delete; // no copy constructor because of pointers in bdd_branch_node

        void init(const ILP_input &input);
        void init();

        template <typename ITERATOR>
        void set_costs(ITERATOR begin, ITERATOR end);

        template <typename ITERATOR>
        bool check_feasibility(ITERATOR var_begin, ITERATOR var_end) const;
        template <typename ITERATOR>
        double evaluate(ITERATOR var_begin, ITERATOR var_end) const;

        double lower_bound() { return lower_bound_; }
        double compute_lower_bound();
        double lower_bound_backward(const std::size_t var, const std::size_t bdd_index);
        double lower_bound_forward(const std::size_t var, const std::size_t bdd_index);

        double compute_upper_bound(const std::vector<char> &primal_solution) const;

        void min_marginal_averaging_iteration();
        void min_marginal_averaging_forward();
        void min_marginal_averaging_backward();
        void min_marginal_averaging_forward_SRMP();
        void min_marginal_averaging_backward_SRMP();
        void min_marginal_averaging_iteration_SRMP();

        void iteration();

        template <typename STREAM>
        void export_dot(STREAM &s) const;

        void set_options(const bdd_min_marginal_averaging_options o) { options = o; }

    protected:
        void init_costs();

        std::array<double, 2> min_marginal(const std::size_t var, const std::size_t bdd_index) const;
        template <typename ITERATOR>
        static std::array<double, 2> average_marginals(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t nr_marginals_to_distribute = std::numeric_limits<std::size_t>::max());
        template <typename ITERATOR>
        std::pair<std::array<double, 2>, bool> average_marginals_forward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const;
        template <typename ITERATOR>
        std::pair<std::array<double, 2>, bool> average_marginals_backward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const;
        void set_marginal(const std::size_t var, const std::size_t bdd_index, const std::array<double, 2> marginals, const std::array<double, 2> min_marginals);
        void set_marginal_forward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double, 2> marginals, const std::array<double, 2> min_marginals, const bool default_avg);
        void set_marginal_backward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double, 2> marginals, const std::array<double, 2> min_marginals, const bool default_avg);

        std::vector<double> costs_;
        double lower_bound_ = -std::numeric_limits<double>::infinity();
        bdd_min_marginal_averaging_options options;
    };


   bdd_min_marginal_averaging_options::bdd_min_marginal_averaging_options(int argc, char** argv)
        {
            TCLAP::CmdLine cmd("Command line parser for bdd min marginal averaging options", ' ', " ");
            TCLAP::ValueArg<std::string> averaging_arg("a","averaging","averaging type",false,"classic","{classic|SRMP}");
            TCLAP::ValueArg<std::string> order_arg("o","order","variable order",false,"input","{input|bfs|cuthill|mindegree}");
            cmd.add(averaging_arg);
            cmd.add(order_arg);

            cmd.parse(argc, argv);

            if(averaging_arg.getValue() == "classic")
                averaging_type = bdd_min_marginal_averaging_options::averaging_type::classic;
            else if(averaging_arg.getValue() == "SRMP")
                averaging_type = bdd_min_marginal_averaging_options::averaging_type::SRMP;
            else
                throw std::runtime_error("averaging type not recognized");

            if(order_arg.getValue() == "input")
                variable_order = bdd_min_marginal_averaging_options::variable_order::input;
            else if(order_arg.getValue() == "bfs")
                variable_order = bdd_min_marginal_averaging_options::variable_order::bfs;
            else if(order_arg.getValue() == "cuthill")
                variable_order = bdd_min_marginal_averaging_options::variable_order::cuthill;
            else if(order_arg.getValue() == "mindegree")
                variable_order = bdd_min_marginal_averaging_options::variable_order::mindegree;
            else
                throw std::runtime_error("variable order not recognized");
        }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init_costs()
    {
        for (size_t var = 0; var < this->nr_variables(); var++)
        for (size_t bdd_index = 0; bdd_index < this->nr_bdds(var); bdd_index++)
        {
            auto & bdd_var = this->bdd_variables_(var, bdd_index);
            for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                this->bdd_branch_nodes_[node_index].variable_cost = & bdd_var.cost;
        }
        costs_.resize(this->nr_variables(), std::numeric_limits<double>::infinity());
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init(const ILP_input& input)
    {
        bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init(input);
        init_costs();
        set_costs(input.objective().begin(), input.objective().end());
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init()
    {
        bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init();
        init_costs();
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    std::array<double,2> bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal(const std::size_t var, const std::size_t bdd_index) const
    {
        assert(var < this->nr_variables());
        assert(bdd_index < this->nr_bdds(var));
        std::array<double,2> m = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        const auto& bdd_var = this->bdd_variables_(var,bdd_index);
        for(std::size_t bdd_node_index=bdd_var.first_node_index; bdd_node_index<bdd_var.last_node_index; ++bdd_node_index) {
            //std::cout << "min marginal for var = " << var << ", bdd_index = " << bdd_index << ", bdd_node_index = " << bdd_node_index << "\n";
            const auto& bdd = this->bdd_branch_nodes_[bdd_node_index];
            const auto [m0,m1] = bdd.min_marginal();
            m[0] = std::min(m[0], m0);
            m[1] = std::min(m[1], m1);
        }
        assert(std::isfinite(m[0]));
        assert(std::isfinite(m[1]));
        return m;
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::set_marginal(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals)
    {
        assert(var < this->nr_variables());
        assert(bdd_index < this->nr_bdds(var));
        assert(min_marginals == min_marginal(var,bdd_index));
        assert(std::isfinite(marginals[0]) && std::isfinite(marginals[1]));
        assert(std::isfinite(min_marginals[0]) && std::isfinite(min_marginals[1]));

        auto& bdd_var = this->bdd_variables_(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        assert(std::isfinite(marginal_diff));
        assert(std::isfinite(marginal_diff_target));
        bdd_var.cost += -marginal_diff + marginal_diff_target; 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::set_marginal_forward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals, const bool default_avg)
    {
        auto& bdd_var = this->bdd_variables_(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        assert(std::isfinite(marginal_diff));
        assert(std::isfinite(marginal_diff_target));
        if (default_avg)
        {
            bdd_var.cost += -marginal_diff + marginal_diff_target;
        }
        else if(this->last_variable_of_bdd(var, bdd_index)) {
            bdd_var.cost -= marginal_diff;
        } else {
            assert(std::isfinite(marginal_diff_target));
            bdd_var.cost += -marginal_diff + marginal_diff_target; 
        }

    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::set_marginal_backward_SRMP(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals, const bool default_avg)
    {
        auto& bdd_var = this->bdd_variables_(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        assert(std::isfinite(marginal_diff));
        assert(std::isfinite(marginal_diff_target));
        if (default_avg)
        {
            bdd_var.cost += -marginal_diff + marginal_diff_target;
        }
        else if(this->first_variable_of_bdd(var, bdd_index)) {
            bdd_var.cost -= marginal_diff;
        } else {
            assert(std::isfinite(marginal_diff_target));
            bdd_var.cost += -marginal_diff + marginal_diff_target; 
        }
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    double bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::compute_lower_bound()
    {
        double lb = 0.0;
        for(long int var=this->nr_variables()-1; var>=0; --var)
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index)
                lb += lower_bound_backward(var,bdd_index);
        this->lower_bound_ = lb;
        return lb;
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    double bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::lower_bound_backward(const std::size_t var, const std::size_t bdd_index)
    {
        const auto& bdd_var = this->bdd_variables_(var,bdd_index);
        if(this->first_variable_of_bdd(var, bdd_index)) {
            assert(bdd_var.nr_bdd_nodes() == 1);
            const auto& bdd = this->get_bdd_branch_node(var, bdd_index, 0);
            return bdd.m;
        } else {
            return 0.0;
        }
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    double bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::lower_bound_forward(const std::size_t var, const std::size_t bdd_index)
    {
        const auto& bdd_var = this->bdd_variables_(var,bdd_index);
        if(this->last_variable_of_bdd(var, bdd_index)) 
        {
            double lb = std::numeric_limits<double>::infinity();
            const std::size_t first_node_index = bdd_var.first_node_index;
            const std::size_t last_node_index = bdd_var.last_node_index;
            for (std::size_t i = first_node_index; i < last_node_index; ++i)
            {
                const auto& node = this->bdd_branch_nodes_[i];
                assert(node.low_outgoing == node.terminal_1() || node.high_outgoing == node.terminal_1());
                const double node_lb = [&]() {
                    double lb = std::numeric_limits<double>::infinity();
                    if(node.low_outgoing == node.terminal_1())
                        lb = std::min(lb, node.m);
                    if(node.high_outgoing == node.terminal_1())
                        lb = std::min(lb, node.m + *node.variable_cost);
                    return lb;
                }();
                lb = std::min(lb, node_lb); 
            }
            return lb;
        }
        else
        {
            return 0.0;
        }
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    double bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::compute_upper_bound(const std::vector<char> & primal_solution) const
    {
        if (primal_solution.size() != this->nr_variables())
            return std::numeric_limits<double>::infinity();
        if (!check_feasibility(primal_solution.begin(), primal_solution.end()))
            return std::numeric_limits<double>::infinity();
        else
            return evaluate(primal_solution.begin(), primal_solution.end());
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal_averaging_iteration()
    {
        const auto begin_time = std::chrono::steady_clock::now();
        min_marginal_averaging_forward();
        const auto after_forward = std::chrono::steady_clock::now();
        // std::cout << "forward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(after_forward - begin_time).count() << " ms, " << std::flush;
        const auto before_backward = std::chrono::steady_clock::now();
        min_marginal_averaging_backward();
        const auto end_time = std::chrono::steady_clock::now();
        // std::cout << "backward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - before_backward).count() << " ms, " << std::flush;
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename ITERATOR>
        std::array<double,2> bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::average_marginals(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t nr_marginals_to_distribute)
        {
            std::array<double,2> average_marginal = {0.0, 0.0};
            for(auto m_iter=marginals_begin; m_iter!=marginals_end; ++m_iter) {
                average_marginal[0] += (*m_iter)[0];
                average_marginal[1] += (*m_iter)[1];
            }
            if(nr_marginals_to_distribute == std::numeric_limits<std::size_t>::max()) 
            {
                const double no_marginals = std::distance(marginals_begin, marginals_end);
                average_marginal[0] /= no_marginals;
                average_marginal[1] /= no_marginals;
            } 
            else
            {
                average_marginal[0] /= nr_marginals_to_distribute;
                average_marginal[1] /= nr_marginals_to_distribute; 
            }

            assert(std::isfinite(average_marginal[0]));
            assert(std::isfinite(average_marginal[1]));
            return average_marginal;
        }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename ITERATOR>
        std::pair<std::array<double,2>, bool> bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::average_marginals_forward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const
        {
            assert(this->nr_bdds(var) == std::distance(marginals_begin, marginals_end));
            std::array<double,2> average_marginal = {0.0, 0.0};
            std::size_t nr_averaged_marginals = 0;
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                    average_marginal[0] += (*(marginals_begin+bdd_index))[0];
                    average_marginal[1] += (*(marginals_begin+bdd_index))[1];
                if(!this->last_variable_of_bdd(var, bdd_index))
                    ++nr_averaged_marginals;
            }
            // if no BDD satisfies forward condition, resort to averaging over all BDDs
            bool default_avg = false;
            if (nr_averaged_marginals == 0)
            {
                nr_averaged_marginals = this->nr_bdds(var);
                default_avg = true;
            }

            average_marginal[0] /= double(nr_averaged_marginals);
            average_marginal[1] /= double(nr_averaged_marginals);

            return std::make_pair(average_marginal, default_avg);
        }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename ITERATOR>
        std::pair<std::array<double,2>, bool> bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::average_marginals_backward_SRMP(ITERATOR marginals_begin, ITERATOR marginals_end, const std::size_t var) const
        {
            assert(this->nr_bdds(var) == std::distance(marginals_begin, marginals_end));
            std::array<double,2> average_marginal = {0.0, 0.0};
            std::size_t nr_averaged_marginals = 0;
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                    average_marginal[0] += (*(marginals_begin+bdd_index))[0];
                    average_marginal[1] += (*(marginals_begin+bdd_index))[1];
                if(!this->first_variable_of_bdd(var, bdd_index))
                    ++nr_averaged_marginals;
            }
            // if no BDD satisfies forward condition, resort to averaging over all BDDs
            bool default_avg = false;
            if (nr_averaged_marginals == 0)
            {
                nr_averaged_marginals = this->nr_bdds(var);
                default_avg = true;
            }

            average_marginal[0] /= double(nr_averaged_marginals);
            average_marginal[1] /= double(nr_averaged_marginals);

            assert(std::isfinite(average_marginal[0]));
            assert(std::isfinite(average_marginal[1]));
            return std::make_pair(average_marginal, default_avg);
        }

    // min marginal averaging
    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal_averaging_forward()
    {
        std::vector<std::array<double,2>> min_marginals;
        for(std::size_t var=0; var<this->nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                this->forward_step(var,bdd_index);
                min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
            } 
        }
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal_averaging_forward_SRMP()
    {
        std::vector<std::array<double,2>> min_marginals;
        for(std::size_t var=0; var<this->nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                this->forward_step(var,bdd_index);
                min_marginals.push_back(min_marginal(var,bdd_index));
            }


            const auto average_marginal = average_marginals_forward_SRMP(min_marginals.begin(), min_marginals.end(), var);
            const std::array<double,2> avg_marg = average_marginal.first;
            const bool default_averaging = average_marginal.second;

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                set_marginal_forward_SRMP(var,bdd_index,avg_marg,min_marginals[bdd_index], default_averaging);
            } 
        }
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal_averaging_backward()
    {
        double lb = 0.0;
        std::vector<std::array<double,2>> min_marginals;

        for(long int var=this->nr_variables()-1; var>=0; --var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
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
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
                this->backward_step(var, bdd_index);
                lb += lower_bound_backward(var,bdd_index);
            }
        }

        lower_bound_ = lb; 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal_averaging_backward_SRMP()
    {
        double lb = 0.0;
        std::vector<std::array<double,2>> min_marginals;
        for(long int var=this->nr_variables()-1; var>=0; --var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            const auto average_marginal = average_marginals_backward_SRMP(min_marginals.begin(), min_marginals.end(), var);
            const std::array<double,2> avg_marg = average_marginal.first;
            const bool default_averaging = average_marginal.second;

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                set_marginal_backward_SRMP(var,bdd_index,avg_marg,min_marginals[bdd_index], default_averaging);
                this->backward_step(var, bdd_index);
                lb += lower_bound_backward(var,bdd_index);
            }
        }

        lower_bound_ = lb; 
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::min_marginal_averaging_iteration_SRMP()
    {
        const auto begin_time = std::chrono::steady_clock::now();
        min_marginal_averaging_forward_SRMP();
        const auto after_forward = std::chrono::steady_clock::now();
        // std::cout << "forward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(after_forward - begin_time).count() << " ms, " << std::flush;
        const auto before_backward = std::chrono::steady_clock::now();
        min_marginal_averaging_backward_SRMP();
        const auto end_time = std::chrono::steady_clock::now();
        // std::cout << "backward " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - before_backward).count() << " ms, " << std::flush;
    }



    //const bdd_branch_instruction& bdd_min_marginal_averaging::get_bdd_branch_instruction(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const
   
    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::iteration()
    {
        if(options.averaging_type == bdd_min_marginal_averaging_options::averaging_type::classic)
            min_marginal_averaging_iteration();
        else if(options.averaging_type == bdd_min_marginal_averaging_options::averaging_type::SRMP)
            min_marginal_averaging_iteration_SRMP();
        else
            assert(false);
    }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename ITERATOR>
        void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::set_costs(ITERATOR begin, ITERATOR end)
        {
            // TODO: remove costs_ array
            std::fill(costs_.begin(), costs_.end(), 0.0);
            assert(std::distance(begin,end) <= this->nr_variables());
            std::copy(begin, end, costs_.begin());
            //std::fill(costs_.begin() + std::distance(begin, end), costs_.end(), 0.0);

            // distribute costs to bdds uniformly
            for(std::size_t v=0; v<this->nr_variables(); ++v) {
                assert(this->bdd_variables_[v].size() > 0);
                for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(v); ++bdd_index) {
                    assert(this->nr_bdds(v) > 0);
                    const double cost = v < std::distance(begin,end) ? *(begin+v)/this->nr_bdds(v) : 0.0; 
                    this->bdd_variables_(v,bdd_index).cost = cost;
                    assert(!std::isnan(this->bdd_variables_(v,bdd_index).cost));
                }
            }
            for(const auto& bdd : this->bdd_branch_nodes_) {
                assert(!std::isnan(*bdd.variable_cost));
            }
            this->backward_run();
        }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename ITERATOR>
        bool bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::check_feasibility(ITERATOR var_begin, ITERATOR var_end) const
        {
            assert(std::distance(var_begin, var_end) == this->nr_variables());

            std::vector<char> bdd_nbranch_node_marks(this->bdd_branch_nodes_.size(), 0);

            std::size_t var = 0;
            for(auto var_iter=var_begin; var_iter!=var_end; ++var_iter, ++var) {
                const bool val = *(var_begin+var);
                for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                    const auto& bdd_var = this->bdd_variables_(var, bdd_index);
                    if(bdd_var.is_first_bdd_variable()) {
                        bdd_nbranch_node_marks[bdd_var.first_node_index] = 1;
                    }
                    for(std::size_t bdd_node_index=bdd_var.first_node_index; bdd_node_index<bdd_var.last_node_index; ++bdd_node_index) {
                        if(bdd_nbranch_node_marks[bdd_node_index] == 1) {
                            const auto& bdd = this->bdd_branch_nodes_[bdd_node_index];
                            const auto* bdd_next_index = [&]() {
                                if(val == false)
                                    return bdd.low_outgoing;
                                else 
                                    return bdd.high_outgoing;
                            }();

                            if(bdd_next_index == BDD_BRANCH_NODE::terminal_0())
                                return false;
                            if(bdd_next_index == BDD_BRANCH_NODE::terminal_1()) {
                            } else { 
                                bdd_nbranch_node_marks[ this->bdd_branch_node_index(bdd_next_index) ] = 1;
                            }
                        }
                    }
                }
            }

            return true;
        }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename ITERATOR>
        double bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::evaluate(ITERATOR var_begin, ITERATOR var_end) const
        {
            assert(std::distance(var_begin, var_end) == this->nr_variables());

            if(!check_feasibility(var_begin, var_end))
                return std::numeric_limits<double>::infinity();

            double cost = 0.0;
            std::size_t var = 0;
            for(auto var_iter=var_begin; var_iter!=var_end; ++var_iter, ++var)
                cost += *var_iter * costs_[var];
            return cost;
        }

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    template<typename STREAM>
        void bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>::export_dot(STREAM& s) const
        {
            return this->bdd_storage_.export_dot(s);
            /*
            s << "digraph bdd_min_marginal_averaging {\n";

            std::vector<char> bdd_node_visited(this->bdd_branch_nodes_.size(), false);
            std::size_t cur_bdd_index = 0;
            for(const auto& bdd : this->bdd_branch_nodes_) {
                const std::size_t i = bdd_branch_node_index(bdd);
                if(bdd_node_visited[i])
                    continue;
                bdd_node_visited[i] = true;
                auto get_node_string = [&](const bdd_branch_node_opt* bdd) -> std::string {
                    if(bdd == bdd_branch_node_opt::terminal_0())
                        return "false";
                    if(bdd == bdd_branch_node_opt::terminal_1())
                        return "true";
                    return std::to_string(this->bdd_branch_node_index(bdd));
                };
                s << i << " -> " << get_node_string(bdd.low_outgoing) << " [label=\"0\"];\n";
                s << i << " -> " << get_node_string(bdd.high_outgoing) << " [label=\"1\"];\n";
            }

            s << "}\n"; 
            */
        }


    ////////////////////////////////////////////////////


    // helper functions
    template <typename BDD_BASE, typename BDD_VARIABLE_FUNC, typename BDD_INDEX_FUNC, typename VAR_ITERATOR>
    std::tuple<std::vector<std::size_t>, two_dim_variable_array<std::size_t>> compute_bdd_mask(
        VAR_ITERATOR variable_begin, VAR_ITERATOR variable_end,
        BDD_BASE& bdd_base, BDD_VARIABLE_FUNC bdd_variable_func, BDD_INDEX_FUNC bdd_index_func
        )
    {
        std::unordered_set<std::array<std::size_t, 2>> affected_bdd_variables_set;
        for (auto variable_it = variable_begin; variable_it != variable_end; ++variable_it)
        {
            const std::size_t var = *variable_it;
            for (std::size_t bdd_index = 0; bdd_index < bdd_base.nr_bdds(var); ++bdd_index)
            {
                affected_bdd_variables_set.insert({var, bdd_index});
                const auto *prev_bdd_variable = bdd_base.get_bdd_variable(var, bdd_index).prev;
                while (prev_bdd_variable != nullptr)
                {
                    const std::size_t var = bdd_variable_func(prev_bdd_variable);
                    const std::size_t bdd_index = bdd_index_func(prev_bdd_variable);
                    affected_bdd_variables_set.insert({var, bdd_index});
                    prev_bdd_variable = prev_bdd_variable->prev;
                }

                const auto *next_bdd_variable = bdd_base.get_bdd_variable(var, bdd_index).next;
                while (next_bdd_variable != nullptr)
                {
                    const std::size_t var = bdd_variable_func(next_bdd_variable);
                    const std::size_t bdd_index = bdd_index_func(next_bdd_variable);
                    affected_bdd_variables_set.insert({var, bdd_index});
                    next_bdd_variable = next_bdd_variable->next;
                }
            }
        }

        std::vector<std::array<std::size_t, 2>> affected_bdd_variables;
        affected_bdd_variables.reserve(affected_bdd_variables_set.size());
        for (const auto x : affected_bdd_variables_set)
            affected_bdd_variables.push_back(x);
        auto bdd_var_compare = [](const auto a, const auto b) { return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end()); };
        std::sort(affected_bdd_variables.begin(), affected_bdd_variables.end(), bdd_var_compare);
        //const auto tmp_it = std::adjacent_find(affected_bdd_variables.begin(), affected_bdd_variables.end(), bdd_var_compare);
        //assert(std::adjacent_find(affected_bdd_variables.begin(), affected_bdd_variables.end(), bdd_var_compare) == affected_bdd_variables.end());
        std::vector<std::size_t> affected_variables;
        std::vector<std::size_t> affected_variables_bdd_index_size;
        for (const auto [var, bdd_index] : affected_bdd_variables)
        {
            if (affected_variables.size() == 0 || affected_variables.back() != var)
            {
                affected_variables.push_back(var);
                affected_variables_bdd_index_size.push_back(0);
            }
            ++affected_variables_bdd_index_size.back();
        }

        two_dim_variable_array<std::size_t> affected_variables_bdd_index(affected_variables_bdd_index_size.begin(), affected_variables_bdd_index_size.end());
        affected_variables.clear();
        std::size_t counter = 0;

        for (const auto [var, bdd_index] : affected_bdd_variables)
        {
            if (affected_variables.size() == 0 || affected_variables.back() != var)
            {
                affected_variables.push_back(var);
                counter = 0;
            }
            const std::size_t var_idx = affected_variables.size() - 1;
            affected_variables_bdd_index(var_idx, counter++) = bdd_index;
        }

        return {affected_variables, affected_variables_bdd_index};
    }
}
