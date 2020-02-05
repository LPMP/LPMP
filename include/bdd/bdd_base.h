#pragma once

#include "two_dimensional_variable_array.hxx"
#include "cuddObj.hh"
#include "ILP_input.h"
#include "bdd_branch_node.h"
#include <cassert>
#include <vector>
#include <unordered_map>
#include <tsl/robin_map.h>
#include <numeric>
#include <chrono> // for now
#include "bdd_storage.h"
#include "tclap/CmdLine.h"

namespace LPMP {

    template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    class bdd_base {
        public:
            bdd_base() {}
            bdd_base(const bdd_base&) = delete; // no copy constructor because of pointers in bdd_branch_node

            template<typename BDD_VARIABLES_ITERATOR>
                void add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end, Cudd& bdd_mgr);

            void init(const ILP_input& input);

            std::size_t nr_variables() const { return bdd_variables.size(); }
            std::size_t nr_bdds() const { return bdd_variables.size()-1; }
            std::size_t nr_bdds(const std::size_t var) const { assert(var<nr_variables()); return bdd_variables[var].size(); }

        private:
            void init_branch_nodes();

            std::size_t bdd_branch_node_index(const BDD_BRANCH_NODE* bdd) const;
            std::size_t bdd_branch_node_index(const BDD_BRANCH_NODE& bdd) const { return bdd_branch_node_index(&bdd); }
            std::size_t variable_index(const BDD_BRANCH_NODE& bdd) const;
            std::size_t variable_index(const BDD_VARIABLE& bdd_var) const;

            std::size_t first_variable_of_bdd(const BDD_VARIABLE& bdd_var) const;
            std::size_t last_variable_of_bdd(const BDD_VARIABLE& bdd_var) const;
            bool first_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const;
            bool last_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const;

            void check_bdd_variable(const bdd_variable& bdd_var, const bool last_variable = false, const bool first_variable = false) const;

            std::vector<BDD_BRANCH_NODE> bdd_branch_nodes;
            two_dim_variable_array<BDD_VARIABLE> bdd_variables;
            bdd_storage bdd_storage_;
    };


    template<typename BDD_VARIABLES_ITERATOR>
        void bdd_base::add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end, Cudd& bdd_mgr)
        {
            bdd_storage_.add_bdd(bdd_mgr, bdd, bdd_vars_begin, bdd_vars_end); 
        }

    void bdd_base::init()
    {
        init_branch_nodes();
        costs_.resize(nr_variables(), std::numeric_limits<double>::infinity());
    }

    void bdd_base::init(const ILP_input& input)
    {
        Cudd bdd_mgr;
        bdd_storage_ = bdd_storage(input, bdd_mgr);
        init();
        set_costs(input.objective().begin(), input.objective().end());
    }

    void bdd_base::init_branch_nodes()
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

        bdd_branch_nodes.resize(bdd_storage_.bdd_nodes().size());
        bdd_variables.resize(nr_bdds_per_variable.begin(), nr_bdds_per_variable.end());

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
            stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, BDD_BRANCH_NODE*>(bdd_storage::bdd_node::terminal_0, bdd_branch_node_terminal_0));
            stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, BDD_BRANCH_NODE*>(bdd_storage::bdd_node::terminal_1, bdd_branch_node_terminal_1)); 
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

                auto& bdd_var = bdd_variables(v, nr_bdds_per_variable[v]);
                if(!variable_covered(v)) {
                    assert(bdd_var.is_initial_state());
                    cover_variable(v);

                    bdd_var.first_branch_instruction = bdd_offset_per_variable[v];
                    bdd_var.last_branch_instruction = bdd_offset_per_variable[v];
                    //std::cout << "bdd level offset for var = " << v << ", bdd index = " << bdd_index << " = " << stored_bdd_node_index << "\n";
                    // TODO: remove, not needed anymore
                    if(next_bdd_var != nullptr) {
                        //assert(next_bdd_var > &bdd_var);
                        bdd_var.next = next_bdd_var;
                        next_bdd_var->prev = &bdd_var;
                    }

                    next_bdd_var = &bdd_var;
                } else {
                    assert(!bdd_var.is_initial_state());
                }

                const std::size_t bdd_branch_nodes_index = bdd_var.last_branch_instruction; 
                bdd_var.last_branch_instruction++;

                BDD_BRANCH_NODE& bdd = bdd_branch_nodes[bdd_branch_nodes_index];
                assert(bdd.is_initial_state());
                //std::cout << "address = " << &bdd << "\n";
                bdd.variable_cost = &bdd_var.variable_cost;

                stored_bdd_node_index_to_bdd_address.insert({stored_bdd_node_index, &bdd});

                assert(stored_bdd_node_index_to_bdd_address.count(stored_bdd.low) > 0);
                BDD_BRANCH_NODE& bdd_low = *(stored_bdd_node_index_to_bdd_address.find(stored_bdd.low)->second);
                bdd.low_outgoing = &bdd_low;
                assert(bdd.low_outgoing != nullptr);
                if(!bdd_low.is_terminal()) {
                    bdd.next_low_incoming = bdd_low.first_low_incoming;
                    bdd_low.first_low_incoming = &bdd;
                }

                assert(stored_bdd_node_index_to_bdd_address.count(stored_bdd.high) > 0);
                BDD_BRANCH_NODE& bdd_high = *(stored_bdd_node_index_to_bdd_address.find(stored_bdd.high)->second);
                bdd.high_outgoing = &bdd_high;
                if(!bdd_high.is_terminal()) {
                    bdd.next_high_incoming = bdd_high.first_high_incoming;
                    bdd_high.first_high_incoming = &bdd;
                }

                //if(!bdd.low_outgoing->is_terminal()) { assert(variable_index(bdd) < variable_index(*bdd.low_outgoing)); }
                //if(!bdd.high_outgoing->is_terminal()) { assert(variable_index(bdd) < variable_index(*bdd.high_outgoing)); }
                check_bdd_branch_node(bdd, v+1 == nr_variables(), v == 0);
                check_bdd_variable(bdd_var, v+1 == nr_variables(), v == 0);
                ++nr_bdd_nodes_per_variable[v]; 
                ++bdd_offset_per_variable[v];
            }

            for(const std::size_t v : variables_covered) {
                ++nr_bdds_per_variable[v];
            }
        }

        for(std::size_t v=0; v<nr_bdds_per_variable.size(); ++v) {
            assert(nr_bdds_per_variable[v] == bdd_variables[v].size());
        }

        for(const auto& bdd : bdd_branch_nodes) {
            check_bdd_branch_node(bdd);
        }

        // TODO: clear bdd_storage
    }

    std::size_t bdd_base::bdd_branch_node_index(const BDD_BRANCH_NODE* bdd) const
    {
        assert(bdd >= &bdd_branch_nodes[0]);
        const std::size_t i = bdd - &bdd_branch_nodes[0];
        assert(i < bdd_branch_nodes.size());
        return i; 
    }

    std::size_t bdd_base::variable_index(const BDD_BRANCH_NODE& bdd) const
    {
        const std::size_t i = bdd_branch_node_index(&bdd);

        std::size_t lb = 0;
        std::size_t ub = nr_variables()-1;
        std::size_t v = nr_variables()/2;
        while (! (i >= bdd_variables(v, 0).first_branch_instruction && i < bdd_variables[v].back().last_branch_instruction))
        {
            if (i > bdd_variables(v, 0).first_branch_instruction)
                lb = v+1;
            else
                ub = v-1;
            v = (lb+ub)/2;
        }
        return v;
    }

    std::size_t bdd_base::variable_index(const BDD_VARIABLE& bdd_var) const
    {
        assert(&bdd_var >= &bdd_variables(0,0));
        assert(&bdd_var <= &bdd_variables.back().back());
        std::size_t lb = 0;
        std::size_t ub = nr_variables()-1;
        std::size_t v = nr_variables()/2;
        while(! (&bdd_var >= &bdd_variables(v,0) && &bdd_var <= &bdd_variables[v].back()) ) {
            if( &bdd_var > &bdd_variables(v,0) ) {
                lb = v+1;
            } else {
                ub = v-1;
            }
            v = (ub+lb)/2; 
        }
        assert(&bdd_var >= &bdd_variables(v,0) && &bdd_var <= &bdd_variables[v].back());
        return v; 
    }

    std::size_t bdd_base::first_variable_of_bdd(const BDD_VARIABLE& bdd_var) const
    {
        BDD_VARIABLE const* p = &bdd_var;
        while(p->prev != nullptr)
            p = p->prev;
        return variable_index(*p);
    }

    std::size_t bdd_base::last_variable_of_bdd(const BDD_VARIABLE& bdd_var) const
    {
        BDD_VARIABLE const* p = &bdd_var;
        while(p->next != nullptr)
            p = p->next;
        return variable_index(*p); 
    }

    bool bdd_base::last_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const
    {
        const BDD_VARIABLE& bdd_var = bdd_variables(var, bdd_index);
        return bdd_var.next == nullptr; 
    }

    bool bdd_base::first_variable_of_bdd(const std::size_t var, const std::size_t bdd_index) const
    {
        const BDD_VARIABLE& bdd_var = bdd_variables(var, bdd_index);
        return bdd_var.prev == nullptr; 
    }

    void bdd_base::check_bdd_variable(const bdd_variable& bdd_var, const bool last_variable, const bool first_variable) const
    {
        // go over all branch instructions and check whether they point to all the same variable cost
        for(std::size_t bdd_node_index = bdd_var.first_branch_instruction; bdd_node_index < bdd_var.last_branch_instruction; ++bdd_node_index) {
            assert(&bdd_var.variable_cost == bdd_branch_instructions[bdd_node_index].variable_cost);
        }

        if(last_variable) {
            assert(bdd_var.next == nullptr);
        }
        if(first_variable) {
            assert(bdd_var.prev == nullptr);
        }
    }
}
