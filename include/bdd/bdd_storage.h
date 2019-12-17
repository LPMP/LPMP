#pragma once

#include "cuddObj.hh"
#include "ILP_input.h"
#include "convert_pb_to_bdd.h"
#include <vector>
#include <unordered_map>
#include <numeric>

namespace LPMP {

    // used for storing all participating BDDs temporarily before exporting them to the bdd solver.
    // TODO: add subsume function, which merges BDDs which have supports that are covered.
    // TODO: subsume dangling BDDs into their endpoints
    class bdd_storage {
        public:

        struct bdd_node {
            constexpr static std::size_t terminal_0 = std::numeric_limits<std::size_t>::max()-1;
            constexpr static std::size_t terminal_1 = std::numeric_limits<std::size_t>::max();

            std::size_t low;
            std::size_t high;
            std::size_t variable;
        };

        bdd_storage() {};
        bdd_storage(const ILP_input& input, Cudd& bdd_mgr);

        void check_node_valid(const bdd_node bdd) const
        {
            assert(bdd.low == bdd_node::terminal_0 || bdd.low == bdd_node::terminal_1 || bdd.low < bdd_nodes_.size());
            assert(bdd.high == bdd_node::terminal_0 || bdd.high == bdd_node::terminal_1 || bdd.high < bdd_nodes_.size());
            assert(bdd.variable < nr_variables());
            if(bdd.low != bdd_node::terminal_0 && bdd.low != bdd_node::terminal_1) {
                assert(bdd.variable < bdd_nodes_[bdd.low].variable);
            }
            if(bdd.high != bdd_node::terminal_0 && bdd.high != bdd_node::terminal_1) {
                assert(bdd.variable < bdd_nodes_[bdd.high].variable);
            }
            //assert(bdd.high != bdd.low); this can be so in ou formulation, but not in ordinary BDDs
        }

        template<typename BDD_VARIABLES_ITERATOR>
            void add_bdd(Cudd& bdd_mgr, BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end);

        template<typename STREAM>
            void export_dot(STREAM& s) const;

        std::size_t nr_bdds() const { return bdd_delimiters().size()-1; }
        std::size_t nr_variables() const { return nr_variables_; }
        std::size_t nr_bdd_nodes(const std::size_t bdd_nr) const { assert(bdd_nr < nr_bdds()); return bdd_delimiters_[bdd_nr+1] - bdd_delimiters_[bdd_nr]; }

        const std::vector<bdd_node>& bdd_nodes() const { return bdd_nodes_; }
        const std::vector<std::size_t>& bdd_delimiters() const { return bdd_delimiters_; }

        std::size_t first_bdd_node(const std::size_t bdd_nr) const;
        std::size_t last_bdd_node(const std::size_t bdd_nr) const;

        private:
        void check_bdd_node(const bdd_node bdd) const;
        std::vector<bdd_node> bdd_nodes_;
        std::vector<std::size_t> bdd_delimiters_ = {0};
        std::vector<std::size_t> nr_bdd_nodes_per_variable;
        std::vector<std::size_t> nr_bdds_per_variable;
        std::size_t nr_variables_ = 0;
    };


    template<typename BDD_VARIABLES_ITERATOR>
        void bdd_storage::add_bdd(Cudd& bdd_mgr, BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end)
        {
            assert(std::is_sorted(bdd_vars_begin, bdd_vars_end));
            assert(std::distance(bdd_vars_begin, bdd_vars_end) > 0);
            nr_variables_ = std::max(nr_variables_, *(bdd_vars_end-1)+1);

            ADD add = bdd.Add();
            std::unordered_map<DdNode*, std::size_t> node_to_index;
            DdNode* node = nullptr;
            DdGen* bdd_node_iter = Cudd_FirstNode(bdd_mgr.getManager(), add.getNode(), &node);
            node_to_index.insert({node, bdd_nodes_.size()});

            auto is_terminal_node = [&](DdNode* node) -> bool {
                if(node == Cudd_ReadZero(bdd_mgr.getManager()) || node == Cudd_ReadOne(bdd_mgr.getManager()))
                    return true;
                return false;
            }; 

            auto get_node_index = [&](DdNode* node) -> std::size_t {
                if(node == Cudd_ReadZero(bdd_mgr.getManager())) {
                    return bdd_node::terminal_0;
                } else if(node == Cudd_ReadOne(bdd_mgr.getManager())) {
                    return bdd_node::terminal_1;
                } else {
                    assert(node_to_index.count(node) > 0);
                    return node_to_index.find(node)->second;
                }
            };

            // node indices of chain pointing to terminal_1
            constexpr static std::size_t pointer_to_terminal_1_not_set = std::numeric_limits<std::size_t>::max()-2;
            std::vector<std::size_t> var_to_bdd_node_terminal_1(std::distance(bdd_vars_begin, bdd_vars_end), pointer_to_terminal_1_not_set);
            var_to_bdd_node_terminal_1.back() = bdd_node::terminal_1;

            // add intermediate nodes if they do not exist, return index in bdd_nodes_ of node to point to

            // TODO: possibly share intermediate nodes not only for terminal node!
            auto add_intermediate_nodes = [&](DdNode* start, DdNode* end) -> std::size_t {

                const std::size_t start_var = Cudd_NodeReadIndex(start); 
                assert(start_var < std::distance(bdd_vars_begin, bdd_vars_end));

                if(!is_terminal_node(end)) {
                    const std::size_t end_var = Cudd_NodeReadIndex(end);
                    std::size_t last_index = get_node_index(end);
                    for(std::size_t i = end_var-1; i != start_var; --i) {
                        assert(i>0);
                        const std::size_t v_intermed = *(bdd_vars_begin + i);
                        bdd_nodes_.push_back({last_index, last_index, v_intermed});
                        last_index = bdd_nodes_.size()-1;
                    }
                    return last_index; 

                } else if(get_node_index(end) == bdd_node::terminal_1) {

                    if(var_to_bdd_node_terminal_1[start_var] == pointer_to_terminal_1_not_set) {
                        for(int i = int(std::distance(bdd_vars_begin, bdd_vars_end))-2; i >= int(start_var); --i) {
                            assert(i >= 0 && i < var_to_bdd_node_terminal_1.size());
                            if(var_to_bdd_node_terminal_1[i] == pointer_to_terminal_1_not_set) {
                                const std::size_t v_intermed = *(bdd_vars_begin + i+1);
                                bdd_nodes_.push_back({var_to_bdd_node_terminal_1[i+1], var_to_bdd_node_terminal_1[i+1], v_intermed}); 
                                check_node_valid(bdd_nodes_.back());
                                var_to_bdd_node_terminal_1[i] = bdd_nodes_.size()-1;
                            }
                        }
                    }
                    return var_to_bdd_node_terminal_1[start_var];
                    
                } else if(get_node_index(end) == bdd_node::terminal_0) {
                    return get_node_index(end);
                } else {
                    assert(false);
                    throw std::runtime_error("invalid node");
                }
            };

            auto add_node = [&](DdNode* node) {
                assert(Cudd_IsNonConstant(node));
                DdNode* high = Cudd_T(node);
                DdNode* low = Cudd_E(node);
                //const std::size_t high_index = get_node_index(high);
                //const std::size_t low_index = get_node_index(low);
                assert(Cudd_NodeReadIndex(node) < std::distance(bdd_vars_begin, bdd_vars_end));
                const std::size_t variable = *(bdd_vars_begin + Cudd_NodeReadIndex(node));
                //std::cout << "bdd var = " << variable << "\n";

                // for now: TODO: add dummy edges
                const std::size_t low_index = add_intermediate_nodes(node, low);
                const std::size_t high_index = add_intermediate_nodes(node, high);
                /*
                if(high_index != bdd_node::terminal_0 && high_index != bdd_node::terminal_1) {
                    const std::size_t bdd_high_var = Cudd_NodeReadIndex(high); 
                    const std::size_t bdd_cur_var = Cudd_NodeReadIndex(node);
                    std::size_t last_index = high_index;
                    for(std::size_t i = bdd_high_var-1; i != bdd_cur_var; --i) {
                        assert(i>0);
                        const std::size_t v_intermed = *(bdd_vars_begin + i);
                        bdd_nodes_.push_back({last_index, bdd_node::terminal_0, v_intermed});
                        last_index = bdd_nodes_.size()-1;

                    }
                    //assert(Cudd_NodeReadIndex(high) == Cudd_NodeReadIndex(node) + 1);
                    assert(variable < bdd_nodes_[high_index].variable);
                }
                if(low_index != bdd_node::terminal_0 && low_index != bdd_node::terminal_1) {
                    assert(Cudd_NodeReadIndex(low) == Cudd_NodeReadIndex(node) + 1);
                    assert(variable < bdd_nodes_[low_index].variable);
                }
                */

                assert(node_to_index.count(node) == 0);
                node_to_index.insert({node, bdd_nodes_.size()});
                bdd_nodes_.push_back(bdd_node{high_index, low_index, variable});
                check_node_valid(bdd_nodes_.back());
            };

            assert(Cudd_IsConstant(node));
            const std::size_t nr_bdd_nodes_begin = bdd_nodes_.size();

            while(Cudd_NextNode(bdd_node_iter, &node)) {
                if(!is_terminal_node(node)) {
                    add_node(node);
                }
            }
            const std::size_t nr_bdd_nodes_end = bdd_nodes_.size();
            std::vector<BDD> bdds = {bdd};
            //assert(nr_bdd_nodes == Cudd_DagSize(bdd.getNode()));
            //std::cout << "nr bdd nodes = " << nr_bdd_nodes << "\n";
            bdd_delimiters_.push_back(bdd_delimiters_.back() + nr_bdd_nodes_end - nr_bdd_nodes_begin);
        }

    bdd_storage::bdd_storage(const ILP_input& input, Cudd& bdd_mgr)
    { 
        std::vector<int> coefficients;
        std::vector<std::size_t> variables;
        bdd_converter converter(bdd_mgr);

        for(const auto& constraint : input.constraints()) {
            coefficients.clear();
            variables.clear();
            for(const auto e : constraint.variables) {
                coefficients.push_back(e.coefficient);
                variables.push_back(e.var);
            }
            assert(std::is_sorted(variables.begin(), variables.end()));
            
            BDD bdd = converter.convert_to_bdd(coefficients, constraint.ineq, constraint.right_hand_side);
            add_bdd(bdd_mgr, bdd, variables.begin(), variables.end());
        }
    }

    std::size_t bdd_storage::first_bdd_node(const std::size_t bdd_nr) const
    {
        assert(bdd_nr < nr_bdds());
        return bdd_nodes_[bdd_delimiters_[bdd_nr]].variable;
    }

    std::size_t bdd_storage::last_bdd_node(const std::size_t bdd_nr) const
    {
        assert(bdd_nr < nr_bdds());
        std::size_t max_node = 0;
        for(std::size_t i=bdd_delimiters_[bdd_nr]; i<bdd_delimiters_[bdd_nr+1]; ++i)
            max_node = std::max(max_node, bdd_nodes_[i].variable);
        return max_node;
    }

    template<typename STREAM>
        void bdd_storage::export_dot(STREAM& s) const
        {
            s << "digraph bdd_base {\n";

            auto get_node_string = [&](const std::size_t i) -> std::string {
                if(i == bdd_node::terminal_0)
                    return "false";
                if(i == bdd_node::terminal_1)
                    return "true";
                return std::to_string(i);
            };

            for(std::size_t i=0; i<bdd_nodes_.size(); ++i) {
                s << i << " -> " << get_node_string(bdd_nodes_[i].low) << " [label=\"0\"];\n";
                s << i << " -> " << get_node_string(bdd_nodes_[i].high) << " [label=\"1\"];\n";
            }
            s << "}\n"; 
        }

}
