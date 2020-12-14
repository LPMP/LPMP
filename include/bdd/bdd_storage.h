#pragma once

#include "bdd.h"
#include "ILP_input.h"
#include "convert_pb_to_bdd.h"
#include "hash_helper.hxx"
#include "bdd_preprocessor.h"
#include "bdd_collection.h"
#include "tclap/CmdLine.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <numeric>

namespace LPMP {

    // used for storing all participating BDDs temporarily before exporting them to the bdd solver.
    class bdd_storage {
        public:

        struct bdd_node {
            constexpr static std::size_t terminal_0 = std::numeric_limits<std::size_t>::max()-1;
            constexpr static std::size_t terminal_1 = std::numeric_limits<std::size_t>::max();
            bool low_is_terminal() const { return low == terminal_0 || low == terminal_1; }
            bool high_is_terminal() const { return low == terminal_0 || low == terminal_1; }

            std::size_t low;
            std::size_t high;
            std::size_t variable;
        };

        bdd_storage(TCLAP::CmdLine& cmd);
        void init(const ILP_input& input);

        template <typename BDD_VARIABLES_ITERATOR>
        void add_bdd(BDD::bdd_mgr& bdd_mgr, BDD::node_ref bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end);

        void add_bdd(BDD::bdd_collection_entry bdd);

        template <typename STREAM>
        void export_dot(STREAM &s) const;

        std::size_t nr_bdds() const { return bdd_delimiters().size()-1; }
        std::size_t nr_variables() const { return nr_variables_; }
        std::size_t nr_bdd_nodes(const std::size_t bdd_nr) const { assert(bdd_nr < nr_bdds()); return bdd_delimiters_[bdd_nr+1] - bdd_delimiters_[bdd_nr]; }

        const std::vector<bdd_node>& bdd_nodes() const { return bdd_nodes_; }
        const std::vector<std::size_t>& bdd_delimiters() const { return bdd_delimiters_; }

        // TODO: rename to ..._variable
        std::size_t first_bdd_node(const std::size_t bdd_nr) const;
        std::size_t last_bdd_node(const std::size_t bdd_nr) const;

        // return all edges with endpoints being variables that are consecutive in some BDD
        std::vector<std::array<size_t,2>> dependency_graph() const;

    private:
        void check_node_valid(const bdd_node bdd) const;

        template<typename BDD_NODE_TYPE, typename BDD_GETTER, typename VARIABLE_GETTER, typename NEXT_BDD_NODE, typename BDD_VARIABLE_ITERATOR>
            void add_bdd_impl(
                    const size_t nr_bdds,
                    BDD_GETTER bdd_getter,
                    VARIABLE_GETTER& get_var,
                    NEXT_BDD_NODE& next_bdd_node,
                    BDD_VARIABLE_ITERATOR bdd_vars_begin, BDD_VARIABLE_ITERATOR bdd_vars_end
                    );

        void check_bdd_node(const bdd_node bdd) const;

        std::vector<bdd_node> bdd_nodes_;
        std::vector<std::size_t> bdd_delimiters_ = {0};
        std::vector<std::size_t> nr_bdd_nodes_per_variable;
        std::vector<std::size_t> nr_bdds_per_variable;
        std::size_t nr_variables_ = 0;

        // for BDD decomposition //
        std::vector<size_t> interval_boundaries;
        std::vector<size_t> intervals;
        void compute_intervals(const size_t nr_intervals);
        size_t interval(const size_t variable) const;
        size_t nr_intervals() const;
        std::tuple<two_dim_variable_array<bdd_storage::bdd_node>, two_dim_variable_array<size_t>> split_bdd_nodes(const size_t nr_intervals);

        TCLAP::MultiArg<std::string> preprocessing_arg;
    };


    bdd_storage::bdd_storage(TCLAP::CmdLine& cmd)
        : preprocessing_arg("","bdd_preprocessing","preprocess BDDs",false,"{none|bridge|subsumption|subsumption_except_one|contiguous_overlap|partial_contiguous_overlap|cliques}", cmd)
    {}

    template<typename BDD_VARIABLES_ITERATOR>
        void bdd_storage::add_bdd(BDD::bdd_mgr& bdd_mgr, BDD::node_ref bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end)
        {
            assert(std::is_sorted(bdd_vars_begin, bdd_vars_end));
            assert(std::distance(bdd_vars_begin, bdd_vars_end) > 0);

            std::vector<BDD::node_ref> bdd_nodes = bdd.nodes_postorder();

            auto get_bdd = [&](const size_t bdd_nr) {
                assert(bdd_nr < bdd_nodes.size());
                return bdd_nodes[bdd_nr];
            };

            auto get_variable = [&](BDD::node_ref& bdd) { 
                const size_t i = bdd.variable();
                return i;
                assert(i < std::distance(bdd_vars_begin, bdd_vars_end));
                return *(bdd_vars_begin + i); 
            };

            auto get_next_node = [&](BDD::node_ref& bdd) -> BDD::node_ref { 
                const ptrdiff_t p = std::distance(&bdd_nodes[0], &bdd);
                assert(p >= 0 && p < bdd_nodes.size());
                return bdd_nodes[p+1];
            };

            add_bdd_impl<BDD::node_ref>(
                    bdd_nodes.size(),
                    get_bdd,
                    get_variable,
                    get_next_node,
                    bdd_vars_begin, bdd_vars_end
                    );

            return;
        }

    void bdd_storage::add_bdd(BDD::bdd_collection_entry bdd)
    {
        const auto vars = bdd.variables();
        std::unordered_map<size_t,size_t> rebase_to_iota;
        for(size_t i=0; i<vars.size(); ++i)
            rebase_to_iota.insert({vars[i], i});
        bdd.rebase(rebase_to_iota);

        auto get_node = [&](const size_t i) {
            const size_t j = bdd.nr_nodes() - 3 - i;
            return bdd[j]; 
        };
        auto get_next_node = [&](BDD::bdd_collection_node node) { return node.next_postorder(); };
        auto get_variable = [&](BDD::bdd_collection_node node) { return node.variable(); };

        add_bdd_impl<BDD::bdd_collection_node>(
                bdd.nr_nodes() - 2, // do not count in top- and botsink
                get_node,
                get_variable,
                get_next_node,
                vars.begin(), vars.end()
                );

        bdd.rebase(vars.begin(), vars.end());
    }

    template<typename BDD_NODE_TYPE, typename BDD_GETTER, typename VARIABLE_GETTER, typename NEXT_BDD_NODE, typename BDD_VARIABLE_ITERATOR>
        void bdd_storage::add_bdd_impl(
                const size_t nr_bdds,
                BDD_GETTER bdd_getter,
                VARIABLE_GETTER& get_var, 
                NEXT_BDD_NODE& next_bdd_node,
                BDD_VARIABLE_ITERATOR bdd_vars_begin, BDD_VARIABLE_ITERATOR bdd_vars_end)
        {
            std::unordered_map<BDD_NODE_TYPE, size_t> node_to_index;

            auto get_node_index = [&](BDD_NODE_TYPE node) -> std::size_t {
                if(node.is_botsink()) {
                    return bdd_node::terminal_0;
                } else if(node.is_topsink()) {
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

            auto add_intermediate_nodes = [&](BDD_NODE_TYPE start, BDD_NODE_TYPE end) -> std::size_t {

                const std::size_t start_var = get_var(start);

                if(!end.is_terminal()) {
                    const size_t end_var = get_var(end);
                    size_t last_index = get_node_index(end);
                    for(std::size_t i = end_var-1; i != start_var; --i) {
                        assert(i>0);
                        const std::size_t v_intermed = *(bdd_vars_begin + i);
                        bdd_nodes_.push_back({last_index, last_index, v_intermed});
                        last_index = bdd_nodes_.size()-1;
                    }
                    return last_index; 

                } else if(get_node_index(end) == bdd_node::terminal_1) {

                    if(var_to_bdd_node_terminal_1[start_var] == pointer_to_terminal_1_not_set) {
                        for(std::ptrdiff_t i = std::ptrdiff_t(std::distance(bdd_vars_begin, bdd_vars_end))-2; i >= std::ptrdiff_t(start_var); --i) {
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

            const size_t nr_bdd_nodes_begin = bdd_nodes_.size();

            for(size_t i=0; i<nr_bdds; ++i)
            {
                auto node = bdd_getter(i);
                const size_t variable = get_var(node);
                nr_variables_ = std::max(*(bdd_vars_begin+variable)+1, nr_variables_);

                const size_t low_index = add_intermediate_nodes(node, node.low());
                const size_t high_index = add_intermediate_nodes(node, node.high());

                assert(node_to_index.count(node) == 0);
                node_to_index.insert({node, bdd_nodes_.size()});
                bdd_nodes_.push_back(bdd_node{high_index, low_index, *(bdd_vars_begin+variable)});
                check_node_valid(bdd_nodes_.back()); 
            } 

            const size_t nr_bdd_nodes_end = bdd_nodes_.size();
            bdd_delimiters_.push_back(bdd_delimiters_.back() + nr_bdd_nodes_end - nr_bdd_nodes_begin);
        }


    void bdd_storage::init(const ILP_input& input)
    {
        BDD::bdd_mgr bdd_mgr;

        bdd_preprocessor bdd_pre;
        const bool preprocess = preprocessing_arg.getValue().size() > 0;

        // first transform linear inequalities into BDDs
        std::vector<int> coefficients;
        std::vector<std::size_t> variables;
        std::vector<BDD::node_ref> bdds;
        bdd_converter converter(bdd_mgr);

        for(const auto& constraint : input.constraints()) {
            coefficients.clear();
            variables.clear();
            for(const auto e : constraint.variables) {
                coefficients.push_back(e.coefficient);
                variables.push_back(e.var);
            }
            assert(std::is_sorted(variables.begin(), variables.end()));
            
            BDD::node_ref bdd = converter.convert_to_bdd(coefficients, constraint.ineq, constraint.right_hand_side);
            if(preprocess)
                bdd_pre.add_bdd(bdd, variables.begin(), variables.end());
            else
                add_bdd(bdd_mgr, bdd, variables.begin(), variables.end());
        }

        // second, preprocess BDDs
        if(preprocessing_arg.getValue().size() > 0)
        {
            for(const std::string& preprocessing : preprocessing_arg.getValue())
            {
                if(preprocessing == "bridge")
                    bdd_pre.set_coalesce_bridge();
                else if(preprocessing == "subsumption")
                    bdd_pre.set_coalesce_subsumption();
                else if(preprocessing == "contiguous_overlap")
                    bdd_pre.set_coalesce_contiguous_overlap();
                else if(preprocessing == "subsumption_except_one")
                    bdd_pre.set_coalesce_subsumption_except_one();
                else if(preprocessing == "partial_contiguous_overlap")
                    bdd_pre.set_coalesce_partial_contiguous_overlap();
                else if(preprocessing == "cliques")
                    bdd_pre.set_coalesce_cliques();
                else
                    throw std::runtime_error("bdd preprocessing argument " + preprocessing + " not recognized.");
            }
            bdd_pre.coalesce_bdd_collection();

            for(size_t bdd_nr=0; bdd_nr<bdd_pre.get_bdd_collection().nr_bdds(); ++bdd_nr)
                add_bdd(bdd_pre.get_bdd_collection()[bdd_nr]);
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

    std::vector<std::array<size_t,2>> bdd_storage::dependency_graph() const
    {
        std::unordered_set<std::array<size_t,2>> edges;
        std::unordered_set<size_t> cur_vars;
        std::vector<size_t> cur_vars_sorted;
        for(size_t bdd_nr=0; bdd_nr<nr_bdds(); ++bdd_nr)
        {
            cur_vars.clear();
            for(size_t i=bdd_delimiters_[bdd_nr]; i<bdd_delimiters_[bdd_nr+1]; ++i)
                cur_vars.insert(bdd_nodes_[i].variable);
            cur_vars_sorted.clear();
            for(const size_t v : cur_vars)
                cur_vars_sorted.push_back(v);
            std::sort(cur_vars_sorted.begin(), cur_vars_sorted.end());
            for(size_t i=0; i+1<cur_vars_sorted.size(); ++i)
                edges.insert({cur_vars_sorted[i], cur_vars_sorted[i+1]});
        }

        return std::vector<std::array<size_t,2>>(edges.begin(), edges.end());
    }

    template<typename STREAM>
        void bdd_storage::export_dot(STREAM& s) const
        {
            s << "digraph bdd_min_marginal_averaging {\n";

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

    inline void bdd_storage::check_node_valid(const bdd_node bdd) const
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

    ///////////////////////////
    // for BDD decomposition //
    ///////////////////////////

    void bdd_storage::compute_intervals(const size_t nr_intervals)
    {
        // first approach: Just partition variables equidistantly.
        // TODO: Take into account number of BDD nodes in each interval
        
        assert(nr_intervals > 1);
        interval_boundaries.clear();
        interval_boundaries.reserve(nr_intervals+1);
        for(size_t interval=0; interval<nr_intervals; ++interval)
            interval_boundaries.push_back(std::round(double(interval*this->nr_variables())/double(nr_intervals)));

        interval_boundaries.push_back(this->nr_variables()); 
        assert(interval_boundaries.size() == nr_intervals+1);

        intervals.clear();
        intervals.reserve(this->nr_variables());
        for(size_t interval = 0; interval+1<interval_boundaries.size(); ++interval)
            for(size_t var=interval_boundaries[interval]; var<interval_boundaries[interval+1]; ++var)
                intervals.push_back(interval);
    }

    size_t bdd_storage::interval(const size_t variable) const
    {
        assert(variable < this->nr_variables());
        assert(interval_boundaries.size() > 2 || interval_boundaries.size() == 0);
        if(interval_boundaries.size() == 0)
            return 0;
        assert(intervals.size() == this->nr_variables());
        return intervals[variable];
    }

    size_t bdd_storage::nr_intervals() const
    {
        assert(interval_boundaries.size() > 2 || interval_boundaries.size() == 0);
        if(interval_boundaries.size() == 0)
            return 1;
        return interval_boundaries.size()-1;
    }

    // take bdd_nodes_ and return two_dim_variable_array<bdd_node> bdd_nodes_split_, two_dim_variable_array<size_t> bdd_delimiters_split_
    std::tuple<two_dim_variable_array<bdd_storage::bdd_node>, two_dim_variable_array<size_t>> bdd_storage::split_bdd_nodes(const size_t nr_intervals)
    {
        compute_intervals(nr_intervals);

        std::vector<size_t> nr_bdd_nodes_per_interval(nr_intervals(), 0);
        std::vector<size_t> nr_bdds_per_inteval(nr_intervals(), 1);
        std::unordered_set<size_t> active_intervals;

        for(size_t bdd_counter=0; bdd_counter<bdd_delimiters_.size()-1; ++bdd_counter)
        {
            const size_t last_bdd_interval = [&]() {
                size_t last_bdd_interval = 0;
                for(auto bdd_node_counter=bdd_delimiters_[bdd_counter]; bdd_node_counter<bdd_delimiters_[bdd_counter+1]; ++bdd_node_counter)
                {
                    const auto& bdd = bdd_nodes_[bdd_node_counter];
                    last_bdd_interval = std::max(interval(bdd.variable), last_bdd_interval);
                }
                return last_bdd_interval;
            }();
            
            // there are the following cases:
            // (i) All bdd nodes are non-terminal and in the same interval
            // (ii) all bdd nodes are non-terminal but in different intervals
            // (iii) one arc is bottom and the other one is in the same interval
            // (iv) one arc is bottom and the other one is in the next interval
            // (v) At least one arc is top -> bdd node is in last interval
            
            // count in which intervals bdd has nodes
            active_intervals.clear();
            for(auto bdd_node_counter=bdd_delimiters_[bdd_counter]; bdd_node_counter<bdd_delimiters_[bdd_counter+1]; ++bdd_node_counter)
            {
                const auto& bdd = bdd_nodes_[bdd_node_counter];
                active_intervals.insert(interval(bdd.variable));
            }
            for(const size_t i : active_intervals)
                ++nr_bdds_per_inteval[i];

            // count number of bdd nodes per interval
            for(auto bdd_node_counter=bdd_delimiters_[bdd_counter]; bdd_node_counter<bdd_delimiters_[bdd_counter+1]; ++bdd_node_counter)
            {
                const auto& bdd = bdd_nodes_[bdd_node_counter];
                // first check if node is split node. If not, increase node count in correct interval. If split node, increment node count in both intervals straddled.
                // case (i)
                if(!bdd.low_is_terminal() && !bdd.high_is_terminal())
                {
                    const bdd_node& low = bdd_nodes_[bdd.low];
                    const bdd_node& high = bdd_nodes_[bdd.high];
                    assert(low.variable == high.variable);
                    if(interval(bdd.variable) == interval(low.variable)) // case (i)
                    {
                        ++nr_bdd_nodes_per_interval[interval(bdd.variable)]; 
                    }
                    else // case (ii)
                    {
                        ++nr_bdd_nodes_per_interval[interval(bdd.variable)];
                        ++nr_bdd_nodes_per_interval[interval(bdd_nodes_[bdd.low].variable)]; // bdd nodes pointed  to by low and high will be counted in next interval again when visiting those 
                    }
                }
                else if(bdd.low_is_terminal() && bdd.high_is_terminal()) // case (v)
                {
                    assert(bdd.low == bdd_node::terminal_1 || bdd.high == bdd_node::terminal_1);
                    ++nr_bdd_nodes_per_interval[interval(bdd.variable)]; 
                }
                else if(bdd.low == bdd_node::terminal_0 || bdd.high == bdd_node::terminal_0)
                {
                    if(bdd.low == bdd_node::terminal_0)
                    {
                        const size_t high_var = bdd_nodes_[bdd.high].variable;
                        if(interval(bdd.variable) == interval(high_var)) // case (iii)
                        {
                            ++nr_bdd_nodes_per_interval[interval(bdd.variable)];
                        }
                        else // case (iv)
                        {
                            // TODO: not necessarily += 2, possibly += 1 if node is shared!
                            ++nr_bdd_nodes_per_interval[interval(bdd.variable)];
                            ++nr_bdd_nodes_per_interval[interval(high_var)];
                        }
                    }
                    else
                    {
                        assert(bdd.high == bdd_node::terminal_0);
                        const size_t low_var = bdd_nodes_[bdd.low].variable;
                        if(interval(bdd.variable) == interval(low_var)) // case (iii)
                        {
                            ++nr_bdd_nodes_per_interval[interval(bdd.variable)];
                        }
                        else // case (iv)
                        {
                            // TODO: not necessarily += 2, possibly += 1 if node is shared!
                            ++nr_bdd_nodes_per_interval[interval(bdd.variable)];
                            ++nr_bdd_nodes_per_interval[interval(low_var)];
                        }
                    } 
                }
                else
                {
                    assert(false); // We should have covered all cases
                }
            }
        }

        two_dim_variable_array<bdd_node> split_bdd_nodes(nr_bdd_nodes_per_interval.begin(), nr_bdd_nodes_per_interval.end());
        std::fill(nr_bdd_nodes_per_interval.begin(), nr_bdd_nodes_per_interval.end(), 0);

        two_dim_variable_array<size_t> split_bdd_delimiters(nr_bdds_per_inteval.begin(), nr_bdds_per_inteval.end());
        std::fill(nr_bdds_per_inteval.begin(), nr_bdds_per_inteval.end(), 1);
        for(size_t i=0; i<split_bdd_delimiters.size(); ++i)
            split_bdd_delimiters(i,0) = 0;

        // fill split bdd nodes, record duplicated arcs
        struct split_node { size_t interval; size_t offset; };
        std::vector<std::array<split_node,2>> duplicated_nodes;
        std::unordered_map<std::array<size_t,2>,size_t> split_bdd_node_indices; // bdd index in bdd_nodes_, interval
        for(size_t bdd_counter=0; bdd_counter<bdd_delimiters_.size()-1; ++bdd_counter)
        {
            split_bdd_node_indices.clear();
            for(auto bdd_node_counter=bdd_delimiters_[bdd_counter]; bdd_node_counter<bdd_delimiters_[bdd_counter+1]; ++bdd_node_counter)
            { 
                const bdd_node& bdd = bdd_nodes_[bdd_node_counter];
                // case (i) & (ii)
                if(!bdd.low_is_terminal() && !bdd.high_is_terminal())
                {
                    const size_t i = interval(bdd.variable);
                    const bdd_node& low = bdd_nodes_[bdd.low];
                    const bdd_node& high = bdd_nodes_[bdd.high];
                    assert(low.variable == high.variable);

                    if(interval(bdd.variable) == interval(low.variable)) // case (i)
                    {
                        assert(split_bdd_node_indices.count({bdd.low, i}) > 0);
                        const size_t low_idx = split_bdd_node_indices.find({bdd.low, i})->second;
                        assert(split_bdd_node_indices.count({bdd.high, i}) > 0);
                        const size_t high_idx = split_bdd_node_indices.find({bdd.high, i})->second;

                        split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, low_idx, high_idx};
                        split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i]));
                        ++nr_bdd_nodes_per_interval[i];
                    }
                    else // case (ii)
                    {
                        // in interval i, low and high arcs should point to topsink
                        split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, bdd_node::terminal_1, bdd_node::terminal_1};
                        split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i]));
                        ++nr_bdd_nodes_per_interval[i];

                        // in next interval
                        const size_t next_i = interval(bdd_nodes_[bdd.low].variable);
                        assert(i < next_i);
                        const size_t next_lo_idx = split_bdd_node_indices.find({bdd.low, next_i})->second;
                        const size_t next_hi_idx = split_bdd_node_indices.find({bdd.high, next_i})->second;
                        split_bdd_nodes(next_i, nr_bdd_nodes_per_interval[next_i]) = {bdd.variable, next_lo_idx, next_hi_idx};
                        split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, next_i}, nr_bdd_nodes_per_interval[next_i]));
                        ++nr_bdd_nodes_per_interval[next_i]; 

                        duplicated_nodes.push_back({split_node{i, nr_bdd_nodes_per_interval[i]-1}, split_node{next_i, nr_bdd_nodes_per_interval[next_i]-1}});
                    }
                }
                else if(bdd.low_is_terminal() && bdd.high_is_terminal()) // case (v)
                {
                    assert(bdd.low == bdd_node::terminal_1 || bdd.high == bdd_node::terminal_1);
                    const size_t i = interval(bdd.variable);
                    split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, bdd.low, bdd.high};
                    split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i])); 
                    ++nr_bdd_nodes_per_interval[i]; 
                }
                else if(bdd.low == bdd_node::terminal_0 || bdd.high == bdd_node::terminal_0)
                {
                    const size_t i = interval(bdd.variable);
                    if(bdd.low == bdd_node::terminal_0)
                    {
                        const size_t high_var = bdd_nodes_[bdd.high].variable;
                        if(i == interval(high_var)) // case (iii)
                        {
                            assert(split_bdd_node_indices.count({bdd.high,i}) > 0);
                            const size_t high_idx = split_bdd_node_indices.find({bdd.high,i})->second;
                            split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, bdd_node::terminal_0, high_idx};
                            split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i])); 
                            ++nr_bdd_nodes_per_interval[i]; 
                        }
                        else // case (iv)
                        {
                            split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, bdd_node::terminal_0, bdd_node::terminal_1};
                            split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i])); 
                            ++nr_bdd_nodes_per_interval[i]; 

                            const size_t next_i = interval(high_var);
                            assert(split_bdd_node_indices.count({bdd.high, next_i}) > 0);
                            const size_t next_high_idx = split_bdd_node_indices.find({bdd.high, next_i})->second;
                            split_bdd_nodes(next_i, nr_bdd_nodes_per_interval[next_i]) = {bdd.variable, bdd_node::terminal_0, next_high_idx};
                            ++nr_bdd_nodes_per_interval[next_i]; 

                            duplicated_nodes.push_back({split_node{i, nr_bdd_nodes_per_interval[i]-1}, split_node{next_i, nr_bdd_nodes_per_interval[next_i]-1}});
                        }
                    }
                    else
                    {
                        assert(bdd.high == bdd_node::terminal_0);
                        const size_t low_var = bdd_nodes_[bdd.low].variable;
                        if(interval(bdd.variable) == interval(low_var)) // case (iii)
                        {
                            assert(split_bdd_node_indices.count({bdd.low,i}) > 0);
                            const size_t low_idx = split_bdd_node_indices.find({bdd.high,i})->second;
                            split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, low_idx, bdd_node::terminal_0};
                            split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i])); 
                            ++nr_bdd_nodes_per_interval[i]; 
                        }
                        else // case (iv)
                        {
                            split_bdd_nodes(i, nr_bdd_nodes_per_interval[i]) = {bdd.variable, bdd_node::terminal_1, bdd_node::terminal_0};
                            split_bdd_node_indices.insert(std::make_pair(std::array<size_t,2>{bdd_node_counter, i}, nr_bdd_nodes_per_interval[i])); 
                            ++nr_bdd_nodes_per_interval[i]; 

                            const size_t next_i = interval(low_var);
                            assert(split_bdd_node_indices.count({bdd.low, next_i}) > 0);
                            const size_t next_low_idx = split_bdd_node_indices.find({bdd.low, next_i})->second;
                            split_bdd_nodes(next_i, nr_bdd_nodes_per_interval[next_i]) = {bdd.variable, next_low_idx, bdd_node::terminal_0};
                            ++nr_bdd_nodes_per_interval[next_i]; 

                            duplicated_nodes.push_back({split_node{i, nr_bdd_nodes_per_interval[i]-1}, split_node{next_i, nr_bdd_nodes_per_interval[next_i]-1}});
                        }
                    } 
                }
            }
            // go over each affected interval and set new bdd delimiter values
            active_intervals.clear();
            for(auto bdd_node_counter=bdd_delimiters_[bdd_counter]; bdd_node_counter<bdd_delimiters_[bdd_counter+1]; ++bdd_node_counter)
            {
                const auto& bdd = bdd_nodes_[bdd_node_counter];
                active_intervals.insert(interval(bdd.variable));
            }
            for(const size_t i : active_intervals)
                split_bdd_delimiters(i, nr_bdds_per_inteval[i]++) = nr_bdd_nodes_per_interval[i];
        }

        return {split_bdd_nodes, split_bdd_delimiters};
    }

} // namespace LPMP
