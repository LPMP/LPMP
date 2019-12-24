#include "two_dimensional_variable_array.hxx"
#include "cuddObj.hh"
#include "ILP_input.h"
#include "convert_pb_to_bdd.h"
#include "bdd_branch_instruction.h"
#include <cassert>
#include <vector>
#include <unordered_map>
#include <numeric>
#include <chrono> // for now
#include "bdd_storage.h"

namespace LPMP {

    //class bdd_branch_instruction;
    //void check_bdd_branch_instruction(const bdd_branch_instruction& bdd, const bool last_variable = false, const bool first_variable = false);

    /*
    class bdd_branch_instruction {
        public:
            bdd_branch_instruction* low_outgoing = nullptr;
            bdd_branch_instruction* high_outgoing = nullptr;
            bdd_branch_instruction* first_low_incoming = nullptr;
            bdd_branch_instruction* first_high_incoming = nullptr;
            bdd_branch_instruction* next_low_incoming = nullptr;
            bdd_branch_instruction* next_high_incoming = nullptr; 

            double* cumulative_sum = nullptr;
            double* variable_cost = nullptr;
            double m = 0.0;

            // From C++20
            //friend bool operator==(const bdd_branch_instruction&, const bdd_branch_instruction&) = default;
            friend bool operator==(const bdd_branch_instruction& x, const bdd_branch_instruction& y);

            void backward_step();
            void forward_step();

            double cost_from_first() const;
            double cost_from_terminal() const;

            bool is_terminal() const;
            bool is_first() const { return first_low_incoming == nullptr && first_high_incoming == nullptr; }

            bool is_initial_state() const { return *this == bdd_branch_instruction{}; }

            std::array<double,2> min_marginal() const;
    };

    static bdd_branch_instruction* const bdd_branch_instruction_terminal_0 = static_cast<bdd_branch_instruction*>(nullptr)+1;;
    static bdd_branch_instruction* const bdd_branch_instruction_terminal_1 = static_cast<bdd_branch_instruction*>(nullptr)+2;

    bool bdd_branch_instruction::is_terminal() const { return this == bdd_branch_instruction_terminal_0 || this == bdd_branch_instruction_terminal_1; }
    */

    class bdd_base {
        public:
            struct bdd_branch_instruction_level {
                std::size_t first_branch_instruction = std::numeric_limits<std::size_t>::max();
                std::size_t last_branch_instruction = std::numeric_limits<std::size_t>::max();
                // TODO: remove cumulative sum
                double cumulative_sum = std::numeric_limits<double>::infinity();
                double variable_cost = std::numeric_limits<double>::infinity(); 
                // TODO: prev and next pointers possibly not needed
                bdd_branch_instruction_level* prev = nullptr;
                bdd_branch_instruction_level* next = nullptr;

                std::size_t no_bdd_nodes() const { return last_branch_instruction - first_branch_instruction; }
                bool first_bdd_variable() const { return prev == nullptr; }
                bool is_initial_state() const { return *this == bdd_branch_instruction_level{}; }
                friend bool operator==(const bdd_branch_instruction_level&, const bdd_branch_instruction_level&);
            };

            bdd_base(Cudd& bdd_mgr) : bdd_mgr_(bdd_mgr) {}

            bdd_base(const bdd_base&) = delete; // no copy constructor because of pointers in bdd_branch_instruction

            template<typename BDD_VARIABLES_ITERATOR>
                void add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end);

            void init(const ILP_input& input);

            std::size_t nr_variables() const { return bdd_branch_instruction_levels.size(); }
            std::size_t nr_bdds(const std::size_t var) const { assert(var<nr_variables()); return bdd_branch_instruction_levels[var].size(); }

            // after all BDDs are added, initialize them so that optimization can start
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

            double lower_bound_backward_run();
            double lower_bound_backward(const std::size_t var, const std::size_t bdd_index);

            double min_marginal_averaging_iteration();
            void min_marginal_averaging_forward();
            double min_marginal_averaging_backward();

            double min_marginal_anisotropic_diffusion();
            void min_marginal_anisotropic_diffusion_forward();
            double min_marginal_anisotropic_diffusion_backward();

            const bdd_branch_instruction& get_bdd_branch_instruction(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const;

            template<typename STREAM>
                void export_dot(STREAM& s) const;

        private:
            void init_branch_instructions();
            std::array<double,2> min_marginal(const std::size_t var, const std::size_t bdd_index) const;
            template<typename ITERATOR>
                static std::array<double,2> average_marginals(ITERATOR marginals_begin, ITERATOR marginals_end);
            void set_marginal(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals);

            //void check_bdd_branch_instruction(const bdd_branch_instruction& bdd, const bool last_variable = false, const bool first_variable = false) const;
            void check_bdd_branch_instruction_level(const bdd_branch_instruction_level& bdd, const bool last_variable = false, const bool first_variable = false) const;

            std::size_t bdd_branch_instruction_index(const bdd_branch_instruction* bdd) const;
            std::size_t bdd_branch_instruction_index(const bdd_branch_instruction& bdd) const { return bdd_branch_instruction_index(&bdd); }
            std::size_t bdd_variable(const bdd_branch_instruction& bdd) const;
            std::size_t bdd_level_variable(const bdd_branch_instruction_level& bdd_level) const;

            std::size_t first_variable_of_bdd(const bdd_branch_instruction_level& bdd_level) const;
            std::size_t last_variable_of_bdd(const bdd_branch_instruction_level& bdd_level) const;

            std::vector<bdd_branch_instruction> bdd_branch_instructions;
            two_dim_variable_array<bdd_branch_instruction_level> bdd_branch_instruction_levels;
            Cudd& bdd_mgr_;
            bdd_storage bdd_storage_;
            std::vector<double> costs_; 
    };

    /*
    bool operator==(const bdd_branch_instruction& x, const bdd_branch_instruction& y)
    {
        const bool equal = (x.low_outgoing == y.low_outgoing &&
            x.high_outgoing == y.high_outgoing &&
            x.first_low_incoming == y.first_low_incoming &&
            x.first_high_incoming == y.first_high_incoming &&
            x.next_low_incoming == y.next_low_incoming &&
            x.next_high_incoming == y.next_high_incoming &&
            x.cumulative_sum == y.cumulative_sum &&
            x.variable_cost == y.variable_cost &&
            x.m == y.m);
        return equal;
    }
    */

    bool operator==(const bdd_base::bdd_branch_instruction_level& x, const bdd_base::bdd_branch_instruction_level& y)
    {
        return (x.first_branch_instruction == y.first_branch_instruction &&
            x.last_branch_instruction == y.last_branch_instruction &&
            x.cumulative_sum == y.cumulative_sum &&
            x.variable_cost == y.variable_cost &&
            x.prev == y.prev &&
            x.next == y.next); 
    }

    /*
    void bdd_branch_instruction::forward_step()
    {
        check_bdd_branch_instruction(*this);

        if(is_first()) {
            m = 0.0;
            //std::cout << "forward step m for " << this << " = " << m << "\n";
            return;
        }

        m = std::numeric_limits<double>::infinity();

        // iterate over all incoming low edges 
        {
            bdd_branch_instruction* cur = first_low_incoming;
            while(cur != nullptr) {
                //m = std::min(m, cur->m + *cumulative_sum - *(cur->cumulative_sum));
                m = std::min(m, cur->m);
                cur = cur->next_low_incoming;
            }
        }

        // iterate over all incoming high edges 
        {
            bdd_branch_instruction* cur = first_high_incoming;
            while(cur != nullptr) {
                //m = std::min(m, cur->m + *variable_cost + *cumulative_sum - *(cur->cumulative_sum));
                m = std::min(m, cur->m + *(cur->variable_cost));
                cur = cur->next_high_incoming;
            }
        }

        //std::cout << "forward step m for " << this << " = " << m << "\n";
        assert(std::isfinite(m));
        assert(std::abs(m - cost_from_first()) <= 1e-8);

        check_bdd_branch_instruction(*this);
    }

    void bdd_branch_instruction::backward_step()
    {
        check_bdd_branch_instruction(*this);

        // low edge
        const double low_cost = [&]() {
            if(low_outgoing == bdd_branch_instruction_terminal_0) {
                return std::numeric_limits<double>::infinity();
            } else if(low_outgoing == bdd_branch_instruction_terminal_1) {
                //return *cumulative_sum;
                return 0.0;
            } else {
                //return low_outgoing->m + *cumulative_sum - *(low_outgoing->cumulative_sum) - std::min(*low_outgoing->variable_cost,0.0);
                return low_outgoing->m;
            }
        }();

        // high edge
        const double high_cost = [&]() {
            if(high_outgoing == bdd_branch_instruction_terminal_0) {
                return std::numeric_limits<double>::infinity(); 
            } else if(high_outgoing == bdd_branch_instruction_terminal_1) {
                //return *cumulative_sum + *variable_cost; 
                return *variable_cost; 
            } else {
                //return high_outgoing->m + *variable_cost + *cumulative_sum - *(high_outgoing->cumulative_sum) - std::min(*high_outgoing->variable_cost,0.0);
                return high_outgoing->m + *variable_cost;
            }
        }();

        //std::cout << "variable cost = " << *variable_cost << ", low cost = " << low_cost << ", high cost = " << high_cost << ", cumulative sum = " << *cumulative_sum << "\n"; 

        assert(!std::isnan(low_cost));
        assert(!std::isnan(high_cost));
        assert(std::isfinite(std::min(low_cost,high_cost)));
        m = std::min(low_cost, high_cost);
        //std::cout << "backward step m for " << this << ", low outgoing = " << low_outgoing << ", high outgoing = " << high_outgoing << " = " << m << "\n";

        check_bdd_branch_instruction(*this);
        assert(std::abs(m - cost_from_terminal()) <= 1e-8);
    }

    double bdd_branch_instruction::cost_from_first() const
    {
        // TODO: only works if no bdd nodes skips variables
        double c = std::numeric_limits<double>::infinity();
        if(is_first())
            return 0.0;
        
        // iterate over all incoming low edges 
        {
            bdd_branch_instruction* cur = first_low_incoming;
            while(cur != nullptr) {
                c = std::min(c, cur->cost_from_first());
                cur = cur->next_low_incoming;
            }
        }

        // iterate over all incoming high edges 
        {
            bdd_branch_instruction* cur = first_high_incoming;
            while(cur != nullptr) {
                c = std::min(c, cur->cost_from_first() + *cur->variable_cost);
                cur = cur->next_high_incoming;
            }
        }

        return c;
    }

    double bdd_branch_instruction::cost_from_terminal() const
    {
        // TODO: only works if no bdd nodes skips variables
        // low edge
        const double low_cost = [&]() {
            if(low_outgoing == bdd_branch_instruction_terminal_0) {
                return std::numeric_limits<double>::infinity();
            } else if(low_outgoing == bdd_branch_instruction_terminal_1) {
                return 0.0;
            } else {
                return low_outgoing->cost_from_terminal();;
            }
        }();

        // high edge
        const double high_cost = [&]() {
            if(high_outgoing == bdd_branch_instruction_terminal_0) {
                return std::numeric_limits<double>::infinity(); 
            } else if(high_outgoing == bdd_branch_instruction_terminal_1) {
                return *variable_cost; 
            } else {
                return high_outgoing->cost_from_terminal() + *variable_cost;
            }
        }();

        return std::min(low_cost, high_cost); 
    }

    std::array<double,2> bdd_branch_instruction::min_marginal() const
    {
        check_bdd_branch_instruction(*this);

        //std::cout << "in min_marginal() for " << this << ", m = " << m << "\n";
        assert(std::abs(m - cost_from_first()) <= 1e-8);
        if(!low_outgoing->is_terminal()) {
            assert(std::abs(low_outgoing->m - low_outgoing->cost_from_terminal()) <= 1e-8);
        }
        if(!high_outgoing->is_terminal()) {
            assert(std::abs(high_outgoing->m - high_outgoing->cost_from_terminal()) <= 1e-8);
        }

        const double m0 = [&]() {
            if(low_outgoing == bdd_branch_instruction_terminal_0)
                return std::numeric_limits<double>::infinity();
            if(low_outgoing == bdd_branch_instruction_terminal_1)
                return this->m;
            return this->m + this->low_outgoing->m;
        }();

        const double m1 = [&]() {
            if(high_outgoing == bdd_branch_instruction_terminal_0)
                return std::numeric_limits<double>::infinity();
            if(high_outgoing == bdd_branch_instruction_terminal_1)
                return this->m + *this->variable_cost;
            return this->m + *this->variable_cost + this->high_outgoing->m;
        }();

        assert(std::isfinite(std::min(m0,m1)));

        return {m0,m1};
    }
    */

    void bdd_base::forward_step(const std::size_t var, const std::size_t bdd_index)
    {
        assert(var < bdd_branch_instruction_levels.size());
        assert(bdd_index < bdd_branch_instruction_levels[var].size());

        auto& bdd_level = bdd_branch_instruction_levels(var,bdd_index);
        assert(var != 0 || bdd_level.prev == nullptr);

        // update cumulative sum
        if(bdd_level.prev != nullptr) {
            bdd_level.cumulative_sum = bdd_level.prev->cumulative_sum + std::min(0.0, bdd_level.prev->variable_cost);
        } else {
            bdd_level.cumulative_sum = 0.0;
        }
        bdd_level.cumulative_sum = 0.0; // TODO: for now
        //std::cout << "forward cumulative sum for var " << var << " = " << bdd_level.cumulative_sum << "\n";

        // iterate over all bdd nodes and make forward step
        const std::size_t first_branch_instruction = bdd_level.first_branch_instruction;
        const std::size_t last_branch_instruction = bdd_level.last_branch_instruction;
        for(std::size_t i=first_branch_instruction; i<last_branch_instruction; ++i) {
            //std::cout << "forward step for var = " << var << ", bdd_index = " << bdd_index << ", bdd_node_index = " << i << "\n";
            bdd_branch_instructions[i].forward_step();
        }
    }

    void bdd_base::forward_step(const std::size_t var)
    {
        assert(var < bdd_branch_instruction_levels.size());
        for(std::size_t bdd_index = 0; bdd_index<bdd_branch_instruction_levels[var].size(); ++bdd_index)
            forward_step(var, bdd_index);
    }

    void bdd_base::backward_step(const std::size_t var, const std::size_t bdd_index)
    {
        assert(var < bdd_branch_instruction_levels.size());
        assert(bdd_index < bdd_branch_instruction_levels[var].size());

        auto& bdd_level = bdd_branch_instruction_levels(var,bdd_index);
        assert(var+1 != nr_variables() || bdd_level.next == nullptr);

        // TODO: can be removed
        // update cumulative sum
        //if(bdd_level.next != nullptr) {
        //    bdd_level.cumulative_sum = bdd_level.next->cumulative_sum + std::min(0.0, bdd_level.next->variable_cost);
        //} else {
        //    bdd_level.cumulative_sum = 0.0;
        //}
        //bdd_level.cumulative_sum = 0.0; // TODO: for now
        //std::cout << "backward cumulative sum for var " << var << " = " << bdd_level.cumulative_sum << "\n";

        // iterate over all bdd nodes and make forward step
        const std::size_t first_branch_instruction = bdd_level.first_branch_instruction;
        const std::size_t last_branch_instruction = bdd_level.last_branch_instruction;
        //std::cout << "no bdd nodes for var " << var << " = " << last_branch_instruction - first_branch_instruction << "\n";
        //std::cout << "bdd branch instruction offset = " << first_branch_instruction << "\n";
        for(std::size_t i=first_branch_instruction; i<last_branch_instruction; ++i) {
            check_bdd_branch_instruction(bdd_branch_instructions[i], var+1 == nr_variables(), var == 0);
            bdd_branch_instructions[i].backward_step();
        }
    }

    void bdd_base::backward_step(const std::size_t var)
    {
        assert(var < bdd_branch_instruction_levels.size());
        for(std::size_t bdd_index = 0; bdd_index<bdd_branch_instruction_levels[var].size(); ++bdd_index)
            backward_step(var, bdd_index);
    }

    template<typename BDD_VARIABLES_ITERATOR>
        void bdd_base::add_bdd(BDD& bdd, BDD_VARIABLES_ITERATOR bdd_vars_begin, BDD_VARIABLES_ITERATOR bdd_vars_end)
        {
            bdd_storage_.add_bdd(bdd_mgr_, bdd, bdd_vars_begin, bdd_vars_end); 
        }

    void bdd_base::init()
    {
        init_branch_instructions();
        costs_.resize(nr_variables(), std::numeric_limits<double>::infinity());
        std::cout << "nr bdds = " << bdd_storage_.nr_bdds() << "\n";
    }

    void bdd_base::init(const ILP_input& input)
    {
        bdd_storage_ = bdd_storage(input, bdd_mgr_);
        init();
        set_costs(input.objective().begin(), input.objective().end());
    }

    void bdd_base::init_branch_instructions()
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
        
        std::vector<std::size_t> bdd_offset_per_variable = {0};
        std::partial_sum(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end()-1, std::back_inserter(bdd_offset_per_variable));
        assert(bdd_offset_per_variable.size() == bdd_storage_.nr_variables());
        assert(bdd_offset_per_variable.back() + nr_bdd_nodes_per_variable.back() == bdd_storage_.bdd_nodes().size());

        bdd_branch_instructions.resize(bdd_storage_.bdd_nodes().size());
        bdd_branch_instruction_levels.resize(nr_bdds_per_variable.begin(), nr_bdds_per_variable.end());

        for(std::size_t v=0; v<nr_bdds_per_variable.size(); ++v) {
            //std::cout << "v = " << v << ", nr bdd nodes per var = " << nr_bdd_nodes_per_variable[v] << ", offset = " << bdd_offset_per_variable[v] << "\n";
        }

        // fill branch instructions into datastructures and set pointers
        std::fill(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end(), 0); // counter for bdd instruction per variable
        std::fill(nr_bdds_per_variable.begin(), nr_bdds_per_variable.end(), 0); // counter for bdd per variable

        // TODO: use tsl::unordered_map
        std::unordered_map<std::size_t, bdd_branch_instruction*> stored_bdd_node_index_to_bdd_address;
        //stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, bdd_branch_instruction*>(bdd_storage::bdd_node::terminal_0, bdd_branch_instruction_terminal_0));
        //stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, bdd_branch_instruction*>(bdd_storage::bdd_node::terminal_1, bdd_branch_instruction_terminal_1)); 

        std::size_t c = 0; // check if everything is read contiguously
        for(std::size_t bdd_index=0; bdd_index<bdd_storage_.nr_bdds(); ++bdd_index) {
            uncover_variables(); 
            stored_bdd_node_index_to_bdd_address.clear();
            stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, bdd_branch_instruction*>(bdd_storage::bdd_node::terminal_0, bdd_branch_instruction_terminal_0));
            stored_bdd_node_index_to_bdd_address.insert(std::pair<std::size_t, bdd_branch_instruction*>(bdd_storage::bdd_node::terminal_1, bdd_branch_instruction_terminal_1)); 
            bdd_branch_instruction_level* next_bdd_level = nullptr;

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

                auto& bdd_level = bdd_branch_instruction_levels(v, nr_bdds_per_variable[v]);
                if(!variable_covered(v)) {
                    assert(bdd_level.is_initial_state());
                    cover_variable(v);

                    bdd_level.first_branch_instruction = bdd_offset_per_variable[v];
                    bdd_level.last_branch_instruction = bdd_offset_per_variable[v];
                    //std::cout << "bdd level offset for var = " << v << ", bdd index = " << bdd_index << " = " << stored_bdd_node_index << "\n";
                    // TODO: remove, not needed anymore
                    if(next_bdd_level != nullptr) {
                        //assert(next_bdd_level > &bdd_level);
                        bdd_level.next = next_bdd_level;
                        next_bdd_level->prev = &bdd_level;
                    }

                    next_bdd_level = &bdd_level;
                } else {
                    assert(!bdd_level.is_initial_state());
                }

                const std::size_t bdd_branch_instructions_index = bdd_level.last_branch_instruction; 
                bdd_level.last_branch_instruction++;

                bdd_branch_instruction& bdd = bdd_branch_instructions[bdd_branch_instructions_index];
                assert(bdd.is_initial_state());
                //std::cout << "address = " << &bdd << "\n";
                bdd.variable_cost = &bdd_level.variable_cost;
                bdd.cumulative_sum = &bdd_level.cumulative_sum;

                stored_bdd_node_index_to_bdd_address.insert({stored_bdd_node_index, &bdd});

                assert(stored_bdd_node_index_to_bdd_address.count(stored_bdd.low) > 0);
                bdd_branch_instruction& bdd_low = *(stored_bdd_node_index_to_bdd_address.find(stored_bdd.low)->second);
                bdd.low_outgoing = &bdd_low;
                assert(bdd.low_outgoing != nullptr);
                if(!bdd_low.is_terminal()) {
                    bdd.next_low_incoming = bdd_low.first_low_incoming;
                    bdd_low.first_low_incoming = &bdd;
                }

                assert(stored_bdd_node_index_to_bdd_address.count(stored_bdd.high) > 0);
                bdd_branch_instruction& bdd_high = *(stored_bdd_node_index_to_bdd_address.find(stored_bdd.high)->second);
                bdd.high_outgoing = &bdd_high;
                if(!bdd_high.is_terminal()) {
                    bdd.next_high_incoming = bdd_high.first_high_incoming;
                    bdd_high.first_high_incoming = &bdd;
                }

                //if(!bdd.low_outgoing->is_terminal()) { assert(bdd_variable(bdd) < bdd_variable(*bdd.low_outgoing)); }
                //if(!bdd.high_outgoing->is_terminal()) { assert(bdd_variable(bdd) < bdd_variable(*bdd.high_outgoing)); }
                check_bdd_branch_instruction(bdd, v+1 == nr_variables(), v == 0);
                check_bdd_branch_instruction_level(bdd_level, v+1 == nr_variables(), v == 0);
                ++nr_bdd_nodes_per_variable[v]; 
                ++bdd_offset_per_variable[v];
            }

            for(const std::size_t v : variables_covered) {
                ++nr_bdds_per_variable[v];
            }
        }

        for(std::size_t v=0; v<nr_bdds_per_variable.size(); ++v) {
            assert(nr_bdds_per_variable[v] == bdd_branch_instruction_levels[v].size());
        }

        for(const auto& bdd : bdd_branch_instructions) {
            check_bdd_branch_instruction(bdd);
        }

        // TODO: clear bdd_storage
    }

    std::array<double,2> bdd_base::min_marginal(const std::size_t var, const std::size_t bdd_index) const
    {
        assert(var < nr_variables());
        assert(bdd_index < nr_bdds(var));
        std::array<double,2> m = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        const auto& bdd_level = bdd_branch_instruction_levels(var,bdd_index);
        for(std::size_t bdd_node_index=bdd_level.first_branch_instruction; bdd_node_index<bdd_level.last_branch_instruction; ++bdd_node_index) {
            //std::cout << "min marginal for var = " << var << ", bdd_index = " << bdd_index << ", bdd_node_index = " << bdd_node_index << "\n";
            const auto& bdd = bdd_branch_instructions[bdd_node_index];
            const auto [m0,m1] = bdd.min_marginal();
            m[0] = std::min(m[0], m0);
            m[1] = std::min(m[1], m1);
        }
        assert(std::isfinite(m[0]));
        assert(std::isfinite(m[1]));
        return m;
    }

    void bdd_base::set_marginal(const std::size_t var, const std::size_t bdd_index, const std::array<double,2> marginals, const std::array<double,2> min_marginals)
    {
        assert(var < nr_variables());
        assert(bdd_index < nr_bdds(var));
        assert(min_marginals == min_marginal(var,bdd_index));
        assert(std::isfinite(marginals[0]) && std::isfinite(marginals[1]));
        assert(std::isfinite(min_marginals[0]) && std::isfinite(min_marginals[1]));

        auto& bdd_level = bdd_branch_instruction_levels(var,bdd_index);
        const double marginal_diff = min_marginals[1] - min_marginals[0];
        const double marginal_diff_target = marginals[1] - marginals[0];
        bdd_level.variable_cost += -marginal_diff + marginal_diff_target; 
    }

    void bdd_base::backward_run()
    {
        for(long int var=nr_variables()-1; var>=0; --var) {
            backward_step(var); 
        } 
    }

    double bdd_base::lower_bound_backward_run()
    {
        backward_run();
        double lb = 0.0;
        for(long int var=nr_variables()-1; var>=0; --var)
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index)
                lb += lower_bound_backward(var,bdd_index);
        return lb;
    }

    double bdd_base::lower_bound_backward(const std::size_t var, const std::size_t bdd_index)
    {
        const auto& bdd_level = bdd_branch_instruction_levels(var,bdd_index);
        if(bdd_level.prev == nullptr) {
            assert(bdd_level.no_bdd_nodes() == 1);
            const auto& bdd = get_bdd_branch_instruction(var, bdd_index, 0);
            return bdd.m;
        } else {
            return 0.0;
        }
    }

    void bdd_base::forward_run()
    {
        for(std::size_t var=0; var<nr_variables(); ++var) {
            forward_step(var); 
        }
    }

    double bdd_base::min_marginal_averaging_iteration()
    {
        const auto begin_time = std::chrono::steady_clock::now();
        min_marginal_averaging_forward();
        const auto after_forward = std::chrono::steady_clock::now();
        std::cout << "forward pass took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(after_forward - begin_time).count() << " milliseconds\n";
        const auto before_backward = std::chrono::steady_clock::now();
        const double lb = min_marginal_averaging_backward();
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "backward pass took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - before_backward).count() << " milliseconds\n";
        return lb;
    }

    template<typename ITERATOR>
        std::array<double,2> bdd_base::average_marginals(ITERATOR marginals_begin, ITERATOR marginals_end)
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

    // min marginal averaging
    void bdd_base::min_marginal_averaging_forward()
    {
        std::vector<std::array<double,2>> min_marginals;
        for(std::size_t var=0; var<nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                forward_step(var,bdd_index);
                min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            if(min_marginals.size() > 1) {
                //std::cout << "min marginals in forward pass: ";
                //for(const auto m : min_marginals)
                    //std::cout << "(" << m[0] << "," << m[1] << "), ";
                //std::cout << "\n";
            }

            //std::cout << "min marginal for variable " << var << ": ";
            //for(const auto& m : min_marginals)
            //    std::cout << m[0] << "," << m[1] << "; ";
            //std::cout << "\n";

            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());
            //std::cout << "averaged marginals = " << average_marginal[0] << "," << average_marginal[1] << "\n";

            // set marginals in each bdd so min marginals match each other
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
            }

        }
    }

    double bdd_base::min_marginal_averaging_backward()
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

        return lb;
    }

    double bdd_base::min_marginal_anisotropic_diffusion()
    {
        two_dim_variable_array<std::size_t> first_bdd_var(bdd_branch_instruction_levels);
        for(std::size_t v=0; v<nr_variables(); ++v) {
            for(std::size_t bdd_index=0; bdd_index<bdd_branch_instruction_levels[v].size(); ++bdd_index) {
                const auto& bdd_level = bdd_branch_instruction_levels(v,bdd_index);
                if(bdd_level.prev == nullptr)
                    first_bdd_var(v,bdd_index) = v;
                else {
                    const auto& prev_bdd_level = *bdd_level.prev;

                }
            }
        }   

        two_dim_variable_array<std::size_t> last_bdd_var(bdd_branch_instruction_levels);

        min_marginal_anisotropic_diffusion_forward();
        return min_marginal_anisotropic_diffusion_backward();
    }

    void bdd_base::min_marginal_anisotropic_diffusion_forward()
    {
        /*
        std::vector<std::array<double,2>> min_marginals;
        std::vector<char> last_variable;
        for(std::size_t var=0; var<nr_variables(); ++var) {

            // collect min marginals
            min_marginals.clear();
            last_variable.clear();
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                forward_step(var,bdd_index);
                min_marginals.push_back(min_marginal(var,bdd_index)); 
                last_variable.push_back(last_bdd_variable(bdd_index,var));
            }

            // set min-marginals of bdds where current var is last one to zero, for other average it
            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());
            }

            for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                if(last_variable[bdd_index]) {
                    const double min_val = std::min(min_marginals[bdd_index][0], min_marginals[bdd_index][1]);
                    set_marginal(var,bdd_index,{min_val, min_val}, min_marginals[bdd_index]);
                } else {
                    set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]); 
                }
            }

        }
        */
    }

    double bdd_base::min_marginal_anisotropic_diffusion_backward()
    {
        return 0.0;
    }

    const bdd_branch_instruction& bdd_base::get_bdd_branch_instruction(const std::size_t var, const std::size_t bdd_index, const std::size_t bdd_node_index) const
    {
        assert(var < nr_variables());
        assert(bdd_index < bdd_branch_instruction_levels[var].size());
        assert(bdd_node_index < bdd_branch_instruction_levels(var,bdd_index).no_bdd_nodes());
        return bdd_branch_instructions[ bdd_branch_instruction_levels(var,bdd_index).first_branch_instruction + bdd_node_index ];
    }

    template<typename ITERATOR>
        void bdd_base::set_costs(ITERATOR begin, ITERATOR end)
        {
            assert(std::distance(begin,end) <= nr_variables());
            std::copy(begin, end, costs_.begin());
            std::fill(costs_.begin() + std::distance(begin, end), costs_.end(), 0.0);

            // distribute costs to bdds uniformly
            for(std::size_t v=0; v<nr_variables(); ++v) {
                assert(bdd_branch_instruction_levels[v].size() > 0);
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(v); ++bdd_index) {
                    bdd_branch_instruction_levels(v,bdd_index).variable_cost = *(begin+v)/nr_bdds(v);
                }
            }
            backward_run();
        }

    template<typename ITERATOR>
        bool bdd_base::check_feasibility(ITERATOR var_begin, ITERATOR var_end) const
        {
            assert(std::distance(var_begin, var_end) == nr_variables());

            std::vector<char> bdd_branch_instruction_marks(bdd_branch_instructions.size(), 0);

            std::size_t var = 0;
            for(auto var_iter=var_begin; var_iter!=var_end; ++var_iter, ++var) {
                const bool val = *(var_begin+var);
                for(std::size_t bdd_index=0; bdd_index<nr_bdds(var); ++bdd_index) {
                    const auto& bdd_level = bdd_branch_instruction_levels(var, bdd_index);
                    if(bdd_level.first_bdd_variable()) {
                        bdd_branch_instruction_marks[bdd_level.first_branch_instruction] = 1;
                    }
                    for(std::size_t bdd_node_index=bdd_level.first_branch_instruction; bdd_node_index<bdd_level.last_branch_instruction; ++bdd_node_index) {
                        if(bdd_branch_instruction_marks[bdd_node_index] == 1) {
                            const auto& bdd = bdd_branch_instructions[bdd_node_index];
                            const auto* bdd_next_index = [&]() {
                                if(val == false)
                                    return bdd.low_outgoing;
                                else 
                                    return bdd.high_outgoing;
                            }();

                            if(bdd_next_index == bdd_branch_instruction_terminal_0)
                                return false;
                            if(bdd_next_index == bdd_branch_instruction_terminal_1) {
                            } else { 
                                bdd_branch_instruction_marks[ bdd_branch_instruction_index(bdd_next_index) ] = 1;
                            }
                        }
                    }
                }
            }

            return true;
        }

    template<typename ITERATOR>
        double bdd_base::evaluate(ITERATOR var_begin, ITERATOR var_end) const
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

    /*
    void check_bdd_branch_instruction(const bdd_branch_instruction& bdd, const bool last_variable, const bool first_variable)
    {
        assert(bdd.low_outgoing != nullptr);
        assert(bdd.low_outgoing == bdd_branch_instruction_terminal_0 || bdd.low_outgoing == bdd_branch_instruction_terminal_1 || bdd.low_outgoing > &bdd);
        assert(bdd.high_outgoing != nullptr);
        assert(bdd.high_outgoing == bdd_branch_instruction_terminal_0 || bdd.high_outgoing == bdd_branch_instruction_terminal_1 || bdd.high_outgoing > &bdd);
        assert(bdd.cumulative_sum != nullptr);
        assert(bdd.variable_cost != nullptr);
        //assert(bdd.low_outgoing != bdd.high_outgoing); // this need not hold true for our BDDs.
        assert(std::isfinite(bdd.m));
        if(last_variable) {
            assert(bdd.low_outgoing == bdd_branch_instruction_terminal_0 || bdd.low_outgoing == bdd_branch_instruction_terminal_1);
            assert(bdd.high_outgoing == bdd_branch_instruction_terminal_0 || bdd.high_outgoing == bdd_branch_instruction_terminal_1);
        }
        if(first_variable) {
            assert(bdd.first_low_incoming == nullptr);
            assert(bdd.first_high_incoming == nullptr);
        } 
        if(bdd.first_low_incoming != nullptr) {
            bdd_branch_instruction* cur = bdd.first_low_incoming;
            while(cur != nullptr) {
                assert(cur < &bdd);
                assert(cur->low_outgoing == &bdd);
                cur = cur->next_low_incoming;
            }
        }
        if(bdd.first_high_incoming != nullptr) {
            bdd_branch_instruction* cur = bdd.first_high_incoming;
            while(cur != nullptr) {
                assert(cur < &bdd);
                assert(cur->high_outgoing == &bdd);
                cur = cur->next_high_incoming;
            }
        }
    }
    */

    void bdd_base::check_bdd_branch_instruction_level(const bdd_branch_instruction_level& bdd_level, const bool last_variable, const bool first_variable) const
    {
        // go over all branch instructions and check whether they point to all the same variable cost and cumulative sum
        for(std::size_t bdd_node_index = bdd_level.first_branch_instruction; bdd_node_index < bdd_level.last_branch_instruction; ++bdd_node_index) {
            assert(&bdd_level.variable_cost == bdd_branch_instructions[bdd_node_index].variable_cost);
            assert(&bdd_level.cumulative_sum == bdd_branch_instructions[bdd_node_index].cumulative_sum);
        }

        if(last_variable) {
            assert(bdd_level.next == nullptr);
        }
        if(first_variable) {
            assert(bdd_level.prev == nullptr);
        }
    }

    std::size_t bdd_base::bdd_branch_instruction_index(const bdd_branch_instruction* bdd) const
    {
        assert(bdd >= &bdd_branch_instructions[0]);
        const std::size_t i = bdd - &bdd_branch_instructions[0];
        assert(i < bdd_branch_instructions.size());
        return i; 
    }

    std::size_t bdd_base::bdd_variable(const bdd_branch_instruction& bdd) const
    {
        // TODO: recursive search
        const std::size_t i = bdd_branch_instruction_index(&bdd);
        for(std::size_t v=0; v<nr_variables(); ++v) {
            for(std::size_t bdd_index=0; bdd_index<bdd_branch_instruction_levels[v].size(); ++bdd_index)
            if(i >= bdd_branch_instruction_levels(v,bdd_index).first_branch_instruction &&
                    i < bdd_branch_instruction_levels(v,bdd_index).last_branch_instruction)
                return v;
        }
        throw std::runtime_error("Could not obtain bdd variable");
    }

    std::size_t bdd_base::bdd_level_variable(const bdd_branch_instruction_level& bdd_level) const
    {
        assert(&bdd_level >= &bdd_branch_instruction_levels(0,0));
        assert(&bdd_level <= &bdd_branch_instruction_levels.back().back());
        std::size_t lb = 0;
        std::size_t ub = nr_variables()-1;
        std::size_t v = nr_variables()/2;
        while(! (&bdd_level >= &bdd_branch_instruction_levels(v,0) && &bdd_level <= &bdd_branch_instruction_levels[v].back()) ) {
            if( &bdd_level > &bdd_branch_instruction_levels(v,0) ) {
                lb = v+1;
            } else {
                ub = v-1;
            }
            v = (ub+lb)/2; 
        }
        assert(&bdd_level >= &bdd_branch_instruction_levels(v,0) && &bdd_level <= &bdd_branch_instruction_levels[v].back());
        return v; 
    }

    std::size_t bdd_base::first_variable_of_bdd(const bdd_branch_instruction_level& bdd_level) const
    {
        bdd_branch_instruction_level const* p = &bdd_level;
        while(p->prev != nullptr)
            p = p->prev;
        return bdd_level_variable(*p);
    }

    std::size_t bdd_base::last_variable_of_bdd(const bdd_branch_instruction_level& bdd_level) const
    {
        bdd_branch_instruction_level const* p = &bdd_level;
        while(p->next != nullptr)
            p = p->next;
        return bdd_level_variable(*p); 
    }

    template<typename STREAM>
        void bdd_base::export_dot(STREAM& s) const
        {
            return bdd_storage_.export_dot(s);
            s << "digraph bdd_base {\n";

            std::vector<char> bdd_node_visited(bdd_branch_instructions.size(), false);
            std::size_t cur_bdd_index = 0;
            for(const auto& bdd : bdd_branch_instructions) {
                const std::size_t i = bdd_branch_instruction_index(bdd);
                if(bdd_node_visited[i])
                    continue;
                bdd_node_visited[i] = true;
                auto get_node_string = [&](const bdd_branch_instruction* bdd) -> std::string {
                    if(bdd == bdd_branch_instruction_terminal_0)
                        return "false";
                    if(bdd == bdd_branch_instruction_terminal_1)
                        return "true";
                    return std::to_string(bdd_branch_instruction_index(bdd));
                };
                s << i << " -> " << get_node_string(bdd.low_outgoing) << " [label=\"0\"];\n";
                s << i << " -> " << get_node_string(bdd.high_outgoing) << " [label=\"1\"];\n";
            }

            s << "}\n"; 
        }
}
