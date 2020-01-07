#pragma once
#include <array>

namespace LPMP {

    class bdd_branch_instruction;
    void check_bdd_branch_instruction(const bdd_branch_instruction& bdd, const bool last_variable = false, const bool first_variable = false);

    class bdd_branch_instruction {
        public:
            bdd_branch_instruction* low_outgoing = nullptr;
            bdd_branch_instruction* high_outgoing = nullptr;
            bdd_branch_instruction* first_low_incoming = nullptr;
            bdd_branch_instruction* first_high_incoming = nullptr;
            bdd_branch_instruction* next_low_incoming = nullptr;
            bdd_branch_instruction* next_high_incoming = nullptr; 

            double* variable_cost = nullptr;
            double m = 0.0; // intermediate value of shortest path from either terminal or first node (depending on algorithm state)

            // From C++20
            //friend bool operator==(const bdd_branch_instruction&, const bdd_branch_instruction&) = default;
            friend bool operator==(const bdd_branch_instruction& x, const bdd_branch_instruction& y);

            void backward_step();
            void forward_step();

            // Debug functions for checking correctness of forward and backward step
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

    bool operator==(const bdd_branch_instruction& x, const bdd_branch_instruction& y)
    {
        const bool equal = (x.low_outgoing == y.low_outgoing &&
            x.high_outgoing == y.high_outgoing &&
            x.first_low_incoming == y.first_low_incoming &&
            x.first_high_incoming == y.first_high_incoming &&
            x.next_low_incoming == y.next_low_incoming &&
            x.next_high_incoming == y.next_high_incoming &&
            x.variable_cost == y.variable_cost &&
            x.m == y.m);
        return equal;
    }

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

    void check_bdd_branch_instruction(const bdd_branch_instruction& bdd, const bool last_variable, const bool first_variable)
    {
        assert(bdd.low_outgoing != nullptr);
        assert(bdd.low_outgoing == bdd_branch_instruction_terminal_0 || bdd.low_outgoing == bdd_branch_instruction_terminal_1 || bdd.low_outgoing > &bdd);
        assert(bdd.high_outgoing != nullptr);
        assert(bdd.high_outgoing == bdd_branch_instruction_terminal_0 || bdd.high_outgoing == bdd_branch_instruction_terminal_1 || bdd.high_outgoing > &bdd);
        //assert(bdd.cumulative_sum != nullptr); // TODO: remove.
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

}

