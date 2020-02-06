#pragma once
#include <array>

namespace LPMP {

    class bdd_branch_node;
    void check_bdd_branch_node(const bdd_branch_node& bdd, const bool last_variable = false, const bool first_variable = false);

    class bdd_variable {
        std::size_t first_node_index = std::numeric_limits<std::size_t>::max();
        std::size_t last_node_index = std::numeric_limits<std::size_t>::max();

        bdd_variable* prev = nullptr;
        bdd_variable* next = nullptr;

        std::size_t nr_feasible_low_arcs;
        std::size_t nr_feasible_high_arcs;

        std::size_t nr_bdd_nodes() const { return last_node_index - first_node_index; }
        bool is_first_bdd_variable() const { return prev == nullptr; }
        bool is_initial_state() const { return *this == bdd_variable{}; }
        friend bool operator==(const bdd_variable&, const bdd_variable&);
    };

    bool operator==(const bdd_variable& x, const bdd_variable& y)
    {
        return (x.first_node_index == y.first_node_index &&
            x.last_node_index == y.last_node_index &&
            x.prev == y.prev &&
            x.next == y.next); 
    }

    class bdd_branch_node {
        public:
            bdd_branch_node* low_outgoing = nullptr;
            bdd_branch_node* high_outgoing = nullptr;
            bdd_branch_node* first_low_incoming = nullptr;
            bdd_branch_node* first_high_incoming = nullptr;
            bdd_branch_node* next_low_incoming = nullptr;
            bdd_branch_node* next_high_incoming = nullptr;

            bdd_branch_node* prev_low_incoming = nullptr;
            bdd_branch_node* prev_high_incoming = nullptr;

            bdd_variable* bdd_var;

            // From C++20
            friend bool operator==(const bdd_branch_node& x, const bdd_branch_node& y);

            bool is_terminal() const;
            bool is_first() const { return first_low_incoming == nullptr && first_high_incoming == nullptr; }

            bool is_initial_state() const { return *this == bdd_branch_node{}; }
    };

    static bdd_branch_node* const bdd_branch_node_terminal_0 = static_cast<bdd_branch_node*>(nullptr)+1;
    static bdd_branch_node* const bdd_branch_node_terminal_1 = static_cast<bdd_branch_node*>(nullptr)+2;

    bool bdd_branch_node::is_terminal() const { return this == bdd_branch_node_terminal_0 || this == bdd_branch_node_terminal_1; }

    bool operator==(const bdd_branch_node& x, const bdd_branch_node& y)
    {
        const bool equal = (x.low_outgoing == y.low_outgoing &&
            x.high_outgoing == y.high_outgoing &&
            x.first_low_incoming == y.first_low_incoming &&
            x.first_high_incoming == y.first_high_incoming &&
            x.prev_low_incoming == y.prev_low_incoming &&
            x.prev_high_incoming == y.prev_high_incoming &&
            x.next_low_incoming == y.next_low_incoming &&
            x.next_high_incoming == y.next_high_incoming &&
            x.bdd_var == y.bdd_var);
        return equal;
    }

    
    void check_bdd_branch_node(const bdd_branch_node& bdd, const bool last_variable, const bool first_variable)
    {
#ifdef NDEBUG
        return;
#endif
        assert(bdd.low_outgoing != nullptr);
        assert(bdd.low_outgoing == bdd_branch_node_terminal_0 || bdd.low_outgoing == bdd_branch_node_terminal_1 || bdd.low_outgoing > &bdd);
        assert(bdd.high_outgoing != nullptr);
        assert(bdd.high_outgoing == bdd_branch_node_terminal_0 || bdd.high_outgoing == bdd_branch_node_terminal_1 || bdd.high_outgoing > &bdd);
        //assert(bdd.cumulative_sum != nullptr); // TODO: remove.
        assert(bdd.variable_cost != nullptr);
        //assert(bdd.low_outgoing != bdd.high_outgoing); // this need not hold true for our BDDs.
        assert(std::isfinite(bdd.m));
        if(last_variable) {
            assert(bdd.low_outgoing == bdd_branch_node_terminal_0 || bdd.low_outgoing == bdd_branch_node_terminal_1);
            assert(bdd.high_outgoing == bdd_branch_node_terminal_0 || bdd.high_outgoing == bdd_branch_node_terminal_1);
        }
        if(first_variable) {
            assert(bdd.first_low_incoming == nullptr);
            assert(bdd.first_high_incoming == nullptr);
        } 
        if(bdd.first_low_incoming != nullptr) {
            bdd_branch_node* cur = bdd.first_low_incoming;
            while(cur != nullptr) {
                assert(cur < &bdd);
                assert(cur->low_outgoing == &bdd);
                cur = cur->next_low_incoming;
            }
        }
        if(bdd.first_high_incoming != nullptr) {
            bdd_branch_node* cur = bdd.first_high_incoming;
            while(cur != nullptr) {
                assert(cur < &bdd);
                assert(cur->high_outgoing == &bdd);
                cur = cur->next_high_incoming;
            }
        }
    }

}

