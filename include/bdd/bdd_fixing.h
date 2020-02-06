#pragma once

#include "bdd_base.h"
#include "bdd_branch_node_fixing.h"

namespace LPMP {

    class bdd_fixing : public bdd_base<bdd_variable, bdd_branch_node> {
        public:
            bool fix_variables(std::vector<size_t> & indices, std::vector<char> & values);

            void init(const ILP_input& input);

        private:
        	bool fix_variable(size_t var, char value);
            bool remove_all_incoming_arcs(bdd_branch_node & bdd_node);
            void remove_all_outgoing_arcs(bdd_branch_node & bdd_node);

        	std::vector<char> primal_solution_ = std::vector<char>(nr_variables(), 2);
    };

    void bdd_fixing::init(const ILP_input& input)
    {
    	bdd_base<bdd_variable, bdd_branch_node>::init(input);

    	for (size_t var = 0; var < nr_variables(); var++)
    	{
    		for (size_t bdd_index = 0; bdd_index < nr_bdds(var); bdd_index++)
    		{
    			auto & bdd_var = bdd_variables_(var, bdd_index);
    			bdd_var.nr_feasible_low_arcs = 0;
    			bdd_var.nr_feasible_high_arcs = 0;
    			for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
    			{
    				auto & bdd_node = bdd_branch_nodes_[node_index];
    				if (bdd_node.low_outgoing != bdd_branch_node_terminal_0)
    					bdd_var.nr_feasible_low_arcs++;
    				if (bdd_node.high_outgoing != bdd_branch_node_terminal_0)
    					bdd_var.nr_feasible_high_arcs++;

    				bdd_node.bdd_var = & bdd_var;

    				auto * low_incoming = bdd_node.first_low_incoming;
    				while (low_incoming != nullptr && low_incoming->next_low_incoming != nullptr)
    				{
    					low_incoming->next_low_incoming->prev_low_incoming = low_incoming;
    					low_incoming = low_incoming.next_low_incoming;
    				}
    				auto * high_incoming = bdd_node.first_high_incoming;
    				while (high_incoming != nullptr && high_incoming->next_high_incoming != nullptr)
    				{
    					high_incoming->next_high_incoming->prev_high_incoming = high_incoming;
    					high_incoming = high_incoming.next_high_incoming;
    				}
    			}
    		}
    	}
    }

    bool bdd_fixing::fix_variable(std::size_t var, char value)
    {
        assert(0 <= value && value <= 1);
        assert(primal_solution_.size() == nr_variables());

        // check if variable is already fixed
        if (primal_solution_[var] == value)
            return true;
        else if (primal_solution_[var] < 2)
            return false;

        // mark variable as fixed
        primal_solution_[var] = value;

        for (size_t bdd_index = 0; bdd_index < nr_bdds(var); bdd_index++)
        {
            auto & bdd_var = bdd_variables_(var, bdd_index);
            for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
            {
                auto & bdd_node = bdd_branch_nodes_[node_index];

                if (value)
                {
                    // remove incoming arc from child node
                    if (!bdd_node.low_outgoing->is_terminal())
                    {
                        if (bdd_node.prev_low_incoming != nullptr)
                            bdd_node.prev_low_incoming->next_low_incoming = bdd_node.next_low_incoming;
                        else
                            bdd_node.low_outgoing->first_low_incoming = bdd_node.next_low_incoming;
                        // restructure child if it is unreachable
                        if (bdd_node.low_outgoing->first_low_incoming == nullptr && bdd_node.low_outgoing->first_high_incoming == nullptr)
                            remove_all_outgoing_arcs(*bdd_node.low_outgoing);
                    }
                    // turn outgoing arc towards 0-terminal
                    bdd_node.low_outgoing = bdd_branch_node_terminal_0;
                    bdd_node.bdd_var->nr_feasible_low_arcs--;
                }
                else
                {
                    if (!bdd_node.high_outgoing->is_terminal())
                    {
                        if (bdd_node.prev_high_incoming != nullptr)
                            bdd_node.prev_high_incoming->next_high_incoming = bdd_node.next_high_incoming;
                        else
                            bdd_node.high_outgoing->first_high_incoming = bdd_node.next_high_incoming;
                        if (bdd_node.high_outgoing->first_low_incoming == nullptr && bdd_node.high_outgoing->first_high_incoming == nullptr)
                            remove_all_outgoing_arcs(*bdd_node.high_outgoing);
                    }
                    bdd_node.high_outgoing = bdd_branch_node_terminal_0;
                    bdd_node.bdd_var->nr_feasible_high_arcs--;
                }

                // restructure parents if node is dead-end
                if (bdd_node.low_outgoing == bdd_branch_node_terminal_0 && bdd_node.high_outgoing == bdd_branch_node_terminal_0)
                {
                    if (!remove_all_incoming_arcs(bdd_node))
                        return false;
                }
            }

            // check if other variables in BDD are now restricted
            auto * cur = bdd_var.prev;
            while (cur != nullptr)
            {
                if (primal_solution_[variable_index(*cur)] < 2)
                {
                    cur = cur->prev;
                    continue;
                }
                if (cur->bdd_var->nr_feasible_low_arcs == 0)
                {
                    if (!fix_variable(variable_index(*cur), 1))
                        return false;
                }
                if (cur->bdd_var->nr_feasible_high_arcs == 0)
                {
                    if (!fix_variable(variable_index(*cur), 0))
                        return false;
                }
                cur = cur->prev;
            }
            cur = bdd_var.next;
            while (cur != nullptr)
            {
                if (primal_solution_[variable_index(*cur)] < 2)
                {
                    cur = cur->next;
                    continue;
                }
                if (cur->bdd_var->nr_feasible_low_arcs == 0)
                {
                    if (!fix_variable(variable_index(*cur), 1))
                        return false;
                }
                if (cur->bdd_var->nr_feasible_high_arcs == 0)
                {
                    if (!fix_variable(variable_index(*cur), 0))
                        return false;
                }
                cur = cur->next;
            }
        }
        return true;
    }

    bool bdd_fixing::remove_all_incoming_arcs(bdd_branch_node & bdd_node)
    {
        if (bdd_node.is_initial_state())
            return false;
        // low arcs
        {
            auto * cur = bdd_node.first_low_incoming;
            while (cur != nullptr)
            {
                cur->low_outgoing = bdd_branch_node_terminal_0;
                cur->bdd_var->nr_feasible_low_arcs--;
                // recursive call if parent is dead-end
                if (cur->high_outgoing == bdd_branch_node_terminal_0)
                {
                    if (!remove_all_incoming_arcs(*cur))
                        return false;
                }
                cur = cur->next_low_incoming;
            }
            bdd_node.first_low_incoming = nullptr;
        }
        // high arcs
        {
            auto * cur = bdd_node.first_high_incoming;
            while (cur != nullptr)
            {
                cur->high_outgoing = bdd_branch_node_terminal_0;
                cur->bdd_var->nr_feasible_high_arcs--;
                if (cur->low_outgoing == bdd_branch_node_terminal_0)
                {
                    if (!remove_all_incoming_arcs(*cur))
                        return false;
                }
                cur = cur->next_high_incoming;
            }
            bdd_node.first_high_incoming = nullptr;        
        }
        return true;
    }

    void bdd_fixing::remove_all_outgoing_arcs(bdd_branch_node & bdd_node)
    {
        // low arc
        if (!bdd_node.low_outgoing->is_terminal())
        {
            if (bdd_node.prev_low_incoming != nullptr)
                bdd_node.prev_low_incoming->next_low_incoming = bdd_node.next_low_incoming;
            else
                bdd_node.low_outgoing->first_low_incoming = bdd_node.next_low_incoming;
            // recursive call if child node is unreachable
            if (bdd_node.low_outgoing->first_low_incoming == nullptr && bdd_node.low_outgoing->first_high_incoming == nullptr)
                remove_all_outgoing_arcs(*bdd_node.low_outgoing);
        }
        bdd_node.low_outgoing = bdd_branch_node_terminal_0;
        bdd_node.bdd_var->nr_feasible_low_arcs--;
        // high arc
        if (!bdd_node.high_outgoing->is_terminal())
        {
            if (bdd_node.prev_high_incoming != nullptr)
                bdd_node.prev_high_incoming->next_high_incoming = bdd_node.next_high_incoming;
            else
                bdd_node.high_outgoing->first_high_incoming = bdd_node.next_high_incoming;
            if (bdd_node.high_outgoing->first_low_incoming == nullptr && bdd_node.high_outgoing->first_high_incoming == nullptr)
                remove_all_outgoing_arcs(*bdd_node.high_outgoing);
        }
        bdd_node.high_outgoing = bdd_branch_node_terminal_0;
        bdd_node.bdd_var->nr_feasible_high_arcs--;
    }

    bool bdd_fixing::fix_variables(const std::vector<size_t> & indices, const std::vector<char> & values)
    {
        assert(primal_solution_.size() == nr_variables());
        assert(indices.size() == values.size());

        std::fill(primal_solution_.begin(), primal_solution_.end(), 2);

        for (size_t var : indices)
        {
            if (!fix_variable(var, values[var]))
                return false;
        }
        return true;
    }
}
