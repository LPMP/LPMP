#pragma once


#include "ILP_input.h"

#include "bdd_variable.h"
#include "bdd_branch_node.h"
#include "bdd_min_marginal_averaging_smoothed.h"

#include <cassert>
#include <vector>
#include <stack>

namespace LPMP {


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

    class bdd_mma_fixing : public bdd_mma_base<bdd_variable_fix, bdd_branch_node_fix> {
        public:
            bool fix_variables();

            bool fix_variables(const std::vector<size_t> & indices, const std::vector<char> & values);
            bool fix_variable(const size_t var, const char value);
            bool is_fixed(const size_t var) const;

            std::vector<double> total_min_marginals();
            std::vector<double> search_space_reduction_coeffs();

            void min_marginal_averaging_forward();
            void min_marginal_averaging_backward();
            void min_marginal_averaging_iteration();

            void init(const ILP_input& input);
            void init();

            void revert_changes(const size_t target_log_size);

            void init_primal_solution() { primal_solution_.resize(nr_variables(), 2); }
            const std::vector<char> & primal_solution() const { return primal_solution_; }
            double compute_upper_bound();
            const size_t log_size() const { return log_.size(); }

        private:
            void init_pointers();

            bool remove_all_incoming_arcs(bdd_branch_node_fix & bdd_node);
            void remove_all_outgoing_arcs(bdd_branch_node_fix & bdd_node);
            void remove_outgoing_low_arc(bdd_branch_node_fix & bdd_node);
            void remove_outgoing_high_arc(bdd_branch_node_fix & bdd_node);

            std::vector<char> primal_solution_;
            std::stack<log_entry, std::deque<log_entry>> log_;

            // JUST FOR INSPECTION
            ILP_input input_;
    };

    double bdd_mma_fixing::compute_upper_bound()
    {
        return bdd_mma_base<bdd_variable_fix, bdd_branch_node_fix>::compute_upper_bound(primal_solution());
    }

    void bdd_mma_fixing::init_pointers()
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
        init_primal_solution(); 
    }

    void bdd_mma_fixing::init(const ILP_input& input)
    {
        bdd_mma_base<bdd_variable_fix, bdd_branch_node_fix>::init(input);
        init_pointers();
        input_ = input;
    }

    void bdd_mma_fixing::init()
    {
        bdd_mma_base<bdd_variable_fix, bdd_branch_node_fix>::init();
        init_pointers();
    }

    bool bdd_mma_fixing::fix_variable(const std::size_t var, const char value)
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
                auto & bdd_node = bdd_branch_nodes_[node_index];

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

    bool bdd_mma_fixing::remove_all_incoming_arcs(bdd_branch_node_fix & bdd_node)
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

    void bdd_mma_fixing::remove_all_outgoing_arcs(bdd_branch_node_fix & bdd_node)
    {
        remove_outgoing_low_arc(bdd_node);
        remove_outgoing_high_arc(bdd_node);
    }

    void bdd_mma_fixing::remove_outgoing_low_arc(bdd_branch_node_fix & bdd_node)
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

    void bdd_mma_fixing::remove_outgoing_high_arc(bdd_branch_node_fix & bdd_node)
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

    bool bdd_mma_fixing::fix_variables(const std::vector<size_t> & variables, const std::vector<char> & values)
    {
        assert(variables.size() == values.size());

        init_primal_solution();

        struct VarFix
        {
            VarFix(const size_t log_size, const size_t index, const char val)
            : log_size_(log_size), index_(index), val_(val) {}
            
            const size_t log_size_;
            const size_t index_;
            const char val_;
        };

        std::stack<VarFix, std::deque<VarFix>> variable_fixes;
        variable_fixes.emplace(log_.size(), 0, 1-values[0]);
        variable_fixes.emplace(log_.size(), 0, values[0]);

        size_t nfixes = 0;
        size_t max_fixes = nr_variables();
        // size_t max_fixes = std::numeric_limits<size_t>::max();
        std::cout << "Search tree node budget: " << max_fixes << std::endl;
        std::cout << "Expanded: " << std::endl;

        while (!variable_fixes.empty())
        {
            nfixes++;
            std::cout << "\r" << nfixes << std::flush;
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

    bool bdd_mma_fixing::is_fixed(const size_t var) const
    {
        assert(primal_solution_.size() == nr_variables());
        assert(var < primal_solution_.size());
        return primal_solution_[var] < 2;
    }

    void bdd_mma_fixing::revert_changes(const size_t target_log_size)
    {
        while (log_.size() > target_log_size)
        {
            log_.top().restore();
            log_.pop();
        }
    }

    std::vector<double> bdd_mma_fixing::total_min_marginals()
    {
        std::vector<double> total_min_marginals;
        for(std::size_t var=0; var<this->nr_variables(); ++var)
        {
            double total_min_marg = 0;
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index)
            {
                this->forward_step(var,bdd_index);
                if (is_fixed(var))
                    continue;
                std::array<double,2> min_marg = min_marginal(var,bdd_index);
                total_min_marg += (min_marg[1] - min_marg[0]);
            }
            total_min_marginals.push_back(total_min_marg);
            // std::cout << input_.get_var_name(var) << " : " << total_min_marginals[var] << std::endl;
        }
        return total_min_marginals;
    }

    std::vector<double> bdd_mma_fixing::search_space_reduction_coeffs()
    {
        std::vector<double> coeffs;
        // solution count backward run
        for (ptrdiff_t var = this->nr_variables()-1; var >= 0; --var)
        {
            for (size_t bdd_index=0; bdd_index<this->nr_bdds(var); bdd_index++)
            {
                auto & bdd_var = bdd_variables_(var, bdd_index);
                for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                {
                    auto & bdd_node = bdd_branch_nodes_[node_index];
                    bdd_node.count_backward_step();
                }
            }
        }
        // solution count forward run
        for (size_t var = 0; var < this->nr_variables(); var++)
        {
            double coeff = 0;
            for (size_t bdd_index=0; bdd_index<this->nr_bdds(var); bdd_index++)
            {
                auto & bdd_var = bdd_variables_(var, bdd_index);
                for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                {
                    auto & bdd_node = bdd_branch_nodes_[node_index];
                    bdd_node.count_forward_step();
                    coeff += bdd_node.count_high() - bdd_node.count_low();
                }
            }
            coeffs.push_back(coeff);
        }
        return coeffs;
    }

    bool bdd_mma_fixing::fix_variables()
    {
        std::vector<double> reduction_coeffs = search_space_reduction_coeffs();
        // additional MMA iteration increases chance of finding a feasible solution
        this->backward_run();
        min_marginal_averaging_iteration();
        std::vector<double> total_min_marginals = this->total_min_marginals();
        std::vector<size_t> variables;
        for (size_t i = 0; i < nr_variables(); i++)
            variables.push_back(i);

        const double eps = -std::numeric_limits<double>::epsilon();

        auto sign = [](const double val) -> double
        {
            if (val < 0)
                return -1.0;
            else if (val > 0)
                return 1.0;
            else
                return 0.0;
        };

        auto order = [&](const size_t a, const size_t b)
        {
            return sign(reduction_coeffs[a]) * total_min_marginals[a] > sign(reduction_coeffs[b]) * total_min_marginals[b];
        };
        std::sort(variables.begin(), variables.end(), order);

        std::vector<char> values;
        for (size_t i = 0; i < variables.size(); i++)
        {
            const char val = (total_min_marginals[variables[i]] < eps) ? 1 : 0;
            values.push_back(val);
        }

        return fix_variables(variables, values);
    }

    // bool bdd_mma_fixing::fix_variables()
    // {
    //     init_primal_solution();
    //     std::vector<double> total_min_marginals = this->total_min_marginals();

    //     while (true)
    //     {
    //         // nfixes++;
    //         // std::cout << "\rExpanded " << nfixes << " out of " << max_fixes << " search tree nodes.." << std::flush;
    //         // if (nfixes > max_fixes)
    //         //     return false;

    //         // double prev_lb = lower_bound_;
    //         // double min_progress = 1e-04;
    //         // int max_iter = 10;
    //         // for (int iter = 0; iter < max_iter; iter++)
    //         // {
    //         //     min_marginal_averaging_iteration();
    //         //     if (std::abs((lower_bound_-prev_lb)/prev_lb) < min_progress)
    //         //     {
    //         //         std::cout << "iters = " << iter+1 << ", " << std::flush;
    //         //         break;
    //         //     }
    //         //     prev_lb = lower_bound_;
    //         //     if (iter+1 == max_iter)
    //         //         std::cout << "iters = " << max_iter << ", " << std::flush;
    //         // }
    //         // std::cout << "lower bound = " << lower_bound_ << ". " << std::flush;

    //         // total_min_marginals = this->total_min_marginals();

    //         const double eps = 1e-12;

    //         double min_score = std::numeric_limits<double>::infinity();
    //         double max_score = -min_score;
    //         double prev_min_score = min_score;
    //         size_t min_var;
    //         size_t max_var;
    //         size_t nconflicts = 0;
    //         for (size_t var = 0; var < nr_variables(); var++)
    //         {
    //             if (is_fixed(var))
    //             {
    //                 if (primal_solution_[var] == 1 && total_min_marginals[var] > eps ||
    //                     primal_solution_[var] == 0 && total_min_marginals[var] < -eps)
    //                     nconflicts++;
    //                 continue;
    //             }

    //             if (total_min_marginals[var] < min_score)
    //             {
    //                 min_score = total_min_marginals[var];
    //                 min_var = var;
    //             }
    //             if (total_min_marginals[var] > max_score)
    //             {
    //                 max_score = total_min_marginals[var];
    //                 max_var = var;
    //             }
    //         }

    //         if (min_score == -std::numeric_limits<double>::infinity())
    //             std::cout << "Min score is -inf" << std::endl;

    //         if (min_score == std::numeric_limits<double>::infinity())
    //             return true;

    //         const size_t lsize = log_size();
    //         const double old_lb = lower_bound_;

    //         const bool small_enough = min_score < -eps;
            
    //         const char val = small_enough ? 1 : 0;
    //         const size_t best_var = small_enough ? min_var : max_var;
    //         const double best_score = small_enough ? min_score : max_score;

    //         bool feasible = fix_variable(best_var, val);
    //         this->backward_run();

    //         std::cout << "Marginal conflicts: " << nconflicts << ". Fixed var " << best_var << " (" << input_.get_var_name(best_var) << ")" << " to " << (int) val << ". score = " << best_score << std::endl;

    //         // TODO enable backtracking
    //         if (!feasible)
    //             return false;

    //         // double prev_lb = lower_bound_;
    //         // double min_progress = 1e-04;
    //         // int max_iter = 10;
    //         // for (int iter = 0; iter < max_iter; iter++)
    //         // {
    //         //     min_marginal_averaging_iteration();
    //         //     if (std::abs((lower_bound_-prev_lb)/prev_lb) < min_progress)
    //         //     {
    //         //         std::cout << "iters = " << iter+1 << ", " << std::flush;
    //         //         break;
    //         //     }
    //         //     prev_lb = lower_bound_;
    //         //     if (iter+1 == max_iter)
    //         //         std::cout << "iters = " << max_iter << ", " << std::flush;
    //         // }
    //         // std::cout << "lower bound = " << lower_bound_ << ". " << std::flush;

    //         // if (std::abs(lower_bound_-old_lb) > std::abs(5*min_score))
    //         // {
    //         //     revert_changes(lsize);
    //         //     bool feasible = fix_variable(best_var, 1-val);
    //         //     if (!feasible)
    //         //         return false;
    //         //     this->backward_run();

    //         //     std::cout << "Changed var " << best_var << " (" << input_.get_var_name(best_var) << ")" << " to " << (int) 1-val << std::endl;

    //         //     double prev_lb = lower_bound_;
    //         //     double min_progress = 1e-04;
    //         //     int max_iter = 10;
    //         //     for (int iter = 0; iter < max_iter; iter++)
    //         //     {
    //         //         min_marginal_averaging_iteration();
    //         //         if (std::abs((lower_bound_-prev_lb)/prev_lb) < min_progress)
    //         //         {
    //         //             std::cout << "iters = " << iter+1 << ", " << std::flush;
    //         //             break;
    //         //         }
    //         //         prev_lb = lower_bound_;
    //         //         if (iter+1 == max_iter)
    //         //             std::cout << "iters = " << max_iter << ", " << std::flush;
    //         //     }
    //         //     std::cout << "lower bound = " << lower_bound_ << ". " << std::flush;
    //         // }
    //     }
    // }

    void bdd_mma_fixing::min_marginal_averaging_forward()
    {
        std::vector<std::array<double,2>> min_marginals;
        for(std::size_t var=0; var<this->nr_variables(); ++var) {

            min_marginals.clear();
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                this->forward_step(var,bdd_index);
                if (!is_fixed(var))
                    min_marginals.push_back(min_marginal(var,bdd_index)); 
            }

            if (is_fixed(var))
                continue;

            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());

            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
            } 
        }
    }

    void bdd_mma_fixing::min_marginal_averaging_backward()
    {
        double lb = 0.0;
        std::vector<std::array<double,2>> min_marginals;

        for(long int var=this->nr_variables()-1; var>=0; --var) {

            min_marginals.clear();
            if (!is_fixed(var))
            {
                for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                    min_marginals.push_back(min_marginal(var,bdd_index)); 
                }
            }
            const std::array<double,2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());
            
            for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(var); ++bdd_index) {
                if (!is_fixed(var))
                    set_marginal(var,bdd_index,average_marginal,min_marginals[bdd_index]);
                this->backward_step(var, bdd_index);
                lb += lower_bound_backward(var,bdd_index);
            }
        }
        lower_bound_ = lb; 
    }

    void bdd_mma_fixing::min_marginal_averaging_iteration()
    {
        this->min_marginal_averaging_forward();
        this->min_marginal_averaging_backward();
    }
}
