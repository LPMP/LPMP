#pragma once

#include "bdd_min_marginal_averaging.h"
#include <array>
#include <vector>

namespace LPMP {

struct bdd_branch_node_exp_sum_entry
{
    std::array<double, 2> sum;
    std::array<double, 2> max;
};

template<typename DERIVED>
class bdd_branch_node_opt_smoothed_base : public bdd_branch_node_opt_base<DERIVED>
{
public:
    // below two are provided by base
    //double *variable_cost = nullptr;
    //double m = 0.0;

    double current_max = 0.0; // intermediate maximum value in the exp sum, used for stabilizing log-sum-exp computation. Also referred to as streamed log sum exp.

    // From C++20
    friend bool operator==(const bdd_branch_node_opt_smoothed_base &x, const bdd_branch_node_opt_smoothed_base &y);

    void smooth_backward_step();
    void smooth_forward_step();

    // Debug functions for checking correctness of forward and backward step
    double smooth_cost_from_first() const;
    double smooth_cost_from_terminal() const;

    bdd_branch_node_exp_sum_entry exp_sums() const;
};

class bdd_branch_node_opt_smoothed : public bdd_branch_node_opt_smoothed_base<bdd_branch_node_opt_smoothed>
{};

    //template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
    //class bdd_mma_base : public bdd_base<BDD_VARIABLE, BDD_BRANCH_NODE>, public bdd_solver_interface

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
class bdd_min_marginal_averaging_smoothed_base : public bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>
{
public:
    bdd_min_marginal_averaging_smoothed_base() {}
    bdd_min_marginal_averaging_smoothed_base(const bdd_min_marginal_averaging_smoothed_base &) = delete; // no copy constructor because of pointers in bdd_branch_node

    void init(const ILP_input &input);
    void init();

    template <typename ITERATOR>
    void set_costs(ITERATOR begin, ITERATOR end); // copy from bdd_min_marginal_averaging

    // TODO: delete, for now only for interface
    double lower_bound() { return -std::numeric_limits<double>::infinity(); } // TODO: not implemented yet
    void iteration() {}

    double smooth_lower_bound() { return -std::numeric_limits<double>::infinity(); } // TODO: not implemented yet
    double compute_smooth_lower_bound();
    double compute_smooth_lower_bound_forward();

    void smooth_forward_run();
    void smooth_backward_run();

    void smooth_averaging_pass_forward();
    void smooth_averaging_pass_backward();
    void smooth_iteration();

    template <typename STREAM>
    void export_dot(STREAM &s) const { return this->bdd_storage_.export_dot(s); }

private:
    void init_costs(); // copy from bdd_min_marginal_averaging
    double smooth_lower_bound_backward(const std::size_t var, const std::size_t bdd_index);
    double smooth_lower_bound_forward(const std::size_t var, const std::size_t bdd_index);

    void smooth_forward_step(const std::size_t var, const std::size_t bdd_index);
    void smooth_backward_step(const std::size_t var, const std::size_t bdd_index);

    void smooth_forward_step(const std::size_t var);
    void smooth_backward_step(const std::size_t var);

    bdd_branch_node_exp_sum_entry exp_sums(const std::size_t var, const std::size_t bdd_index) const;
    template <typename ITERATOR>
    static double average_exp_sums(ITERATOR exp_sums_begin, ITERATOR exp_sums_end);

    void update_Lagrange_multiplier(const std::size_t var, const std::size_t bdd_index, const bdd_branch_node_exp_sum_entry exp_sums, const double average_exp_sums);

    double lower_bound_ = -std::numeric_limits<double>::infinity();
    double cost_scaling_ = 1.0;
};

class bdd_min_marginal_averaging_smoothed : public bdd_min_marginal_averaging_smoothed_base<bdd_variable_mma, bdd_branch_node_opt_smoothed>
{};

////////////////////
// Implementation //
////////////////////

// TODO: C++20: make default operator==
bool operator==(const bdd_branch_node_exp_sum_entry &x, const bdd_branch_node_exp_sum_entry &y) 
{ 
    return x.sum == y.sum && x.max == y.max;
} 

template<typename DERIVED>
void bdd_branch_node_opt_smoothed_base<DERIVED>::smooth_forward_step()
{
    check_bdd_branch_node(*this);

    if (this->is_first()) {
        this->m = 1.0; // == exp(0);
        current_max = 0.0;
        return;
    }

    this->m = 0.0; 
    current_max = -std::numeric_limits<double>::infinity();

    // iterate over all incoming low edges
    {
        for(bdd_branch_node_opt_smoothed *cur = this->first_low_incoming; cur != nullptr; cur = cur->next_low_incoming)
        {
            if(cur->current_max < current_max)
            {
                this->m += std::exp(cur->current_max - current_max) * cur->m;
            }
            else
            {
                this->m *= std::exp(current_max - cur->current_max);
                this->m += cur->m;
                current_max = cur->current_max;
            }
            assert(std::isfinite(this->m));
        }
    }

    // iterate over all incoming high edges
    {
        for(bdd_branch_node_opt_smoothed *cur = this->first_high_incoming; cur != nullptr; cur = cur->next_high_incoming)
        {
            if(cur->current_max -*cur->variable_cost < current_max)
            {
                this->m += std::exp((cur->current_max -*cur->variable_cost) - current_max) * cur->m;
            } 
            else
            {
                this->m *= std::exp(current_max - (cur->current_max - *cur->variable_cost));
                this->m += cur->m; //1.0;
                current_max = cur->current_max - *cur->variable_cost;
            }
            //m += std::exp(-*(cur->variable_cost)) * cur->m;
            assert(std::isfinite(this->m));
        }
    }

    assert(std::isfinite(this->m));
    assert(this->m > 0.0);
    assert(this->m < 10000.0);
    check_bdd_branch_node(*this);
}

template<typename DERIVED>
void bdd_branch_node_opt_smoothed_base<DERIVED>::smooth_backward_step()
{
    check_bdd_branch_node(*this);

    // low edge
    const auto [low_cost, low_max] = [&]() -> std::array<double,2> {
        if (this->low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        else if (this->low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return {std::exp(0.0), 0.0};
        else
            return {this->low_outgoing->m, this->low_outgoing->current_max};
    }();

    // high edge
    const auto [high_cost, high_max] = [&]() -> std::array<double,2> {
        //if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
        //    return {0.0, -std::numeric_limits<double>::infinity()};
        //else if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
        //    return {std::exp(-*variable_cost - low_max)), -*variable_cost};
        //else
        //    return {std::exp(-*variable_cost) * high_outgoing->m, -*variable_cost + high_outgoing->current_max};
        if (this->high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        else if (this->high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return {std::exp(0.0), 0.0};
        else
            return {this->high_outgoing->m, this->high_outgoing->current_max};
    }();

    assert(std::isfinite(low_cost));
    assert(std::isfinite(high_cost));
    this->m = low_cost;
    current_max = low_max;
    //current_max = 0;
    //m += std::exp(-*variable_cost)*high_cost;
    //return;

    if (std::isfinite(high_max))
    {
        if (high_max - *this->variable_cost < low_max)
        {
            this->m += std::exp((high_max - *this->variable_cost) - current_max) * high_cost;
            //m += std::exp(-*variable_cost - current_max) * high_cost;
        }
        else
        {
            this->m *= std::exp(current_max - (high_max - *this->variable_cost));
            this->m += high_cost;//1.0;
            current_max = high_max - *this->variable_cost;
        }
    }
    //m = low_cost + high_cost;

    if (!std::isfinite(this->m) || !(std::abs(this->m) < 10000.0))
    {
        std::cout << low_cost << "," << low_max << "\n";
        std::cout << high_cost << "," << high_max << "\n";
    }
    assert(std::isfinite(this->m));
    assert(std::abs(this->m) < 10000.0);
    assert(std::isfinite(current_max));

    check_bdd_branch_node(*this);
    //assert(std::abs(m - cost_from_terminal()) <= 1e-8);
}

template<typename DERIVED>
double bdd_branch_node_opt_smoothed_base<DERIVED>::smooth_cost_from_first() const
{
    double c = 0.0;

    if (this->is_first())
        return 0.0;

    // iterate over all incoming low edges
    for (bdd_branch_node_opt_smoothed *cur = this->first_low_incoming; cur != nullptr; cur = cur->next_low_incoming)
        c += cur->smooth_cost_from_first();

    // iterate over all incoming high edges
    for (bdd_branch_node_opt_smoothed *cur = this->first_high_incoming; cur != nullptr; cur = cur->next_high_incoming)
        c += std::exp(-*cur->variable_cost) * cur->smooth_cost_from_first(); // ??

    return c;
}

template<typename DERIVED>
double bdd_branch_node_opt_smoothed_base<DERIVED>::smooth_cost_from_terminal() const
{
    // TODO: only works if no bdd nodes skips variables
    // low edge
    const double low_cost = [&]() {
        if (this->low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return 0.0;
        else if (this->low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return 1.0;
        else
            return this->low_outgoing->smooth_cost_from_terminal();
    }();

    // high edge
    const double high_cost = [&]() {
        if (this->high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return 0.0;
        else if (this->high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return std::exp(-*this->variable_cost);
        else
            return this->high_outgoing->smooth_cost_from_terminal() + std::exp(-*this->variable_cost);
    }();

    return low_cost + high_cost;
}

template<typename DERIVED>
bdd_branch_node_exp_sum_entry bdd_branch_node_opt_smoothed_base<DERIVED>::exp_sums() const
{
    check_bdd_branch_node(*this);

    // assert(std::abs(m - cost_from_first()) <= 1e-8);
    if (!bdd_branch_node_opt_smoothed::is_terminal(this->low_outgoing))
    {
        //assert(std::abs(low_outgoing->m - low_outgoing->cost_from_terminal()) <= 1e-8);
    }
    if (!bdd_branch_node_opt_smoothed::is_terminal(this->high_outgoing))
    {
        //assert(std::abs(high_outgoing->m - high_outgoing->cost_from_terminal()) <= 1e-8);
    }

    bdd_branch_node_exp_sum_entry e;

    std::tie(e.sum[0], e.max[0]) = [&]() -> std::tuple<double, double> {
        if (this->low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        if (this->low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return {this->m, current_max};
        else
            return {this->m * this->low_outgoing->m, current_max + this->low_outgoing->current_max};
    }();

    std::tie(e.sum[1], e.max[1]) = [&]() -> std::tuple<double, double> {
        if (this->high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        if (this->high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
        {
            //const double new_max = std::max(this->current_max, this->current_max - *this->variable_cost);
            //return {this->m * std::exp(-*this->variable_cost + this->current_max - new_max), new_max};
            return {this->m, this->current_max - *this->variable_cost};
        }
        else
        {
            //const double new_max = std::max({this->current_max, this->current_max - *this->variable_cost, this->high_outgoing->current_max});
            //return {this->m * std::exp(-*this->variable_cost + this->current_max + this->high_outgoing->current_max - new_max) * this->high_outgoing->m, new_max};
            return {this->m * this->high_outgoing->m, this->current_max - *this->variable_cost + this->high_outgoing->current_max};
        }
    }();

    assert(std::isfinite(e.sum[0]));
    assert(std::isfinite(e.sum[1]));
    assert(e.sum[0] >= 0.0);
    assert(e.sum[1] >= 0.0);
    if(!(e.sum[0] > 0 || e.sum[1] > 0))
        std::cout << e.sum[0] << ";" << e.sum[1] << "\n";
    assert(e.sum[0] > 0 || e.sum[1] > 0);
    if(this->low_outgoing != bdd_branch_node_opt_smoothed::terminal_0() && this->low_outgoing != bdd_branch_node_opt_smoothed::terminal_1())
    {
        assert(e.sum[0] > 0);
        assert(std::isfinite(e.max[0]));
    }
    if (this->high_outgoing != bdd_branch_node_opt_smoothed::terminal_0() && this->high_outgoing != bdd_branch_node_opt_smoothed::terminal_1())
    {
        assert(e.sum[1] > 0);
        assert(std::isfinite(e.max[1]));
    }

    assert(std::abs(e.sum[0]) < 10000.0);
    assert(std::abs(e.sum[1]) < 10000.0);

    //std::cout << e.sum[0] << "," << e.sum[1] << " ; " << e.max[0] << "," << e.max[1] << "\n";
    return e;
}

/*
std::array<double, 2> bdd_branch_node_opt_smoothed::min_marginal_debug() const
{
    check_bdd_branch_node(*this);

    const double m_debug = cost_from_first();

    const double m0 = [&]() {
        if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return std::numeric_limits<double>::infinity();
        if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return m_debug;
        return m_debug + this->low_outgoing->cost_from_terminal();
    }();

    const double m1 = [&]() {
        if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return std::numeric_limits<double>::infinity();
        if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return m_debug + *this->variable_cost;
        return m_debug + *this->variable_cost + this->high_outgoing->cost_from_terminal();
    }();

    assert(std::isfinite(std::min(m0, m1)));

    return {m0, m1};
}
*/

/////////////////////////////////////////
// bdd_min_marginal_averaging_smoothed //
/////////////////////////////////////////

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init_costs()
{
    for (size_t var = 0; var < this->nr_variables(); var++)
        for (size_t bdd_index = 0; bdd_index < this->nr_bdds(var); bdd_index++)
        {
            auto &bdd_var = this->bdd_variables_(var, bdd_index);
            for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                this->bdd_branch_nodes_[node_index].variable_cost = &bdd_var.cost;
        }
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init(const ILP_input &input)
{
    bdd_base<bdd_variable_mma, bdd_branch_node_opt_smoothed>::init(input);
    init_costs();
    set_costs(input.objective().begin(), input.objective().end());
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::init()
{
    bdd_base<bdd_variable_mma, bdd_branch_node_opt_smoothed>::init();
    init_costs();
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
template <typename ITERATOR>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::set_costs(ITERATOR begin, ITERATOR end)
{
    const double max_entry = *std::max_element(begin, end);
    const double min_entry = *std::min_element(begin, end);
    assert(std::abs(max_entry) > 0.0 || std::abs(min_entry) > 0.0);
    cost_scaling_ = 1.0; //std::max(std::abs(max_entry), std::abs(min_entry));
    // distribute costs to bdds uniformly
    for (std::size_t v = 0; v < this->nr_variables(); ++v)
    {
        assert(this->bdd_variables_[v].size() > 0);
        for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(v); ++bdd_index)
        {
            const double cost = v < std::distance(begin, end) ? *(begin + v) / this->nr_bdds(v) : 0.0;
            this->bdd_variables_(v, bdd_index).cost = 1.0/cost_scaling_ * cost;
            assert(!std::isnan(this->bdd_variables_(v, bdd_index).cost));
        }
    }
    for (const auto &bdd : this->bdd_branch_nodes_)
    {
        assert(!std::isnan(*bdd.variable_cost));
    }
    smooth_backward_run();
    const double init_lb = compute_smooth_lower_bound();
    std::cout << "initial lb = " << init_lb << "\n";
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
double bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::compute_smooth_lower_bound()
{
    double lb = 0.0;
    for (std::ptrdiff_t var = this->nr_variables() - 1; var >= 0; --var)
        for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
        {
            lb += smooth_lower_bound_backward(var, bdd_index);
            //std::cout << "intermediate lb =  " << lb << "\n";
        }
    return lb;
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
double bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::compute_smooth_lower_bound_forward()
{
    double lb = 0.0;
    for (std::size_t var = 0; var < this->nr_variables(); ++var)
        for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
            lb += smooth_lower_bound_forward(var, bdd_index);
    return lb;
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
double bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_lower_bound_forward(const std::size_t var, const std::size_t bdd_index)
{
    const auto &bdd_var = this->bdd_variables_(var, bdd_index);
    if (this->last_variable_of_bdd(var, bdd_index))
    {
        double lb = 0;
        double cur_max = -std::numeric_limits<double>::infinity();
        const std::size_t first_node_index = bdd_var.first_node_index;
        const std::size_t last_node_index = bdd_var.last_node_index;
        for (std::size_t i = first_node_index; i < last_node_index; ++i)
        {
            const auto &node = this->bdd_branch_nodes_[i];
            assert(node.low_outgoing == node.terminal_1() || node.high_outgoing == node.terminal_1());
            if(node.low_outgoing == node.terminal_1())
            {
                if(node.current_max < cur_max)
                {
                    lb += std::exp(node.current_max - cur_max) * node.m;
                }
                else
                {
                    lb *= std::exp(cur_max - node.current_max);
                    lb += node.m; //1.0;
                    cur_max = node.current_max; 
                } 
            }

            if(node.high_outgoing == node.terminal_1())
            {
                if(node.current_max - *node.variable_cost < cur_max)
                {
                    lb += std::exp(node.current_max - *node.variable_cost - cur_max) * node.m; 
                }
                else
                {
                    lb *= std::exp(cur_max - (node.current_max - *node.variable_cost));
                    lb += node.m;
                    cur_max = node.current_max - *node.variable_cost;
                }
            }
        }
        return -cost_scaling_ * (std::log(lb) + cur_max);
    }
    else
    {
        return 0.0;
    }
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
double bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_lower_bound_backward(const std::size_t var, const std::size_t bdd_index)
{
    const auto &bdd_var = this->bdd_variables_(var, bdd_index);
    if (this->first_variable_of_bdd(var, bdd_index))
    {
        assert(bdd_var.nr_bdd_nodes() == 1);
        const auto &bdd = this->get_bdd_branch_node(var, bdd_index, 0);
        //std::cout << "bdd lb = " << - std::log(bdd.m) << "\n";
        //std::cout << "bdd cur max = " << bdd.current_max << "\n";
        assert(std::isfinite(std::log(bdd.m)));
        return -cost_scaling_ * (std::log(bdd.m) + bdd.current_max);
    }
    else
    {
        return 0.0;
    }
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_forward_step(const std::size_t var, const std::size_t bdd_index)
{
    assert(var < this->bdd_variables_.size());
    assert(bdd_index < this->bdd_variables_[var].size());

    auto &bdd_var = this->bdd_variables_(var, bdd_index);
    assert(var != 0 || bdd_var.prev == nullptr);

    // iterate over all bdd nodes and make forward step
    const std::size_t first_node_index = bdd_var.first_node_index;
    const std::size_t last_node_index = bdd_var.last_node_index;
    for (std::size_t i = first_node_index; i < last_node_index; ++i)
        this->bdd_branch_nodes_[i].smooth_forward_step();
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_backward_step(const std::size_t var, const std::size_t bdd_index)
{
    assert(var < this->bdd_variables_.size());
    assert(bdd_index < this->bdd_variables_[var].size());

    auto &bdd_var = this->bdd_variables_(var, bdd_index);
    assert(var + 1 != this->nr_variables() || bdd_var.next == nullptr);
    // iterate over all bdd nodes and make forward step
    const std::ptrdiff_t first_node_index = bdd_var.first_node_index;
    const std::ptrdiff_t last_node_index = bdd_var.last_node_index;
    for (std::ptrdiff_t i = last_node_index - 1; i >= first_node_index; --i)
    {
        check_bdd_branch_node(this->bdd_branch_nodes_[i], var + 1 == this->nr_variables(), var == 0);
        this->bdd_branch_nodes_[i].smooth_backward_step();
    }
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_forward_step(const std::size_t var)
{
    assert(var < this->bdd_variables_.size());
    for (std::size_t bdd_index = 0; bdd_index < this->bdd_variables_[var].size(); ++bdd_index)
        smooth_forward_step(var, bdd_index);
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_backward_step(const std::size_t var)
{
    assert(var < this->bdd_variables_.size());
    for (std::ptrdiff_t bdd_index = this->bdd_variables_[var].size()-1; bdd_index >= 0; --bdd_index)
        smooth_backward_step(var, bdd_index);
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_forward_run()
{
    for (std::size_t var = 0; var < this->nr_variables(); ++var)
        smooth_forward_step(var);
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_backward_run()
{
    for (std::ptrdiff_t var = this->nr_variables() - 1; var >= 0; --var)
        smooth_backward_step(var);
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_iteration()
{
    smooth_averaging_pass_forward();
    smooth_averaging_pass_backward();
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
bdd_branch_node_exp_sum_entry bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::exp_sums(const std::size_t var, const std::size_t bdd_index) const
{
    assert(var < this->nr_variables());
    assert(bdd_index < this->nr_bdds(var));
    bdd_branch_node_exp_sum_entry e;
    e.sum = {0.0, 0.0};
    e.max = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    std::array<double, 2> s = {0.0, 0.0};
    std::array<double, 2> current_max = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    const auto &bdd_var = this->bdd_variables_(var, bdd_index);
    for (std::size_t bdd_node_index = bdd_var.first_node_index; bdd_node_index < bdd_var.last_node_index; ++bdd_node_index)
    {
        const auto &bdd = this->bdd_branch_nodes_[bdd_node_index];
        const auto bdd_exp_sums = bdd.exp_sums();
        //std::cout << "var" << var << ", bdd index = " << bdd_index << ": " << "\n";
        //std::cout << bdd_exp_sums.sum[0] << "," << bdd_exp_sums.sum[1] << ";" << bdd_exp_sums.max[0] << "," << bdd_exp_sums.max[1] << "\n";

        if (bdd_exp_sums.sum[0] > 0)
        {
            if (current_max[0] > bdd_exp_sums.max[0])
            {
                s[0] += bdd_exp_sums.sum[0] * std::exp(bdd_exp_sums.max[0] - current_max[0]);
            }
            else
            {
                s[0] *= std::exp(current_max[0] - bdd_exp_sums.max[0]);
                s[0] += bdd_exp_sums.sum[0];
                current_max[0] = bdd_exp_sums.max[0];
            }
        }

        if (bdd_exp_sums.sum[1] > 0)
        {
            if (current_max[1] > bdd_exp_sums.max[1])
            {
                s[1] += bdd_exp_sums.sum[1] * std::exp(bdd_exp_sums.max[1] - current_max[1]);
            }
            else
            {
                s[1] *= std::exp(current_max[1] - bdd_exp_sums.max[1]);
                s[1] += bdd_exp_sums.sum[1];
                current_max[1] = bdd_exp_sums.max[1];
            }
        }
        assert(std::isfinite(s[0]));
        assert(std::isfinite(s[1]));
    }
    //std::cout << s[0] << "," << s[1] << ";" << current_max[0] << "," << current_max[1] << "\n";
    assert(s[0] > 0);
    assert(s[1] > 0);
    return {s, current_max};
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
template <typename ITERATOR>
double bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::average_exp_sums(ITERATOR exp_sums_begin, ITERATOR exp_sums_end)
{
    const double n = std::distance(exp_sums_begin, exp_sums_end);
    double average = 0.0;
    for(auto it=exp_sums_begin; it!=exp_sums_end; ++it)
    {
        const double sum_0 = (*it).sum[0];
        const double sum_1 = (*it).sum[1];
        assert(std::abs(sum_0) <= 100000.0);
        assert(std::abs(sum_1) <= 100000.0);
        assert(std::abs((*it).max[0]) <= 100000.0);
        assert(std::abs((*it).max[1]) <= 100000.0);
        //average += 1.0 / double(n) * (std::log(sum_1 / sum_0) - ((*it).max[1] - (*it).max[0]));
        average -= 1.0 / double(n) * (std::log(sum_1) + (*it).max[1] - (std::log(sum_0) + (*it).max[0]));
        assert(std::isfinite(average));
    }

    return average;
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::update_Lagrange_multiplier(const std::size_t var, const std::size_t bdd_index, const bdd_branch_node_exp_sum_entry exp_sums, const double average_exp_sums)
{
    assert(var < this->nr_variables());
    assert(bdd_index < this->nr_bdds(var));
    assert(exp_sums == this->exp_sums(var, bdd_index));
    assert(std::isfinite(exp_sums.sum[0]) && std::isfinite(exp_sums.sum[1]));
    assert(std::isfinite(average_exp_sums));

    auto &bdd_var = this->bdd_variables_(var, bdd_index);
    const double diff = (std::log(exp_sums.sum[1]) + exp_sums.max[1] - (std::log(exp_sums.sum[0]) + exp_sums.max[0])) + average_exp_sums;
    assert(std::isfinite(diff));
    bdd_var.cost += diff;
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_averaging_pass_forward()
{
    std::vector<bdd_branch_node_exp_sum_entry> exp_sums;
    for (std::size_t var = 0; var < this->nr_variables(); ++var)
    {
        if (this->nr_bdds(var) > 1)
        {
            // collect min marginals
            exp_sums.clear();
            for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
            {
                smooth_forward_step(var, bdd_index);
                exp_sums.push_back(this->exp_sums(var, bdd_index));
            }

            const double average_exp_sums = this->average_exp_sums(exp_sums.begin(), exp_sums.end());

            // set marginals in each bdd so min marginals match each other
            for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
                update_Lagrange_multiplier(var, bdd_index, exp_sums[bdd_index], average_exp_sums);
        }
        else
            smooth_forward_step(var, 0);
    }
}

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::smooth_averaging_pass_backward()
{
    double lb = 0.0;
    std::vector<bdd_branch_node_exp_sum_entry> exp_sums;

    for (long int var = this->nr_variables() - 1; var >= 0; --var)
    {
        if(this->nr_bdds(var) > 1)
        {
            exp_sums.clear();
            for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
                exp_sums.push_back(this->exp_sums(var, bdd_index));

            const double average_exp_sums = this->average_exp_sums(exp_sums.begin(), exp_sums.end());

            for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
            {
                update_Lagrange_multiplier(var, bdd_index, exp_sums[bdd_index], average_exp_sums);
                smooth_backward_step(var, bdd_index);
                lb += smooth_lower_bound_backward(var, bdd_index);
            }
        }
        else
        {
            smooth_backward_step(var, 0);
            lb += smooth_lower_bound_backward(var, 0);
        }
    }

    lower_bound_ = lb;
}

} // namespace LPMP