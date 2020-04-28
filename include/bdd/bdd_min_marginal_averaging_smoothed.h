#pragma once

#include "bdd_min_marginal_averaging.h"
#include <array>
#include <vector>

namespace LPMP {


/////////////////////////////////////////
// bdd_min_marginal_averaging_smoothed //
/////////////////////////////////////////

template<typename BDD_VARIABLE, typename BDD_BRANCH_NODE>
class bdd_min_marginal_averaging_smoothed_base : public bdd_mma_base<BDD_VARIABLE, BDD_BRANCH_NODE>
{
public:
    bdd_min_marginal_averaging_smoothed_base() {}
    bdd_min_marginal_averaging_smoothed_base(const bdd_min_marginal_averaging_smoothed_base &) = delete; // no copy constructor because of pointers in bdd_branch_node

    double smooth_lower_bound() { return -std::numeric_limits<double>::infinity(); } // TODO: not implemented yet
    double compute_smooth_lower_bound();
    double compute_smooth_lower_bound_forward();

    void smooth_forward_run();
    void smooth_backward_run();

    void smooth_averaging_pass_forward();
    void smooth_averaging_pass_backward();
    void smooth_iteration();

    void set_cost_scaling(const double scale);

private:
    //void init_costs(); // copy from bdd_min_marginal_averaging
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
void bdd_min_marginal_averaging_smoothed_base<BDD_VARIABLE, BDD_BRANCH_NODE>::set_cost_scaling(const double scale)
{
    assert(scale > 0.0);

    // rescale Lagrange multipliers
    for (std::size_t var = 0; var < this->nr_variables(); ++var)
    {
        for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
        {
            auto &bdd_var = this->bdd_variables_(var, bdd_index);
            bdd_var.cost *= cost_scaling_ / scale;
        }
    }

    cost_scaling_ = scale;
    // TODO: backward run must be performed
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