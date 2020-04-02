#pragma once

#include "bdd_min_marginal_averaging.h"

namespace LPMP {

class bdd_branch_node_opt_smoothed : public bdd_branch_node<bdd_branch_node_opt_smoothed>
{
public:
    double *variable_cost = nullptr;
    double m = 0.0; // intermediate value of exp-sum from either terminal or first node (depending on algorithm state)
    double current_max = 0.0; // intermediate maximum value in the exp sum, used for stabilizing log-sum-exp computation. Also referred to as streamed log sum exp.

    // From C++20
    friend bool operator==(const bdd_branch_node_opt_smoothed &x, const bdd_branch_node_opt_smoothed &y);

    void backward_step();
    void forward_step();

    // Debug functions for checking correctness of forward and backward step
    double cost_from_first() const;
    double cost_from_terminal() const;

    struct exp_sum_entry
    {
        std::array<double,2> sum;
        std::array<double,2> max;
    };
    exp_sum_entry exp_sums() const;
};

// TODO: C++20: make default operator==
bool operator==(const bdd_branch_node_opt_smoothed::exp_sum_entry &x, const bdd_branch_node_opt_smoothed::exp_sum_entry &y) { return x.sum == y.sum && x.max == y.max; } 

class bdd_min_marginal_averaging_smoothed : public bdd_base<bdd_variable_mma, bdd_branch_node_opt_smoothed>, public bdd_solver_interface
{
public:
    bdd_min_marginal_averaging_smoothed() {}
    bdd_min_marginal_averaging_smoothed(const bdd_min_marginal_averaging_smoothed &) = delete; // no copy constructor because of pointers in bdd_branch_node

    void init(const ILP_input &input);
    void init();

    template <typename ITERATOR>
    void set_costs(ITERATOR begin, ITERATOR end); // copy from bdd_min_marginal_averaging

    double lower_bound() { return -std::numeric_limits<double>::infinity(); } // TODO: not implemented yet
    double compute_lower_bound();
    double compute_lower_bound_forward();

    void averaging_pass_forward();
    void averaging_pass_backward();
    void iteration();

    template <typename STREAM>
    void export_dot(STREAM &s) const { return this->bdd_storage_.export_dot(s); }

private:
    void init_costs(); // copy from bdd_min_marginal_averaging
    double lower_bound_backward(const std::size_t var, const std::size_t bdd_index);
    double lower_bound_forward(const std::size_t var, const std::size_t bdd_index);

    bdd_branch_node_opt_smoothed::exp_sum_entry exp_sums(const std::size_t var, const std::size_t bdd_index) const;
    template <typename ITERATOR>
    static double average_exp_sums(ITERATOR exp_sums_begin, ITERATOR exp_sums_end);

    void update_Lagrange_multiplier(const std::size_t var, const std::size_t bdd_index, const bdd_branch_node_opt_smoothed::exp_sum_entry exp_sums, const double average_exp_sums);

    double lower_bound_ = -std::numeric_limits<double>::infinity();
    double cost_scaling_ = 1.0;
};

////////////////////
// Implementation //
////////////////////

void bdd_branch_node_opt_smoothed::forward_step()
{
    check_bdd_branch_node(*this);

    if (is_first()) {
        m = 1.0; // == exp(0);
        current_max = 0.0;
        return;
    }

    m = 0.0; 
    current_max = -std::numeric_limits<double>::infinity();

    // iterate over all incoming low edges
    {
        // TODO: replace??
        //for(bdd_branch_node_opt_smoothed *cur = first_low_incoming; cur != nullptr; cur = cur->next_low_incoming) {}
        bdd_branch_node_opt_smoothed *cur = first_low_incoming;
        while (cur != nullptr)
        {
            if(cur->current_max < current_max)
            {
                m += std::exp(cur->current_max - current_max) * cur->m;
            }
            else
            {
                m *= std::exp(current_max - cur->current_max);
                m += cur->m;
                current_max = cur->current_max;
            }
            assert(std::isfinite(m));
            cur = cur->next_low_incoming;
        }
    }

    // iterate over all incoming high edges
    {
        // TODO: replace??
        //for(bdd_branch_node_opt_smoothed *cur = first_high_incoming; cur != nullptr; cur = cur->next_high_incoming) {}
        bdd_branch_node_opt_smoothed *cur = first_high_incoming;
        while (cur != nullptr)
        {
            if(cur->current_max -*cur->variable_cost < current_max)
            {
                m += std::exp((cur->current_max -*cur->variable_cost) - current_max) * cur->m;
            } 
            else
            {
                m *= std::exp(current_max - (cur->current_max - *cur->variable_cost));
                m += cur->m; //1.0;
                current_max = cur->current_max - *cur->variable_cost;
            }
            //m += std::exp(-*(cur->variable_cost)) * cur->m;
            assert(std::isfinite(m));
            cur = cur->next_high_incoming;
        }
    }

    assert(std::isfinite(m));
    assert(m > 0.0);
    assert(m < 10000.0);
    check_bdd_branch_node(*this);
}

void bdd_branch_node_opt_smoothed::backward_step()
{
    check_bdd_branch_node(*this);

    // low edge
    const auto [low_cost, low_max] = [&]() -> std::array<double,2> {
        if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        else if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return {std::exp(0.0), 0.0};
        else
            return {low_outgoing->m, low_outgoing->current_max};
    }();

    // high edge
    const auto [high_cost, high_max] = [&]() -> std::array<double,2> {
        //if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
        //    return {0.0, -std::numeric_limits<double>::infinity()};
        //else if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
        //    return {std::exp(-*variable_cost - low_max)), -*variable_cost};
        //else
        //    return {std::exp(-*variable_cost) * high_outgoing->m, -*variable_cost + high_outgoing->current_max};
        if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        else if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return {std::exp(0.0), 0.0};
        else
            return {high_outgoing->m, high_outgoing->current_max};
    }();

    assert(std::isfinite(low_cost));
    assert(std::isfinite(high_cost));
    m = low_cost;
    current_max = low_max;
    //current_max = 0;
    //m += std::exp(-*variable_cost)*high_cost;
    //return;

    if (std::isfinite(high_max))
    {
        if (high_max - *variable_cost < low_max)
        {
            m += std::exp((high_max - *variable_cost) - current_max) * high_cost;
            //m += std::exp(-*variable_cost - current_max) * high_cost;
        }
        else
        {
            m *= std::exp(current_max - (high_max - *variable_cost));
            m += high_cost;//1.0;
            current_max = high_max - *variable_cost;
        }
    }
    //m = low_cost + high_cost;

    if (!std::isfinite(m) || !(std::abs(m) < 10000.0))
    {
        std::cout << low_cost << "," << low_max << "\n";
        std::cout << high_cost << "," << high_max << "\n";
    }
    assert(std::isfinite(m));
    assert(std::abs(m) < 10000.0);
    assert(std::isfinite(current_max));

    check_bdd_branch_node(*this);
    //assert(std::abs(m - cost_from_terminal()) <= 1e-8);
}

double bdd_branch_node_opt_smoothed::cost_from_first() const
{
    double c = 0.0;

    if (is_first())
        return 0.0;

    // iterate over all incoming low edges
    {
        bdd_branch_node_opt_smoothed *cur = first_low_incoming;
        while (cur != nullptr)
        {
            c += cur->cost_from_first();
            cur = cur->next_low_incoming;
        }
    }

    // iterate over all incoming high edges
    {
        bdd_branch_node_opt_smoothed *cur = first_high_incoming;
        while (cur != nullptr)
        {
            c += std::exp(-*cur->variable_cost) * cur->cost_from_first(); // ??
            cur = cur->next_high_incoming;
        }
    }

    return c;
}

double bdd_branch_node_opt_smoothed::cost_from_terminal() const
{
    // TODO: only works if no bdd nodes skips variables
    // low edge
    const double low_cost = [&]() {
        if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return 0.0;
        else if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return 1.0;
        else
            return low_outgoing->cost_from_terminal();
    }();

    // high edge
    const double high_cost = [&]() {
        if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return 0.0;
        else if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return std::exp(-*variable_cost);
        else
            return high_outgoing->cost_from_terminal() + std::exp(-*variable_cost);
    }();

    return low_cost + high_cost;
}

bdd_branch_node_opt_smoothed::exp_sum_entry bdd_branch_node_opt_smoothed::exp_sums() const
{
    check_bdd_branch_node(*this);

    // assert(std::abs(m - cost_from_first()) <= 1e-8);
    if (!bdd_branch_node_opt_smoothed::is_terminal(low_outgoing))
    {
        //assert(std::abs(low_outgoing->m - low_outgoing->cost_from_terminal()) <= 1e-8);
    }
    if (!bdd_branch_node_opt_smoothed::is_terminal(high_outgoing))
    {
        //assert(std::abs(high_outgoing->m - high_outgoing->cost_from_terminal()) <= 1e-8);
    }

    bdd_branch_node_opt_smoothed::exp_sum_entry e;

    std::tie(e.sum[0], e.max[0]) = [&]() -> std::tuple<double, double> {
        if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        if (low_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
            return {this->m, current_max};
        else
            return {this->m * this->low_outgoing->m, current_max + this->low_outgoing->current_max};
    }();

    std::tie(e.sum[1], e.max[1]) = [&]() -> std::tuple<double, double> {
        if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_0())
            return {0.0, -std::numeric_limits<double>::infinity()};
        if (high_outgoing == bdd_branch_node_opt_smoothed::terminal_1())
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
    if(low_outgoing != bdd_branch_node_opt_smoothed::terminal_0() && low_outgoing != bdd_branch_node_opt_smoothed::terminal_1())
    {
        assert(e.sum[0] > 0);
        assert(std::isfinite(e.max[0]));
    }
    if (high_outgoing != bdd_branch_node_opt_smoothed::terminal_0() && high_outgoing != bdd_branch_node_opt_smoothed::terminal_1())
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

void bdd_min_marginal_averaging_smoothed::init_costs()
{
    for (size_t var = 0; var < nr_variables(); var++)
        for (size_t bdd_index = 0; bdd_index < nr_bdds(var); bdd_index++)
        {
            auto &bdd_var = bdd_variables_(var, bdd_index);
            for (size_t node_index = bdd_var.first_node_index; node_index < bdd_var.last_node_index; node_index++)
                bdd_branch_nodes_[node_index].variable_cost = &bdd_var.cost;
        }
}

void bdd_min_marginal_averaging_smoothed::init(const ILP_input &input)
{
    bdd_base<bdd_variable_mma, bdd_branch_node_opt_smoothed>::init(input);
    init_costs();
    set_costs(input.objective().begin(), input.objective().end());
}

void bdd_min_marginal_averaging_smoothed::init()
{
    bdd_base<bdd_variable_mma, bdd_branch_node_opt_smoothed>::init();
    init_costs();
}

template <typename ITERATOR>
void bdd_min_marginal_averaging_smoothed::set_costs(ITERATOR begin, ITERATOR end)
{
    const double max_entry = *std::max_element(begin, end);
    const double min_entry = *std::min_element(begin, end);
    assert(std::abs(max_entry) > 0.0 || std::abs(min_entry) > 0.0);
    cost_scaling_ = 1.0; //std::max(std::abs(max_entry), std::abs(min_entry));
    // distribute costs to bdds uniformly
    for (std::size_t v = 0; v < nr_variables(); ++v)
    {
        assert(bdd_variables_[v].size() > 0);
        for (std::size_t bdd_index = 0; bdd_index < nr_bdds(v); ++bdd_index)
        {
            const double cost = v < std::distance(begin, end) ? *(begin + v) / nr_bdds(v) : 0.0;
            bdd_variables_(v, bdd_index).cost = 1.0/cost_scaling_ * cost;
            assert(!std::isnan(bdd_variables_(v, bdd_index).cost));
        }
    }
    for (const auto &bdd : bdd_branch_nodes_)
    {
        assert(!std::isnan(*bdd.variable_cost));
    }
    backward_run();
    const double init_lb = compute_lower_bound();
    std::cout << "initial lb = " << init_lb << "\n";
}

double bdd_min_marginal_averaging_smoothed::compute_lower_bound()
{
    double lb = 0.0;
    for (std::ptrdiff_t var = nr_variables() - 1; var >= 0; --var)
        for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
        {
            lb += lower_bound_backward(var, bdd_index);
            //std::cout << "intermediate lb =  " << lb << "\n";
        }
    return lb;
}

double bdd_min_marginal_averaging_smoothed::compute_lower_bound_forward()
{
    double lb = 0.0;
    for(std::size_t var=0; var<nr_variables(); ++var)
        for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
            lb += lower_bound_forward(var, bdd_index);
    return lb;
}

double bdd_min_marginal_averaging_smoothed::lower_bound_forward(const std::size_t var, const std::size_t bdd_index)
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

double bdd_min_marginal_averaging_smoothed::lower_bound_backward(const std::size_t var, const std::size_t bdd_index)
{
    const auto &bdd_var = this->bdd_variables_(var, bdd_index);
    if (this->first_variable_of_bdd(var, bdd_index))
    {
        assert(bdd_var.nr_bdd_nodes() == 1);
        const auto &bdd = get_bdd_branch_node(var, bdd_index, 0);
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

void bdd_min_marginal_averaging_smoothed::iteration()
{
    averaging_pass_forward();
    averaging_pass_backward();
}

bdd_branch_node_opt_smoothed::exp_sum_entry bdd_min_marginal_averaging_smoothed::exp_sums(const std::size_t var, const std::size_t bdd_index) const
{
    assert(var < nr_variables());
    assert(bdd_index < nr_bdds(var));
    bdd_branch_node_opt_smoothed::exp_sum_entry e;
    e.sum = {0.0, 0.0};
    e.max = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    std::array<double, 2> s = {0.0, 0.0};
    std::array<double, 2> current_max = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    const auto &bdd_var = bdd_variables_(var, bdd_index);
    for (std::size_t bdd_node_index = bdd_var.first_node_index; bdd_node_index < bdd_var.last_node_index; ++bdd_node_index)
    {
        const auto &bdd = bdd_branch_nodes_[bdd_node_index];
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

template <typename ITERATOR>
double bdd_min_marginal_averaging_smoothed::average_exp_sums(ITERATOR exp_sums_begin, ITERATOR exp_sums_end)
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

void bdd_min_marginal_averaging_smoothed::update_Lagrange_multiplier(const std::size_t var, const std::size_t bdd_index, const bdd_branch_node_opt_smoothed::exp_sum_entry exp_sums, const double average_exp_sums)
{
    assert(var < nr_variables());
    assert(bdd_index < nr_bdds(var));
    assert(exp_sums == this->exp_sums(var, bdd_index));
    assert(std::isfinite(exp_sums.sum[0]) && std::isfinite(exp_sums.sum[1]));
    assert(std::isfinite(average_exp_sums));

    auto &bdd_var = bdd_variables_(var, bdd_index);
    const double diff = (std::log(exp_sums.sum[1]) + exp_sums.max[1] - (std::log(exp_sums.sum[0]) + exp_sums.max[0])) + average_exp_sums;
    assert(std::isfinite(diff));
    bdd_var.cost += diff;
}

void bdd_min_marginal_averaging_smoothed::averaging_pass_forward()
{
    std::vector<bdd_branch_node_opt_smoothed::exp_sum_entry> exp_sums;
    for (std::size_t var = 0; var < nr_variables(); ++var)
    {
        if(nr_bdds(var) > 1) {
            // collect min marginals
            exp_sums.clear();
            for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
            {
                forward_step(var, bdd_index);
                exp_sums.push_back(this->exp_sums(var, bdd_index));
            }

            const double average_exp_sums = this->average_exp_sums(exp_sums.begin(), exp_sums.end());

            // set marginals in each bdd so min marginals match each other
            for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
                update_Lagrange_multiplier(var, bdd_index, exp_sums[bdd_index], average_exp_sums);
        }
        else
            forward_step(var, 0);
    }
}

void bdd_min_marginal_averaging_smoothed::averaging_pass_backward()
{
    double lb = 0.0;
    std::vector<bdd_branch_node_opt_smoothed::exp_sum_entry> exp_sums;

    for (long int var = nr_variables() - 1; var >= 0; --var)
    {
        if(nr_bdds(var) > 1)
        {
            exp_sums.clear();
            for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
                exp_sums.push_back(this->exp_sums(var, bdd_index));

            const double average_exp_sums = this->average_exp_sums(exp_sums.begin(), exp_sums.end());

            for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
            {
                update_Lagrange_multiplier(var, bdd_index, exp_sums[bdd_index], average_exp_sums);
                backward_step(var, bdd_index);
                lb += lower_bound_backward(var, bdd_index);
            }
        }
        else
        {
            backward_step(var, 0);
            lb += lower_bound_backward(var, 0);
        }
    }

    lower_bound_ = lb;
}

} // namespace LPMP