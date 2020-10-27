#pragma once

#include <vector>
#include <tuple>
#include <chrono>
#include "bdd_variable.h"
#include "bdd_min_marginal_averaging.h"
#include "two_dimensional_variable_array.hxx"

namespace LPMP
{

// do message passing on subsets of relevant variables only
class bdd_min_marginal_averaging_restricted : public bdd_min_marginal_averaging
{
public:
    template <typename ITERATOR>
    double variable_score(ITERATOR marginals_begin, ITERATOR marginals_end, const double th) const;

    template <typename VAR_ITERATOR, typename BDD_MASK>
    void min_marginal_averaging_forward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask);

    template <typename VAR_ITERATOR, typename BDD_MASK>
    double min_marginal_averaging_backward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask);

    void init(const ILP_input &input);

    std::vector<std::size_t> full_forward_iteration();

    std::array<double, 2> push_to_masked_backward_iteration(const std::vector<std::size_t> affected_variables, const two_dim_variable_array<std::size_t> affected_bdd_mask);
    void iteration();

private:
    template <typename ITERATOR>
    std::tuple<std::vector<std::size_t>, two_dim_variable_array<std::size_t>> compute_bdd_mask(ITERATOR variable_begin, ITERATOR variable_end) const;

    std::size_t variable(const bdd_variable_mma *bdd_var) const;
    std::size_t bdd_index(const bdd_variable_mma *bdd_var) const;

    struct bdd_variable_offset_inverse
    {
        std::size_t var;
        std::size_t bdd_index;
    };
    std::vector<bdd_variable_offset_inverse> bdd_variable_offset_inverse_;

    // how many inconsistent variables may exist at most so that restricted iterations are tried out
    constexpr static double restricted_iter_nr_inconsistent_vars_th = 0.02;
    constexpr static std::size_t max_nr_restricted_iter = 50;
    double prev_iteration_lb_restricted = -std::numeric_limits<double>::infinity();
    double initial_lb_restricted = -std::numeric_limits<double>::infinity();
    double initial_lb = -std::numeric_limits<double>::infinity();

    double full_iteration_lb_increase = -std::numeric_limits<double>::infinity();
    double full_iteration_lb_increase_per_sec = -std::numeric_limits<double>::infinity();
};

void bdd_min_marginal_averaging_restricted::init(const ILP_input &input)
{
    bdd_min_marginal_averaging::init(input);

    bdd_variable_offset_inverse_.reserve(this->bdd_variables_.data().size());
    for (std::size_t var = 0; var < this->bdd_variables_.size(); ++var)
        for (std::size_t bdd_index = 0; bdd_index < this->bdd_variables_[var].size(); ++bdd_index)
            bdd_variable_offset_inverse_.push_back({var, bdd_index});
}

template <typename ITERATOR>
double bdd_min_marginal_averaging_restricted::variable_score(ITERATOR marginals_begin, ITERATOR marginals_end, const double th) const
{
    double largest_positive = -std::numeric_limits<double>::infinity();
    double smallest_negative = std::numeric_limits<double>::infinity();
    std::size_t zero_min_marginal_diff;
    for (auto marginals_it = marginals_begin; marginals_it != marginals_end; ++marginals_it)
    {
        const double marginal_diff = (*marginals_it)[1] - (*marginals_it)[0];
        if (marginal_diff > th)
            largest_positive = std::max(largest_positive, marginal_diff);
        else if (marginal_diff < -th)
            smallest_negative = std::min(smallest_negative, marginal_diff);
        else
            zero_min_marginal_diff++;
    }

    if (std::isfinite(largest_positive) && std::isfinite(smallest_negative))
    {
        const double lower_bound_gain = std::min(largest_positive, -smallest_negative);
        return lower_bound_gain;
    }
    else if (zero_min_marginal_diff == std::distance(marginals_begin, marginals_end))
    {
        return std::numeric_limits<double>::infinity();
    }
    return 0.0;
}

template <typename VAR_ITERATOR, typename BDD_MASK>
void bdd_min_marginal_averaging_restricted::min_marginal_averaging_forward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask)
{
    assert(std::is_sorted(var_begin, var_end));
    std::vector<std::array<double, 2>> min_marginals;
    for (std::size_t i = 0; i < std::distance(var_begin, var_end); ++i)
    {
        const std::size_t var = *(var_begin + i);
        //std::cout << "variable = " << var << "; ";

        // collect min marginals
        min_marginals.clear();
        for (std::size_t j = 0; j < bdd_mask[i].size(); ++j)
        {
            const std::size_t bdd_index = bdd_mask(i, j);
            forward_step(var, bdd_index);
            min_marginals.push_back(min_marginal(var, bdd_index));
        }

        // set marginals in each bdd so min marginals match each other
        if (min_marginals.size() > 1)
        {
            const std::array<double, 2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());
            for (std::size_t j = 0; j < bdd_mask[i].size(); ++j)
            {
                const std::size_t bdd_index = bdd_mask(i, j);
                set_marginal(var, bdd_index, average_marginal, min_marginals[j]);
            }
        }
    }
}

template <typename VAR_ITERATOR, typename BDD_MASK>
double bdd_min_marginal_averaging_restricted::min_marginal_averaging_backward_restricted(VAR_ITERATOR var_begin, VAR_ITERATOR var_end, BDD_MASK bdd_mask)
{
    //assert(std::is_sorted(var_begin, var_end, std::greater_equal<std::size_t>()));
    assert(std::is_sorted(var_begin, var_end));
    std::vector<std::array<double, 2>> min_marginals;
    double lb = 0;

    for (std::ptrdiff_t i = std::distance(var_begin, var_end) - 1; i >= 0; --i)
    {
        const std::size_t var = *(var_begin + i);

        // collect min marginals
        min_marginals.clear();
        for (std::size_t j = 0; j < bdd_mask[i].size(); ++j)
        {
            const std::size_t bdd_index = bdd_mask(i, j);
            min_marginals.push_back(min_marginal(var, bdd_index));
        }

        const std::array<double, 2> average_marginal = min_marginals.size() > 0 ? average_marginals(min_marginals.begin(), min_marginals.end()) : std::array<double, 2>{0.0, 0.0};

        // set marginals in each bdd so min marginals match each other
        for (std::size_t j = 0; j < bdd_mask[i].size(); ++j)
        {
            const std::size_t bdd_index = bdd_mask(i, j);
            set_marginal(var, bdd_index, average_marginal, min_marginals[j]);
            backward_step(var, bdd_index);
            lb += lower_bound_backward(var, bdd_index);
        }
    }
    return lb;
}

std::size_t bdd_min_marginal_averaging_restricted::variable(const bdd_variable_mma *bdd_var) const
{
    const std::size_t bdd_var_offset = std::distance(&this->bdd_variables_.data()[0], bdd_var);
    assert(&this->bdd_variables_.data()[0] <= bdd_var);
    assert(bdd_var_offset < bdd_variable_offset_inverse_.size());
    return bdd_variable_offset_inverse_[bdd_var_offset].var;
}

std::size_t bdd_min_marginal_averaging_restricted::bdd_index(const bdd_variable_mma *bdd_var) const
{
    const std::size_t bdd_var_offset = std::distance(&this->bdd_variables_.data()[0], bdd_var);
    assert(&this->bdd_variables_.data()[0] <= bdd_var);
    assert(bdd_var_offset < bdd_variable_offset_inverse_.size());
    return bdd_variable_offset_inverse_[bdd_var_offset].bdd_index;
}

// given variables, mark all associated bdds.
// first return argument: All variables one needs to iterate over.
// second argument: All relevant bdd indices for respective variable.
template <typename ITERATOR>
std::tuple<std::vector<std::size_t>, two_dim_variable_array<std::size_t>> bdd_min_marginal_averaging_restricted::compute_bdd_mask(ITERATOR variable_begin, ITERATOR variable_end) const
{
    auto bdd_variable_func = [&](const bdd_variable_mma *bdd_var) { return this->variable(bdd_var); };
    auto bdd_index_func = [&](const bdd_variable_mma *bdd_var) { return this->bdd_index(bdd_var); };
    return LPMP::compute_bdd_mask(variable_begin, variable_end, *this, bdd_variable_func, bdd_index_func);
}

// return variables with inconsistent marginals
std::vector<std::size_t> bdd_min_marginal_averaging_restricted::full_forward_iteration()
{
    constexpr double th = 1e-4;
    std::vector<std::size_t> inconsistent_vars;

    std::vector<std::array<double, 2>> min_marginals;
    for (std::size_t var = 0; var < nr_variables(); ++var)
    {
        min_marginals.clear();
        for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
        {
            forward_step(var, bdd_index);
            min_marginals.push_back(this->min_marginal(var, bdd_index));
        }

        const double s = variable_score(min_marginals.begin(), min_marginals.end(), th);
        if (s >= th)
            inconsistent_vars.push_back(var);

        const std::array<double, 2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());

        for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
            set_marginal(var, bdd_index, average_marginal, min_marginals[bdd_index]);
    }

    return inconsistent_vars;
}

std::array<double, 2> bdd_min_marginal_averaging_restricted::push_to_masked_backward_iteration(const std::vector<std::size_t> affected_variables, const two_dim_variable_array<std::size_t> affected_bdd_mask)
{
    assert(std::is_sorted(affected_variables.begin(), affected_variables.end()));
    assert(affected_variables.size() == affected_bdd_mask.size());
    std::ptrdiff_t affected_variables_counter = affected_variables.size() - 1;
    std::vector<std::array<double, 2>> min_marginals;

    // backward pass
    double affected_bdd_lb = 0.0;
    double non_affected_bdd_lb = 0.0;
    for (std::ptrdiff_t var = nr_variables() - 1; var >= 0; --var)
    {
        min_marginals.clear();
        for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
            min_marginals.push_back(min_marginal(var, bdd_index));

        if (affected_variables_counter >= 0 && var == affected_variables[affected_variables_counter])
        {
            for (const auto x : affected_bdd_mask[affected_variables_counter])
                //std::cout << x << ",";
                //std::cout << "\n";
                assert(std::is_sorted(affected_bdd_mask[affected_variables_counter].begin(), affected_bdd_mask[affected_variables_counter].end()));
            assert(std::adjacent_find(affected_bdd_mask[affected_variables_counter].begin(), affected_bdd_mask[affected_variables_counter].end()) == affected_bdd_mask[affected_variables_counter].end());
            assert(affected_bdd_mask[affected_variables_counter].back() < nr_bdds(var));
            const std::array<double, 2> average_marginal = this->average_marginals(min_marginals.begin(), min_marginals.end(), affected_bdd_mask[affected_variables_counter].size());
            //const std::array<double, 2> average_marginal = this->average_marginals(min_marginals.begin(), min_marginals.end());
            // push marginals into affected bdds
            std::size_t bdd_idx_counter = 0;
            for (std::size_t bdd_index = 0; bdd_index < this->nr_bdds(var); ++bdd_index)
            {
                min_marginals.push_back(min_marginal(var, bdd_index));
                // affected bdd
                if (bdd_idx_counter < affected_bdd_mask[affected_variables_counter].size() && affected_bdd_mask(affected_variables_counter, bdd_idx_counter) == bdd_index)
                {
                    set_marginal(var, bdd_index, average_marginal, min_marginals[bdd_index]);
                    ++bdd_idx_counter;
                    backward_step(var, bdd_index);
                    affected_bdd_lb += lower_bound_backward(var, bdd_index);
                }
                else // unaffected bdd
                {
                    set_marginal(var, bdd_index, {0.0, 0.0}, min_marginals[bdd_index]);
                    //set_marginal(var, bdd_index, average_marginal, min_marginals[bdd_index]);
                    backward_step(var, bdd_index);
                    non_affected_bdd_lb += lower_bound_backward(var, bdd_index);
                }
            }
            assert(bdd_idx_counter == affected_bdd_mask[affected_variables_counter].size());
            --affected_variables_counter;
        }
        else
        {
            // perform normal min-marginal averaging
            const std::array<double, 2> average_marginal = average_marginals(min_marginals.begin(), min_marginals.end());
            for (std::size_t bdd_index = 0; bdd_index < nr_bdds(var); ++bdd_index)
            {
                set_marginal(var, bdd_index, average_marginal, min_marginals[bdd_index]);
                backward_step(var, bdd_index);
                non_affected_bdd_lb += lower_bound_backward(var, bdd_index);
            }
        }
    }
    return {affected_bdd_lb, non_affected_bdd_lb};
}

void bdd_min_marginal_averaging_restricted::iteration()
{
    auto increase_rate = [](const auto begin_time, const auto end_time, const double begin_lb, const double end_lb) {
        assert(end_lb >= begin_lb - 1e-6);
        const std::chrono::duration<double> seconds_passed = end_time - begin_time;
        assert(seconds_passed.count() >= 0.0);
        return (end_lb - begin_lb) / seconds_passed.count();
    };

    const double initial_lb = this->lower_bound_;
    const auto begin_time = std::chrono::system_clock::now();
    auto inconsistent_vars = full_forward_iteration();
    if (double(inconsistent_vars.size()) / double(nr_variables()) < restricted_iter_nr_inconsistent_vars_th)
    {
        auto [affected_variables, affected_bdd_mask] = compute_bdd_mask(inconsistent_vars.begin(), inconsistent_vars.end());
        const auto [affected_bdd_lb, non_affected_bdd_lb] = push_to_masked_backward_iteration(affected_variables, affected_bdd_mask);
        const auto end_full_backward_iter = std::chrono::system_clock::now();
        const double increase_rate_full_backward = increase_rate(begin_time, end_full_backward_iter, initial_lb, affected_bdd_lb + non_affected_bdd_lb);
        assert(std::abs(affected_bdd_lb + non_affected_bdd_lb - this->compute_lower_bound()) <= 1e-6);
        this->lower_bound_ = affected_bdd_lb + non_affected_bdd_lb;
        while (double(affected_bdd_mask.size()) / double(nr_variables()) < restricted_iter_nr_inconsistent_vars_th)
        {
            std::cout << "#affected variables = " << affected_variables.size() << "\n";
            const auto begin_restricted_iteration_time = std::chrono::system_clock::now();
            double prev_non_affected_bdd_lb = non_affected_bdd_lb;

            for (std::size_t i = 0; i < max_nr_restricted_iter; ++i)
            {
                const auto begin_restricted_time_current_iter = std::chrono::system_clock::now();
                min_marginal_averaging_forward_restricted(affected_variables.begin(), affected_variables.end(), affected_bdd_mask);
                const double lb_after_restricted_iter = min_marginal_averaging_backward_restricted(affected_variables.begin(), affected_variables.end(), affected_bdd_mask);
                const auto end_restricted_iteration_time = std::chrono::system_clock::now();
                assert(lb_after_restricted_iter >= non_affected_bdd_lb - 1e-6);
                const double increase_rate_restricted = increase_rate(begin_restricted_iteration_time, end_restricted_iteration_time, prev_non_affected_bdd_lb, lb_after_restricted_iter);
                if (increase_rate_restricted <= increase_rate_full_backward)
                {
                    std::cout << "#restricted iterations = " << i << "\n";
                    break;
                }
                prev_non_affected_bdd_lb = lb_after_restricted_iter;
            }
            std::tie(affected_variables, affected_bdd_mask) = compute_bdd_mask(affected_variables.begin(), affected_variables.end());
        }
    }
    else
    {
        bdd_min_marginal_averaging_restricted::min_marginal_averaging_backward();
        const double lb = this->lower_bound();
        const auto end_time = std::chrono::system_clock::now();
        const std::chrono::duration<double> seconds_passed = end_time - begin_time;
        full_iteration_lb_increase = lb - initial_lb;
        full_iteration_lb_increase_per_sec = full_iteration_lb_increase / seconds_passed.count();
        this->lower_bound_ = lb;
    }
}

} // namespace LPMP