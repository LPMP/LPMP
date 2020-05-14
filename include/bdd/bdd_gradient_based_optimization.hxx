#pragma once

#include "bdd_anisotropic_diffusion.h"
#include "bdd_variable.h"
#include "bdd_min_marginal_averaging_smoothed.h"
#include <numeric>
#include <vector>

namespace LPMP
{

class bdd_gradient_opt : public bdd_base_consecutive<bdd_branch_node_opt_smoothed, bdd_variable_mma>
{
public:
    std::vector<double> gradient();
    void apply_and_project(const std::vector<double>& gradient);
    void smooth_forward_step(const std::size_t bdd_index, const std::size_t j);
    std::vector<double> projection_free_gradient(const std::vector<double> &a) const;
    std::vector<double> transpose(const std::vector<double> &b) const;
    template <typename VEC>
    void apply_delta(const VEC &delta_begin);
    std::size_t projection_free_gradient_size() const { return gradient_size() - this->nr_variables(); }
    double compute_lower_bound();

private:
    std::size_t gradient_size() const;
    double gradient(const std::size_t bdd_index, const std::size_t j) const;
    bdd_branch_node_exp_sum_entry exp_sums(const std::size_t bdd_index, const std::size_t j) const;
};

std::size_t bdd_gradient_opt::gradient_size() const
{
    std::size_t s = 0;
    for(std::size_t bdd_idx= 0; bdd_idx<this->nr_bdds(); ++bdd_idx)
        s += this->nr_variables(bdd_idx);
    return s;
}

std::vector<double> bdd_gradient_opt::gradient()
{
    std::vector<double> g;
    g.reserve(gradient_size());

    for (std::ptrdiff_t i = this->bdd_branch_nodes_.size()-1; i>=0; --i)
        bdd_branch_nodes_[i].smooth_backward_step();

    double f = 0.0;
    for(std::size_t bdd_idx= 0; bdd_idx<this->nr_bdds(); ++bdd_idx)
    {
        const auto &bdd = this->bdd_branch_nodes_[bdd_variable_delimiters(bdd_idx,0).bdd_branch_node_offset];
        assert(std::isfinite(std::log(bdd.m)));
        const double cost_scaling_ = 1.0;
        f += -cost_scaling_ * (std::log(bdd.m) + bdd.current_max); 
    }

    for(std::size_t bdd_idx= 0; bdd_idx<this->nr_bdds(); ++bdd_idx)
    {
        for (std::size_t j = 0; j < nr_variables(bdd_idx); ++j)
        {
            smooth_forward_step(bdd_idx, j);
            g.push_back(gradient(bdd_idx, j));
        }
    }

    return g;
}

// map gradient onto constraint free space.
std::vector<double> bdd_gradient_opt::projection_free_gradient(const std::vector<double>& a) const
{
    assert(a.size() == gradient_size());
    std::vector<std::size_t> bdd_index_counter;
    bdd_index_counter.reserve(this->nr_bdds());
    for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(); ++bdd_index)
        bdd_index_counter.push_back(this->nr_variables(bdd_index));
    std::vector<std::size_t> bdd_index_offset;
    bdd_index_offset.reserve(this->nr_bdds());
    bdd_index_offset.push_back(0);
    std::partial_sum(bdd_index_counter.begin(), bdd_index_counter.end() - 1, std::back_inserter(bdd_index_offset));
    std::fill(bdd_index_counter.begin(), bdd_index_counter.end(), 0);

    std::vector<double> b;
    b.reserve(gradient_size() - this->nr_variables());
    assert(this->nr_bdds() == this->bdd_variable_delimiters.size());

    for (std::size_t i = 0; i < this->nr_variables(); ++i)
    {
        const std::size_t last_bdd_index = this->Lagrange_multipliers(i, this->Lagrange_multipliers[i].size()-1).bdd_nr;
        const double last_gradient = a[bdd_index_offset[last_bdd_index] + bdd_index_counter[last_bdd_index]++];
        for(std::size_t j=0; j<this->Lagrange_multipliers[i].size()-1; ++j)
        {
            const std::size_t bdd_index = this->Lagrange_multipliers(i, j).bdd_nr;
            const double gradient = a[bdd_index_offset[bdd_index] + bdd_index_counter[bdd_index]++];
            b.push_back(gradient - last_gradient);
        }
    }

    assert(b.size() == gradient_size() - this->nr_variables());
    return b;

/*
    std::vector<std::size_t> nr_bdd_nodes_per_variable;
    nr_bdd_nodes_per_variable.reserve(this->nr_variables());
    for(const auto& x : this->Lagrange_multipliers)
        nr_bdd_nodes_per_variable.push_back(x.size()-1);

    std::vector<std::size_t> bdd_offset_per_variable;
    bdd_offset_per_variable.reserve(bdd_storage_.nr_variables());
    bdd_offset_per_variable.push_back(0);
    std::partial_sum(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end() - 1, std::back_inserter(bdd_offset_per_variable));
    std::fill(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end(), 0);

    std::size_t idx = 0;
    for(std::size_t i=0; i<this->nr_bdds(); ++i)
    {
        assert(this->nr_variables(i) == this->bdd_variable_delimiters[i].size());
        for(std::size_t j=0; j<this->nr_variables(i); ++j)
        {
            const std::size_t var = this->bdd_variable_delimiters(i,j).variable;
            if(nr_bdd_nodes_per_variable[var] + 1 < this->Lagrange_multipliers[var].size())
            {
                b[bdd_offset_per_variable[var] + nr_bdd_node_per_variable[var]++] += a[idx++];
            }
            else
            {
                for(std::size_t k=0; k+1<nr_bdd_nodes_per_variable[var]; ++k)
                {
                    b[bdd_offset_per_variable[var] + k] -= a[idx];
                }
                ++idx;
            }
        } 
    }
    assert(a.size() == idx);

    return b;
    */
}

// apply change of Lagrange multipliers to costs
template<typename VEC>
void bdd_gradient_opt::apply_delta(const VEC& delta)
{
    assert(delta.size() == gradient_size() - this->nr_variables());
    for(std::size_t i=0; i<delta.size(); ++i)
    {
        assert(std::isfinite(delta[i]));
    }
    std::vector<std::size_t> bdd_index_counter;
    bdd_index_counter.reserve(this->nr_bdds());
    for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(); ++bdd_index)
        bdd_index_counter.push_back(this->nr_variables(bdd_index));
    std::vector<std::size_t> bdd_index_offset;
    bdd_index_offset.reserve(this->nr_bdds());
    bdd_index_offset.push_back(0);
    std::partial_sum(bdd_index_counter.begin(), bdd_index_counter.end() - 1, std::back_inserter(bdd_index_offset));
    std::fill(bdd_index_counter.begin(), bdd_index_counter.end(), 0);

    std::vector<double> a(gradient_size(), 0);
    std::size_t idx = 0;
    for (std::size_t i = 0; i < this->nr_variables(); ++i)
    {
        const std::size_t last_bdd_index = this->Lagrange_multipliers[i].size() - 1;
        for (std::size_t j = 0; j < last_bdd_index; ++j)
        {
            const std::size_t bdd_index = this->Lagrange_multipliers(i, j).bdd_nr;
            this->Lagrange_multipliers(i,j).bdd_variable.cost += delta[idx];
            this->Lagrange_multipliers(i, last_bdd_index).bdd_variable.cost -= delta[idx++];
        }
    }
    assert(idx == delta.size());
}

// transpose of above
std::vector<double> bdd_gradient_opt::transpose(const std::vector<double>& b) const
{
    assert(b.size() == gradient_size() - this->nr_variables());

    std::vector<std::size_t> bdd_index_counter;
    bdd_index_counter.reserve(this->nr_bdds());
    for(std::size_t bdd_index=0; bdd_index<this->nr_bdds(); ++bdd_index)
        bdd_index_counter.push_back(this->nr_variables(bdd_index));
    std::vector<std::size_t> bdd_index_offset;
    bdd_index_offset.reserve(this->nr_bdds());
    bdd_index_offset.push_back(0);
    std::partial_sum(bdd_index_counter.begin(), bdd_index_counter.end() - 1, std::back_inserter(bdd_index_offset));
    std::fill(bdd_index_counter.begin(), bdd_index_counter.end(), 0);

    std::vector<double> a(gradient_size(), 0);
    std::size_t idx = 0;
    for (std::size_t i = 0; i < this->nr_variables(); ++i)
    {
        const std::size_t last_bdd_index = this->Lagrange_multipliers(i, this->Lagrange_multipliers[i].size()-1).bdd_nr;
        for(std::size_t j=0; j<this->Lagrange_multipliers[i].size()-1; ++j)
        {
            const std::size_t bdd_index = this->Lagrange_multipliers(i, j).bdd_nr;
            a[bdd_index_offset[bdd_index] + bdd_index_counter[bdd_index]++] = b[idx];
            a[bdd_index_offset[last_bdd_index] + bdd_index_counter[last_bdd_index]] -= b[idx++];
        }
        bdd_index_counter[last_bdd_index]++;
    }

    return a;

/*
    std::vector<std::size_t> nr_bdd_nodes_per_variable;
    nr_bdd_nodes_per_variable.reserve(this->nr_variables());
    for(const auto& x : this->Lagrange_multipliers)
        nr_bdd_nodes_per_variable(x.size()-1);

    std::vector<std::size_t> bdd_offset_per_variable;
    bdd_offset_per_variable.reserve(bdd_storage_.nr_variables());
    bdd_offset_per_variable.push_back(0);
    std::partial_sum(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end() - 1, std::back_inserter(bdd_offset_per_variable));
    std::fill(nr_bdd_nodes_per_variable.begin(), nr_bdd_nodes_per_variable.end(), 0);

    std::vector<double> a;
    a.reserve(gradient_size());
    for(std::size_t i=0; i<this->nr_bdds(); ++i)
    {
        assert(this->nr_variables(i) == this->bdd_variable_delimiters_[i].size()));
        for(std::size_t j=0; j<this->nr_variables(i); ++j)
        {
            const std::size_t var = bdd_variable_delimiters_(i,j).variable;
            if(nr_bdd_nodes_per_variable[var] + 1 < this->Lagrange_multipliers_[var].size())
            {
                a.push_back(b[bdd_offset_per_variable[var] + nr_bdd_nodes_per_variable[var]++]);
            }
            else
            {
                const double x = std::accumulate(b.begin() + bdd_offset_per_variable[var], b.begin() + bdd_offset_per_variable[var] + this->Lagrange_multipliers_[var].size()-1), 0.0);
                a.push_back(-x);
            } 
        }
    }
    assert(a.size() == gradient_size());
    return a;
    */
}

void bdd_gradient_opt::smooth_forward_step(const std::size_t bdd_index, const std::size_t j)
{
    assert(bdd_index < this->nr_bdds());
    assert(j < this->nr_variables(bdd_index));
    const auto [first_bdd_node, last_bdd_node] = this->bdd_branch_node_range(bdd_index, j);
    for(std::size_t bdd_branch_node_idx = first_bdd_node; bdd_branch_node_idx < last_bdd_node; ++bdd_branch_node_idx)
        this->bdd_branch_nodes_[bdd_branch_node_idx].smooth_forward_step();
}

bdd_branch_node_exp_sum_entry bdd_gradient_opt::exp_sums(const std::size_t bdd_index, const std::size_t j) const
{
    assert(bdd_index < this->nr_bdds());
    assert(j < this->nr_variables(bdd_index));
    const auto [first_bdd_node, last_bdd_node] = this->bdd_branch_node_range(bdd_index, j);
    auto first_bdd_branch_node_it = this->bdd_branch_nodes_.begin() + first_bdd_node;
    auto last_bdd_branch_node_it = this->bdd_branch_nodes_.begin() + last_bdd_node;
    return bdd_branch_node_opt_smoothed::exp_sums(first_bdd_branch_node_it, last_bdd_branch_node_it);
}

double bdd_gradient_opt::gradient(const std::size_t bdd_index, const std::size_t j) const
{
    const auto es = exp_sums(bdd_index, j);
    //const double g = es.sum[1]*exp(es.max[1]) / (es.sum[1]*exp(es.max[1]) + es.sum[0]*exp(es.max[0]));
    // computes above without overflow
    const double g = [&]() {
        if(es.max[1] > es.max[0])
            return es.sum[1] / (es.sum[1] + es.sum[0] * std::exp(es.max[0] - es.max[1]));
            else
            return es.sum[1]*std::exp(es.max[1] - es.max[0]) / (es.sum[1]*std::exp(es.max[1] - es.max[0]) + es.sum[0]);
    }();
    assert(std::isfinite(g));
    return g;
}

void bdd_gradient_opt::apply_and_project(const std::vector<double> &gradient)
{
    assert(gradient.size() == gradient_size());
    std::vector<double> orig_costs(this->Lagrange_multipliers.size(), 0.0);
    std::vector<double> updated_costs(orig_costs.size(), 0.0);

    std::size_t c = 0;
    for(std::size_t bdd_idx= 0; bdd_idx<this->nr_bdds(); ++bdd_idx)
    {
        for (std::size_t j = 0; j < nr_variables(bdd_idx); ++j)
        {
            const std::size_t var = this->bdd_variable_delimiters(bdd_idx, j).variable;
            const std::size_t offset = this->bdd_variable_delimiters(bdd_idx, j).Lagrange_multipliers_index;
            orig_costs[var] += this->Lagrange_multipliers(var,offset).bdd_variable.cost;
            this->Lagrange_multipliers(var, offset).bdd_variable.cost += +10.0 * gradient[c++];
            updated_costs[var] += this->Lagrange_multipliers(var, offset).bdd_variable.cost;
        }
    }
    assert(c == gradient.size());

    for(std::size_t bdd_idx= 0; bdd_idx<this->nr_bdds(); ++bdd_idx)
    {
        for (std::size_t j = 0; j < nr_variables(bdd_idx); ++j)
        {
            const std::size_t var = this->bdd_variable_delimiters(bdd_idx, j).variable;
            const std::size_t offset = this->bdd_variable_delimiters(bdd_idx, j).Lagrange_multipliers_index;
            this->Lagrange_multipliers(var, offset).bdd_variable.cost -= (updated_costs[var] - orig_costs[var])/this->Lagrange_multipliers[var].size();
        }
    }
}

double bdd_gradient_opt::compute_lower_bound()
{
    for (std::ptrdiff_t i = this->bdd_branch_nodes_.size()-1; i>=0; --i)
        bdd_branch_nodes_[i].smooth_backward_step();

    double lb = 0.0;
    for(std::size_t bdd_idx= 0; bdd_idx<this->nr_bdds(); ++bdd_idx)
    {
        const auto &bdd = this->bdd_branch_nodes_[bdd_variable_delimiters(bdd_idx,0).bdd_branch_node_offset];
        assert(std::isfinite(std::log(bdd.m)));
        const double cost_scaling_ = 1.0;
        lb += -cost_scaling_ * (std::log(bdd.m) + bdd.current_max); 
    }

    return lb; 
}

}