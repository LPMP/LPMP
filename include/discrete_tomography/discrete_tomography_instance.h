#pragma once

#include "mrf/mrf_input.h"
#include <vector>
#include <algorithm>

namespace LPMP {

class discrete_tomography_instance {
public:
    mrf_input mrf;

    std::vector<std::vector<std::size_t>> projection_variables;
    std::vector<std::vector<double>> projection_costs;

    void propagate_projection_costs();

    template<typename STREAM>
    void write_to_lp(STREAM& s) const;

    private:
    bool is_sum_constraint(const std::size_t i) const;
    std::size_t sum_constraint_value(const std::size_t i) const;
};

inline void discrete_tomography_instance::propagate_projection_costs()
{
    // propagate all constraints that constrain variables to be zero
    for(std::size_t p=0; p<projection_costs.size(); ++p)
        if (is_sum_constraint(p) && sum_constraint_value(p) == 0)
            for (const std::size_t i : projection_variables[p])
                for (std::size_t l = 0; l < mrf.cardinality(i); ++l)
                    mrf.unaries(i, l) = std::numeric_limits<double>::infinity();
}

inline bool discrete_tomography_instance::is_sum_constraint(const std::size_t i) const
{
    assert(i < projection_costs.size());
    if(std::count(projection_costs[i].begin(), projection_costs[i].end(), std::numeric_limits<double>::infinity()) + 1 == projection_costs[i].size())
        return true;
    return false;
}

inline std::size_t discrete_tomography_instance::sum_constraint_value(const std::size_t i) const
{
    assert(is_sum_constraint(i));
    const auto pos = std::find_if(projection_costs[i].begin(), projection_costs[i].end(), [](const double c) { return c != std::numeric_limits<double>::infinity(); });
    assert(*pos == 0.0);
    return std::distance(projection_costs[i].begin(), pos);
}

template <typename STREAM>
void discrete_tomography_instance::write_to_lp(STREAM &s) const
{
    s << "Minimize\n";
    mrf.write_to_lp_objective(s);

    s << "Subject To\n";
    mrf.write_to_lp_constraints(s);

    // projection constraints
    assert(projection_costs.size() == projection_variables.size());
    for(std::size_t p=0; p<projection_variables.size(); ++p)
    {
        if(is_sum_constraint(p))
        {
            bool variable_printed = false;
            for (const std::size_t i : projection_variables[p])
            {
                for (std::size_t l = 1; l < mrf.cardinality(i); ++l)
                {
                    if(mrf.unary_variable_active(i, l))
                    {
                        if(variable_printed)
                            s << " + ";
                        else
                            variable_printed = true;
                        s << l << " " << mrf.unary_variable_identifier(i, l);
                    }
                }
            }
            s << " = " << sum_constraint_value(p) << "\n";
        }
        else
        {
            throw std::runtime_error("export of general projection constraints not implemented.");
        }
    }

    s << "Bounds\nBinaries\n";
    mrf.write_to_lp_variables(s);

    s << "End\n";

}

}