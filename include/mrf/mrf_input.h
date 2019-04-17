#pragma once

#include <vector>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"

namespace LPMP {

   struct mrf_input {
       two_dim_variable_array<double> unaries;
       three_dimensional_variable_array<double> pairwise_values;
       std::vector<std::array<std::size_t,2>> pairwise_indices;
       // higher order potentials not supported currently

       std::size_t no_variables() const { return unaries.size(); }
       std::size_t cardinality(const std::size_t i) const { assert(i<no_variables()); return unaries[i].size(); }
       std::size_t no_pairwise_factors() const { return pairwise_indices.size(); }
       auto get_unary(const std::size_t i) const { assert(i<no_variables()); return unaries[i]; }
       std::array<std::size_t,2> get_pairwise_variables(const std::size_t i) const { assert(i<no_pairwise_factors()); return pairwise_indices[i]; }
       auto get_pairwise_potential(const std::size_t i) const { assert(i<no_pairwise_factors()); return pairwise_values[i]; }
   };

} // namespace LPMP
