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

       template<typename STREAM>
           void write(STREAM& s) const;
   };


   template<typename STREAM>
       void mrf_input::write(STREAM& s) const
       {
           auto unary_var = [](const std::size_t i, const std::size_t l) -> std::string {
               return "x_" + std::to_string(i) + "(" + std::to_string(l) + ")";
           };

           auto pairwise_var = [](const std::array<std::size_t,2> vars, const std::array<std::size_t,2> labels) {
               return "x_{" + std::to_string(vars[0]) + "," + std::to_string(vars[1]) + "}(" + std::to_string(labels[0]) + "," + std::to_string(labels[1]) + ")";
           };

           // objective
           s << "Minimize\n";
           for(std::size_t i=0; i<unaries.size(); ++i) {
               for(std::size_t l=0; l<unaries[i].size(); ++l) {
                   s << (unaries(i,l) >= 0 ? "+ " : "- ") << std::abs(unaries(i,l)) << " " << unary_var(i,l) << "\n";
               }
           }

           assert(pairwise_indices.size() == pairwise_values.size());
           for(std::size_t pairwise_idx=0; pairwise_idx<pairwise_indices.size(); ++pairwise_idx) {
               const auto [i,j] = pairwise_indices[pairwise_idx];
               assert(pairwise_values.dim2(pairwise_idx) == unaries[i].size());
               assert(pairwise_values.dim3(pairwise_idx) == unaries[j].size());
               for(std::size_t l_i=0; l_i<pairwise_values.dim2(pairwise_idx); ++l_i) {
                   for(std::size_t l_j=0; l_j<pairwise_values.dim3(pairwise_idx); ++l_j) {
                       const double val = pairwise_values(pairwise_idx, l_i, l_j); 
                       s << (val >= 0.0 ? "+ " : "- ") << std::abs(val) << " " << pairwise_var({i,j},{l_i,l_j}) << "\n";
                   }
               }
           }

           // constraints
           s << "Subject To\n";

           // simplex constraints
           for(std::size_t i=0; i<unaries.size(); ++i) {
               for(std::size_t l=0; l<unaries[i].size(); ++l) {
                   s << unary_var(i,l);
                   if(l+1 != unaries[i].size())
                       s << " + ";
               }
               s << " = 1\n";
           }

           for(std::size_t pairwise_idx=0; pairwise_idx<pairwise_indices.size(); ++pairwise_idx) {
               const auto [i,j] = pairwise_indices[pairwise_idx];
               assert(i < j);
               for(std::size_t l_i=0; l_i<pairwise_values.dim2(pairwise_idx); ++l_i) {
                   for(std::size_t l_j=0; l_j<pairwise_values.dim3(pairwise_idx); ++l_j) {
                       s << pairwise_var({i,j},{l_i,l_j});
                       if(l_i+1 < unaries[i].size() || l_j+1 < unaries[i].size())
                           s << " + ";
                   }
               }
               s << " = 1\n";
           }
           
           // marginalization constraints
           for(std::size_t pairwise_idx=0; pairwise_idx<pairwise_indices.size(); ++pairwise_idx) {
               const auto [i,j] = pairwise_indices[pairwise_idx];
               assert(i < j);
               for(std::size_t l_i=0; l_i<unaries[i].size(); ++l_i) {
                   s << unary_var(i,l_i);
                   for(std::size_t l_j=0; l_j<unaries[j].size(); ++l_j) {
                       s << " - " << pairwise_var({i,j},{l_i,l_j});
                   }
                   s << " = 0\n";
               }

               for(std::size_t l_j=0; l_j<unaries[j].size(); ++l_j) {
                   s << unary_var(j,l_j);
                   for(std::size_t l_i=0; l_i<unaries[i].size(); ++l_i) {
                       s << " - " << pairwise_var({i,j},{l_i,l_j});
                   }
                   s << " = 0\n";
               }
           }

           s << "End\n";
       }

} // namespace LPMP
