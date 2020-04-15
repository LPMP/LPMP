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

       template <typename STREAM>
       void write_to_lp_objective(STREAM &s) const;

       template <typename STREAM>
       void write_to_lp_constraints(STREAM &s) const;

       template <typename STREAM>
       void write_to_lp_variables(STREAM &s) const;

       template <typename STREAM>
       void write(STREAM &s) const;

       std::string unary_variable_identifier(const std::size_t i, const std::size_t l) const;
       std::string pairwise_variable_identifier(const std::array<std::size_t,2> vars, const std::array<std::size_t,2> labels) const;

       bool unary_variable_active(const std::size_t i, const std::size_t j) const;
       bool pairwise_variable_active(const std::size_t pairwise_index, const std::array<std::size_t, 2> labels) const;
   };

   inline std::string mrf_input::unary_variable_identifier(const std::size_t i, const std::size_t l) const
   {
       assert(i < no_variables());
       assert(l < cardinality(i));
       return "x_" + std::to_string(i) + "(" + std::to_string(l) + ")";
   };

   inline bool mrf_input::unary_variable_active(const std::size_t i, const std::size_t j) const
   {
       assert(i < no_variables());
       assert(j < cardinality(i));
       return unaries(i,j) < std::numeric_limits<double>::infinity();
   }

   inline std::string mrf_input::pairwise_variable_identifier(const std::array<std::size_t, 2> vars, const std::array<std::size_t, 2> labels) const
   {
       assert(vars[0] < no_variables());
       assert(vars[1] < no_variables());
       assert(labels[0] < cardinality(vars[0]));
       assert(labels[1] < cardinality(vars[1]));

       return "x_{" + std::to_string(vars[0]) + "," + std::to_string(vars[1]) + "}(" + std::to_string(labels[0]) + "," + std::to_string(labels[1]) + ")";
   };

   inline bool mrf_input::pairwise_variable_active(const std::size_t pairwise_index, const std::array<std::size_t, 2> labels) const
   {
       assert(pairwise_index < no_pairwise_factors());
       const auto [i,j] = get_pairwise_variables(pairwise_index);
       assert(labels[0] < cardinality(i));
       assert(labels[1] < cardinality(j));

       return unary_variable_active(i, labels[0]) && unary_variable_active(j, labels[1]) && get_pairwise_potential(i)(labels[0], labels[1]) < std::numeric_limits<double>::infinity();
   }

   template <typename STREAM>
   void mrf_input::write_to_lp_objective(STREAM &s) const
   {
       for (std::size_t i = 0; i < unaries.size(); ++i)
           for (std::size_t l = 0; l < unaries[i].size(); ++l)
               if (unary_variable_active(i, l))
                   s << (unaries(i, l) >= 0 ? "+ " : "- ") << std::abs(unaries(i, l)) << " " << unary_variable_identifier(i, l) << "\n";

       assert(pairwise_indices.size() == pairwise_values.size());
       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(pairwise_values.dim2(pairwise_idx) == unaries[i].size());
           assert(pairwise_values.dim3(pairwise_idx) == unaries[j].size());
           for (std::size_t l_i = 0; l_i < pairwise_values.dim2(pairwise_idx); ++l_i)
           {
               for (std::size_t l_j = 0; l_j < pairwise_values.dim3(pairwise_idx); ++l_j)
               {
                   if(pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                   {
                       const double val = pairwise_values(pairwise_idx, l_i, l_j);
                       s << (val >= 0.0 ? "+ " : "- ") << std::abs(val) << " " << pairwise_variable_identifier({i, j}, {l_i, l_j}) << "\n";
                   }
               }
           }
       }
   }

   template <typename STREAM>
   void mrf_input::write_to_lp_constraints(STREAM &s) const
   {
       // simplex constraints
       for (std::size_t i = 0; i < unaries.size(); ++i)
       {
           bool variable_printed = false;
           for (std::size_t l = 0; l < unaries[i].size(); ++l)
           {
               if (unary_variable_active(i, l))
               {
                   if (variable_printed)
                       s << " + ";
                   else
                       variable_printed = true;
                   s << unary_variable_identifier(i, l);
               }
               }
           s << " = 1\n";
       }

       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(i < j);
           bool variable_printed = false;
           for (std::size_t l_i = 0; l_i < pairwise_values.dim2(pairwise_idx); ++l_i)
           {
               for (std::size_t l_j = 0; l_j < pairwise_values.dim3(pairwise_idx); ++l_j)
               {
                   if(pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                   {
                       if (variable_printed)
                           s << " + ";
                       else
                           variable_printed = true;
                       s << pairwise_variable_identifier({i, j}, {l_i, l_j});
                   }
               }
           }
           s << " = 1\n";
       }

       // marginalization constraints
       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(i < j);
           for (std::size_t l_i = 0; l_i < unaries[i].size(); ++l_i)
           {
               if(unary_variable_active(i, l_i))
               {
                   s << unary_variable_identifier(i, l_i);
                   for (std::size_t l_j = 0; l_j < unaries[j].size(); ++l_j)
                       if (pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                           s << " - " << pairwise_variable_identifier({i, j}, {l_i, l_j});
                   s << " = 0\n";
               }
           }

           for (std::size_t l_j = 0; l_j < unaries[j].size(); ++l_j)
           {
               if(unary_variable_active(j, l_j))
               {
               s << unary_variable_identifier(j, l_j);
               for (std::size_t l_i = 0; l_i < unaries[i].size(); ++l_i)
                   if (pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                       s << " - " << pairwise_variable_identifier({i, j}, {l_i, l_j});
               s << " = 0\n";
               }
           }
       }
   }

   template <typename STREAM>
   void mrf_input::write_to_lp_variables(STREAM &s) const
   {
       for (std::size_t i = 0; i < unaries.size(); ++i)
           for (std::size_t l = 0; l < unaries[i].size(); ++l)
               if (unary_variable_active(i, l))
                   s << unary_variable_identifier(i, l) << "\n";

       assert(pairwise_indices.size() == pairwise_values.size());
       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(pairwise_values.dim2(pairwise_idx) == unaries[i].size());
           assert(pairwise_values.dim3(pairwise_idx) == unaries[j].size());
           for (std::size_t l_i = 0; l_i < pairwise_values.dim2(pairwise_idx); ++l_i)
               for (std::size_t l_j = 0; l_j < pairwise_values.dim3(pairwise_idx); ++l_j)
                   if (pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                       s << pairwise_variable_identifier({i, j}, {l_i, l_j}) << "\n";
       } 
   }


   template <typename STREAM>
   void mrf_input::write(STREAM &s) const
   {
       s << "Minimize\n";
       write_to_lp_objective(s);

       s << "Subject To\n";
       write_to_lp_constraints(s);

       s << "Bounds\nBinaries\n";
       write_to_lp_variables(s);

       s << "End\n";
   }

} // namespace LPMP
