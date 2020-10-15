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
       const auto get_pairwise_potential(const std::size_t i) const { assert(i<no_pairwise_factors()); return pairwise_values[i]; }
       auto get_pairwise_potential(const std::size_t i) { assert(i<no_pairwise_factors()); return pairwise_values[i]; }

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
       std::string Potts_pairwise_variable_identifier(const std::size_t i, const std::size_t j) const;

       bool unary_variable_active(const std::size_t i, const std::size_t j) const;
       bool unary_variable_active(const std::size_t i) const;
       std::size_t forced_label(const std::size_t i) const;
       bool pairwise_variable_active(const std::size_t pairwise_index, const std::array<std::size_t, 2> labels) const;
       bool pairwise_variable_active(const std::size_t pairwise_index) const;

       void propagate();

       bool is_Potts(const std::size_t pairwise_idx) const;
       double Potts_strength(const std::size_t pairwise_index) const;
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

   inline bool mrf_input::unary_variable_active(const std::size_t i) const
   {
       std::size_t nr_active_labels = 0;
       for (std::size_t l = 0; l < cardinality(i); ++l)
       {
           if (unary_variable_active(i, l))
               ++nr_active_labels;
       }
       assert(nr_active_labels <= cardinality(i)); // or < ?
       assert(nr_active_labels > 0);
       return nr_active_labels > 1;
   }

   inline std::size_t mrf_input::forced_label(const std::size_t i) const
   {
       assert(!unary_variable_active(i));
       for (std::size_t l = 0; l < cardinality(i); ++l)
           if (unary_variable_active(i, l))
               return l;
       throw std::runtime_error("no label active");
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

   inline bool mrf_input::pairwise_variable_active(const std::size_t pairwise_index) const
   {
       assert(pairwise_index < no_pairwise_factors());
       const auto [i,j] = get_pairwise_variables(pairwise_index);
       return unary_variable_active(i) && unary_variable_active(j);
   }

   inline std::string mrf_input::Potts_pairwise_variable_identifier(const std::size_t i, const std::size_t j) const
   {
       assert(i < no_variables());
       assert(j < no_variables());
       return "y_" + std::to_string(i) + "_" + std::to_string(j);
   }

   inline void mrf_input::propagate()
   {
       // if exactly one unary variable is active, incorporate pairwise costs into unary ones
       for(std::size_t pairwise_idx = 0; pairwise_idx < no_pairwise_factors(); ++pairwise_idx)
       {
           const auto [i, j] = get_pairwise_variables(pairwise_idx);
           auto pot = get_pairwise_potential(pairwise_idx);
           if(unary_variable_active(i) && !unary_variable_active(j))
           {
               const std::size_t l_j = forced_label(j);
               for(std::size_t l_i = 0; l_i < cardinality(i); ++l_i)
               {
                   const double delta = pot(l_i, l_j);
                   unaries(l_i, l_i) += delta;
                   for (std::size_t ll_j = 0; ll_j < cardinality(j); ++ll_j)
                       pot(l_i, ll_j) -= delta;
               } 
           }
           if(!unary_variable_active(i) && unary_variable_active(j))
           {
               const std::size_t l_i = forced_label(i);
               for(std::size_t l_j = 0; l_j < cardinality(j); ++l_j)
               {
                   const double delta = pot(l_i, l_j);
                   unaries(l_i, l_i) += delta;
                   for (std::size_t ll_i = 0; ll_i < cardinality(i); ++ll_i)
                       pot(ll_i, l_j) -= delta;
               } 
           }
       }
   }

   inline bool mrf_input::is_Potts(const std::size_t pairwise_idx) const
   {
       assert(pairwise_idx < no_pairwise_factors());

       const auto [i, j] = get_pairwise_variables(pairwise_idx);

       for(std::size_t l=0; l<cardinality(i); ++l)
           if (!unary_variable_active(i, l))
               return false;

       for(std::size_t l=0; l<cardinality(j); ++l)
           if (!unary_variable_active(j, l))
               return false;

       if (cardinality(i) != cardinality(j))
           return false;
       auto pot = get_pairwise_potential(pairwise_idx);
       assert(cardinality(i) > 1);
       const double diagonal = pot(0, 0);
       if(diagonal != 0.0)
           return false;
       const double off_diagonal = pot(0, 1);
       if(off_diagonal < 0.0)
           return false;
       for (std::size_t l_i = 0; l_i < cardinality(i); ++l_i)
           for (std::size_t l_j = 0; l_j < cardinality(j); ++l_j)
               if (l_i != l_j && pot(l_i, l_j) != off_diagonal)
                   return false;
       for (std::size_t l = 0; l < cardinality(i); ++l)
           if (pot(l, l) != diagonal)
               return false;
       return true;
   }

   inline double mrf_input::Potts_strength(const std::size_t pairwise_idx) const
   {
       assert(is_Potts(pairwise_idx));
       auto pot = get_pairwise_potential(pairwise_idx);
       return pot(0,1);
   }

   template <typename STREAM>
   void mrf_input::write_to_lp_objective(STREAM &s) const
   {
       for (std::size_t i = 0; i < unaries.size(); ++i)
           if (unary_variable_active(i))
               for (std::size_t l = 0; l < unaries[i].size(); ++l)
                   if (unary_variable_active(i, l))
                       s << (unaries(i, l) >= 0 ? "+ " : "- ") << std::abs(unaries(i, l)) << " " << unary_variable_identifier(i, l) << "\n";

       assert(pairwise_indices.size() == pairwise_values.size());
       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           const auto [i, j] = pairwise_indices[pairwise_idx];
           if (!pairwise_variable_active(pairwise_idx))
               continue;
           if (is_Potts(pairwise_idx))
           {
               s << " + " << Potts_strength(pairwise_idx) << " " << Potts_pairwise_variable_identifier(i, j) << "\n";
           }
           else
           {
               assert(pairwise_values.dim2(pairwise_idx) == unaries[i].size());
               assert(pairwise_values.dim3(pairwise_idx) == unaries[j].size());
               for (std::size_t l_i = 0; l_i < pairwise_values.dim2(pairwise_idx); ++l_i)
               {
                   for (std::size_t l_j = 0; l_j < pairwise_values.dim3(pairwise_idx); ++l_j)
                   {
                       if (pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                       {
                           const double val = pairwise_values(pairwise_idx, l_i, l_j);
                           s << (val >= 0.0 ? "+ " : "- ") << std::abs(val) << " " << pairwise_variable_identifier({i, j}, {l_i, l_j}) << "\n";
                       }
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
           if (!unary_variable_active(i))
               continue;
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
           if (!pairwise_variable_active(pairwise_idx))
               continue;
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(i < j);
           if(is_Potts(pairwise_idx))
           {
               // no constraints since only single variable is here. Difference constraints are below
           }
           else
           {
               bool variable_printed = false;
               for (std::size_t l_i = 0; l_i < pairwise_values.dim2(pairwise_idx); ++l_i)
               {
                   for (std::size_t l_j = 0; l_j < pairwise_values.dim3(pairwise_idx); ++l_j)
                   {
                       if (pairwise_variable_active(pairwise_idx, {l_i, l_j}))
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
       }

       // marginalization constraints
       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           if(!pairwise_variable_active(pairwise_idx))
               continue;
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(i < j);
           if(is_Potts(pairwise_idx))
           {
               for (std::size_t l = 0; l < unaries[i].size(); ++l)
               {
                   s << Potts_pairwise_variable_identifier(i, j) << " - " << unary_variable_identifier(i, l) << " + " << unary_variable_identifier(j, l) << " >= 0\n";
                   s << Potts_pairwise_variable_identifier(i, j) << " - " << unary_variable_identifier(j, l) << " + " << unary_variable_identifier(i, l) << " >= 0\n";
               }
           }
           else
           {
               for (std::size_t l_i = 0; l_i < unaries[i].size(); ++l_i)
               {
                   if (unary_variable_active(i, l_i))
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
                   if (unary_variable_active(j, l_j))
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
   }

   template <typename STREAM>
   void mrf_input::write_to_lp_variables(STREAM &s) const
   {
       for (std::size_t i = 0; i < unaries.size(); ++i)
           if (unary_variable_active(i))
               for (std::size_t l = 0; l < unaries[i].size(); ++l)
                   if (unary_variable_active(i, l))
                       s << unary_variable_identifier(i, l) << "\n";

       assert(pairwise_indices.size() == pairwise_values.size());
       for (std::size_t pairwise_idx = 0; pairwise_idx < pairwise_indices.size(); ++pairwise_idx)
       {
           if(!pairwise_variable_active(pairwise_idx))
               continue;
           const auto [i, j] = pairwise_indices[pairwise_idx];
           assert(pairwise_values.dim2(pairwise_idx) == unaries[i].size());
           assert(pairwise_values.dim3(pairwise_idx) == unaries[j].size());
           if(is_Potts(pairwise_idx))
           {
               s << Potts_pairwise_variable_identifier(i, j) << "\n";
           }
           else
           {
               for (std::size_t l_i = 0; l_i < pairwise_values.dim2(pairwise_idx); ++l_i)
                   for (std::size_t l_j = 0; l_j < pairwise_values.dim3(pairwise_idx); ++l_j)
                       if (pairwise_variable_active(pairwise_idx, {l_i, l_j}))
                           s << pairwise_variable_identifier({i, j}, {l_i, l_j}) << "\n";
           }
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
