#pragma once

#include "mrf_input.h"
#include "binary_MRF_instance.hxx"

#include <string>
#include <cassert>

namespace LPMP {

// file format described in http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php
namespace binary_MRF_uai_input {

   binary_MRF_instance parse_file(const std::string& filename);
   binary_MRF_instance parse_string(const std::string& input);

}

namespace mrf_uai_input {

   mrf_input parse_file(const std::string& filename);
   mrf_input parse_string(const std::string& input);

/*


   template<typename MRF_CONSTRUCTOR>
      void build_mrf(MRF_CONSTRUCTOR& mrf, const MrfInput& input)
      {
         assert(input.number_of_cliques_ == input.clique_scopes_.size());
         assert(input.number_of_cliques_ == input.function_tables_.size());

         // only unary and pairwise potentials supported right now
         for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            assert(input.clique_scopes_[i].size() < 3);
         }

         // first input the unaries, as pairwise potentials need them to be able to link to them
         // add unary factors with cost zero for each variables. There are models where unaries are not explicitly added.
         for(INDEX i=0; i<input.number_of_variables_; ++i) {
            const INDEX noLabels = input.cardinality_[i];
            mrf.add_unary_factor(i,std::vector<REAL>(noLabels,0.0));
         }

         REAL initial_lb = 0.0;
         for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            if(input.clique_scopes_[i].size() == 1) {
               const INDEX var = input.clique_scopes_[i][0];
               //std::cout << "unary potential for variable " << var << ":\n";
               auto* f = mrf.get_unary_factor(var);
               assert(input.function_tables_[i].size() == input.cardinality_[var]);
               initial_lb += *std::min_element(input.function_tables_[i].begin(), input.function_tables_[i].end());
               for(INDEX x=0; x<input.function_tables_[i].size(); ++x) {
                  //std::cout << input.function_tables_[i][x] << " ";
                  assert( (*f->get_factor())[x] == 0.0);
                  (*f->get_factor())[x] = input.function_tables_[i][x];
               }
               //std::cout << "\n";
            }
         }

         //std::cout << "initial lower bound unaries = " << initial_lb << "\n"; 

         // now the pairwise potentials. 
         for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            if(input.clique_scopes_[i].size() == 2) {
               const INDEX var1 = input.clique_scopes_[i][0];
               const INDEX var2 = input.clique_scopes_[i][1];
               const INDEX dim1 = mrf.get_number_of_labels(var1);
               const INDEX dim2 = mrf.get_number_of_labels(var2);
               assert(var1<var2 && var2 < input.number_of_variables_);
               assert(input.function_tables_[i].size() == input.cardinality_[var1]*input.cardinality_[var2]);
               assert(input.function_tables_[i].size() == dim1*dim2);
               matrix<REAL> pairwise_cost(dim1,dim2);
               initial_lb += *std::min_element(input.function_tables_[i].begin(), input.function_tables_[i].end());
               //std::cout << "pairwise potential on (" << var1 << "," << var2 << "):\n";
               for(INDEX l1=0; l1<dim1; ++l1) {
                  for(INDEX l2=0; l2<dim2; ++l2) {
                     pairwise_cost(l1,l2) = input.function_tables_[i][l1*dim2 + l2];
               //      std::cout << input.function_tables_[i][l1*dim2 + l2] << " ";
                  }
               //   std::cout << "\n";
               }
               //std::cout << pairwise_cost;
               mrf.add_pairwise_factor(var1,var2,pairwise_cost); // or do we have to transpose the values?
            }
         }
      }


   template<typename SOLVER, INDEX PROBLEM_CONSTRUCTOR_NO>
   bool ParseString(const std::string& instance, SOLVER& s)
   {
      std::cout << "parsing string\n";
      MrfInput input;
      //bool read_suc = pegtl::parse_string<grammar, action>(instance,"",input);
      bool read_suc = pegtl::parse<grammar, action>(instance,"",input);
      if(read_suc) {
         auto& mrf_constructor = s.template GetProblemConstructor<PROBLEM_CONSTRUCTOR_NO>();
         build_mrf(mrf_constructor, input);
      }
      return read_suc;
   }

   template<typename SOLVER, INDEX PROBLEM_CONSTRUCTOR_NO>
   bool ParseStringDD(const std::string& instance, SOLVER& s)
   {
      const bool success = ParseString<SOLVER,PROBLEM_CONSTRUCTOR_NO>(instance, s);
      auto& mrf = s.template GetProblemConstructor<0>();
      // decompose mrf into trees automatically.
      auto trees = mrf.compute_forest_cover();
      for(auto& tree : trees) {
          s.GetLP().add_tree(tree);
      }

      return success;
   }

   template<typename SOLVER>
   bool ParseProblem(const std::string& filename, SOLVER& s)
   {
      std::cout << "parsing " << filename << "\n";
      pegtl::file_parser problem(filename);
      MrfInput input;
      bool read_suc = problem.parse< grammar, action >(input);
      if(read_suc) {
         auto& mrf_constructor = s.template GetProblemConstructor<0>();
         build_mrf(mrf_constructor, input);
      }
      return read_suc;
   }

   template<typename SOLVER>
   bool ParseProblemDD(const std::string& filename, SOLVER& s)
   {
     const bool success = ParseProblem(filename, s);
     auto& mrf = s.template GetProblemConstructor<0>();
     // decompose mrf into trees automatically.
     auto trees = mrf.compute_forest_cover();
     for(auto& tree : trees) {
       s.GetLP().add_tree(tree);
     }

     return success;
   }
   */

} // namespace mrf_uai_input

} // namespace LPMP
