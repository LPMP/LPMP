#ifndef LPMP_GRAPH_MATCHING_INPUT_H
#define LPMP_GRAPH_MATCHING_INPUT_H

#include <vector>
#include <string>
#include <iostream>
#include "matching_problem_input.h"

namespace LPMP {

// grammar for reading in files in the format of the Dual Decomposition algorithm of Torresani, Kolmogorov and Rother
namespace TorresaniEtAlInput {

   graph_matching_input parse_file(const std::string& filename);

   graph_matching_input parse_string(const std::string& input);

   template<typename SOLVER>
   bool parse_problem(const std::string& filename, SOLVER& s)
   {
       const auto gm_input = parse_file(filename);
       std::cout << "parsing " << filename << "\n";
       auto& gm_constructor = s.template GetProblemConstructor<0>();

       gm_constructor.read_input(gm_input);
       gm_constructor.construct();

       return true;
   }

   template<typename SOLVER>
   bool ParseProblemGM_trees(const std::string& filename, SOLVER& s)
   {
      /*
      auto input = parse_file(filename);
      constexpr PairwiseConstruction pc = FmcConstruction(typename SOLVER::FMC{});

      if(pc == PairwiseConstruction::Left) {
         auto& mrf_left = s.template GetProblemConstructor<0>();
         mrf_left.SetGraph(input.leftGraph_);
         construct_left_mrf(input, mrf_left, 1.0, 1.0);
         auto trees_left = mrf_left.compute_forest_cover();
         for(auto& tree : trees_left) {
            s.GetLP().add_tree(tree);
         }
      } else if(pc == PairwiseConstruction::Right) {
         auto& mrf_right = s.template GetProblemConstructor<1>();
         mrf_right.SetGraph(input.rightGraph_);
         construct_right_mrf(input, mrf_right, 1.0, 1.0);
         auto trees_right = mrf_right.compute_forest_cover();
         for(auto& tree : trees_right) {
            s.GetLP().add_tree(tree);
         }
      } else {
         assert(false);
      }

       */
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemMCF_trees(const std::string& filename, SOLVER& s)
   {
       /*
      using FMC = typename SOLVER::FMC;
      auto input = ParseFile(filename);
      construct_gm( s, input );
      construct_mp( s, input );
      factor_tree<FMC> mcf_tree;
      construct_mcf( s, input, &mcf_tree );
      s.GetLP().add_tree(mcf_tree);

      auto& mrf_left = s.template GetProblemConstructor<0>();
      auto trees_left = mrf_left.compute_forest_cover();
      for(auto& tree : trees_left) {
         s.GetLP().add_tree(tree);
      }

      auto& mrf_right = s.template GetProblemConstructor<1>();
      auto trees_right = mrf_right.compute_forest_cover();
      for(auto& tree : trees_right) {
         s.GetLP().add_tree(tree);
      }
      */
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemLocalSubproblems_trees(const std::string& filename, SOLVER& s)
   {
       /*
      auto input = ParseFile(filename);
      construct_gm( s, input );

      using FMC = typename SOLVER::FMC;
      factor_tree<FMC> mcf_tree;
      construct_mcf( s, input, &mcf_tree );
      s.GetLP().add_tree(mcf_tree);

      auto& mrf_left = s.template GetProblemConstructor<0>();
      auto trees_left = mrf_left.compute_forest_cover();
      for(auto& tree : trees_left) {
         s.GetLP().add_tree(tree);
      }

      auto& mrf_right = s.template GetProblemConstructor<1>();
      auto trees_right = mrf_right.compute_forest_cover();
      for(auto& tree : trees_right) {
         s.GetLP().add_tree(tree);
      }

      auto& local_subproblem_constructor_left = s.template GetProblemConstructor<2>();
      auto local_trees_left = local_subproblem_constructor_left.add_local_subproblems(mrf_left);
      for(auto& tree : local_trees_left) { s.GetLP().add_tree(tree); }

      auto& local_subproblem_constructor_right = s.template GetProblemConstructor<3>();
      auto local_trees_right = local_subproblem_constructor_right.add_local_subproblems(mrf_right);
      for(auto& tree : local_trees_right) { s.GetLP().add_tree(tree); }

*/
      return true; 
   }


} // namespace TorresaniEtAlInput

} // namespace LPMP

#endif // LPMP_GRAPH_MATCHING_INPUT_H

