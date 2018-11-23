#ifndef LPMP_MULTIGRAPH_MATCHING_SYNCHRONIZATION_H
#define LPMP_MULTIGRAPH_MATCHING_SYNCHRONIZATION_H

#include "matching_problem_input.h"
#include <eigen3/Eigen/Eigenvalues>
#include <vector>
#include "vector.hxx"
#include <iostream> // TODO: remove

namespace LPMP {

   void synchronize_multigraph_matching(const multigraph_matching_input::graph_size& mgm_size, multigraph_matching_input::labeling& labeling, const double roundin_th = 0.1);
   void synchronize_multigraph_matching(const multigraph_matching_input::graph_size& mgm_size, multigraph_matching_input::labeling& labeling, const matrix<int>& allowed_matchings, const double roundin_th = 0.1);

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_SYNCHRONIZATION_H
