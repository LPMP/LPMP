#include "multigraph_matching/multigraph_matching_synchronization.h"
#include <Eigen/Eigenvalues>
#include <vector>
#include "union_find.hxx"

namespace LPMP {

   Eigen::MatrixXd compute_eigenvector_matrix(const multigraph_matching_input::graph_size& mgm_size, const multigraph_matching_input::labeling& labeling)
   {
      // first build matching matrix
      Eigen::MatrixXd W(mgm_size.total_no_nodes(), mgm_size.total_no_nodes());
      W.setZero();
      for(std::size_t i=0; i<W.cols(); ++i)
         W(i,i) = 1.0;

      for(const auto& gm : labeling) {
         for(std::size_t i=0; i<gm.labeling.size(); ++i) {
            const std::size_t left_index = mgm_size.node_no(gm.left_graph_no, i);
            if(gm.labeling[i] != std::numeric_limits<std::size_t>::max()) {
               const std::size_t right_index = mgm_size.node_no(gm.right_graph_no, gm.labeling[i]);
               W(left_index, right_index) = 1;
               W(right_index, left_index) = 1;
            }
         }
      }

      //std::cout << W<< "\n";

      // number of eigenvectors to take into account
      const std::size_t k = [&]() {
         std::size_t max_no_nodes = 0;
         for(std::size_t p=0; p<mgm_size.no_graphs(); ++p)
            max_no_nodes = std::max(max_no_nodes, mgm_size.no_nodes(p));
         return max_no_nodes;
      }();

      Eigen::SelfAdjointEigenSolver<decltype(W)> U_tmp(W*W.transpose());
      Eigen::MatrixXd U = U_tmp.eigenvectors();
      auto largest_k_eigenvectors = U.block(0,U.rows()-k,U.cols(),k);
      //std::cout << "eigenvectors:\n" << U << "\n";
      //std::cout << "k largest eigenvectors:\n" << largest_k_eigenvectors << "\n";

      // normalize U such that each column has l2 norm of 1
      Eigen::VectorXd n = Eigen::VectorXd::Constant(largest_k_eigenvectors.rows(),1.0).array() / largest_k_eigenvectors.rowwise().norm().array();
      largest_k_eigenvectors = n.asDiagonal() * largest_k_eigenvectors;
      // test if each of k largest eigenvectors of U has l2-norm 1
      //for(std::size_t i=0; i<largest_k_eigenvectors.rows(); ++i) {
      //   assert(std::abs( largest_k_eigenvectors.row(i).norm() - 1.0) < 1e-8);
      //}

      U = largest_k_eigenvectors * largest_k_eigenvectors.transpose();
      assert(U.cols() == mgm_size.total_no_nodes() && U.rows() == mgm_size.total_no_nodes());
      //std::cout << "eigenvalues = " << U_tmp.eigenvalues().transpose() << "\n";
      //std::cout << U << "\n";

      return U; 
   }

   void transform_universe_matching_to_multigraph_matching(const multigraph_matching_input::graph_size& mgm_size, const std::vector<std::size_t>& Y, multigraph_matching_input::labeling& labeling)
   {
      const std::size_t no_universe_objects = *std::max_element(Y.begin(), Y.end()) + 1;
      union_find uf(Y.size());
      for(std::size_t c=0; c<no_universe_objects; ++c) {
         std::size_t prev_i = std::numeric_limits<std::size_t>::max();
         for(std::size_t i=0; i<Y.size(); ++i) {
            if(Y[i] == c) {
               if(prev_i != std::numeric_limits<std::size_t>::max()) {
                  uf.merge(i, prev_i);
               }
               prev_i = i;
            } 
         }
      }

      for(auto& gm : labeling) {
         std::fill(gm.labeling.begin(), gm.labeling.end(), std::numeric_limits<std::size_t>::max());
         for(std::size_t i=0; i<mgm_size.no_nodes(gm.left_graph_no); ++i) {
            for(std::size_t j=0; j<mgm_size.no_nodes(gm.right_graph_no); ++j) {
               assert(mgm_size.node_no(gm.left_graph_no,i) != mgm_size.node_no(gm.right_graph_no,j));
               if(uf.connected(mgm_size.node_no(gm.left_graph_no,i), mgm_size.node_no(gm.right_graph_no,j))) {
                  gm.labeling[i] = j;
               }
            }
         }
      }
   }

   static constexpr std::size_t node_not_taken = std::numeric_limits<std::size_t>::max();
   void synchronize_multigraph_matching(const multigraph_matching_input::graph_size& mgm_size, multigraph_matching_input::labeling& labeling, const double rounding_th)
   {
      const Eigen::MatrixXd U = compute_eigenvector_matrix(mgm_size, labeling);

      std::vector<std::size_t> Y(mgm_size.total_no_nodes(), node_not_taken);

      std::size_t c = 0; // current object in universe which we want to assign to
      for(std::size_t i=0; i<mgm_size.total_no_nodes(); ++i) {
         if(Y[i] == node_not_taken) {
            Y[i] = c;
            const auto [no_graph, no_node] = mgm_size.graph_node_no(i);
            for(std::size_t q=0; q<mgm_size.no_graphs(); ++q) {
               if(q != no_graph) {
                  std::size_t max_j_idx = std::numeric_limits<std::size_t>::max();
                  double max_j_val = -std::numeric_limits<double>::max();
                  for(std::size_t j=0; j<mgm_size.no_nodes(q); ++j) {
                     const std::size_t idx = mgm_size.node_no(q,j);
                     if(Y[idx] == node_not_taken && max_j_val < U(i,idx)) {
                        max_j_idx = j;
                        max_j_val = U(i,idx);
                     }
                  } 
                  if(max_j_val > rounding_th) {
                     Y[mgm_size.node_no(q, max_j_idx)] = c;
                  } 
               } 
            } 
            ++c;
         } 
      }

      transform_universe_matching_to_multigraph_matching(mgm_size, Y, labeling);
      //labeling.write_primal_matching(std::cout);
      assert(labeling.check_primal_consistency());
   }

   void synchronize_multigraph_matching(const multigraph_matching_input::graph_size& mgm_size, multigraph_matching_input::labeling& labeling, const matrix<int>& allowed_matchings, const double rounding_th)
   {
      assert(mgm_size.total_no_nodes() == allowed_matchings.dim1());
      assert(mgm_size.total_no_nodes() == allowed_matchings.dim2());

      const Eigen::MatrixXd U = compute_eigenvector_matrix(mgm_size, labeling);

      std::vector<std::size_t> Y(mgm_size.total_no_nodes(), node_not_taken); // node to universe matching
      std::size_t c = 0; // current object in universe which we want to assign to
      std::vector<std::size_t> current_cluster_elements;
      for(std::size_t i=0; i<mgm_size.total_no_nodes(); ++i) {
         if(Y[i] == node_not_taken) {
            current_cluster_elements.clear();
            Y[i] = c;
            current_cluster_elements.push_back(i);
            const auto [no_graph, no_node] = mgm_size.graph_node_no(i);
            for(std::size_t q=0; q<mgm_size.no_graphs(); ++q) {
               if(q != no_graph) {
                  std::size_t max_j_idx = std::numeric_limits<std::size_t>::max();
                  double max_j_val = -std::numeric_limits<double>::max();
                  for(std::size_t j=0; j<mgm_size.no_nodes(q); ++j) {
                     const std::size_t idx = mgm_size.node_no(q,j);
                     if(Y[idx] == node_not_taken && max_j_val < U(i,idx)) {
                        bool feasible = true;
                        for(const std::size_t prev_i : current_cluster_elements)
                           if(!allowed_matchings(prev_i, mgm_size.node_no(q, j)))
                              feasible = false;
                        if(feasible) {
                           max_j_idx = j;
                           max_j_val = U(i,idx);
                        }
                     }
                  } 
                  if(max_j_val > rounding_th) {
                     Y[mgm_size.node_no(q, max_j_idx)] = c;
                     current_cluster_elements.push_back(mgm_size.node_no(q, max_j_idx));
                  } 
               } 
            }
            ++c;
         }
      }

      transform_universe_matching_to_multigraph_matching(mgm_size, Y, labeling);
      if(!labeling.check_primal_consistency())
         throw std::runtime_error("labeling after rounding not consistent.");

      // check whether labeling is feasible
      for(const auto& gm : labeling) {
         for(std::size_t i=0; i<gm.labeling.size(); ++i) {
            if(gm.labeling[i] < std::numeric_limits<std::size_t>::max()) {
               if(!allowed_matchings(mgm_size.node_no(gm.left_graph_no, i), mgm_size.node_no(gm.right_graph_no, gm.labeling[i]))) {
                  std::cout << "(" << gm.left_graph_no << "," << i << " -> (" << gm.right_graph_no << "," << gm.labeling[i] << ")\n";
               }
               assert(allowed_matchings(mgm_size.node_no(gm.left_graph_no, i), mgm_size.node_no(gm.right_graph_no, gm.labeling[i])));
            }
         }
      }
   }

} // namespace LPMP
