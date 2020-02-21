#pragma once

#include "permutation.hxx"
#include <Eigen/OrderingMethods>
#include <vector>
#include <limits>

namespace LPMP {

    template<typename ADJACENCY_GRAPH>
    permutation minimum_degree_ordering(const ADJACENCY_GRAPH& adj)
    {
        Eigen::AMDOrdering<int> ordering;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm(adj.size());
        Eigen::SparseMatrix<double> A(adj.size(), adj.size()); 
        std::vector< Eigen::Triplet<double> > adjacency_list;
        for(std::size_t i=0; i<adj.size(); ++i) {
            for(const std::size_t j : adj[i]) {
                assert(j <= std::numeric_limits<int>::max());
                adjacency_list.push_back({i,j,1.0});
                adjacency_list.push_back({i,j,1.0});
            }
        }
        A.setFromTriplets(adjacency_list.begin(), adjacency_list.end());
        ordering(A, perm);
        permutation o;
        for(std::size_t i=0; i<perm.indices().size(); ++i) {
            o.push_back(perm.indices()[i]);
        }
        return o; 
    }


} 
