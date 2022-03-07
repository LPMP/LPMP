/*
 *  Clean and efficient reimplementation of separation for k-ary cycle inequalities, as proposed by David Sontag, Do Kook Choe and Yitao Li in UAI 2012 paper.
 */

#ifndef LPMP_CYCLE_INEQUALITIES_HXX
#define LPMP_CYCLE_INEQUALITIES_HXX

#include <limits>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
//#include <list>
//#include <map>
//#include <queue>

#include "config.hxx"
#include "vector.hxx"
#include "graph.hxx"
#include "union_find.hxx"

namespace LPMP {

   struct triplet_candidate {
      triplet_candidate(const std::size_t i1, const std::size_t i2, const std::size_t i3, const double c)
         : cost(c)
      {
         i = std::min({i1,i2,i3});
         j = std::max(std::min(i1,i2), std::min(std::max(i1,i2),i3));
         k = std::max({i1,i2,i3}); 
         assert(i < j && j < k);
      }

      bool operator<(const triplet_candidate& o) const { return std::abs(cost) > std::abs(o.cost); }
      bool operator==(const triplet_candidate& o) const { return i == o.i && j == o.j && k == o.k; }

      double cost;
      std::size_t i,j,k; 
   };

template<typename MRF_CONSTRUCTOR>
class triplet_search 
{
public:
   typedef std::vector<std::vector<std::pair<std::size_t, double> > > adj_type;

   triplet_search(const MRF_CONSTRUCTOR& mrf, const double epsilon = eps)
      : gm_(mrf),
      eps_(epsilon)
   {}
   ~triplet_search() {};


   std::vector<triplet_candidate> search()
   {
      // Initialize adjacency list (filled in later) TODO: only do this when needed
      std::vector<std::vector<std::size_t> > adjacency_list(gm_.get_number_of_variables());

      // Construct adjacency list for the graph
      // Iterate over all of the edges (we do this by looking at the edge intersection sets)
      for(size_t factorId=0; factorId<gm_.get_number_of_pairwise_factors(); factorId++) {
         auto vars = gm_.get_pairwise_variables(factorId);
         // Get the two nodes i & j
         const size_t i=std::get<0>(vars);
         const size_t j=std::get<1>(vars);
         assert(i<j);
         adjacency_list[i].push_back(j);
         adjacency_list[j].push_back(i);
      }

      // Sort the adjacency list, for fast intersections later
//#pragma omp parallel for schedule(guided)
      for(std::size_t i=0; i < adjacency_list.size(); i++) {
         std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
      }

      std::vector<triplet_candidate> triplet_candidates;

      // Iterate over all of the edge intersection sets
//#pragma omp parallel 
      {
         std::vector<std::size_t> commonNodes(gm_.get_number_of_variables()-1);
         std::vector<triplet_candidate> triplet_candidates_local;
//#pragma omp for schedule(guided)
         for(size_t factorId=0; factorId<gm_.get_number_of_pairwise_factors(); factorId++) {
            auto vars = gm_.get_pairwise_variables(factorId);
            const std::size_t i=std::get<0>(vars);
            const std::size_t j=std::get<1>(vars);
            const auto& factor_ij = *gm_.get_pairwise_factor(i,j)->get_factor();
            const double lb_ij = factor_ij.LowerBound();

            // Now find all neighbors of both i and j to see where the triangles are
            auto intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());
            assert(adjacency_list[i].size() <= commonNodes.size() && adjacency_list[j].size() <= commonNodes.size());

            for(auto n=commonNodes.begin(); n != intersects_iter_end; ++n) {
               std::size_t k = *n;

               // Since a triplet shows up multiple times, we only consider it for the case when i<j<k 
               if(!(j<k))
                  continue;

               const auto& factor_ik = *gm_.get_pairwise_factor(i,k)->get_factor();
               const auto& factor_jk = *gm_.get_pairwise_factor(j,k)->get_factor();
               const double boundIndep = lb_ij + factor_ik.LowerBound() + factor_jk.LowerBound();
               const double boundCycle = minimizeTriangle(factor_ij, factor_ik, factor_jk);

               const double bound = boundCycle - boundIndep; 
               assert(bound >=  - eps);
               if(bound > eps_) {
                  triplet_candidate t(i,j,k, bound);
                  triplet_candidates_local.push_back(t);
               }
            }
         }
//#pragma omp critical
         {
            triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_local.begin(), triplet_candidates_local.end());
         }
      }

      std::sort(triplet_candidates.begin(), triplet_candidates.end());

      return std::move(triplet_candidates);
   }

protected:
   template<typename PAIRWISE_REPAM>
   double minimizeTriangle(const PAIRWISE_REPAM& factor_ij, const PAIRWISE_REPAM& factor_ik, const PAIRWISE_REPAM& factor_jk) const
   {
      double max_val = std::numeric_limits<double>::infinity();

      // Fix value of the first variable
      const std::size_t dim1 = factor_ij.dim1();
      const std::size_t dim2 = factor_ij.dim2();
      const std::size_t dim3 = factor_ik.dim2();
      for(std::size_t i1=0; i1<dim1; ++i1) {
         for(std::size_t i2=0; i2<dim2; ++i2) {
            for(std::size_t i3=0; i3<dim3; ++i3) {
               max_val = std::min(max_val, factor_ij(i1,i2) + factor_ik(i1,i3) + factor_jk(i2,i3));
            }
         }
      }
      return max_val;
   }

   const MRF_CONSTRUCTOR& gm_;
   const double eps_;
};

// search for violated k-ary cycles either with k-projection graph or with expanded projection graph (indicated by EXTENDED template flag)
template<typename MRF_CONSTRUCTOR, bool EXTENDED=false>
class k_ary_cycle_inequalities_search
{
public:
   k_ary_cycle_inequalities_search(const MRF_CONSTRUCTOR& mrf, const double epsilon = eps) 
      : gm_(mrf),
      eps_(epsilon)
   {}

   std::vector<triplet_candidate> search(const std::size_t max_triplets = std::numeric_limits<std::size_t>::max())
   {
      if(!EXTENDED) {
         construct_projection_graph();
      } else {
         auto partitions = compute_partitions(); 
         construct_projection_graph(std::move(partitions));
      }
      return find_cycles(max_triplets);
   }

   template<typename PAIRWISE_REPAM>
      static matrix<double> row_minima(const PAIRWISE_REPAM& f);
   template<typename PAIRWISE_REPAM>
      static matrix<double> column_minima(const PAIRWISE_REPAM& f);
   template<typename PAIRWISE_REPAM>
      static matrix<double> principal_minima(const PAIRWISE_REPAM& f, const matrix<double>& _column_minima);

protected:

   template<typename PAIRWISE_REPAM> 
   std::pair<std::vector<bool>,std::vector<bool>> compute_partitions(const PAIRWISE_REPAM& f);
   std::vector<std::vector<std::vector<bool>>> compute_partitions(); 

   void construct_projection_graph(const std::vector<std::vector<std::vector<bool>>> partitions = {});

   template<typename PAIRWISE_REPAM, typename V>
   double compute_projection_weight_singleton_2(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::size_t x2, const V& row_minima);
   template<typename PAIRWISE_REPAM>
   double compute_projection_weight_singleton_2_enumerate(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::size_t x2);
   template<typename PAIRWISE_REPAM, typename V>
   double compute_projection_weight_singleton_1(const PAIRWISE_REPAM& f, const std::size_t x1, const std::vector<bool>& part_j, const V& column_minima);
   template<typename PAIRWISE_REPAM>
   double compute_projection_weight_singleton_1_enumerate(const PAIRWISE_REPAM& f, const std::size_t x1, const std::vector<bool>& part_j);
   template<typename PAIRWISE_REPAM>
   double compute_projection_weight_on_partitions(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::vector<bool>& part_j);

   std::vector<triplet_candidate> find_cycles(const std::size_t max_triplets);
   void triangulate(std::vector<triplet_candidate>& triplet_candidates, std::vector<std::size_t> cycle, const double th);

   struct weighted_edge : public std::array<std::size_t,2> 
   { 
       weighted_edge(const std::size_t i, const std::size_t j, const double c) : std::array<std::size_t,2>({i,j}), cost(c) {}
       double cost; 
       bool operator<(const weighted_edge& o) const { return cost > o.cost; } 
   };
   std::vector<weighted_edge> projection_edges_;
   std::vector<std::size_t> proj_graph_offsets_;
   std::vector<std::size_t> proj_graph_to_gm_node_;
   std::vector<weighted_edge> proj_graph_edges_;
   graph<double> proj_graph_;

   const MRF_CONSTRUCTOR& gm_;
   const double eps_;
};


// possibly do the three claculations below in one pass for greater efficiency
template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
matrix<double> k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::row_minima(const PAIRWISE_REPAM& f)
{
   matrix<double> _row_minima(f.dim1(),2);
   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      // find smallest and second smallest entry w.r.t. second index given current first index
      double smallest = std::numeric_limits<double>::infinity();
      double second_smallest = std::numeric_limits<double>::infinity();
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         const double val = f(x1,x2);
         const double min = std::min(smallest, val);
         const double max = std::max(smallest, val);
         smallest = min;
         second_smallest = std::min(max, second_smallest);
      }
      _row_minima(x1,0) = smallest;
      _row_minima(x1,1) = second_smallest;
   }
   return std::move(_row_minima);
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
matrix<double> k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::column_minima(const PAIRWISE_REPAM& f)
{
   matrix<double> _column_minima(f.dim2(),2);
   for(size_t i=0; i<_column_minima.size(); ++i)
       _column_minima[i] = std::numeric_limits<double>::infinity();
   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         const double val = f(x1,x2);
         const double min = std::min(_column_minima(x2,0), val);
         const double max = std::max(_column_minima(x2,0), val);
         _column_minima(x2,0) = min;
         _column_minima(x2,1) = std::min(max, _column_minima(x2,1));
      }
   }
   return std::move(_column_minima);
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
matrix<double> k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::principal_minima(const PAIRWISE_REPAM& f, const matrix<double>& _column_minima)
{
   // possibly can be computed more efficiently. Note that at most 4 values need to be stored
   matrix<double> _principal_minima(f.dim1(), f.dim2());
   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      double smallest = std::numeric_limits<double>::infinity();
      double second_smallest = std::numeric_limits<double>::infinity();
      std::size_t smallest_ind = 0;
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         const double val_xij = f(x1,x2);
         const double x = _column_minima(x2,0) == val_xij ? _column_minima(x2,1) : _column_minima(x2,0);
         if(x<smallest) {
            smallest = x;
            smallest_ind = x2;
         }
      }
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         const double val_xij = f(x1,x2);
         const double x = _column_minima(x2,0) == val_xij ? _column_minima(x2,1) : _column_minima(x2,0);
         if(x < second_smallest && x2 != smallest_ind) {
            second_smallest = x;
         }
      }
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         _principal_minima(x1,x2) = smallest;
      }
      _principal_minima(x1,smallest_ind) = second_smallest; 
   }
   return std::move(_principal_minima);
}

// Given an undirected graph, finds odd-signed cycles.
template<typename MRF_CONSTRUCTOR, bool EXTENDED>
std::vector<triplet_candidate> 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::find_cycles(const std::size_t max_triplets)
{
   double largest_th = 0.0;
   union_find uf(proj_graph_.no_nodes());
   std::size_t e=0;

   auto merge_edge = [&uf](const std::size_t m, const std::size_t n, const double s) {
      assert(s != 0.0);
      if(s < 0.0) {
        uf.merge(2*n,2*m);
        uf.merge(2*n+1,2*m+1);
     } else {
        uf.merge(2*n,2*m+1);
        uf.merge(2*n+1,2*m);
     }
   };

   assert(std::is_sorted(projection_edges_.begin(), projection_edges_.end()));

   for(; e<projection_edges_.size(); ++e) {
      const std::size_t i = projection_edges_[e][0];
      const std::size_t j = projection_edges_[e][1];
      const double s = projection_edges_[e].cost;
      merge_edge(i,j,s);

      if(uf.connected(2*i,2*i+1) || uf.connected(2*j,2*j+1)) {
         largest_th = std::abs(s);
         break;
      } 
   }

   std::vector<triplet_candidate> triplet_candidates;
   std::vector<bool> already_searched(proj_graph_to_gm_node_.size(),false);
   // first update union find datastructure by merging additional edges with cost greater than th
   double th = 0.5*largest_th;
   for(std::size_t iter=0; iter<8 && th>=eps_; ++iter, th*=0.1) {
      // update connectivity information
      for(; e<projection_edges_.size(); ++e) {
         const std::size_t i = projection_edges_[e][0];
         const std::size_t j = projection_edges_[e][1];
         const double s = projection_edges_[e].cost;
         if(std::abs(s) >= th) {
            merge_edge(i,j,s);
         } else {
            break;
         }
      }

      // now actually search for odd signed cycles
#pragma omp parallel
      {
         std::vector<triplet_candidate> triplet_candidates_local;
         bfs_data bfs(proj_graph_);
#pragma omp for schedule(guided)
         for(std::size_t i=0; i<proj_graph_to_gm_node_.size(); ++i) {
            if(!already_searched[i] && uf.thread_safe_connected(2*i, 2*i+1)) {
               already_searched[i] = true;
               double th = std::numeric_limits<double>::infinity();
               auto cycle = bfs.find_path(2*i, 2*i+1, decltype(bfs)::no_mask_op, [&th](const std::size_t i, const std::size_t j, const double w) { th = std::min(th,w); });
               assert(cycle.size() >= 3);
               if(cycle.size() >= 3) {
                  triangulate(triplet_candidates_local, std::move(cycle), th);
               }
            }
         }
#pragma omp critical
         {
            triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_local.begin(), triplet_candidates_local.end()); 
         }
      }
      if(triplet_candidates.size() > max_triplets) {
         break;
      }
   }

   std::sort(triplet_candidates.begin(), triplet_candidates.end());
   if(triplet_candidates.size() > 0) {
      assert(triplet_candidates[0].cost >= triplet_candidates.back().cost);
   }
   triplet_candidates.erase( unique( triplet_candidates.begin(), triplet_candidates.end() ), triplet_candidates.end() );
   return std::move(triplet_candidates);
}

 
template<typename MRF_CONSTRUCTOR, bool EXTENDED>
void
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::triangulate(std::vector<triplet_candidate>& triplet_candidates, std::vector<std::size_t> cycle, const double th)
{
   // halve indices, so we get node indices in projection graph, then project back onto indices of graphical model
   for(std::size_t& i : cycle) {
      i = proj_graph_to_gm_node_[i/2];
   }
   assert(cycle[0] == cycle.back());
   cycle.resize(cycle.size()-1);

   // possibly cycle has a sub-cycle, or even uses an edge twice consecutively, still add it

   assert(cycle.size() >= 3);

   for(std::size_t i=1; i<cycle.size()-1; ++i) {
      const std::size_t j = cycle[i];
      const std::size_t k = cycle[i+1];
      assert(j!=k);
      if(j != cycle[0] && k != cycle[0]) {
         triplet_candidates.push_back(triplet_candidate(cycle[0], j,k, th));
      }
   }
   /*
   auto min_node = std::min_element(cycle.begin(), cycle.end());
   for(auto it=cycle.begin(); it<min_node-2; ++it) {
      triplet_candidates.push_back(triplet_candidate(*it,*(it+1),*min_node, th));
   }
   for(auto it=min_node+1; it<cycle.end()-2; ++it) {
      triplet_candidates.push_back(triplet_candidate(*it,*(it+1),*min_node, th));
   }
   if(min_node != cycle.begin() && min_node != cycle.end()-1) { // when min_node is somewhere in the middle.
      triplet_candidates.push_back(triplet_candidate(*cycle.begin(), *cycle.rbegin(), *min_node, th));
   }
   */
}


template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM> 
std::pair<std::vector<bool>,std::vector<bool>> 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_partitions(const PAIRWISE_REPAM& f)
{
   if(f.dim1() <= 3 && f.dim2() <= 3) {
      return std::make_pair(std::vector<bool>{}, std::vector<bool>{});
   }
   std::vector<std::tuple<std::size_t,std::size_t,double>> sorted_factor(f.dim1()*f.dim2());
   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         sorted_factor[x1*f.dim2() + x2] = std::make_tuple(x1,x2,f(x1,x2));
      }
   }
   std::sort(sorted_factor.begin(), sorted_factor.end(), [](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });
   union_find uf(f.dim1() + f.dim2());
   std::vector<bool> part_i(f.dim1());
   std::vector<bool> part_j(f.dim2());
   for(auto x : sorted_factor) {
      const std::size_t x1 = std::get<0>(x);
      const std::size_t x2 = std::get<1>(x);

      if(!uf.connected(x1, f.dim1() + x2)) {
         // check if merging c1 and c2 would result in one partition
         if(uf.count() > 2) {
            uf.merge(x1, f.dim1() + x2);
         } else {
            assert(uf.count() == 2);
            const std::size_t c1 = uf.find(x1);
            //const std::size_t c2 = uf.find(f.dim1() + x2);
            // record labels of partition c1
            for(std::size_t y1=0; y1<f.dim1(); ++y1) {
               if(uf.find(y1) == c1) {
                  part_i[y1] = true; 
               } else {
                  part_i[y1] = false; 
               }
            }
            for(std::size_t y2=0; y2<f.dim2(); ++y2) {
               if(uf.find(f.dim1() + y2) == c1) {
                  part_j[y2] = true; 
               } else {
                  part_j[y2] = false; 
               }
            }
            break;
         }
      }
   }
   // normalize partitions
   if(part_i[0] == false) {
      part_i.flip();
   }
   if(part_j[0] == false) {
      part_j.flip();
   }
   // if partition consists of a singleton on any side, do not add it
   const std::size_t sum_1 = std::count(part_i.begin(), part_i.end(), true);
   if(sum_1 == 1 || sum_1 == f.dim1()-1) {
      part_i.clear();
   }
   const std::size_t sum_2 = std::count(part_j.begin(), part_j.end(), true);
   if(sum_2 == 1 || sum_2 == f.dim2()-1) {
      part_j.clear();
   }
   if(f.dim1() <= 3) {
      part_i.clear();
   }
   if(f.dim2() <= 3) {
      part_j.clear();
   }

   return std::make_pair(std::move(part_i), std::move(part_j));
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
std::vector<std::vector<std::vector<bool>>> 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_partitions()
{
   std::vector<std::vector<std::vector<bool>>> partitions(gm_.get_number_of_variables());
#pragma omp parallel
   {
      auto partitions_local = partitions;
#pragma omp for schedule(guided) nowait
      for(size_t factorId=0; factorId<gm_.get_number_of_pairwise_factors(); factorId++) {
         auto vars = gm_.get_pairwise_variables(factorId);
         const auto& factor = *gm_.get_pairwise_factor(factorId)->get_factor();
         const size_t i=std::get<0>(vars);
         const size_t j=std::get<1>(vars);
         assert(i<j);
         const auto part_ij = compute_partitions(factor);
         const auto& part_i = std::get<0>(part_ij);
         if(part_i.size() > 0) {
            partitions_local[i].push_back(std::move(part_i));
         }
         const auto& part_j = std::get<1>(part_ij);
         if(part_j.size() > 0) {
            partitions_local[j].push_back(std::move(part_j));
         }
      }
#pragma omp critical
      {
         for(std::size_t i=0; i<partitions.size(); ++i) {
            partitions[i].insert(partitions[i].end(), partitions_local[i].begin(), partitions_local[i].end());
         }
      }
   }

   // remove duplicate partitions
#pragma omp parallel for schedule(guided)
   for(std::size_t i=0; i<partitions.size(); ++i) {
      std::sort(partitions[i].begin(), partitions[i].end());
      partitions[i].erase( unique( partitions[i].begin(), partitions[i].end() ), partitions[i].end() );
   }

   return partitions;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM, typename V>
double 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_singleton_2(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::size_t x2, const V& _row_minima)
{
   // same dimension if part_i[x1] and x2
   assert(f.dim1() == part_i.size());
   assert(f.dim1() == _row_minima.dim1() && 2 == _row_minima.dim2());
   double min_same_part = std::numeric_limits<double>::infinity();
   double min_different_part = std::numeric_limits<double>::infinity();
   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      if(part_i[x1]) {
         min_same_part = std::min(min_same_part, f(x1,x2)); 
      } else {
         min_different_part = std::min(min_different_part, f(x1,x2));
      }
   }

   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      const double val_xij = f(x1,x2);
      const double row_min = _row_minima(x1,0) == val_xij ? _row_minima(x1,1) : _row_minima(x1,0);
      if(part_i[x1]) {
         min_different_part = std::min(min_different_part, row_min);
      } else {
         min_same_part = std::min(min_same_part, row_min);
      }
   }

   assert(min_same_part - min_different_part == compute_projection_weight_singleton_2_enumerate(f, part_i, x2));

   return min_same_part - min_different_part; 
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
double 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_singleton_2_enumerate(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::size_t x2)
{
double min_same_part = std::numeric_limits<double>::infinity();
   double min_different_part = std::numeric_limits<double>::infinity();
   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      for(std::size_t x2_iter=0; x2_iter<f.dim2(); ++x2_iter) {
         const double val_xij = f(x1,x2_iter);
         if(x2 == x2_iter) {
            min_same_part = std::min(min_same_part, val_xij);
         } else {
            min_different_part = std::min(min_different_part, val_xij);
         }
      }
   }
   assert(min_same_part == min_same_part && min_different_part == min_different_part);

   return min_same_part - min_different_part; 
}


template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM, typename V>
double 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_singleton_1(const PAIRWISE_REPAM& f, const std::size_t x1, const std::vector<bool>& part_j, const V& _column_minima)
{
   assert(f.dim2() == part_j.size());
   assert(f.dim2() == _column_minima.dim1() && 2 == _column_minima.dim2());
   double min_same_part = std::numeric_limits<double>::infinity();
   double min_different_part = std::numeric_limits<double>::infinity();
   for(std::size_t x2=0; x2<f.dim2(); ++x2) {
      if(part_j[x2]) {
         min_same_part = std::min(min_same_part, f(x1,x2)); 
      } else {
         min_different_part = std::min(min_different_part, f(x1,x2));
      }
   }

   for(std::size_t x2=0; x2<f.dim2(); ++x2) {
      const double val_xij = f(x1,x2);
      const double column_min = _column_minima(x2,0) == val_xij ? _column_minima(x2,1) : _column_minima(x2,0);
      if(part_j[x2]) {
         min_different_part = std::min(min_different_part, column_min);
      } else {
         min_same_part = std::min(min_same_part, column_min);
      }
   }

   assert(min_same_part - min_different_part == compute_projection_weight_singleton_1_enumerate(f,x1,part_j));
   return min_same_part - min_different_part;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
double 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_singleton_1_enumerate(const PAIRWISE_REPAM& f, const std::size_t x1, const std::vector<bool>& part_j)
{
   double min_same_part = std::numeric_limits<double>::infinity();
   double min_different_part = std::numeric_limits<double>::infinity();
   for(std::size_t x1_iter=0; x1_iter<f.dim1(); ++x1_iter) {
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         const double val_xij = f(x1_iter,x2);
         if(x1 == x1_iter) {
            min_same_part = std::min(min_same_part, val_xij);
         } else {
            min_different_part = std::min(min_different_part, val_xij);
         }
      }
   }
   return min_same_part - min_different_part;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
double 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_on_partitions(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::vector<bool>& part_j)
{
   assert(f.dim1() == part_i.size());
   assert(f.dim2() == part_j.size());
   double min_same_part = std::numeric_limits<double>::infinity();
   double min_different_part = std::numeric_limits<double>::infinity();

   for(std::size_t x1=0; x1<f.dim1(); ++x1) {
      for(std::size_t x2=0; x2<f.dim2(); ++x2) {
         const double val_xij = f(x1,x2);
         if(part_i[x1] == part_j[x2]) {
            min_same_part = std::min(min_same_part, val_xij);
         } else {
            min_different_part = std::min(min_different_part, val_xij); 
         }
      }
   }
   return min_same_part - min_different_part;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
void 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::construct_projection_graph(const std::vector<std::vector<std::vector<bool>>> partitions)
{
   proj_graph_offsets_ = std::vector<std::size_t>(gm_.get_number_of_variables());
   proj_graph_offsets_[0] = 0;
   for(std::size_t i=1; i<gm_.get_number_of_variables(); i++) {
      proj_graph_offsets_[i] = gm_.get_number_of_labels(i-1);
      if(EXTENDED) { 
         proj_graph_offsets_[i] += partitions[i-1].size();
      }
   }

   std::partial_sum(proj_graph_offsets_.begin(), proj_graph_offsets_.end(), proj_graph_offsets_.begin());
   std::size_t proj_graph_nodes = proj_graph_offsets_.back() + gm_.get_number_of_labels(gm_.get_number_of_variables()-1);
   if(EXTENDED) {
      proj_graph_nodes += partitions.back().size();
   }

   proj_graph_to_gm_node_ = std::vector<std::size_t>(proj_graph_nodes);
   {
      std::size_t c=0;
      for(std::size_t i=0; i<gm_.get_number_of_variables(); i++) {
         std::size_t no_proj_graph_nodes_for_label = gm_.get_number_of_labels(i);
         if(EXTENDED) {
            no_proj_graph_nodes_for_label += partitions[i].size();
         }
         for(std::size_t x=0; x<no_proj_graph_nodes_for_label; ++x) {
            proj_graph_to_gm_node_[c] = i;
            ++c;
         }
      }
      assert(c == proj_graph_nodes);
   }

   projection_edges_.clear();

   auto add_to_projection_edges = [this](auto& projection_edges, const std::size_t n, const std::size_t m, const double val) {
      if(std::abs(val) >= eps_ && !std::isnan(val)) {          
         projection_edges.push_back(weighted_edge(m,n,val));
      } 
   };

#pragma omp parallel
   {
      auto projection_edges_local = projection_edges_;
#pragma omp for schedule(guided)
      for(size_t factorId=0; factorId<gm_.get_number_of_pairwise_factors(); factorId++) {
         // Get the two nodes i & j and the edge intersection set. Put in right order.
         const std::size_t i = std::get<0>(gm_.get_pairwise_variables(factorId));
         const std::size_t j = std::get<1>(gm_.get_pairwise_variables(factorId));

         // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
         if(gm_.get_number_of_labels(i) <= 1 || gm_.get_number_of_labels(j) <= 1)
            continue;

         // For each of their singleton states efficiently compute edge weights
         const auto& factor_ij = *gm_.get_pairwise_factor(i,j)->get_factor(); // better retrieve by factor id

         const auto row_min = row_minima(factor_ij);
         const auto col_min = column_minima(factor_ij);
         const auto principal_min = principal_minima(factor_ij, col_min);

         assert(i<j);
         for(std::size_t xi=0; xi<factor_ij.dim1(); xi++) {
            const std::size_t m = proj_graph_offsets_[i] + xi;
            assert(i == proj_graph_to_gm_node_[m]);

            for(std::size_t xj=0; xj<factor_ij.dim2(); xj++) {
               const std::size_t n = proj_graph_offsets_[j] + xj;
               assert(j == proj_graph_to_gm_node_[n]);

               const double val_xij = factor_ij(xi,xj);

               const double val_not_xi = row_min(xi,0) == val_xij ? row_min(xi,1) : row_min(xi,0);
               const double val_not_xj = col_min(xj,0) == val_xij ? col_min(xj,1) : col_min(xj,0);

               // val_s < 0 means same projection < different projection, > 0 the opposite
               // Hence we search for a cycle with an odd number of entries > 0      
               const double cost_projection_same = std::min(val_xij, principal_min(xi,xj));
               const double cost_projection_different = std::min(val_not_xi, val_not_xj);
               const double val_s = cost_projection_same - cost_projection_different;

               // TODO: use threshold here, to make next stage faster
               add_to_projection_edges(projection_edges_local,n,m,val_s);
            }
         }

         if(EXTENDED) {
            // add edge weights between each general projection and each singleton state 
            for(std::size_t x1=0; x1<factor_ij.dim1(); ++x1) {
               const std::size_t m = proj_graph_offsets_[i] + x1;
               for(std::size_t p2=0; p2<partitions[j].size(); ++p2) {
                  const std::size_t n = proj_graph_offsets_[j] + gm_.get_number_of_labels(j) + p2;
                  const double val_s = compute_projection_weight_singleton_1(factor_ij, x1, partitions[j][p2], col_min);
                  add_to_projection_edges(projection_edges_local,n,m,val_s);
               }
            }
            for(std::size_t x2=0; x2<factor_ij.dim2(); ++x2) {
               const std::size_t m = proj_graph_offsets_[j] + x2;
               for(std::size_t p1=0; p1<partitions[i].size(); ++p1) {
                  const std::size_t n = proj_graph_offsets_[i] + gm_.get_number_of_labels(i) + p1;
                  const double val_s = compute_projection_weight_singleton_2(factor_ij, partitions[i][p1], x2, row_min);
                  add_to_projection_edges(projection_edges_local,n,m,val_s);
               }
            }

            // compute edge weights between general projections
            for(std::size_t p1=0; p1<partitions[i].size(); ++p1) {
               const std::size_t n = proj_graph_offsets_[i] + gm_.get_number_of_labels(i) + p1;
               for(std::size_t p2=0; p2<partitions[j].size(); ++p2) {
                  const std::size_t m = proj_graph_offsets_[j] + gm_.get_number_of_labels(j) + p2;
                  const double val_s = compute_projection_weight_on_partitions(factor_ij, partitions[i][p1], partitions[j][p2]);
                  add_to_projection_edges(projection_edges_local,n,m,val_s);
               }
            }
         }
      }
#pragma omp critical
      {
         projection_edges_.insert(projection_edges_.end(), projection_edges_local.begin(), projection_edges_local.end());
      }
   }

    std::vector<weighted_edge> doubled_edges;
    doubled_edges.reserve(2*projection_edges_.size());

    for(const auto edge : projection_edges_) {
        const auto cost = edge.cost;
        const auto n = edge[0];
        const auto m = edge[1];
        if(cost < 0.0) {
            doubled_edges.push_back(weighted_edge(2*n, 2*m, -cost));
            doubled_edges.push_back(weighted_edge(2*n+1, 2*m+1, -cost));
        } else {
            doubled_edges.push_back(weighted_edge(2*n, 2*m+1, cost));
            doubled_edges.push_back(weighted_edge(2*n+1, 2*m, cost));
        }
    }
   proj_graph_ = graph<double>(doubled_edges.begin(), doubled_edges.end(), [](const auto& weighted_graph) { return weighted_graph.cost; });

   std::sort(projection_edges_.begin(), projection_edges_.end());
}
} // end namespace LPMP

#endif // LPMP_CYCLE_INEQUALITIES_HXX

