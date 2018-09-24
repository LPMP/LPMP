#ifndef LPMP_SEQUENTIAL_EDGE_FIXATION_HXX
#define LPMP_SEQUENTIAL_EDGE_FIXATION_HXX

#include "union_find.hxx"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "multicut_text_input.h"

namespace LPMP {

// aka mutex watershed (Wolf & Pape et al, ECCV 2018, when edges are sorted), but more efficient implementation
class multicut_sequential_edge_fixation {
public:
   multicut_sequential_edge_fixation(const std::size_t no_nodes)
      : uf(no_nodes),
      constraint_connectivity(no_nodes)
   {}

   template<typename EDGE_ITERATOR, typename EDGE_FUNC, typename EDGE_SORT_FUNC>
   multicut_sequential_edge_fixation(const std::size_t no_nodes, EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end,
         EDGE_FUNC&& edge_func, EDGE_SORT_FUNC&& edge_sort_func, const double chunk_size_multiplier = 2)
   : multicut_sequential_edge_fixation(no_nodes)
   {
      assert(chunk_size_multiplier > 0.0);
      // TODO: measure whether this is faster than just sorting all edges at once
      while(!multicut_determined() && edge_begin != edge_end) {
         const std::size_t no_edges = std::distance(edge_begin, edge_end);
         const std::size_t no_considered_edges = std::min(std::size_t(std::ceil(chunk_size_multiplier*double(no_nodes))), no_edges);
         auto edge_middle = edge_begin + no_considered_edges;
         std::nth_element(edge_begin, edge_middle, edge_end, edge_sort_func);
         std::sort(edge_begin, edge_middle, edge_sort_func);
         add_constraints(edge_begin, edge_middle, edge_func);
         edge_begin = edge_middle;
      } 
   }

   // EDGE_FUNC takes as input an edge and outputs whether edge shall be connected
   template<typename EDGE_ITERATOR, typename EDGE_FUNC>
   multicut_sequential_edge_fixation(const std::size_t no_nodes, EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_FUNC&& edge_func)
   : multicut_sequential_edge_fixation(no_nodes)
   {
      add_constraints(edge_begin, edge_end, edge_func);
   }

   static bool edge_contractive(const multicut_instance::weighted_edge& e) { return e.cost > 0.0; }
   static bool edge_sort(const multicut_instance::weighted_edge& e1, const multicut_instance::weighted_edge& e2) { return std::abs(e1.cost) > std::abs(e2.cost); }

   multicut_sequential_edge_fixation(multicut_instance& input)
   : multicut_sequential_edge_fixation(input.no_nodes, input.edges.begin(), input.edges.end(), edge_contractive, edge_sort)
   {}

   template<typename EDGE_ITERATOR, typename EDGE_FUNC>
   bool add_constraints(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_FUNC&& edge_func)
   {
      bool changed = false;
      for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
         if(multicut_determined()) break;
         const std::size_t i = (*edge_it)[0];
         assert(i < no_nodes());
         const std::size_t j = (*edge_it)[1];
         assert(j < no_nodes());
         if(!constraint_connectivity.connected(i,j)) {
            changed = true;
            constraint_connectivity.merge(i,j);
            if( edge_func(*edge_it) ) {
               uf.merge(i,j);
            }
         } 
      } 
      return changed;
   } 

   bool connected(const std::size_t i, const std::size_t j) const
   {
      assert(multicut_determined());
      return uf.connected(i,j);
   } 

   bool multicut_determined() const
   {
      return constraint_connectivity.count() == 1;
   }

   std::size_t no_nodes() const
   {
      return uf.size();
   }

   template<typename STREAM, typename EDGE_ITERATOR>
   void write_edge_labeling(STREAM& s, EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end) const
   {
      assert(multicut_determined());
      for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
         const std::size_t i = (*edge_it)[0];
         assert(i < no_nodes());
         const std::size_t j = (*edge_it)[1];
         assert(j < no_nodes());
         s << i << " " << j << " ";
         s << (connected(i,j) ? 0 : 1) << "\n";
      } 
   }

   template<typename STREAM>
   void write_node_labeling(STREAM& s) const
   {
      assert(multicut_determined());
      const auto ids = uf.get_contiguous_ids();

      for(std::size_t i=0; i<no_nodes(); ++i)
         s << i << ": " << ids[uf.find(i)];
   }

private:
   mutable union_find uf;
   union_find constraint_connectivity; 
};

} // namespace LPMP

#endif // LPMP_SEQUENTIAL_EDGE_FIXATION_HXX
