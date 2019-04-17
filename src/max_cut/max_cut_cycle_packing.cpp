#include "max_cut/max_cut_cycle_packing.h"
#include "graph.hxx"
#include "cut_base/cut_base_apply_packing.hxx"
#include "sequence_compression.h"
#include <iostream>
#include <array>
#include <vector>

namespace LPMP {

constexpr double tolerance = 1e-8; 
struct weighted_edge : public std::array<std::size_t,2> { double cost; };

cycle_packing max_cut_cycle_packing_impl(const max_cut_instance input, const bool record_cycles)
{
   double lower_bound = 0.0;

   std::vector<weighted_edge> bipartite_edges;
   for(const auto& e : input.edges()) {
      if(e.cost[0] < -tolerance) {
         bipartite_edges.push_back({e[0], e[1] + input.no_nodes(), -e.cost[0]});
         bipartite_edges.push_back({e[0] + input.no_nodes(), e[1], -e.cost[0]});
      } else if(e.cost[0] > tolerance) {
         bipartite_edges.push_back(weighted_edge{e[0], e[1], e.cost[0]});
         bipartite_edges.push_back(weighted_edge{e[0] + input.no_nodes(), e[1] + input.no_nodes(), e.cost[0]}); 
      }
      lower_bound += std::min(0.0, e.cost[0]);
   }

   std::cout << "max-cut cycle packing\n";
   std::cout << "initial lower bound = " << lower_bound << "\n";

   graph<double> g(bipartite_edges.begin(), bipartite_edges.end(), [](const weighted_edge& e) { return e.cost; });
   bfs_data<graph<double>> bfs(g);

   union_find uf(2*input.no_nodes());
   auto compute_connectivity = [&]() {
        uf.reset();
        g.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost) {
              if(cost >= tolerance) 
                 uf.merge(i,j);
        });
   };
   auto connected = [&](const std::size_t i, const std::size_t j) -> bool { return uf.connected(i,j); };
   auto cycles_present = [&]() {
      for(std::size_t i=0; i<input.no_nodes(); ++i)
         if(connected(i,i+input.no_nodes()))
            return true;
      return false;
   };

   cycle_packing cp;
   sequence_compression sc(input.no_nodes()); // for detecting subcycles
   // iteratively pack cycles of given length
   std::array<std::size_t,11> cycle_lengths = {1,2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
   for(const std::size_t cycle_length : cycle_lengths) {
      std::cout << "find cycles of length " << cycle_length << "\n";

      compute_connectivity();
      if(!cycles_present())
         break;

      for(std::size_t i=0; i<input.no_nodes(); ++i) {

          if(!connected(i, i+input.no_nodes()))
              continue;

         auto mask_small_edges = [cycle_length](const std::size_t i, const std::size_t j, const double cost, const std::size_t distance) { 
            if(cost <= tolerance) return false;
            if(distance >= cycle_length) return false;
            return true;
         };

         while(1) {
            double cycle_cap = std::numeric_limits<double>::infinity();
            auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const double cost) { cycle_cap = std::min(cycle_cap, cost); };
            auto cycle = bfs.find_path(i, i + input.no_nodes(), mask_small_edges, cycle_capacity); // TODO: possibly do not reallocate for every search but give as argument to find_path
            if(cycle.size() == 0) 
               break;

            lower_bound += cycle_cap;

            sc.reset();
            for(std::size_t c=1; c<cycle.size(); ++c) {
                const std::size_t node = cycle[c] % input.no_nodes();
                if(sc.index_added(node)) { // find first position of repeated
                    const std::size_t first_occurence = std::find_if(cycle.begin()+1, cycle.begin()+c-1, [=](const std::size_t x) { return x % input.no_nodes() == node; }) - cycle.begin();
                    std::copy(cycle.begin() + first_occurence, cycle.begin() + c+1, cycle.begin());
                    cycle.resize(c - first_occurence);
                    assert(cycle.size() > 0);
                    assert(cycle[0] % input.no_nodes() == cycle.back() % input.no_nodes());
                    break;
                }
                sc.add_index(node); 
            }

            // subtract minimum weight and remove edges of negligible weight
            for(std::size_t c=1; c<cycle.size(); ++c) {
               const std::size_t u = cycle[c-1];
               const std::size_t cu = u/input.no_nodes();
               const std::size_t v = cycle[c];
               const std::size_t cv = v/input.no_nodes();
               
               if(cu + cv == 0) {
                  g.edge(u,v) -= cycle_cap;
                  g.edge(v,u) -= cycle_cap;
                  g.edge(u + input.no_nodes(),v + input.no_nodes()) -= cycle_cap;
                  g.edge(v + input.no_nodes(),u + input.no_nodes()) -= cycle_cap;
               } else if(cu + cv == 1) {
                  const std::size_t uc1 = u % input.no_nodes();
                  const std::size_t uc2 = uc1 + input.no_nodes();
                  const std::size_t vc1 = v % input.no_nodes();
                  const std::size_t vc2 = vc1 + input.no_nodes();
                  g.edge(uc1, vc2) -= cycle_cap;
                  g.edge(vc2, uc1) -= cycle_cap;
                  g.edge(uc2, vc1) -= cycle_cap;
                  g.edge(vc1, uc2) -= cycle_cap;
               }  else if(cu + cv == 2) {
                  g.edge(u,v) -= cycle_cap;
                  g.edge(v,u) -= cycle_cap;
                  g.edge(u % input.no_nodes(),v % input.no_nodes()) -= cycle_cap;
                  g.edge(v % input.no_nodes(),u % input.no_nodes()) -= cycle_cap;
               } else {
                  assert(false);
               }
            }

            // compute cycle on original graph
            cycle.resize(cycle.size()-1);
            for(auto& x : cycle)
                x = x % input.no_nodes();

            if(record_cycles)
                cp.add_cycle(cycle.begin(), cycle.end(), cycle_cap); 
         }
      }
   }

   std::cout << "final lower bound after cycle packing= " << lower_bound << "\n";
   return cp;
}

void max_cut_cycle_packing(const max_cut_instance& input)
{
   max_cut_cycle_packing_impl(input, false); 
}
cycle_packing compute_max_cut_cycle_packing(const max_cut_instance& input)
{
   return max_cut_cycle_packing_impl(input, true);
}

triplet_max_cut_instance pack_max_cut_instance(const max_cut_instance& input, const cycle_packing& cp)
{
    return pack_cut_base_instance<max_cut_edge_triplet_message_0, max_cut_edge_triplet_message_1, max_cut_edge_triplet_message_2 , triplet_max_cut_instance>(input, cp);
}

}
