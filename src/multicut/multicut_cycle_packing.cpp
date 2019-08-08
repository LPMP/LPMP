#include "cut_base/cut_base_instance.hxx"
#include "multicut/multicut_cycle_packing.h"
#include "multicut/multicut_factors_messages.h"
#include "cut_base/cut_base_apply_packing.hxx"
#include "graph.hxx"
#include "union_find.hxx"
#include <array>
#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <cassert>

namespace LPMP {

constexpr static double tolerance = 1e-8;

struct weighted_edge : public std::array<std::size_t,2> { double cost; };

cycle_packing multicut_cycle_packing_impl(const multicut_instance& input, const bool record_cycles)
{
   const auto begin_time = std::chrono::steady_clock::now();

   double lower_bound = 0.0;
   std::vector<weighted_edge> repulsive_edges;
   std::vector<weighted_edge> positive_edges;
   for (const auto& e : input.edges()) {
      if(e.cost < 0.0)
         repulsive_edges.push_back({e[0], e[1], e.cost});
      else if(e.cost > 0.0)
         positive_edges.push_back({e[0], e[1], e.cost});
      lower_bound += std::min(0.0, e.cost[0]);
   }

   std::cout << "cycle packing\n";
   std::cout << "initial lower bound = " << lower_bound << "\n";
   std::cout << "#repulsive edges = " << repulsive_edges.size() << "\n";
   std::cout << "#attractive edges = " << positive_edges.size() << "\n";

   graph<double> pos_edges_graph(positive_edges.begin(), positive_edges.end(), [](const weighted_edge& e) { return e.cost; });
   cycle_packing cp;

   const auto initialization_end_time = std::chrono::steady_clock::now();
   std::cout << "initialization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(initialization_end_time - begin_time).count() << " milliseconds\n";

   // build up connectivity w.r.t. current positive edges
   union_find uf(input.no_nodes());
   auto compute_connectivity = [&]() {
        uf.reset();
        pos_edges_graph.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost) {
              if(cost >= tolerance) uf.merge(i,j);
        });
   };
   auto connected = [&](const std::size_t i, const std::size_t j) -> bool { return uf.connected(i,j); };

   // iteratively pack cycles of given length
   const std::array<std::size_t,11> cycle_lengths = {1,2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
   for(const std::size_t cycle_length : cycle_lengths) {
      std::cout << "find cycles of length " << cycle_length << ", with #repulsive edges = " << repulsive_edges.size() << " remaining, lower bound = " << lower_bound << "\n";
      //compute_connectivity();

      // shuffling can give great speed-up if edges with similar indices are spatially close in the graph
      std::random_shuffle(repulsive_edges.begin(), repulsive_edges.end());

      std::size_t progress = 0;

      bfs_data<graph<double>> bfs(pos_edges_graph);

      for(std::size_t i=0; i<repulsive_edges.size(); ++i) {
         auto& re = repulsive_edges[i];

         // periodically update component labeling for speed-up
         progress++;
         if (progress > 0.05 * repulsive_edges.size()) {
            compute_connectivity();
            progress = 0;
         }

         // check if conflicted cycle exists and repulsive edge has positive weight
         if(-re.cost <= tolerance || !connected(re[0], re[1]) || std::max(re[0], re[1]) >= pos_edges_graph.no_nodes()) {
            const std::size_t imax = repulsive_edges.size() - 1;
            repulsive_edges[i] = repulsive_edges[imax];
            repulsive_edges.resize(imax);
            i--;
            continue;
         }

         // balance short cycles as long as available and positive weight left
         auto mask_small_edges = [cycle_length](const std::size_t i, const std::size_t j, const double cost, const std::size_t distance) {
            if(cost <= tolerance) return false;
            if(distance >= cycle_length) return false;
            return true;
         };
         std::vector<std::size_t> cycle;

         while(-re.cost > tolerance) {
            double cycle_cap = std::numeric_limits<double>::infinity();
            auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const double cost) { cycle_cap = std::min(cycle_cap, cost); };
            std::vector<std::size_t> cycle = bfs.find_path(re[0], re[1], mask_small_edges, cycle_capacity); // possibly do not reallocate for every search but give as argument to find_path
            cycle_cap = std::min(cycle_cap, -re.cost);
            if(cycle.size() > 0 && cycle_cap >= tolerance) {

               if(record_cycles)
                  cp.add_cycle(cycle.begin(), cycle.end(), cycle_cap);

               // subtract minimum weight and remove edges of negligible weight
               assert(-re.cost >= cycle_cap && re.cost <= 0.0);
               re.cost += cycle_cap;
               assert(re.cost <= 0.0);
               for(std::size_t c=1; c<cycle.size(); ++c) {
                  pos_edges_graph.edge(cycle[c-1], cycle[c]) -= cycle_cap;
                  pos_edges_graph.edge(cycle[c], cycle[c-1]) -= cycle_cap;

                  assert(pos_edges_graph.edge(cycle[c-1], cycle[c]) >= 0.0);
                  assert(pos_edges_graph.edge(cycle[c], cycle[c-1]) >= 0.0);
                  assert(pos_edges_graph.edge(cycle[c-1], cycle[c]) == pos_edges_graph.edge(cycle[c], cycle[c-1]));
               }

               lower_bound += cycle_cap;

               // remove repulsive edge if weight is negligible
               if (-re.cost <= tolerance) {
                  auto imax = repulsive_edges.size() - 1;
                  repulsive_edges[i] = repulsive_edges[imax];
                  repulsive_edges.resize(imax);
                  i--;
                  break;
               }
            } else {
               break;
            }
         }
         }

         // terminate if no conflicted cycle remains
         if (repulsive_edges.empty())
            break;
      }

      const auto end_time = std::chrono::steady_clock::now();
      std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

      std::cout << "final lower bound = " << lower_bound << "\n";
      std::cout << "#repulsive edges = " << repulsive_edges.size() << "\n";

      return cp;
}

void multicut_cycle_packing(const multicut_instance& input)
{
   multicut_cycle_packing_impl(input, false);
}
cycle_packing compute_multicut_cycle_packing(const multicut_instance& input)
{
   return multicut_cycle_packing_impl(input, true);
}

// apply the cycle packing
triplet_multicut_instance pack_multicut_instance(const multicut_instance& input, const cycle_packing& cp)
{
    return pack_cut_base_instance<multicut_edge_triplet_message_0, multicut_edge_triplet_message_1, multicut_edge_triplet_message_2 , triplet_multicut_instance>(input, cp);
}

}
