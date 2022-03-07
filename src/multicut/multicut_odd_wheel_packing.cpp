#include "multicut/multicut_odd_wheel_packing.h"
#include "bipartite_graph_helper.hxx"
#include "union_find.hxx"
#include "cut_base/cut_base_apply_packing.hxx"
#include <vector>
#include <unordered_set>

namespace LPMP {

constexpr static double tolerance = 1e-8;

struct weighted_edge : public std::array<std::size_t,2> { double cost; };

std::array<std::size_t,2> incident_edges_indices(const std::size_t center_node_index)
{
   if(center_node_index == 0)
      return {0,1};
   else if(center_node_index == 1)
      return {0,2};
   else {
       assert(center_node_index == 2);
       return {1,2};
   }
} 

// compute difference between minimal costs such that exactly one edge incident to center node is cut against cost when when zero or two incident to it are cut
double compute_triangle_th(const std::size_t center_node_index, const multicut_triplet_factor& t)
{
   assert(center_node_index <= 2);

   size_t first_incident_edge, second_incident_edge;
   const auto edges = incident_edges_indices(center_node_index);
   std::tie(first_incident_edge, second_incident_edge) = std::make_tuple(edges[0], edges[1]);
   //const auto [first_incident_edge, second_incident_edge] = incident_edges_indices(center_node_index);
   double min_exactly_one_incident_cut = std::numeric_limits<double>::infinity();
   double min_other_cases = 0.0;
   auto update_costs = [first_incident_edge, second_incident_edge, &min_exactly_one_incident_cut, &min_other_cases](const std::bitset<3> labeling, const double cost) {
      if(std::size_t(labeling[first_incident_edge]) + std::size_t(labeling[second_incident_edge]) == 1)
         min_exactly_one_incident_cut = std::min(cost, min_exactly_one_incident_cut);
      else
         min_other_cases = std::min(cost, min_other_cases); 
   };
   t.for_each_labeling(update_costs);
   return - min_exactly_one_incident_cut + min_other_cases;
}

void reparametrize_triplet(multicut_triplet_factor& t, const std::size_t center_node_index, const double weight)
{
   assert(center_node_index <= 2);
   assert(weight >= 0.0);

   size_t first_incident_edge, second_incident_edge;
   const auto edges = incident_edges_indices(center_node_index);
   std::tie(first_incident_edge, second_incident_edge) = std::make_tuple(edges[0], edges[1]);
   //const auto [first_incident_edge, second_incident_edge] = incident_edges_indices(center_node_index);
   auto update_costs = [&](const std::bitset<3> labeling, double& cost) {
      if(std::size_t(labeling[first_incident_edge]) + std::size_t(labeling[second_incident_edge]) == 1)
         cost -= weight;
   };
   t.for_each_labeling(update_costs);
}

odd_wheel_packing multicut_odd_wheel_packing_impl(const triplet_multicut_instance& input, const bool record_odd_wheels)
{
   odd_wheel_packing owp;

   // prepare triangles
   std::vector<std::size_t> no_incident_triangles(input.no_nodes(), 0);
   for(const auto& t : input.triplets()) {
      no_incident_triangles[t[0]]++;
      no_incident_triangles[t[1]]++;
      no_incident_triangles[t[2]]++; 
   }
   std::vector<multicut_triplet_factor> triplets;
   triplets.reserve(input.triplets().size());
   struct odd_wheel_edge : public std::array<std::size_t,2> {
      multicut_triplet_factor* triplet; 
   };
   two_dim_variable_array<odd_wheel_edge> triangle_thresholds(no_incident_triangles.begin(), no_incident_triangles.end());
   std::fill(no_incident_triangles.begin(), no_incident_triangles.end(), 0);
   for(auto& t : input.triplets()) {
      triplets.push_back(t.cost);

      // TODO: possibly only include if multicut triplet factor has cost for given center node >= tolerance
      auto& e12 = triangle_thresholds( t[0], no_incident_triangles[t[0]]++ );
      e12.triplet = &triplets.back();
      e12[0] = t[1];
      e12[1] = t[2];

      auto& e02 = triangle_thresholds( t[1], no_incident_triangles[t[1]]++ );
      e02.triplet = &triplets.back();
      e02[0] = t[0];
      e02[1] = t[2];

      auto& e01 = triangle_thresholds( t[2], no_incident_triangles[t[2]]++ );
      e01.triplet = &triplets.back();
      e01[0] = t[0];
      e01[1] = t[1];
   }
   for(std::size_t i=0; i<triangle_thresholds.size(); ++i) {
      assert(no_incident_triangles[i] == triangle_thresholds[i].size());
      std::sort(triangle_thresholds[i].begin(), triangle_thresholds[i].end(), [](const auto& e1, const auto& e2) { return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end()); });
      for(std::size_t c=1; c<triangle_thresholds[i].size(); ++c) {
         assert(
               triangle_thresholds[i][c-1][0] != triangle_thresholds[i][c][0] || 
               triangle_thresholds[i][c-1][1] != triangle_thresholds[i][c][1]
               );
      }
   }

   double lower_bound = input.lower_bound();
   std::cout << "odd wheel packing\n";
   std::cout << "initial lower bound = " << lower_bound << "\n";
   std::cout << "#triplets = " << input.triplets().size() << "\n";

   struct triplet_edge {
       double cost;
       multicut_triplet_factor* triplet;
       std::size_t center_node_index;
   };
   compressed_bipartite_graph_helper<triplet_edge> bfs_helper(input.no_nodes());

   std::chrono::duration<double> graph_construction_time;
   std::chrono::duration<double> cycle_search_time;

   // TODO: iterate over vertices in random order
   for(std::size_t i=0; i<input.no_nodes(); ++i) {
       auto center_node_index_func = [&](const odd_wheel_edge& e) {
           if(i < e[0]) return 0;
           if(i < e[1]) return 1;
           else return 2;
       };

       auto start = std::chrono::system_clock::now();
       bfs_helper.construct_compressed_bipartite_graph(triangle_thresholds[i].begin(), triangle_thresholds[i].end(), 
               [&](const odd_wheel_edge& e) -> std::array<std::size_t,2> { 
               return {e[0],e[1]};
               },
               [&](const odd_wheel_edge& e) { 
               const std::size_t center_node_index = center_node_index_func(e);
               const double cost = compute_triangle_th(center_node_index, *e.triplet);
               return triplet_edge{cost, e.triplet, center_node_index};
               }
               );
       graph_construction_time += std::chrono::system_clock::now()- start;

       start = std::chrono::system_clock::now();
      // TODO: regularly recompute union find?
      // TODO: interchange cycle length and iterating over all nodes
      const std::array<std::size_t,6> cycle_lengths = {3,4,5,6,7,std::numeric_limits<std::size_t>::max()};
      for(const std::size_t cycle_length : cycle_lengths) {
         for(std::size_t ci=0; ci<bfs_helper.no_compressed_nodes(); ++ci) {
            if(ci + bfs_helper.no_compressed_nodes() < bfs_helper.get_graph().no_nodes() && (true)) {// || uf.connected(ci, ci+no_compressed_nodes))) { // TODO: activate uf again

               // balance short cycles as long as available and positive weight left
               auto mask_small_edges = [cycle_length](const std::size_t i, const std::size_t j, const triplet_edge& e, const std::size_t distance) { 
                  if(e.cost <= tolerance) return false;
                  if(distance >= cycle_length) return false;
                  return true;
               };

               double cycle_cap = std::numeric_limits<double>::infinity();
               auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const triplet_edge& e) { cycle_cap = std::min(cycle_cap, e.cost); };

               // TODO: find multiple cycles connecting current node
               std::vector<std::size_t> cycle = bfs_helper.get_bfs().find_path(ci, ci+bfs_helper.no_compressed_nodes(), mask_small_edges, cycle_capacity); 
               if(cycle.size() == 0)
                   continue;
               assert(cycle_cap >= -1e8);
               //std::cout << "found cycle for node " << i << " with weight " << cycle_cap << "\n";
               assert(cycle.size() == 0 || cycle[0] % bfs_helper.no_compressed_nodes() == cycle.back() % bfs_helper.no_compressed_nodes());

               for(std::size_t c=1; c<cycle.size(); ++c) {
                   assert(bfs_helper.get_graph().edge_present(cycle[c-1], cycle[c]));
               }

               std::transform(cycle.begin(), cycle.end(), cycle.begin(), [&](const std::size_t val) { return val % bfs_helper.no_compressed_nodes(); });
               cycle = find_subcycle(cycle);
               assert(!has_subcycles(cycle));

               lower_bound += cycle_cap; 

               // substract weights from triplets and update underlying triplets
               for(std::size_t c=1; c<cycle.size(); ++c) {
                  const std::size_t ci = cycle[c-1] % bfs_helper.no_compressed_nodes();
                  const std::size_t cj = cycle[c] % bfs_helper.no_compressed_nodes();
                  assert(ci < bfs_helper.no_compressed_nodes() && cj < bfs_helper.no_compressed_nodes());

                  assert(bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes()).cost >= 1e-8);
                  assert(bfs_helper.get_graph().edge(cj+bfs_helper.no_compressed_nodes(), ci).cost >= 1e-8);
                  assert(bfs_helper.get_graph().edge(cj, ci+bfs_helper.no_compressed_nodes()).cost >= 1e-8);
                  assert(bfs_helper.get_graph().edge(ci+bfs_helper.no_compressed_nodes(), cj).cost >= 1e-8);

                  bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes()).cost -= cycle_cap;
                  bfs_helper.get_graph().edge(cj+bfs_helper.no_compressed_nodes(), ci).cost -= cycle_cap;
                  bfs_helper.get_graph().edge(cj, ci+bfs_helper.no_compressed_nodes()).cost -= cycle_cap;
                  bfs_helper.get_graph().edge(ci+bfs_helper.no_compressed_nodes(), cj).cost -= cycle_cap;

                  assert(bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes()).cost >= -1e-8);
                  assert(bfs_helper.get_graph().edge(cj+bfs_helper.no_compressed_nodes(), ci).cost >= -1e-8);
                  assert(bfs_helper.get_graph().edge(cj, ci+bfs_helper.no_compressed_nodes()).cost >= -1e-8);
                  assert(bfs_helper.get_graph().edge(ci+bfs_helper.no_compressed_nodes(), cj).cost >= -1e-8);

                  auto& e = bfs_helper.get_graph().edge(ci, cj+bfs_helper.no_compressed_nodes());
                  reparametrize_triplet( *e.triplet, e.center_node_index, cycle_cap);
               }
               // transform back to original nodes
               cycle.resize(cycle.size()-1); 
               bfs_helper.compressed_path_to_original(cycle);
               if(record_odd_wheels)
                  owp.add_odd_wheel(i, cycle.begin(), cycle.end(), cycle_cap);
            }
         }
      }
      cycle_search_time += std::chrono::system_clock::now()- start;
   }

      std::cout << "time for compressed graph construction = " << graph_construction_time.count() << "\n";
      std::cout << "time for cycle search = " << cycle_search_time.count() << "\n";


   std::cout << "final lower bound = " << lower_bound << "\n";
   return owp;
}

void multicut_odd_wheel_packing(const triplet_multicut_instance& input)
{
   multicut_odd_wheel_packing_impl(input, false); 
}
odd_wheel_packing compute_multicut_odd_wheel_packing(const triplet_multicut_instance& input)
{
   return multicut_odd_wheel_packing_impl(input, true); 
}

quadruplet_multicut_instance pack_multicut_instance(const triplet_multicut_instance& input, const odd_wheel_packing& owp)
{
    return pack_triplet_cut_base_instance<triplet_multicut_instance, quadruplet_multicut_instance, multicut_triplet_quadruplet_message_012, multicut_triplet_quadruplet_message_013, multicut_triplet_quadruplet_message_023, multicut_triplet_quadruplet_message_123>(input, owp);
}


}
