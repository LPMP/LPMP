#pragma once

#include "multicut_triplet_constructor.hxx"
#include "multicut_quadruplet_constructor.hxx"
//#include "multicut_quintuplet_constructor.hxx"

//#include "cut_base/cut_base_constructor.hxx"
//#include "multicut_factors_messages.hxx"
//#include "multicut_instance.h"
//#include "cut_base/cut_base_packing.h"
//#include "multicut_cycle_packing.h"
//#include "multicut_odd_wheel_packing.h"
//#include "union_find.hxx"
//#include "graph.hxx"
//#include <omp.h>
//
//#include "andres/graph/graph.hxx"
//#include "andres/graph/grid-graph.hxx"
//#include "andres/graph/multicut/kernighan-lin.hxx"
//#include "andres/graph/multicut/greedy-additive.hxx"
//#include "andres/graph/multicut-lifted/kernighan-lin.hxx"
//#include "andres/graph/multicut-lifted/greedy-additive.hxx"
//
//namespace LPMP {
//
//template<
//   class FACTOR_MESSAGE_CONNECTION, 
//   typename UNARY_FACTOR, typename TRIPLET_FACTOR,
//   typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2
//   >
//class multicut_constructor : public cut_constructor_base<
//			     multicut_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0,UNARY_TRIPLET_MESSAGE_1,  UNARY_TRIPLET_MESSAGE_2 >, 
//			     FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2> 
//{
//public:
//   using type = multicut_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;
//   using FMC = FACTOR_MESSAGE_CONNECTION;
//   using base_constructor = cut_constructor_base<
//	   multicut_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0,UNARY_TRIPLET_MESSAGE_1,  UNARY_TRIPLET_MESSAGE_2 >, 
//	   FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;
//
//   using base_constructor::base_constructor;
//
//   using weighted_graph = graph<double>;
//
//   // add triplet factor with prespecified cost
//   // edges on triplet are (uv,uw,vw) (lexicographical)
//   typename base_constructor::triplet_factor* add_higher_order_triplet(
//         const INDEX u, const INDEX v, const INDEX w, 
//         const REAL c000, const REAL c011, const REAL c101, const REAL c110, const REAL c111)
//   {
//      auto* f = this->AddTripletFactor(u,v,w);
//      auto* t = f->get_factor();
//      this->lp_->add_to_constant(c000);
//      (*t)[0] = c110 - c000;
//      (*t)[1] = c101 - c000;
//      (*t)[2] = c011 - c000;
//      (*t)[3] = c111 - c000;
//      return f;
//   }
//
//   static REAL triplet_cost(const REAL cost_ij, const REAL cost_ik, const REAL cost_jk)
//   {
//      return std::min({0.0, cost_ij+cost_ik, cost_ij+cost_jk, cost_ik+cost_jk, cost_ij+cost_ik+cost_jk}); 
//   }
//
//   INDEX find_violated_cycles(const INDEX max_triplets_to_add)
//   {
//      std::vector<std::tuple<std::size_t, std::size_t, double, bool> > negative_edges; // endpoints, edge cost, searched positive path with given endpoints?
//
//      struct weighted_edge : public std::array<std::size_t,2> { double cost; double& operator()() { return cost; } };
//      std::vector<weighted_edge> positive_edges;
//
//      // we can speed up compution by skipping path searches for node pairs which lie in different connected components. Connectedness is stored in a union find structure
//      double pos_th = 0.0;
//      for(auto& it : this->unaryFactorsVector_) {
//         const double v = (*it.second->get_factor())[0];
//         const std::size_t i = std::get<0>(it.first);
//         const std::size_t j = std::get<1>(it.first);
//         if(v >= 0.0) {
//            pos_th = std::max(pos_th, v);
//            positive_edges.push_back({i,j,v});
//         } else {
//            negative_edges.push_back(std::make_tuple(i,j,v,false));
//         }
//      }
//
//      if(negative_edges.size() == 0 || negative_edges.size() == this->unaryFactorsVector_.size()) { return 0; }
//
//      weighted_graph posEdgesGraph(positive_edges.begin(), positive_edges.end()); // graph consisting of positive edges
//
//      // do zrobienia: possibly add reparametrization of triplet factors additionally. Do it analogously as in graph matching via send_messages_to_unaries()?
//
//      std::sort(negative_edges.begin(), negative_edges.end(), [](const auto& e1, const auto& e2)->bool { return std::get<2>(e1) < std::get<2>(e2); });
//
//
//         // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
//         // here, longest path is sought after only the edges with positive taken into account
//         // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min
//
//         //MostViolatedPathData mp(posEdgesGraph);
//         //BfsData mp(posEdgesGraph);
//
//         union_find uf(this->noNodes_);
//         INDEX triplets_added = 0;
//         const REAL initial_th = 0.6*std::min(-std::get<2>(negative_edges[0]), pos_th);
//         bool zero_th_iteration = true;
//
//         for(REAL th=initial_th; th>=eps || zero_th_iteration; th*=0.1) {
//            if(th < eps) {
//               if(triplets_added <= 0.01*max_triplets_to_add) {
//                  // we would first like to go on with odd wheels and odd bicycle wheels before actually adding triplets with no guaranteed dual increase.
//                  //std::cout << "additional separation with no guaranteed dual increase, i.e. th = 0\n";
//                  //th = 0.0;
//                  //zero_th_iteration = false;
//               } else {
//                  break;
//               }
//            }
//            // first update union find datastructure
//            for(auto& it : this->unaryFactorsVector_) {
//               const double v = (*it.second->get_factor())[0];
//               if(v >= th) {
//                  const std::size_t i = std::get<0>(it.first);
//                  const std::size_t j = std::get<1>(it.first);
//                  uf.merge(i,j);   
//               }
//            }
//
//            std::vector<typename base_constructor::triplet_candidate> triplet_candidates;
//
//#pragma omp parallel 
//            {
//               std::vector<typename base_constructor::triplet_candidate> triplet_candidates_local;
//               bfs_data<weighted_graph> mp2(posEdgesGraph);
//#pragma omp for schedule(guided) nowait
//               for(std::size_t c=0; c<negative_edges.size(); ++c) {
//                  const std::size_t i = std::get<0>(negative_edges[c]);
//                  const std::size_t j = std::get<1>(negative_edges[c]);
//                  const double v = std::get<2>(negative_edges[c]);
//                  const bool already_used_for_path_search = std::get<3>(negative_edges[c]);
//                  //if(-v <= th) break;
//                  //if(already_used_for_path_search) continue;
//                  if(-v > th && !already_used_for_path_search && uf.thread_safe_connected(i,j)) {
//                     //auto cycle = mp.FindPath(i,j,posEdgesGraph);
//                     auto mask_small_edges = [th](const std::size_t i, const std::size_t j, const double cost, const std::size_t distance) { return cost >= th; };
//                     double cycle_cap = std::numeric_limits<double>::infinity();
//                     auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const double cost) { cycle_cap = std::min(cycle_cap, cost); };
//                     auto cycle = mp2.find_path(i, j, mask_small_edges, cycle_capacity);
//                     const REAL dualIncrease = std::min(-v, cycle_cap);
//                     assert(cycle.size() > 0);
//                     if(cycle.size() > 0) {
//                        this->triangulate_cycle(dualIncrease, cycle.begin(), cycle.end(), triplet_candidates_local);
//                        //triplets_added += AddCycle(std::get<1>(cycle));
//                        //if(triplets_added > max_triplets_to_add) {
//                        //   return triplets_added;
//                        //}
//                     } else {
//                        throw std::runtime_error("No path found although there should be one"); 
//                     }
//                  }
//               }
//#pragma omp critical
//               {
//                  triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_local.begin(), triplet_candidates_local.end()); 
//            }
//         }
//
//         // sort by guaranteed increase in decreasing order
//         std::sort(triplet_candidates.begin(), triplet_candidates.end());
//
//         if(triplet_candidates.size() > 0 && diagnostics()) {
//            std::cout << "best triplet candidate in triplet search has guaranteed dual improvement " << triplet_candidates[0].cost << "\n";
//         }
//
//         for(const auto& triplet_candidate : triplet_candidates) {
//            const INDEX i = triplet_candidate.nodes[0];
//            const INDEX j = triplet_candidate.nodes[1];
//            const INDEX k = triplet_candidate.nodes[2];
//            if(!this->HasTripletFactor(i,j,k)) {
//               this->AddTripletFactor(i,j,k);
//               triplets_added++;
//               if(triplets_added > max_triplets_to_add) {
//                  break;
//               } 
//            }
//         }
//
//         triplet_candidates.clear();
//      }
//
//      return triplets_added;
//   }
//
//   void construct(const multicut_instance& mc)
//   {
//      for(const auto& e : mc.edges()) {
//         this->AddUnaryFactor(e[0], e[1], e.cost);
//      }
//      this->no_original_edges_ = this->unaryFactorsVector_.size();
//
//      // compute cycle packing and add returned cycles to problem formulation. Also reparametrize edges (possibly do not do this?)
//      cycle_packing cp = compute_multicut_cycle_packing(mc);
//      for(std::size_t c=0; c<cp.no_cycles(); ++c) {
//         const auto [cycle_begin, cycle_end] = cp.get_cycle(c);
//         //this->add_cycle(cycle_begin, cycle_end);
//         this->add_cycle(cycle_begin, cycle_end, cp.get_cycle_weight(c));
//      }
//   }
//
//   bool CheckPrimalConsistency() const
//   {
//     if(debug()) {
//       std::cout << "checking primal feasibility for multicut ... ";
//     }
//     union_find uf(this->noNodes_);
//     for(const auto& e : this->unaryFactorsVector_) {
//       auto* f = e.second; 
//       if(f->get_factor()->primal()[0] == false) {
//         // connect components 
//         const INDEX i = std::get<0>(e.first);
//         const INDEX j = std::get<1>(e.first);
//         uf.merge(i,j);
//       }
//     }
//     for(const auto& e : this->unaryFactorsVector_) {
//       auto* f = e.second; 
//       if(f->get_factor()->primal()[0] == true) {
//         const INDEX i = std::get<0>(e.first);
//         const INDEX j = std::get<1>(e.first);
//         // there must not be a path from i1 to i2 consisting of edges with primal value false only
//         if(uf.connected(i,j)) {
//           if(debug()) {
//             std::cout << "solution infeasible: (" << i << "," << j << ") = true, yet there exists a path with false values only\n";
//           }
//           return false;
//         }
//       }
//     }
//
//     if(debug()) {
//       std::cout << "solution feasible\n";
//     }
//     return true; 
//   }
//
//   static std::vector<char> round(std::vector<typename base_constructor::edge> edges)
//   {
//     INDEX no_nodes = 0;
//     for(auto e : edges) {
//       no_nodes = std::max({no_nodes, e.i+1, e.j+1});
//     }
//      andres::graph::Graph<> graph(no_nodes);
//      std::vector<double> edge_values;
//      edge_values.reserve(edges.size());
//      for(const auto e : edges) {
//         graph.insertEdge(e.i, e.j); 
//        edge_values.push_back(e.cost); 
//      }
//
//      std::vector<char> labeling(edges.size(), 0);
//      if(graph.numberOfEdges() > 0) {
//         andres::graph::multicut::greedyAdditiveEdgeContraction(graph, edge_values, labeling);
//         andres::graph::multicut::kernighanLin(graph, edge_values, labeling, labeling);
//      }
//      return labeling; 
//   }
//};
//
//
//
//template<
//   class BASE_CONSTRUCTOR,
//   typename ODD_3_WHEEL_FACTOR,
//   typename TRIPLET_ODD_3_WHEEL_MESSAGE_012, typename TRIPLET_ODD_3_WHEEL_MESSAGE_013, typename TRIPLET_ODD_3_WHEEL_MESSAGE_023, typename TRIPLET_ODD_3_WHEEL_MESSAGE_123
//>
//class multicut_odd_wheel_constructor : public odd_wheel_constructor_base<
//				       multicut_odd_wheel_constructor<BASE_CONSTRUCTOR, ODD_3_WHEEL_FACTOR, TRIPLET_ODD_3_WHEEL_MESSAGE_012, TRIPLET_ODD_3_WHEEL_MESSAGE_013, TRIPLET_ODD_3_WHEEL_MESSAGE_023, TRIPLET_ODD_3_WHEEL_MESSAGE_123>,
//				       BASE_CONSTRUCTOR, ODD_3_WHEEL_FACTOR, TRIPLET_ODD_3_WHEEL_MESSAGE_012, TRIPLET_ODD_3_WHEEL_MESSAGE_013, TRIPLET_ODD_3_WHEEL_MESSAGE_023, TRIPLET_ODD_3_WHEEL_MESSAGE_123 > {
//public:
//   using FMC = typename BASE_CONSTRUCTOR::FMC;
//
//   using base_constructor = odd_wheel_constructor_base<
//				       multicut_odd_wheel_constructor<BASE_CONSTRUCTOR, ODD_3_WHEEL_FACTOR, TRIPLET_ODD_3_WHEEL_MESSAGE_012, TRIPLET_ODD_3_WHEEL_MESSAGE_013, TRIPLET_ODD_3_WHEEL_MESSAGE_023, TRIPLET_ODD_3_WHEEL_MESSAGE_123>,
//				       BASE_CONSTRUCTOR, ODD_3_WHEEL_FACTOR, TRIPLET_ODD_3_WHEEL_MESSAGE_012, TRIPLET_ODD_3_WHEEL_MESSAGE_013, TRIPLET_ODD_3_WHEEL_MESSAGE_023, TRIPLET_ODD_3_WHEEL_MESSAGE_123 >;
//
//   using base_constructor::base_constructor;
//
//   /*
//   quadruplet labeling order:
//   1,1,0,1,0,0
//   1,0,1,0,1,0
//   0,1,1,0,0,1
//   0,0,0,1,1,1
//   0,1,1,1,1,0
//   1,0,1,1,0,1
//   1,1,0,0,1,1
//   0,1,1,1,1,1
//   1,0,1,1,1,1
//   1,1,0,1,1,1
//   1,1,1,0,1,1
//   1,1,1,1,0,1
//   1,1,1,1,1,0
//   1,1,1,1,1,1
//   */ 
//   template<typename ITERATOR>
//   typename base_constructor::odd_3_wheel_factor_container* add_higher_order_quadruplet(
//         const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3,
//         ITERATOR cost_begin, ITERATOR cost_end)
//   {
//      auto* f = this->add_odd_3_wheel_factor(i0,i1,i2,i3);
//      assert(std::distance(cost_begin, cost_end) == 15);
//      this->lp_->add_to_constant(*cost_begin);
//      auto* t = f->get_factor();
//
//      std::size_t i = 0;
//      for(auto it = cost_begin+1; it != cost_end; ++it, ++i) {
//         (*t)[i] = *it - *cost_begin; 
//      }
//
//      return f;
//   }
//
//
//   // node i is center node, (j,k) is cycle edge
//   REAL ComputeTriangleTh(const INDEX i, const INDEX j, const INDEX k)
//   {
//      std::array<INDEX,3> triplet{i,j,k};
//      std::sort(triplet.begin(),triplet.end());
//      assert(this->HasUnaryFactor(triplet[0],triplet[1]) && this->HasUnaryFactor(triplet[0],triplet[2]) && this->HasUnaryFactor(triplet[1],triplet[2]));
//      std::array<REAL,4> cost;
//      std::fill(cost.begin(),cost.end(),0.0);
//      if(this->HasTripletFactor(triplet[0],triplet[1],triplet[2])) {
//         auto* t = this->GetTripletFactor(triplet[0],triplet[1],triplet[2])->get_factor();
//         assert(t->size() == 4);
//         cost[0] += (*t)[0];
//         cost[1] += (*t)[1];
//         cost[2] += (*t)[2];
//         cost[3] += (*t)[3];
//      } 
//      // get cost directly from edge factors
//      const REAL c01 = (*this->GetUnaryFactor(triplet[0], triplet[1])->get_factor())[0];
//      const REAL c02 = (*this->GetUnaryFactor(triplet[0], triplet[2])->get_factor())[0];
//      const REAL c12 = (*this->GetUnaryFactor(triplet[1], triplet[2])->get_factor())[0];
//      cost[0] += c01 + c02;
//      cost[1] += c01 + c12;
//      cost[2] += c02 + c12;
//      cost[3] += c01 + c02 + c12;
//
//      assert(j<k); // if not, below computation is not valid
//      // compute difference between cost such that exactly one edge incident to center node is 1 against cost when when zero or two incident to it are 1
//      if(std::min(j,k) == triplet[0] && std::max(j,k) == triplet[1]) { // jk is first edge
//         return std::min({0.0,cost[3],cost[2]}) - std::min(cost[0],cost[1]);
//      } else if(std::min(j,k) == triplet[0] && std::max(j,k) == triplet[2]) { // jk is second edge
//         return std::min({0.0,cost[3],cost[1]}) - std::min(cost[0],cost[2]);
//      } else { // jk is third edge
//         assert(std::min(j,k) == triplet[1] && std::max(j,k) == triplet[2]);
//         return std::min({0.0,cost[3],cost[0]}) - std::min(cost[1],cost[2]);
//      }
//      assert(false);
//   }
//
//   template<typename ADJ_LIST>
//   void ComputeTriangles(const INDEX i, const ADJ_LIST& adjacencyList, const REAL minTh, 
//         std::unordered_map<INDEX,INDEX>& origToCompressedNode, 
//         std::vector<INDEX>& compressedToOrigNode, 
//         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges)
//   {
//      //std::unordered_map<INDEX,INDEX> origToCompressedNode {tripletByIndices_[i].size()}; // compresses node indices
//      //origToCompressedNode.max_load_factor(0.7);
//      //std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
//      //std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;
//      origToCompressedNode.clear();
//      compressedToOrigNode.clear();
//      compressedEdges.clear();
//
//      // find all triangles ijk
//      std::vector<INDEX> commonNodes(this->noNodes_); // for detecting triangles
//      for(INDEX j : adjacencyList[i]) {
//         auto commonNodesEnd = std::set_intersection(adjacencyList[i].begin(), adjacencyList[i].end(), adjacencyList[j].begin(), adjacencyList[j].end(), commonNodes.begin());
//         for(auto it=commonNodes.begin(); it!=commonNodesEnd; ++it) {
//            const INDEX k = *it; 
//            if(j<k) { // edge is encountered twice
//               const REAL th = ComputeTriangleTh(i,j,k);
//               if(th >= minTh) {
//                  // add j and k to compressed nodes
//                  if(origToCompressedNode.find(j) == origToCompressedNode.end()) {
//                     origToCompressedNode.insert(std::make_pair(j, origToCompressedNode.size()));
//                     compressedToOrigNode.push_back(j);
//                  }
//                  if(origToCompressedNode.find(k) == origToCompressedNode.end()) {
//                     origToCompressedNode.insert(std::make_pair(k, origToCompressedNode.size()));
//                     compressedToOrigNode.push_back(k);
//                  }
//                  const INDEX jc = origToCompressedNode[j];
//                  const INDEX kc = origToCompressedNode[k];
//                  assert(jc != kc);
//                  compressedEdges.push_back(std::make_tuple(jc,kc,th));
//               }
//            }
//         }
//      }
//   }
//
//   // what is the lowest threshold so that an odd cycle exists
//   REAL compute_odd_cycle_threshold(
//         std::unordered_map<INDEX,INDEX>& origToCompressedNode,
//         std::vector<INDEX>& compressedToOrigNode,
//         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges 
//         )
//   {
//      std::sort(compressedEdges.begin(), compressedEdges.end(), [](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });
//
//      const INDEX noCompressedNodes = origToCompressedNode.size();
//      const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
//      union_find uf(noBipartiteCompressedNodes); 
//      // construct bipartite graph based on triangles
//      for(auto& e : compressedEdges) {
//         const INDEX jc = std::get<0>(e);
//         const INDEX kc = std::get<1>(e);
//         uf.merge(jc,noCompressedNodes + kc);
//         uf.merge(noCompressedNodes + jc,kc);
//         if(uf.connected(jc, noCompressedNodes + jc)) {
//            assert(uf.connected(kc, noCompressedNodes + kc));
//            return std::get<2>(e);
//         }
//      }
//      return -std::numeric_limits<REAL>::infinity(); // no constraint found
//   }
//
//   // returns nodes of odd wheel without center node
//   std::vector<INDEX> compute_path_in_bipartite_graph(
//         // original compressed nodes and edges. Double them and search.
//         std::unordered_map<INDEX,INDEX>& origToCompressedNode,
//         std::vector<INDEX>& compressedToOrigNode,
//         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges,
//         const REAL th
//         )
//   {
//      const INDEX noCompressedNodes = origToCompressedNode.size();
//      const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
//      struct weighted_edge : public std::array<std::size_t,2> { double cost; double& operator()() { return cost; } }; // TODO: make one declaration
//      std::vector<weighted_edge> positive_edges;
//      union_find uf(noBipartiteCompressedNodes);
//      for(const auto e : compressedEdges) {
//         const std::size_t i = std::get<0>(e);
//         const std::size_t j = std::get<1>(e);
//         const double cost = std::get<2>(e);
//         if(cost >= th) {
//            positive_edges.push_back({i, noCompressedNodes + j, cost});
//            positive_edges.push_back({j, noCompressedNodes + i, cost});
//            uf.merge(i,noCompressedNodes + j);
//            uf.merge(j, noCompressedNodes + i);
//         }
//      }
//      graph<double> g(positive_edges.begin(), positive_edges.end());
//      bfs_data<graph<double>> mp(g);
//      // now check whether path exists between any given edges on graph
//      for(std::size_t j=0; j<noCompressedNodes; ++j) { // not nice: this has to be original number of nodes and bipartiteNumberOfNodes
//         // find path from node j to node noNodes+j in g
//         if(uf.connected(j,noCompressedNodes+j)) {
//            auto path = mp.find_path(j,noCompressedNodes+j);
//            assert(path.size() > 3);
//            path.resize(path.size()-1); // first and last node coincide
//            for(INDEX k=0; k<path.size(); ++k) { // note: last node is copy of first one
//               //assert(compressedToOrigNode.find(path[k]%noCompressedNodes) != compressedToOrigNode.end());
//               path[k] = compressedToOrigNode[path[k]%noCompressedNodes];
//            }
//
//            if(HasUniqueValues(path)) { // possibly already add the subpath that is unique and do not search for it later. Indicate this with a std::vector<bool>
//               //assert(HasUniqueValues(path)); // if not, a shorter subpath has been found. This subpath will be detected or has been deteced and has been added
//               this->cycle_normal_form(path.begin(), path.end());
//               //cycle_normal_form called unnecesarily in EnforceOddWheel
//               return path;
//            } 
//         }
//      }
//      return std::vector<INDEX>(0);
//   }
//
//   struct odd_3_wheel_candidate {
//      std::array<INDEX,4> nodes;
//      REAL cost;
//   };
//
//   template<typename ITERATOR>
//   void triangulate_odd_wheel(const INDEX i, const REAL cost, ITERATOR path_begin, ITERATOR path_end, std::vector<odd_3_wheel_candidate>& candidates)
//   {
//      assert(std::distance(path_begin, path_end) >= 3);
//      this->cycle_normal_form(path_begin, path_end);
//      const INDEX first_node = *path_begin;
//      for(auto it=path_begin+1; it+1!=path_end; ++it) {
//         std::array<INDEX,4> nodes({i,first_node, *it, *(it+1)});
//         std::sort(nodes.begin(), nodes.end());
//         assert(HasUniqueValues(nodes));
//         candidates.push_back({nodes, cost});
//      } 
//   }
//
//   INDEX find_odd_wheels(const INDEX max_factors_to_add)
//   {
//      // first prepare datastructures for threshold finding and violated constraint search
//
//      // search for all triangles present in the graph, also in places where no triplet factor has been added
//      // Construct adjacency list
//
//      std::vector<INDEX> adjacency_list_count(this->noNodes_,0);
//      // first determine size for adjacency_list
//      // to do: possibly parallelize
//      for(auto& it : this->unaryFactorsVector_) {
//         const INDEX i = std::get<0>(it.first);
//         const INDEX j = std::get<1>(it.first);
//         adjacency_list_count[i]++;
//         adjacency_list_count[j]++; 
//      }
//      two_dim_variable_array<INDEX> adjacency_list(adjacency_list_count);
//      std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
//      for(auto& it : this->unaryFactorsVector_) {
//         const INDEX i = std::get<0>(it.first);
//         const INDEX j = std::get<1>(it.first);
//         assert(i<j);
//         adjacency_list[i][adjacency_list_count[i]] = j;
//         adjacency_list_count[i]++;
//         adjacency_list[j][adjacency_list_count[j]] = i;
//         adjacency_list_count[j]++;
//      }
//
//      // Sort the adjacency list, for fast intersections later
//#pragma omp parallel for schedule(guided)
//      for(int i=0; i < adjacency_list.size(); i++) {
//         std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
//      }
//
//      // find maximum threshold where still some cycle can be added for each node
//      std::vector<odd_3_wheel_candidate> odd_3_wheel_candidates;
//
//#pragma omp parallel 
//      {
//         std::unordered_map<INDEX,INDEX> origToCompressedNode; // compresses node indices
//         std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
//         std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;
//         std::vector<odd_3_wheel_candidate> odd_3_wheel_candidates_local;
//#pragma omp for schedule(guided) nowait 
//         for(INDEX i=0; i<this->noNodes_; ++i) {
//            ComputeTriangles(i, adjacency_list, eps, origToCompressedNode, compressedToOrigNode, compressedEdges);
//            const REAL th = compute_odd_cycle_threshold(origToCompressedNode, compressedToOrigNode, compressedEdges);
//            if(th >= eps) {
//               auto oddWheel = compute_path_in_bipartite_graph(origToCompressedNode, compressedToOrigNode, compressedEdges, th);
//               assert(oddWheel.size() > 0);
//               triangulate_odd_wheel(i, th, oddWheel.begin(), oddWheel.end(), odd_3_wheel_candidates_local);
//            }
//         }
//#pragma omp critical
//         {
//            odd_3_wheel_candidates.insert(odd_3_wheel_candidates.end(), odd_3_wheel_candidates_local.begin(), odd_3_wheel_candidates_local.end()); 
//         }
//      }
//      std::sort(odd_3_wheel_candidates.begin(), odd_3_wheel_candidates.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });
//      INDEX factors_added = 0;
//      for(const auto& odd_3_wheel : odd_3_wheel_candidates) {
//         if(!this->has_odd_3_wheel_factor(odd_3_wheel.nodes[0], odd_3_wheel.nodes[1], odd_3_wheel.nodes[2], odd_3_wheel.nodes[3])) {
//            this->add_odd_3_wheel_factor(odd_3_wheel.nodes[0], odd_3_wheel.nodes[1], odd_3_wheel.nodes[2], odd_3_wheel.nodes[3]);
//            ++factors_added;
//            if(factors_added > max_factors_to_add) {
//               break;
//            }
//         } 
//      }
//      return factors_added;
//   }
//};
//
//template<typename ODD_WHEEL_CONSTRUCTOR, 
//   typename ODD_BICYCLE_3_WHEEL_FACTOR,
//   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123,
//   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124,
//   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134,
//   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234,
//   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234
//   >
//class multicut_odd_bicycle_wheel_constructor : public odd_bicycle_wheel_constructor_base<
//					       multicut_odd_bicycle_wheel_constructor< ODD_WHEEL_CONSTRUCTOR, ODD_BICYCLE_3_WHEEL_FACTOR, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234>,
//					       ODD_WHEEL_CONSTRUCTOR, ODD_BICYCLE_3_WHEEL_FACTOR, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234>
//{
//public:
//   using FMC = typename ODD_WHEEL_CONSTRUCTOR::FMC;
//   using base_constructor = odd_bicycle_wheel_constructor_base<
//	   multicut_odd_bicycle_wheel_constructor< ODD_WHEEL_CONSTRUCTOR, ODD_BICYCLE_3_WHEEL_FACTOR, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234>,
//	   ODD_WHEEL_CONSTRUCTOR, ODD_BICYCLE_3_WHEEL_FACTOR, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234, ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234>;
//
//   using odd_bicycle_3_wheel_factor_container = ODD_BICYCLE_3_WHEEL_FACTOR;
//   using odd_3_wheel_odd_bicycle_3_wheel_message_0123_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123;
//   using odd_3_wheel_odd_bicycle_3_wheel_message_0124_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124;
//   using odd_3_wheel_odd_bicycle_3_wheel_message_0134_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134;
//   using odd_3_wheel_odd_bicycle_3_wheel_message_0234_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234;
//   using odd_3_wheel_odd_bicycle_3_wheel_message_1234_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234;
//
//   template<typename SOLVER>
//   multicut_odd_bicycle_wheel_constructor(SOLVER& s) : base_constructor(s) {}
//
//
//   template<typename ITERATOR>
//   auto* add_higher_order_quintuplet(
//         const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3, const INDEX i4,
//         ITERATOR cost_begin, ITERATOR cost_end)
//   {
//      auto* f = this->add_odd_bicycle_3_wheel(i0,i1,i2,i3);
//      auto* t = f->get_factor();
//      assert(std::distance(cost_begin, cost_end) + 1 == t->size());
//      this->add_to_constant(*cost_begin);
//
//      std::size_t i = 0;
//      for(auto it = cost_begin+1; it != cost_end; ++it, ++i) {
//         (*t)[i] = *it - *cost_begin; 
//      }
//
//      return f;
//   }
//
//   // ij is the axle, uv is the wheel edge
//   // labelings with edge ij and uv cut and (i) iu and jv or (ii) iv and ju not cut must have smaller cost than other configurations
//   REAL compute_edge_cost(const INDEX i, const INDEX j, const INDEX u, const INDEX v)
//   {
//      std::array<INDEX,4> idx{i,j,u,v};
//      std::sort(idx.begin(), idx.end());
//      // can it happen that odd bicycle wheel inequalities can tighten the polytope but odd wheel inequalities are not also violated? Only in this case the below construction makes sense
//      if(!this->has_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])) { // create a fake odd 3 wheel factor and reparamaetrize all underlying triplets into it
//
//	      multicut_odd_3_wheel_factor f;
//
//	      // possibly do not create new messages, but use static functions in these messages directly
//	      if(this->HasTripletFactor(idx[0],idx[1],idx[2])) {
//		      const auto& t = *(this->GetTripletFactor(idx[0], idx[1], idx[2])->get_factor());
//		      multicut_triplet_odd_3_wheel_message_012 m;
//		      m.RepamRight(f,t);
//	      }
//
//	      if(this->HasTripletFactor(idx[0],idx[1],idx[3])) {
//		      const auto& t = *(this->GetTripletFactor(idx[0], idx[1], idx[3])->get_factor());
//		      multicut_triplet_odd_3_wheel_message_013 m;
//		      m.RepamRight(f,t);
//	      }
//
//	      if(this->HasTripletFactor(idx[0],idx[2],idx[3])) {
//		      const auto& t = *(this->GetTripletFactor(idx[0], idx[2], idx[3])->get_factor());
//		      multicut_triplet_odd_3_wheel_message_023 m;
//		      m.RepamRight(f,t);
//	      }
//
//            if(this->HasTripletFactor(idx[1],idx[2],idx[3])) {
//               const auto& t = *(this->GetTripletFactor(idx[1], idx[2], idx[3])->get_factor());
//               multicut_triplet_odd_3_wheel_message_123 m;
//               m.RepamRight(f,t);
//            }
//            return this->compute_edge_cost_from_odd_3_wheel(i,j,u,v, f);
//
//      } else {
//         auto& f = *(this->get_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])->get_factor());
//         return this->compute_edge_cost_from_odd_3_wheel(i,j,u,v, f);
//      }
//   }
//
//   double compute_edge_cost_from_odd_3_wheel(const std::size_t i, const std::size_t j, const std::size_t u, const std::size_t v, const multicut_odd_3_wheel_factor& f)
//   {
//      assert(i < j && u < v);
//      std::array<std::size_t,4> idx{i,j,u,v};
//      std::sort(idx.begin(), idx.end());
//      double min_participating_labelings = std::numeric_limits<double>::infinity();
//      double min_non_participating_labelings;
//      // non participating labelings have always != 4 cut edges
//      min_non_participating_labelings = std::min({0.0, f[0], f[1], f[2], f[3], f[7], f[8], f[9], f[10], f[11], f[12], f[13]});
//      // for participating edges, the axle and wheel edge must be cut
//      if(idx[0] == i && idx[1] == j) { // first and sixth edge
//         assert(idx[2] == u && idx[3] == v);
//         min_participating_labelings = std::min(f[5], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
//      } else if(idx[0] == i && idx[2] == j) { // second and fifth
//         assert(idx[1] == u && idx[3] == v);
//         min_participating_labelings = std::min(f[4], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
//      } else if(idx[0] == i && idx[3] == j) { // third and fourth
//         assert(idx[1] == u && idx[2] == v);
//         min_participating_labelings = std::min(f[4], f[5]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
//      } else if(idx[1] == i && idx[2] == j) { // fourth and third
//         assert(idx[0] == u && idx[3] == v);
//         min_participating_labelings = std::min(f[4], f[5]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
//      } else if(idx[1] == i && idx[3] == j) { // fifth and second
//         assert(idx[0] == u && idx[2] == v);
//         min_participating_labelings = std::min(f[4], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
//      } else if(idx[2] == i && idx[3] == j) { // sixth and first
//         assert(idx[0] == u && idx[1] == v);
//         min_participating_labelings = std::min(f[5], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
//      } else {
//         assert(false);
//      }
//      //assert(false);
//
//      return min_participating_labelings - min_non_participating_labelings; 
//   }
//
//   /*
//      static std::array<REAL,2> odd_3_wheel_xxx_to_do(const INDEX i, const INDEX j, const INDEX u, const INDEX v, const typename base_constructor::odd_3_wheel_factor& f)
//      {
//      REAL min_participating_labelings = std::numeric_limits<REAL>::infinity();
//      REAL min_non_participating_labelings;
//
//      min_non_participating_labelings = std::min({0.0, f[0], f[1], f[2], f[3], f[7], f[8], f[9], f[10], f[11], f[12], f[13]});
//
//      // for participating edges, the axle and wheel edge must be cut
//      if(idx[0] == i && idx[1] == j) { // first and sixth edge
//         assert(idx[2] == u && idx[3] == v);
//         min_participating_labelings = std::min(f[5], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
//      } else if(idx[0] == i && idx[2] == j) { // second and fifth
//         assert(idx[1] == u && idx[3] == v);
//         min_participating_labelings = std::min(f[4], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
//      } else if(idx[0] == i && idx[3] == j) { // third and fourth
//         assert(idx[1] == u && idx[2] == v);
//         min_participating_labelings = std::min(f[4], f[5]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
//      } else if(idx[1] == i && idx[2] == j) { // fourth and third
//         assert(idx[0] == u && idx[3] == v);
//         min_participating_labelings = std::min(f[4], f[5]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
//      } else if(idx[1] == i && idx[3] == j) { // fifth and second
//         assert(idx[0] == u && idx[2] == v);
//         min_participating_labelings = std::min(f[4], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
//      } else if(idx[2] == i && idx[3] == j) { // sixth and first
//         assert(idx[0] == u && idx[1] == v);
//         min_participating_labelings = std::min(f[5], f[6]);
//         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
//      } else {
//         assert(false);
//      }
//
//      return {min_participating_labelings, min_non_participating_labelings};
//   }
//*/
//
//   using triangle_intersection_type = std::tuple<INDEX,INDEX, typename base_constructor::triplet_factor*, typename base_constructor::triplet_factor*>;
//
//   void compute_triangles( // triangle are pyramids, though, search for better name
//         const INDEX i, const INDEX j, const REAL minTh, 
//         const typename base_constructor::triplet_connections& connected_triplets,
//         std::vector<triangle_intersection_type>& common_edges,
//         std::unordered_map<INDEX,INDEX>& origToCompressedNode, 
//         std::vector<INDEX>& compressedToOrigNode, 
//         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges)
//   {
//      assert(i < j);
//      origToCompressedNode.clear();
//      compressedToOrigNode.clear();
//      compressedEdges.clear();
//
//      auto merge = [](const auto a, const auto b) -> triangle_intersection_type { 
//         assert(a.nodes[0] == b.nodes[0] && a.nodes[1] == b.nodes[1]);
//         return std::make_tuple(a.nodes[0], a.nodes[1], a.f, b.f); 
//      };
//
//      // find all triangles ijk
//
//      // find all edges uv such that there exist edge triplets iuv and juv. 
//      // this is done by sorting all triplets which have node i and node j, and intersecting the set
//      auto intersects_iter_end = set_intersection_merge(
//            connected_triplets[i].begin(), connected_triplets[i].end(),
//            connected_triplets[j].begin(), connected_triplets[j].end(),
//            common_edges.begin(), [](const auto& a, const auto& b) { return a.operator<(b); }, merge);
//
//      for(auto n=common_edges.begin(); n != intersects_iter_end; ++n) {
//         const INDEX u = std::get<0>(*n);
//         const INDEX v = std::get<1>(*n);
//         assert(u < v);
//         const auto& iuv = std::get<2>(*n)->get_factor();
//         const auto& juv = std::get<3>(*n)->get_factor();
//
//         const REAL dual_increase = compute_edge_cost(i,j,u,v);
//
//         if(dual_increase >= minTh) { // add edge uv to bipartite graph
//
//            if(origToCompressedNode.find(u) == origToCompressedNode.end()) {
//               origToCompressedNode.insert(std::make_pair(u, origToCompressedNode.size()));
//               compressedToOrigNode.push_back(u);
//            }
//            if(origToCompressedNode.find(v) == origToCompressedNode.end()) {
//               origToCompressedNode.insert(std::make_pair(v, origToCompressedNode.size()));
//               compressedToOrigNode.push_back(v);
//            }
//            const INDEX uc = origToCompressedNode[u];
//            const INDEX vc = origToCompressedNode[v];
//            assert(uc != vc);
//            compressedEdges.push_back(std::make_tuple(uc,vc, dual_increase));
//         }
//      } 
//   }
//
//   struct bicycle_candidate {
//      std::array<INDEX,5> idx;
//      REAL cost;
//   };
//
//   template<typename ITERATOR>
//   void triangulate_odd_bicycle_wheel(const INDEX i, const INDEX j, const REAL cost, ITERATOR path_begin, ITERATOR path_end, std::vector<bicycle_candidate>& candidates)
//   {
//      assert(i < j);
//      assert(std::distance(path_begin, path_end) >= 3);
//      this->cycle_normal_form(path_begin, path_end);
//      const INDEX first_node = *path_begin;
//      for(auto it=path_begin+1; it+1!=path_end; ++it) {
//         std::array<INDEX,5> nodes({i,j,first_node, *it, *(it+1)});
//         std::sort(nodes.begin(), nodes.end());
//         assert(HasUniqueValues(nodes));
//         candidates.push_back({nodes, cost});
//      }
//   }
//
//
//   INDEX find_violated_odd_bicycle_wheels(const INDEX max_factors_to_add)
//   {
//      if(this->number_of_edges() > 2) {
//         // preprocessing: sort triplets for fast intersection later
//         auto connected_triplets = this->compute_connected_triplets();
//
//         std::vector<std::tuple<INDEX,REAL>> threshold(this->unaryFactorsVector_.size()); // edge number and threshold
//         std::vector<bicycle_candidate> odd_bicycle_candidates;
//         // given a cut axle edge (negative cost), find cut wheel edges such that among the four spokes exactly two are cut and two are connected.
//     
//#pragma omp parallel
//         {
//            using intersection_type = std::tuple<INDEX,INDEX, typename base_constructor::triplet_factor*, typename base_constructor::triplet_factor*>;
//            std::vector<intersection_type> common_edges(this->number_of_edges()); // possibly this is a bit large!
//            auto merge = [](const auto a, const auto b) -> intersection_type { 
//               assert(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b));
//               return std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a), std::get<2>(b)); 
//            };
//            std::vector<bicycle_candidate> odd_bicycle_candidates_local;
//
//            std::unordered_map<INDEX,INDEX> origToCompressedNode;
//            std::vector<INDEX> compressedToOrigNode;
//            std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;
//
//#pragma omp for schedule(guided) nowait
//            for(INDEX e=0; e<this->unaryFactorsVector_.size(); ++e) {
//
//               // edge ij will be treated as axle of odd bicycle wheel
//               const INDEX i = std::get<0>(this->unaryFactorsVector_[e])[0];
//               const INDEX j = std::get<0>(this->unaryFactorsVector_[e])[1];
//               const REAL cost_ij = std::get<1>(this->unaryFactorsVector_[e])->get_factor()->operator[](0); 
//               if(cost_ij < -eps) {
//                  //origToCompressedNode.clear();
//                  //compressedToOrigNode.clear();
//                  //compressedEdges.clear(); 
//
//                  compute_triangles(i, j, eps, connected_triplets, common_edges, origToCompressedNode, compressedToOrigNode, compressedEdges); 
//                  const REAL th = this->compute_odd_cycle_threshold(origToCompressedNode, compressedToOrigNode, compressedEdges);
//
//                  if(th > eps) {
//                     auto path = this->compute_path_in_bipartite_graph(origToCompressedNode, compressedToOrigNode, compressedEdges, th);
//                     triangulate_odd_bicycle_wheel(i,j, th, path.begin(), path.end(), odd_bicycle_candidates_local);
//                  }
//               } 
//            }
//#pragma omp critical
//            odd_bicycle_candidates.insert(odd_bicycle_candidates.end(), odd_bicycle_candidates_local.begin(), odd_bicycle_candidates_local.end());
//         }
//
//         std::sort(odd_bicycle_candidates.begin(), odd_bicycle_candidates.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });
//         INDEX no_factors_added = 0;
//         for(INDEX i=0; i<odd_bicycle_candidates.size(); ++i) {
//            if(!this->has_odd_bicycle_3_wheel( odd_bicycle_candidates[i].idx )) {
//               this->add_odd_bicycle_3_wheel( odd_bicycle_candidates[i].idx );
//               ++no_factors_added;
//               if(no_factors_added >= max_factors_to_add) {
//                  break;
//               }
//            } 
//         }
//         return no_factors_added;
//      } else {
//         return 0;
//      }
//   }
//
//private:
//
//   std::unordered_map<std::array<INDEX,5>, odd_bicycle_3_wheel_factor_container*> odd_bicycle_3_wheel_factors_;
//
//   //std::unordered_map<std::vector<std::tuple<INDEX,INDEX, typename base_constructor::odd_3_wheel_factor_container*>>> odd_3_wheel_factor_by_indices_; // if odd 3 wheel factor with indices (i1,i2,i3,i4) exists, then (i1,i2,i3,i4) will be in the hash indexed by all two-subsets of the indices.
//
//};
//
//template<class FACTOR_MESSAGE_CONNECTION, INDEX MULTICUT_CONSTRUCTOR_NO, INDEX MRF_CONSTRUCTOR_NO, INDEX MULTICUT_POTTS_MESSAGE_NO>
//class multiway_cut_constructor
//{
//   public:
//      using FMC = FACTOR_MESSAGE_CONNECTION;
//
//      template<typename SOLVER>
//         multiway_cut_constructor(SOLVER& s) 
//         :
//            mc_constructor(s.template GetProblemConstructor<MULTICUT_CONSTRUCTOR_NO>()),
//            mrf_constructor(s.template GetProblemConstructor<MRF_CONSTRUCTOR_NO>()),
//            lp_(&s.GetLP())
//   {}
//
//      INDEX Tighten(const INDEX no_factors_to_add)
//      {
//         // add multicut edges only after first tightening. Before that, it makes no sense to have them.
//         if(multicut_added_ == false) {
//            multicut_added_ = true;
//            for(INDEX i=0; i<mrf_constructor.GetNumberOfPairwiseFactors(); ++i) {
//               auto vars = mrf_constructor.GetPairwiseVariables(i);
//               const INDEX i1 = std::get<0>(vars);
//               const INDEX i2 = std::get<1>(vars);
//               auto* f = mc_constructor.AddUnaryFactor(i1, i2, 0.0);
//               auto* m = new typename FMC::multicut_edge_potts_message_container(f, mrf_constructor.GetPairwiseFactor(i)); 
//               lp_->AddMessage(m); 
//               lp_->add_factor_relation(mrf_constructor.GetUnaryFactor(i1), f);
//               lp_->add_factor_relation(f, mrf_constructor.GetUnaryFactor(i2));
//            }
//
//            // the multicut constructor will not have done anything (no factors yet). Hence call it again from here
//
//         }
//         return 0;
//
//         //return MULTICUT_CONSTRUCTOR::Tighten(no_factors_to_add);
//      }
//   private:
//
//      using multicut_constructor_type = meta::at_c<typename FMC::ProblemDecompositionList, MULTICUT_CONSTRUCTOR_NO>;
//      multicut_constructor_type& mc_constructor;
//
//      using mrf_constructor_type = meta::at_c<typename FMC::ProblemDecompositionList, MRF_CONSTRUCTOR_NO>;
//      mrf_constructor_type& mrf_constructor;
//
//      LP<FMC>* lp_; // TODO: not needed, base constructors holds lp already
//
//      bool multicut_added_ = false;
//};
//
//template< typename MULTICUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG>
//class lifted_multicut_constructor : public lifted_constructor< 
//				    lifted_multicut_constructor<MULTICUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG>,
//				    MULTICUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG> 
//{
//public:
//   using base_constructor = lifted_constructor<
//	   lifted_multicut_constructor<MULTICUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG>,
//	   MULTICUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG>;
//
//   using base_constructor::base_constructor;
//
//   // use GAEC and Kernighan&Lin algorithm of andres graph package to compute primal solution
//   static std::vector<char> round(std::vector<typename base_constructor::Edge> base_edges, std::vector<typename base_constructor::Edge> lifted_edges, std::vector<REAL> edge_values)
//   {
//      std::vector<char> labeling(edge_values.size(),0);
//      if(edge_values.size() > 0) {
//	      INDEX no_nodes = 0;
//	      for(auto e : lifted_edges) {
//		      no_nodes = std::max({no_nodes, e[0]+1, e[1]+1});
//	      }
//
//	      andres::graph::Graph<> original_graph(no_nodes);
//	      andres::graph::Graph<> lifted_graph(no_nodes);
//
//	      for(const auto& e : base_edges) {
//		      original_graph.insertEdge(e[0], e[1]);
//		      lifted_graph.insertEdge(e[0],e[1]);
//	      }
//	      for(const auto& e : lifted_edges) {
//		      lifted_graph.insertEdge(e[0], e[1]);
//	      }
//
//	      andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(original_graph,lifted_graph,edge_values,labeling);
//	      andres::graph::multicut_lifted::kernighanLin(original_graph,lifted_graph,edge_values,labeling,labeling);
//      }
//      return labeling;
//   }
//};
//
//} // end namespace LPMP
