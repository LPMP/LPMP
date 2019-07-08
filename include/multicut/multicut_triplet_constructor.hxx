#pragma once

#include <vector>
#include <cassert>
#include <future>
#include "cut_base/cut_base_triplet_constructor.hxx"
#include "multicut_instance.h"
#include "multicut_cycle_packing.h"
#include "graph.hxx"
#include "multicut_greedy_additive_edge_contraction.h"
#include "multicut_greedy_edge_fixation.h"

namespace LPMP {

template<
   class FACTOR_MESSAGE_CONNECTION, 
   typename UNARY_FACTOR, typename TRIPLET_FACTOR,
   typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2
   >
class multicut_triplet_constructor : public cut_base_triplet_constructor<
			     //multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0,UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2 >, 
			     FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2> 
{
public:
   using type = multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;
   using FMC = FACTOR_MESSAGE_CONNECTION;
   using base_constructor = cut_base_triplet_constructor<
	   //multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0,UNARY_TRIPLET_MESSAGE_1,  UNARY_TRIPLET_MESSAGE_2 >, 
	   FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;

   using weighted_graph = graph<double>;

   template<typename SOLVER>
       multicut_triplet_constructor(SOLVER& s);

   // add triplet factor with prespecified cost
   // edges on triplet are (uv,uw,vw) (lexicographical)
   TRIPLET_FACTOR* add_higher_order_triplet(
         const std::size_t u, const std::size_t v, const std::size_t w, 
         const double c000, const double c011, const double c101, const double c110, const double c111);

   static double triplet_cost(const double cost_ij, const double cost_ik, const double cost_jk);

   std::size_t find_violated_cycles(const std::size_t max_triplets_to_add);

   void construct(multicut_instance mc);
   bool CheckPrimalConsistency() const;
   std::size_t Tighten(const std::size_t no_constraints);
   //static std::vector<char> round(std::vector<typename base_constructor::edge> edges);
   static multicut_edge_labeling round(const multicut_instance input, const std::string method);
   void ComputePrimal();
   void Begin();
   void End();
protected:
    std::future<multicut_edge_labeling> primal_result_handle_;

    TCLAP::ValueArg<std::string> rounding_method_arg_;
    TCLAP::SwitchArg no_informative_factors_arg_;
    TCLAP::SwitchArg no_tightening_packing_arg_;
};

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    template<typename SOLVER>
multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::multicut_triplet_constructor(SOLVER& s)
    :
        rounding_method_arg_("", "multicutRounding", "method for rounding primal solution", false, "gaec", "{gaec|gef}", s.get_cmd()),
        no_informative_factors_arg_("", "noInformativeFactorReparametrization", "do not make factors informative when rounding and tightening", s.get_cmd(), false),
        no_tightening_packing_arg_("", "noTighteningPacking", "do not pack inequalities after tightening", s.get_cmd(), false),
        base_constructor(s)
{}

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    TRIPLET_FACTOR* multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::add_higher_order_triplet
(
 const std::size_t u, const std::size_t v, const std::size_t w, 
 const double c000, const double c011, const double c101, const double c110, const double c111
 )
{
    auto* f = this->add_triplet_factor(u,v,w);
    auto* t = f->get_factor();
    this->lp_->add_to_constant(c000);
    (*t)[0] = c110 - c000;
    (*t)[1] = c101 - c000;
    (*t)[2] = c011 - c000;
    (*t)[3] = c111 - c000;
    return f;
}

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
double multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::triplet_cost(const double cost_ij, const double cost_ik, const double cost_jk)
{
    return std::min({0.0, cost_ij+cost_ik, cost_ij+cost_jk, cost_ik+cost_jk, cost_ij+cost_ik+cost_jk}); 
}

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
   std::size_t multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::find_violated_cycles(const std::size_t max_triplets_to_add)
   {
      std::vector<std::tuple<std::size_t, std::size_t, double, bool> > negative_edges; // endpoints, edge cost, searched positive path with given endpoints?

      struct weighted_edge : public std::array<std::size_t,2> { double cost; double& operator()() { return cost; } };
      std::vector<weighted_edge> positive_edges;

      // we can speed up compution by skipping path searches for node pairs which lie in different connected components. Connectedness is stored in a union find structure
      double pos_th = 0.0;
      for(auto& it : this->unary_factors_vector_) {
         const double v = (*it.second->get_factor())[0];
         const std::size_t i = std::get<0>(it.first);
         const std::size_t j = std::get<1>(it.first);
         if(v >= 0.0) {
            pos_th = std::max(pos_th, v);
            positive_edges.push_back({i,j,v});
         } else {
            negative_edges.push_back(std::make_tuple(i,j,v,false));
         }
      }

      if(negative_edges.size() == 0 || negative_edges.size() == this->unary_factors_vector_.size()) { return 0; }

      weighted_graph posEdgesGraph(positive_edges.begin(), positive_edges.end()); // graph consisting of positive edges

      // do zrobienia: possibly add reparametrization of triplet factors additionally. Do it analogously as in graph matching via send_messages_to_unaries()?

      std::sort(negative_edges.begin(), negative_edges.end(), [](const auto& e1, const auto& e2)->bool { return std::get<2>(e1) < std::get<2>(e2); });


         // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
         // here, longest path is sought after only the edges with positive taken into account
         // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min

         //MostViolatedPathData mp(posEdgesGraph);
         //BfsData mp(posEdgesGraph);

         union_find uf(this->no_nodes());
         std::size_t triplets_added = 0;
         const double initial_th = 0.6*std::min(-std::get<2>(negative_edges[0]), pos_th);
         bool zero_th_iteration = true;

         for(double th=initial_th; th>=eps || zero_th_iteration; th*=0.1) {
            if(th < eps) {
               if(triplets_added <= 0.01*max_triplets_to_add) {
                  // we would first like to go on with odd wheels and odd bicycle wheels before actually adding triplets with no guaranteed dual increase.
                  //std::cout << "additional separation with no guaranteed dual increase, i.e. th = 0\n";
                  //th = 0.0;
                  //zero_th_iteration = false;
               } else {
                  break;
               }
            }
            // first update union find datastructure
            for(auto& it : this->unary_factors_vector_) {
               const double v = (*it.second->get_factor())[0];
               if(v >= th) {
                  const std::size_t i = std::get<0>(it.first);
                  const std::size_t j = std::get<1>(it.first);
                  uf.merge(i,j);   
               }
            }

            std::vector<typename base_constructor::triplet_candidate> triplet_candidates;

#pragma omp parallel 
            {
               std::vector<typename base_constructor::triplet_candidate> triplet_candidates_local;
               bfs_data<weighted_graph> mp2(posEdgesGraph);
#pragma omp for schedule(guided) nowait
               for(std::size_t c=0; c<negative_edges.size(); ++c) {
                  const std::size_t i = std::get<0>(negative_edges[c]);
                  const std::size_t j = std::get<1>(negative_edges[c]);
                  const double v = std::get<2>(negative_edges[c]);
                  const bool already_used_for_path_search = std::get<3>(negative_edges[c]);
                  //if(-v <= th) break;
                  //if(already_used_for_path_search) continue;
                  if(-v > th && !already_used_for_path_search && uf.thread_safe_connected(i,j)) {
                     //auto cycle = mp.FindPath(i,j,posEdgesGraph);
                     auto mask_small_edges = [th](const std::size_t i, const std::size_t j, const double cost, const std::size_t distance) { return cost >= th; };
                     double cycle_cap = std::numeric_limits<double>::infinity();
                     auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const double cost) { cycle_cap = std::min(cycle_cap, cost); };
                     auto cycle = mp2.find_path(i, j, mask_small_edges, cycle_capacity);
                     const double dualIncrease = std::min(-v, cycle_cap);
                     assert(cycle.size() > 0);
                     if(cycle.size() > 0) {
                        this->triangulate_cycle(dualIncrease, cycle.begin(), cycle.end(), triplet_candidates_local);
                        //triplets_added += AddCycle(std::get<1>(cycle));
                        //if(triplets_added > max_triplets_to_add) {
                        //   return triplets_added;
                        //}
                     } else {
                        throw std::runtime_error("No path found although there should be one"); 
                     }
                  }
               }
#pragma omp critical
               {
                  triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_local.begin(), triplet_candidates_local.end()); 
            }
         }

         // sort by guaranteed increase in decreasing order
         std::sort(triplet_candidates.begin(), triplet_candidates.end());

         if(triplet_candidates.size() > 0 && diagnostics()) {
            std::cout << "best triplet candidate in triplet search has guaranteed dual improvement " << triplet_candidates[0].cost << "\n";
         }

         for(const auto& triplet_candidate : triplet_candidates) {
            const std::size_t i = triplet_candidate.nodes[0];
            const std::size_t j = triplet_candidate.nodes[1];
            const std::size_t k = triplet_candidate.nodes[2];
            if(!this->HasTripletFactor(i,j,k)) {
               this->add_triplet_factor(i,j,k);
               triplets_added++;
               if(triplets_added > max_triplets_to_add) {
                  break;
               } 
            }
         }

         triplet_candidates.clear();
      }

      return triplets_added;
   }

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
   void multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::construct(multicut_instance mc)
   {
       mc.normalize();
       for(const auto& e : mc.edges())
           this->add_edge_factor(e[0], e[1], e.cost);
       this->no_original_edges_ = this->unary_factors_vector_.size();
   }

template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
   bool multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::CheckPrimalConsistency() const
   {
     if(debug()) {
       std::cout << "checking primal feasibility for multicut ... ";
     }
     union_find uf(this->no_nodes());
     for(const auto& e : this->unary_factors_vector_) {
       auto* f = e.second; 
       if(f->get_factor()->primal()[0] == false) {
         // connect components 
         const std::size_t i = std::get<0>(e.first);
         const std::size_t j = std::get<1>(e.first);
         uf.merge(i,j);
       }
     }
     for(const auto& e : this->unary_factors_vector_) {
       auto* f = e.second; 
       if(f->get_factor()->primal()[0] == true) {
         const std::size_t i = std::get<0>(e.first);
         const std::size_t j = std::get<1>(e.first);
         // there must not be a path from i1 to i2 consisting of edges with primal value false only
         if(uf.connected(i,j)) {
           if(debug()) {
             std::cout << "solution infeasible: (" << i << "," << j << ") = true, yet there exists a path with false values only\n";
           }
           return false;
         }
       }
     }

     if(debug()) {
       std::cout << "solution feasible\n";
     }
     return true; 
   }

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    std::size_t multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::Tighten(const std::size_t no_constraints)
{
    if(!no_informative_factors_arg_.isSet())
        this->send_messages_to_edges(); // TODO: check if this is always helpful or only for triplet tightening
    const auto mc = this->template export_edges<multicut_instance>();
    cycle_packing cp = compute_multicut_cycle_packing(mc);
    if(debug())
        std::cout << "Found " << cp.no_cycles() << " cycles through cycle packing\n";
    const std::size_t no_triplets_before = this->number_of_triplets();
    for(std::size_t c=0; c<cp.no_cycles(); ++c) {
        const auto [cycle_begin, cycle_end] = cp.get_cycle(c);
        const double cycle_weight = cp.get_cycle_weight(c);
        this->triangulate_cycle(cycle_begin, cycle_end, cycle_weight, !no_tightening_packing_arg_.isSet());
    } 
    if(debug())
        std::cout << "Added " << this->number_of_triplets() - no_triplets_before << " triplets\n";
    return cp.no_cycles();
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    multicut_edge_labeling multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::round(const multicut_instance input, const std::string method)
{
    if(method == "gaec")
        return greedy_additive_edge_contraction(input);
    else if(method == "gef")
        return multicut_greedy_edge_fixation(input);
    else 
        throw std::runtime_error(std::string("rounding method not recognized: ") + method);
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
    void multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::ComputePrimal()
{
    if(!no_informative_factors_arg_.isSet())
        this->send_messages_to_edges();
    if(primal_result_handle_.valid() && primal_result_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        if(debug()) 
            std::cout << "read in primal multicut solution\n";
        auto multicut_sol = primal_result_handle_.get();
        this->write_labeling_into_factors(multicut_sol); 
    }

    if(!primal_result_handle_.valid() || primal_result_handle_.wait_for(std::chrono::seconds(0)) == std::future_status::deferred) {
        if(debug()) 
            std::cout << "export multicut problem for rounding\n";

        const auto instance = this->template export_edges<multicut_instance>();
        primal_result_handle_ = std::async(std::launch::async, round, std::move(instance), rounding_method_arg_.getValue());
    } 
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
void multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::Begin()
{
    base_constructor::Begin();
    const auto instance = this->template export_edges<multicut_instance>();
    primal_result_handle_ = std::async(std::launch::async, round, instance, rounding_method_arg_.getValue());

    // compute cycle packing and add returned cycles to problem formulation. Also reparametrize edges (possibly do not do this?)
    Tighten(std::numeric_limits<std::size_t>::max());
    //cycle_packing cp = compute_multicut_cycle_packing(instance);
    //for(std::size_t c=0; c<cp.no_cycles(); ++c) {
    //    const auto [cycle_begin, cycle_end] = cp.get_cycle(c);
    //    const double cycle_weight = cp.get_cycle_weight(c);
    //    this->triangulate_cycle(cycle_begin, cycle_end, cycle_weight);
    //}
}

    template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
void multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::End()
{
    base_constructor::End();
    if(primal_result_handle_.valid()) {
        primal_result_handle_.wait();
        if(debug()) 
            std::cout << "read in primal multicut solution\n";
        auto multicut_sol = primal_result_handle_.get();
        this->write_labeling_into_factors(multicut_sol); 
    }
}

//template<class FACTOR_MESSAGE_CONNECTION, typename UNARY_FACTOR, typename TRIPLET_FACTOR, typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2>
//   std::vector<char> multicut_triplet_constructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>::round(std::vector<typename base_constructor::edge> edges)
//   {
//       assert(false);
//       return {};
//     std::size_t no_nodes = 0;
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

} // namespace LPMP
