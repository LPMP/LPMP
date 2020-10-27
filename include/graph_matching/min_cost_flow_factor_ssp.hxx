#pragma once

#include "mcf_ssp.hxx"
#include "graph_matching/matching_problem_input.h"
#include "config.hxx"

namespace LPMP {

// mcf holds reverse copies of edges. arcs are ordered lexicographically by (tail,node).
// this creates problems for some applications, where only one direction of arcs is needed, but it may be advantageous for other applications.
class min_cost_flow_factor {
public:
   struct Edge {
      INDEX start_node, end_node, lower_bound, upper_bound; // we assume that the minimum capacity is 0
      REAL cost;
   };
   
   min_cost_flow_factor(const std::vector<Edge>& edges, const std::vector<SIGNED_INDEX>& supply_demand) 
      : mcf_(supply_demand.size(), edges.size())
   {
      INDEX noNodes_ = supply_demand.size();
      INDEX noEdges_ = edges.size();

      for(auto& e : edges) {
         mcf_.add_edge(e.start_node, e.end_node, e.lower_bound, e.upper_bound, e.cost);
      }
      for(INDEX i=0; i<supply_demand.size(); ++i) {
         mcf_.add_node_excess(i, supply_demand[i]);
      }
      mcf_.order();
      mcf_.solve();
   }

   min_cost_flow_factor(const linear_assignment_problem_input& input, const double scaling = 1.0)
   : mcf_(input.no_mcf_nodes(), input.no_mcf_edges())
   {
      input.initialize_mcf(mcf_);
   }

   REAL EvaluatePrimal() const
   {
      return mcf_.objective();
      // do zrobienia: we must check feasiblity of primal as regards flow constraints. For this, a feasible flow on auxiliary variables must be built.
      // check feasibility
      /*
      std::vector<SIGNED_INDEX> excess(minCostFlow_->GetNodeNum(),0);
      for(INDEX e=0; e<repam.size(); ++e) {
         const INDEX tail = minCostFlow_->GetTailNodeId(e);
         const INDEX head = minCostFlow_->GetHeadNodeId(e);
         excess[tail] += primal[e];
         excess[head] -= primal[e];
      }
      for(INDEX i=0; i<excess.size(); ++i) {
         if(excess[i] != demand_[i]) {
            return std::numeric_limits<REAL>::infinity();
         }
      }
      */
      // to do: this is a very hacky implementation and it need not be correct for anything except assignment problems!
      // first check whether primal belongs to a feasible flow
      /*
      for(INDEX i=0; i<minCostFlow_->GetNodeNum(); ++i) {
         SIGNED_INDEX excess = 0;
         INDEX a_idx = minCostFlow_->StartingArc(i);
         for(INDEX c=0; c<minCostFlow_->NoArcs(i); ++c, ++a_idx) {
            if(primal[a_idx] == unknownState) {
               excess = minCostFlow_->GetDemand(i);
               continue;
            }
            const SIGNED_INDEX sign = minCostFlow_->GetCap(a_idx) == 1 ? 1 : -1;
            excess += sign*SIGNED_INDEX(primal[a_idx]);
         }
         if(excess != minCostFlow_->GetDemand(i)) {
            std::cout << "excess = " << excess << " demand = " << minCostFlow_->GetDemand(i) << "\n";
            return std::numeric_limits<REAL>::infinity();
         }
      }
      REAL cost = 0.0;
      assert(repam.size() == minCostFlow_->GetArcNum());
      for(INDEX e=0; e<repam.size(); e++) {
         if(minCostFlow_->GetCap(e) != 0 && minCostFlow_->GetCap(e) != 1) {
            assert(primal[e] == false || repam[e] == 0.0);
         }
         if(minCostFlow_->GetCap(e) == 1 && repam[e] != 0.0) { // some edges are auxiliary and have cost zero. Primal value is not set for them.
            assert(minCostFlow_->GetCap( minCostFlow_->GetReverseArcId(e) ) == 0);
            assert(primal[e] == primal[ minCostFlow_->GetReverseArcId(e) ]);
            assert(primal[e] != unknownState);
            cost += repam[e]*primal[e]; 
         }
      }
      return 0.5*cost;
      */
   }

   void MaximizePotential()
   { 
      mcf_.solve();
   }

   void MaximizePotentialAndComputePrimal()
   { 
      mcf_.solve();
   }

   REAL LowerBound() const
   {
      auto obj = mcf_.solve();
      //std::cout << "mcf objective = " << obj << " = " << mcf_.objective() << "\n";
      return obj;
   }

   // for the SendMessages update step, compute maximal cost change such that current labeling stays optimal
   // do zrobienia: is possibly an array with different values from which to compute the maximal perturbation
   // do zrobienia: rename active_arcs
   template<class REPAM_ARRAY, typename ACTIVE_EDGES_ARRAY>
   void MaximallyPerturbCosts(const REPAM_ARRAY& repam, const ACTIVE_EDGES_ARRAY& active_edges) const
   {
      assert(false);
      // first read reduced cost of mcf into repam update
      /*
      for(INDEX i=0; i<minCostFlow_->GetArcNum(); ++i) {
         if(minCostFlow_->GetTailNodeId(i) < minCostFlow_->GetHeadNodeId(i)) { // avoid copying twice
            repamUpdateFlow_->SetCost(i,minCostFlow_->GetCost(i));
         }
      }
      // set capacities
      assert(active_edges.size() == minCostFlow_->GetArcNum());
      for(INDEX i=0; i<minCostFlow_->GetArcNum(); ++i) {
         if(active_edges[i]) {
            if(minCostFlow_->GetRCap(i) == 0) {
               repamUpdateFlow_->SetCap(i,-1);
            } else {
               assert(minCostFlow_->GetRCap(i) == 1);
               // beware of overflows!
               repamUpdateFlow_->SetCap(i,std::numeric_limits<SIGNED_INDEX>::max()/16);
            }
         } else {
            repamUpdateFlow_->SetCap(i,std::numeric_limits<SIGNED_INDEX>::max()/16);
         }
      }
      repamUpdateFlow_->Solve();
      */
      /*
      // we assume here that repam is the current one
      assert(active_edges.size() == size());
      const SIGNED_INDEX max_cap = 100000;
      for(INDEX a=0; a<minCostFlow_->no_arcs(); ++a) {
         // do zrobienia: make checks in terms of residual capacity
         const SIGNED_INDEX flow = minCostFlow_->get_flow(a);
         if(active_edges[a] == true) {
            if(flow == 0) { 
               repamUpdateFlow_->set_cap(a,max_cap);
            } else if(flow == 1) {
               repamUpdateFlow_->set_cap(a,-1);
            } else if(flow == -1) { // reverse edge, also treat it
               repamUpdateFlow_->set_cap(a,1);
            } else {
               assert(false);
            }
         } else {
            repamUpdateFlow_->set_cap(a,max_cap);
         }
      }
      // copy cost
      repamUpdateFlow_->copy_cost(minCostFlow_);
      // run
      repamUpdateFlow_->cs2_cost_restart();
      // read out dual variables which will go into reparametrization update
      std::vector<REAL> repam_cost;
      INDEX noActiveEdges = 0;
      for(INDEX e=0; e<active_edges.size(); ++e) {
         noActiveEdges += active_edges[e];
      }
      repam_cost.reserve(noActiveEdges);
      for(INDEX e=0; e<active_edges.size(); ++e) {
         if(active_edges[e] == true) {
            // do zrobienia: why should this be true?
            //const REAL reducedCost = Descale(Scale(minCostFlowSolver_->get_cost(e) + repamUpdateFlow_->.potential(graph_->source(arcs_[e])) - minCostFlowRepamUpdate.potential(graph_->target(arcs_[e])));
            const REAL reducedCost = Descale(repamUpdateFlow_->get_reduced_cost(e));
            if(minCostFlow_->get_flow(e) == 1) {
               repam_cost.push_back( -reducedCost );
            } else if(minCostFlow_->get_flow(e) == 0) {
               repam_cost.push_back( -reducedCost );
            } else if(minCostFlow_->get_flow(e) == -1) {
               repam_cost.push_back( -reducedCost );
            } else {
               assert(false); // flow is not 0/1
            }
         }
      }
      )return repam_cost;
      */

      /*
      GraphType::ArcMap<LONG_SIGNED_INDEX> cost(*graph_);
      for(INDEX e=0; e<active_edges.size(); ++e) {
         const SIGNED_INDEX flow = minCostFlow_->flow(arcs_[e]);
         if(active_edges[e] == true) {
            assert( flow == 0 || flow == 1); 
            if(flow == 1) {
               lower[arcs_[e]] = -max_cap;
               upper[arcs_[e]] = -1;
            } else if(flow == 0) {
               lower[arcs_[e]] = 1;
               upper[arcs_[e]] = max_cap;
            }
         } else {
            if(flow == 0) {
               lower[arcs_[e]] = 0;
               upper[arcs_[e]] = +max_cap;
            } else {
               lower[arcs_[e]] = -max_cap;
               upper[arcs_[e]] = +max_cap;
            }
         }
         assert(lower[arcs_[e]] < upper[arcs_[e]]);
         cost[arcs_[e]] = Scale(repam[e]);
      }
      auto ret = minCostFlowRepamUpdate.lowerMap(lower).upperMap(upper).costMap(cost).run();
      assert(ret == MinCostFlowSolverType::OPTIMAL);
      if(ret != MinCostFlowSolverType::OPTIMAL) {
         throw std::runtime_error("Reparametrization problem could not be solved to optimality");
      }

      std::vector<REAL> repam_cost;
      INDEX noActiveEdges = 0;
      for(INDEX e=0; e<active_edges.size(); ++e) {
         noActiveEdges += active_edges[e];
      }
      repam_cost.reserve(noActiveEdges);
      for(INDEX e=0; e<active_edges.size(); ++e) {
         if(active_edges[e] == true) {
            // do zrobienia: why should this be true?
            const REAL reducedCost = Descale(Scale(repam[e]) + minCostFlowRepamUpdate.potential(graph_->source(arcs_[e])) - minCostFlowRepamUpdate.potential(graph_->target(arcs_[e])));
            if(minCostFlow_->flow(arcs_[e]) == 1) {
               repam_cost.push_back( -reducedCost );
            } else if(minCostFlow_->flow(arcs_[e]) == 0) {
               repam_cost.push_back( -reducedCost );
            } else {
               assert(false); // flow is not 0/1
            }
         }
      }
      // do zrobienia: construct test problem with perturbed cost and see whether solution is still optimal and just at edge
      return repam_cost;
      */
   }

   /*
   INDEX GetNumberOfAuxVariables() const { return minCostFlow_->GetEdgeNum() - no_binary_edges_; }

  void CreateConstraints(LpInterfaceAdapter* lp) const {
     // the first no_binary_edges go into usual variables, the last ones into auxiliary
     for(INDEX e=no_binary_edges_; e<minCostFlow_->GetEdgeNum(); ++e) {
        auto var = lp->GetAuxVariable(e-no_binary_edges_);
        lp->SetVariableBound(var, -REAL(minCostFlow_->GetReverseCap(e)), REAL(minCostFlow_->GetCap(e)));
     }
     // flow conservation constraints
     std::vector<LinExpr> flow_conservation(minCostFlow_->GetNodeNum());
     for(auto& c : flow_conservation) {
        c = lp->CreateLinExpr();
     }
     for(INDEX e=0; e<no_binary_edges_; ++e) {
        auto var = lp->GetVariable(e);
        const INDEX i=minCostFlow_->GetTailNodeId(e);
        const INDEX j=minCostFlow_->GetHeadNodeId(e);
        flow_conservation[i] += var;
        flow_conservation[j] -= var;
     }
     for(INDEX e=no_binary_edges_; e<minCostFlow_->GetEdgeNum(); ++e) {
        auto var = lp->GetAuxVariable(e-no_binary_edges_);
        const INDEX i=minCostFlow_->GetTailNodeId(e);
        const INDEX j=minCostFlow_->GetHeadNodeId(e);
        flow_conservation[i] += var;
        flow_conservation[j] -= var;
     }
     for(INDEX i=0; i<minCostFlow_->GetNodeNum(); ++i) {
        auto rhs = lp->CreateLinExpr();
        rhs += demand_[i]; // do zrobienia: store this directly and not in min cost flow problem
        lp->addLinearEquality(flow_conservation[i],rhs);
     }
  }
  */

  void init_primal() { }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) 
  { 
     assert(false);
  }

  auto export_variables() const { return std::tie(); }

  template<typename SOLVER>
  void construct_constraints(SOLVER& s)
  { 
    assert(false); 
  }

  template<typename SOLVER>
  void convert_primal(SOLVER& s)
  { 
    assert(false); 
  }

   mutable MCF::SSP<long,REAL> mcf_; // LowerBound is a const function, yet we call solve there.


private:
};

// used in graph matching via mcf. 
// do zrobienia: rename
template<INDEX COVERING_FACTOR> // COVERING_FACTOR specifies how often an mcf edge is covered at most. Usually, either one or two
class unary_min_cost_flow_message {
// right factor is MinCostFlowFactor describing an 1:1 assignment as constructed above
// left factors are the left and right simplex corresponding to it. 
// Those are ordered as the left and right nodes
// templatize edgeIndex_ to allow for more compact representation. Possibly hold reference to original edgeId structure instead of full vector
public:
   unary_min_cost_flow_message(const INDEX first_arc, const INDEX no_arcs) 
   : first_arc_(first_arc),
   no_arcs_(no_arcs)
   {}

   template<typename LEFT_POT, typename MSG_ARRAY>
   void send_message_to_right(const LEFT_POT& leftPot, MSG_ARRAY& msg, const REAL omega) 
   {
      assert(leftPot.size() == no_arcs_);

      for(INDEX i=0; i<leftPot.size(); ++i) {
         msg[i] -= omega*(leftPot[i]);
      }
   }

   // currently this function must be present, although SendMessagesToLeft is called
   template<typename RIGHT_POT, typename MSG_ARRAY>
   void send_message_to_left(const RIGHT_POT& r, MSG_ARRAY& msg, const REAL omega) 
   {
      throw std::runtime_error("message not implemented");
   }

   // to do: update this, remove RIGHT_REPAM
   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessagesToLeft_deactivated(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omega_begin)
   {
       throw std::runtime_error("message not implemented");
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      for(INDEX e=0; e<no_arcs_; ++e) {
         const auto flow = r.mcf_.flow(first_arc_+e);
         assert(flow == -1 || flow == 0 || flow == 1);
         if(std::abs(flow) == 1) {
            l.primal() = e;
         }
      }
   }

   template<typename LEFT_FACTOR>
   void RepamLeft(LEFT_FACTOR& left, const REAL msg, const INDEX dim) 
   {
      assert(dim < no_arcs_);
      assert(left.size() == no_arcs_);
      left[dim] += msg;
   }

   template<typename RIGHT_FACTOR>
   void RepamRight(RIGHT_FACTOR& right, const REAL msg, const INDEX dim) 
   {
      assert(dim < no_arcs_);
      //std::cout << "arc = " << first_arc_ + dim << ", sign = " << arc_sign(right.mcf_, first_arc_+dim) << ", msg = " << msg << ", cost of mcf = " << right.mcf_.cost(first_arc_+dim) << "\n";
      if(arc_sign(right.mcf_, first_arc_+dim) == -1.0) { assert(std::abs(msg) < 1e-8); }
      right.mcf_.update_cost(first_arc_ + dim, arc_sign(right.mcf_, first_arc_+dim)*msg);
   }
   template<typename MCF_SOLVER>
   static REAL arc_sign(MCF_SOLVER& mcf, const INDEX arc) 
   {
      // return -1.0 if arc has bounds [-1,0] and +1.0 if arc has bounds [0,1]
      if(mcf.upper_bound(arc) == 0.0) {
         assert(mcf.lower_bound(arc) == 1.0);
         return -1.0;
      } else {
         assert(mcf.upper_bound(arc) == 1.0);
         assert(mcf.lower_bound(arc) == 0.0);
         return 1.0;
      }
   }

   template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void construct_constraints(SOLVER& s, LEFT_FACTOR& l, typename SOLVER::vector unary_variables, RIGHT_FACTOR& r)
   { 
     assert(false);
   }

private:
   const INDEX first_arc_;
   const INDEX no_arcs_;
};

} // end namespace LPMP
