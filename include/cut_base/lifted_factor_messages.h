#ifndef LP_MP_LIFTED_MULTICUT_FACTOR_MESSAGES_H
#define LP_MP_LIFTED_MULTICUT_FACTOR_MESSAGES_H

#include <cassert>

namespace LPMP {

// message between edge of the base graph and lifted multicut factor
class CutEdgeLiftedMulticutFactorMessage {
public:
   CutEdgeLiftedMulticutFactorMessage(const std::size_t i) : i_(i) {}

   constexpr static std::size_t size() { return 1; }

   template<typename RIGHT_FACTOR, typename G2>
   void send_message_to_left(const RIGHT_FACTOR& r, G2& msg, const double omega) 
   {
      msg[0] -= omega*r.CutEdgeBreakpoint(i_);
   }

   template<typename LEFT_FACTOR, typename G3>
   void send_message_to_right(const LEFT_FACTOR& l, G3& msg, const double omega)
   {
      msg[0] -= omega*l[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const double msg, const std::size_t msg_dim)
   {
      assert(msg_dim == 0);
      assert(repamPot.size() == 1);
      repamPot[0] += msg;
   }

   template<typename G>
   void RepamRight(G& repamPot, const double msg, const std::size_t msg_dim)
   {
      assert(msg_dim == 0);
      //const double cutEdgeContribDiff = -std::min(0.0,repamPot[i_]) + std::min(0.0,msg);
      //repamPot.GetFactor()->CutEdgeContrib() += cutEdgeContribDiff;
      
      repamPot.update_cut_edge_contrib(i_, msg);
      return;

      repamPot.CutEdgeContrib() -= std::min(0.0,repamPot[i_]);
      const double prevRepam = repamPot[i_];
      repamPot[i_] += msg; 
      repamPot.CutEdgeContrib() += std::min(0.0,repamPot[i_]);
 
      // update maxCutEdgeVal_, if it needs to be updated.
      if(repamPot[i_] > repamPot.MaxCutEdgeVal()) {
         repamPot.MaxCutEdgeVal() = repamPot[i_];
      } else if(prevRepam == repamPot.MaxCutEdgeVal() && msg < 0.0) {
         // search through all indices to find possibly new maxCutEdgeVal_
         double maxCutEdgeVal = -std::numeric_limits<double>::max();
         for(std::size_t i=0; i<repamPot.NoCutEdges(); ++i) {
            maxCutEdgeVal = std::max(repamPot[i],maxCutEdgeVal);
         }
         repamPot.MaxCutEdgeVal() = maxCutEdgeVal;
      }
   }

   /*
    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void
    ComputeRightFromLeftPrimal(LEFT_FACTOR* l, RIGHT_FACTOR* r)
    {
       assert(right[i_] == unknownState);
       //std::cout << "cut edge index = " << i_ << ", primal value = " << std::size_t(left[0]) << "\n";
       right[i_] = left[0];
    }
    */

   template<typename EXTERNAL_SOLVER, typename EDGE_FACTOR, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, EDGE_FACTOR& l, VECTOR left_vars, LiftedMulticutCutFactor& r, VECTOR right_vars)
   {
	   s.make_equal(left_vars[0], right_vars[i_]); 
   }

private:
   std::size_t i_; // index of cut edge in the multicut factor
};

// message between lifted edge and lifted multicut factor
class LiftedEdgeLiftedMulticutFactorMessage {
public:
   LiftedEdgeLiftedMulticutFactorMessage(const std::size_t i) : i_(i) {}

   constexpr static std::size_t size() { return 1; }

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const double omega) 
   {
      //assert(msg.size() == 1);
      msg[0] -= omega*r.LiftedEdgeBreakpoint(i_);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const double omega)
   {
      //assert(msg.size() == 1);
      msg[0] -= omega*l[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const double msg, const std::size_t msg_dim)
   {
      assert(msg_dim == 0);
      assert(repamPot.size() == 1);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const double msg, const std::size_t msg_dim)
   {
      assert(msg_dim == 0);
      //const double liftedEdgeContribDiff = -std::min(0.0,repamPot[i_]) + std::min(0.0,msg);
      //repamPot.GetFactor()->LiftedEdgeContrib() += liftedEdgeContribDiff;
      //const double liftedEdgeForcedContribDiff = -std::max(0.0,repamPot[i_]) + std::max(0.0,msg);
      //repamPot.GetFactor()->LiftedEdgeForcedContrib() += liftedEdgeForcedContribDiff;
      repamPot.update_lifted_edge_contrib(i_, msg);
      return;

      repamPot.LiftedEdgeContrib() -= std::min(0.0,repamPot[i_]);
      repamPot.LiftedEdgeForcedContrib() -= std::max(0.0,repamPot[i_]);
      repamPot[i_] += msg; 
      repamPot.LiftedEdgeContrib() += std::min(0.0,repamPot[i_]);
      repamPot.LiftedEdgeForcedContrib() += std::max(0.0,repamPot[i_]);
   }

   template<typename EXTERNAL_SOLVER, typename EDGE_FACTOR, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, EDGE_FACTOR& l, VECTOR left_vars, LiftedMulticutCutFactor& r, VECTOR right_vars)
   {
	   s.make_equal(left_vars[0], right_vars[r.NoCutEdges() + i_]); 
   }

private:
   std::size_t i_; // index of lifted edge in the multicut factor. Must be bigger than number of cut edges, smaller than size of cut factor
};

} // end namespace LP_MP

#endif // LP_MP_LIFTED_MULTICUT_FACTOR_MESSAGES_H
