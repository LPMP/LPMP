#ifndef LP_MP_LIFTED_MULTICUT_FACTOR_MESSAGES_HXX
#define LP_MP_LIFTED_MULTICUT_FACTOR_MESSAGES_HXX

#include "LP_MP.h"

namespace LP_MP {

class LiftedMulticutCutFactor : public std::vector<REAL> {
public:
   using repam_storage = std::vector<REAL>; // investigate using std::list with pool allocator. Possibly, the number of cut edges can be stored in a vector, and the lifted edges can be added via a list
   LiftedMulticutCutFactor(const INDEX noCutEdges) 
      : repam_storage(noCutEdges, 0.0),
      noCutEdges_(noCutEdges),
      noLiftedEdges_(0),
      maxCutEdgeVal_(0.0),
      cutEdgeContrib_(0.0),
      liftedEdgeContrib_(0.0),
      liftedEdgeForcedContrib_(0.0),
      primal_(noCutEdges_)
   {
   assert(false);
   }
   //LiftedMulticutCutFactor(const INDEX noCutEdges, const INDEX noLiftedEdges) : noCutEdges_(noCutEdges), noLiftedEdges_(noLiftedEdges) {}

   // the feasible set is: When the ordinary edges are all one (cut), then the lifted edges must be one as well
   // if at least one ordinary edge is zero (not cut), then the lifted edges may be arbitrary
   REAL LowerBound() const 
   {
      assert(noCutEdges_ > 0);
      assert(noLiftedEdges_ > 0);
      return cutEdgeContrib_ + liftedEdgeContrib_ + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeVal_));


      //std::cout << "lifted factor repam: ";
      //for(INDEX i=0; i<noCutEdges_; ++i) {
      //   std::cout << repam[i] << ", ";
      //}
      //std::cout << ";;; ";
      //for(INDEX i=0; i<noLiftedEdges_; ++i) {
      //   std::cout << repam[i+noCutEdges_] << ", ";
      //}
      //std::cout << "\n";
      // do zrobienia: recompute all statistics here from scratch, once they are introduced

      // do zrobienia: start at index 1, initialize with zero value
      REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
      REAL cutEdgeContrib = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         maxCutEdgeVal = std::max(maxCutEdgeVal, (*this)[i]);
         cutEdgeContrib += std::min(0.0,(*this)[i]);
      }

      // do zrobienia: start at index 1, initialize with zero value
      REAL liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      REAL liftedEdgeForcedContrib = 0.0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         liftedEdgeContrib += std::min(0.0,(*this)[i + noCutEdges_]);
         liftedEdgeForcedContrib += std::max(0.0,(*this)[i + noCutEdges_]);
      }

      assert(std::abs(liftedEdgeContrib - liftedEdgeContrib_) < eps);
      assert(std::abs(liftedEdgeForcedContrib - liftedEdgeForcedContrib_) < eps);
      assert(std::abs(maxCutEdgeVal - maxCutEdgeVal_) < eps);
      assert(std::abs(cutEdgeContrib - cutEdgeContrib_) < eps);

      // all one <=> maxCutEdgeVal <= 0
      return cutEdgeContrib + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));

   }

   REAL EvaluatePrimal() const
   {
      return std::numeric_limits<REAL>::infinity();
      INDEX noCutEdgesOne = 0;
      REAL x = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         noCutEdgesOne += primal_[i];
         x += primal_[i]*(*this)[i];
      }
      INDEX noLiftedEdgesOne = 0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         noLiftedEdgesOne += primal_[i + noCutEdges_];
         x += primal_[i+noCutEdges_]*(*this)[i+noCutEdges_];
      }
      if(noCutEdgesOne < noCutEdges_ || noLiftedEdgesOne == noLiftedEdges_) {
         return x;
      } else {
         assert(false);
         //std::cout << "solution infeasible: #cut edges = 1: " << noCutEdgesOne << ", #cut edges = " << noCutEdges_ << ", #lifted edges = 1: " << noLiftedEdgesOne << ", #lifted edges = " << noLiftedEdges_ << "\n\n";
         return std::numeric_limits<REAL>::infinity();
      }
   }

   // do zrobienia: use own reparametrization storage and make it bigger here as well.
   void IncreaseLifted()
   {
      ++noLiftedEdges_;
      this->push_back(0.0); 
      primal_.push_back(true);
   }

   // do zrobienia: make this more efficient by keeping track of needed informations
   // compute by how much we can change the cut edge's cost
   REAL CutEdgeBreakpoint(const INDEX c) const 
   { 
      assert(c < noCutEdges_);

      const INDEX edgeIndex = c;
      if((*this)[edgeIndex] < maxCutEdgeVal_ || liftedEdgeForcedContrib_ <= -maxCutEdgeVal_) { // do zrobienia: take into account std::max(0.0,...)
         return (*this)[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeVal_));
      } else { // we must recompute maxCutEdgeVal_ without the active repam value
         REAL maxCutEdgeValExcl = -std::numeric_limits<REAL>::max();
         for(INDEX i=0; i<noCutEdges_; ++i) {
            if(i != c) {
               maxCutEdgeValExcl = std::max(maxCutEdgeValExcl, (*this)[i]);
            }
         }
         return (*this)[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeValExcl));
      }
      assert(false);

      //const INDEX edgeIndex = c;
      //const REAL zeroAssignment = cutEdgeContrib - std::min(0.0,repam[edgeIndex]) + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
      //const REAL oneAssignment = cutEdgeContrib - std::min(0.0,repam[edgeIndex]) + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeValExcl_[edgeIndex]));
      //return oneAssignment - zeroAssignment;
 

      // do zrobienia: start at index 1, initialize with zero value
      /*
      REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
      REAL cutEdgeContrib = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         if(i != c) {
            maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
            cutEdgeContrib += std::min(0.0,repam[i]);
         }
      }

      // do zrobienia: start at index 1, initialize with 0 value
      REAL liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      REAL liftedEdgeForcedContrib = 0.0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
         liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
      }


      // approach: compute the cost for assignment =1 and assignment -0 of current variable c and let change be difference of those two values
      const REAL zeroAssignment = cutEdgeContrib + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
      const REAL oneAssignment = cutEdgeContrib + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      const REAL diff = repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      assert(std::abs(oneAssignment - zeroAssignment - diff) < eps);
      return diff;
      */ 
      
      assert(false);
   }

   REAL LiftedEdgeBreakpoint(const INDEX c) const
   {
      assert(c < noCutEdges_+noLiftedEdges_);
      assert(c >= noCutEdges_);
      const INDEX edgeIndex = c;

      return (*this)[edgeIndex] + std::min(0.0,maxCutEdgeVal_) + std::max(0.0,std::min(liftedEdgeForcedContrib_ - std::max(0.0,(*this)[edgeIndex]),-maxCutEdgeVal_));

      /*
      REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
      REAL cutEdgeContrib = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
         cutEdgeContrib += std::min(0.0,repam[i]);
      }

      REAL liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      REAL liftedEdgeForcedContrib = 0.0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         if(i + noCutEdges_ != c) {
            liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
            liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
         }
      }


      const REAL zeroAssignment = cutEdgeContrib + liftedEdgeContrib - std::min(0.0,maxCutEdgeVal); // this assignment allows at most noCutEdges_-1 cut edges to be active
      const REAL oneAssignment = cutEdgeContrib + liftedEdgeContrib + repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      const REAL diff = repam[edgeIndex] + std::min(0.0,maxCutEdgeVal) + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      assert(std::abs(oneAssignment - zeroAssignment - diff) < eps);
      return oneAssignment - zeroAssignment;
      */
   }

   // possibly update summary values implicitly when writing to them by going through a proxy object. Too complicated?
   void update_cut_edge_contrib(const INDEX c, const REAL msg)
   {
      assert(c < NoCutEdges());

      CutEdgeContrib() -= std::min(0.0,(*this)[c]);
      const REAL prevRepam = (*this)[c];
      (*this)[c] += msg; 
      CutEdgeContrib() += std::min(0.0,(*this)[c]);
 
      // update maxCutEdgeVal_, if it needs to be updated.
      if((*this)[c] > MaxCutEdgeVal()) {
         MaxCutEdgeVal() = (*this)[c];
      } else if(prevRepam == MaxCutEdgeVal() && msg < 0.0) {
         // search through all indices to find possibly new maxCutEdgeVal_
         REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
         for(INDEX i=0; i<NoCutEdges(); ++i) {
            maxCutEdgeVal = std::max((*this)[i],maxCutEdgeVal);
         }
         MaxCutEdgeVal() = maxCutEdgeVal;
      }

   }

   void update_lifted_edge_contrib(const INDEX c, const REAL msg)
   {
      assert(c < NoLiftedEdges());

      LiftedEdgeContrib() -= std::min(0.0,(*this)[c]);
      LiftedEdgeForcedContrib() -= std::max(0.0,(*this)[c]);
      (*this)[c] += msg; 
      LiftedEdgeContrib() += std::min(0.0,(*this)[c]);
      LiftedEdgeForcedContrib() += std::max(0.0,(*this)[c]);
   }

   // statistics with which evaluation is fast.
   REAL& CutEdgeContrib() { return cutEdgeContrib_; }
   REAL& LiftedEdgeContrib() { return liftedEdgeContrib_; }
   REAL& LiftedEdgeForcedContrib() { return liftedEdgeForcedContrib_; }
   REAL& MaxCutEdgeVal() { return maxCutEdgeVal_; }

   INDEX NoLiftedEdges() const { return noLiftedEdges_; }
   INDEX NoCutEdges() const { return noCutEdges_; }

   void init_primal() {}
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<repam_storage*>(this) ); ar( maxCutEdgeVal_, cutEdgeContrib_, liftedEdgeContrib_, liftedEdgeForcedContrib_ ); } 
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

   auto export_variables() { return std::tie( *static_cast<repam_storage*>(this) ); }

   template<typename EXTERNAL_SOLVER, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, VECTOR vars)
   {
	   auto cut_active = s.max(vars.begin(), vars.begin() + noCutEdges_);
	   for(std::size_t i=0; i<noLiftedEdges_; ++i) {
		   // if lifted edge is cut, at least one cut edge must be on as well
		   s.add_implication( vars[noCutEdges_ + i], cut_active );
	   }
   }

   template<typename EXTERNAL_SOLVER, typename VECTOR>
   void convert_primal(EXTERNAL_SOLVER& s, VECTOR vars)
   {
	   for(INDEX i=0; i<vars.size(); ++i) {
		   primal_[i] = s.solution(vars[i]);
	   } 
   }
private:
   std::vector<char> primal_;
   // edges are arranged as follows: first come the lifted edges, then the original ones.
   // do zrobienia: make const
   INDEX noLiftedEdges_; // number of lifted edges that have endpoints in different components
   INDEX noCutEdges_; // number of cut edges in the original graph. Possibly can be automatically inferred from this->size() - noLiftedEdges_

   // hold these quantities explicitly to avoid having to recompute them after every update -> constant time operations
   REAL maxCutEdgeVal_;
   REAL cutEdgeContrib_;
   REAL liftedEdgeContrib_;
   REAL liftedEdgeForcedContrib_;
};

// message between edge of the base graph and lifted multicut factor
class CutEdgeLiftedMulticutFactorMessage {
public:
   CutEdgeLiftedMulticutFactorMessage(const INDEX i) : i_(i) {}

   constexpr static INDEX size() { return 1; }

   template<typename RIGHT_FACTOR, typename G2>
   void send_message_to_left(const RIGHT_FACTOR& r, G2& msg, const REAL omega) 
   {
      msg[0] -= omega*r.CutEdgeBreakpoint(i_);
   }

   template<typename LEFT_FACTOR, typename G3>
   void send_message_to_right(const LEFT_FACTOR& l, G3& msg, const REAL omega)
   {
      msg[0] -= omega*l[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      assert(repamPot.size() == 1);
      repamPot[0] += msg;
   }

   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      //const REAL cutEdgeContribDiff = -std::min(0.0,repamPot[i_]) + std::min(0.0,msg);
      //repamPot.GetFactor()->CutEdgeContrib() += cutEdgeContribDiff;
      
      repamPot.update_cut_edge_contrib(i_, msg);
      return;

      repamPot.CutEdgeContrib() -= std::min(0.0,repamPot[i_]);
      const REAL prevRepam = repamPot[i_];
      repamPot[i_] += msg; 
      repamPot.CutEdgeContrib() += std::min(0.0,repamPot[i_]);
 
      // update maxCutEdgeVal_, if it needs to be updated.
      if(repamPot[i_] > repamPot.MaxCutEdgeVal()) {
         repamPot.MaxCutEdgeVal() = repamPot[i_];
      } else if(prevRepam == repamPot.MaxCutEdgeVal() && msg < 0.0) {
         // search through all indices to find possibly new maxCutEdgeVal_
         REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
         for(INDEX i=0; i<repamPot.NoCutEdges(); ++i) {
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
       //std::cout << "cut edge index = " << i_ << ", primal value = " << INDEX(left[0]) << "\n";
       right[i_] = left[0];
    }
    */

   template<typename EXTERNAL_SOLVER, typename EDGE_FACTOR, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, EDGE_FACTOR& l, VECTOR left_vars, LiftedMulticutCutFactor& r, VECTOR right_vars)
   {
	   s.make_equal(left_vars[0], right_vars[i_]); 
   }

private:
   INDEX i_; // index of cut edge in the multicut factor
};

// message between lifted edge and lifted multicut factor
class LiftedEdgeLiftedMulticutFactorMessage {
public:
   LiftedEdgeLiftedMulticutFactorMessage(const INDEX i) : i_(i) {}

   constexpr static INDEX size() { return 1; }

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega) 
   {
      //assert(msg.size() == 1);
      msg[0] -= omega*r.LiftedEdgeBreakpoint(i_);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
   {
      //assert(msg.size() == 1);
      msg[0] -= omega*l[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      assert(repamPot.size() == 1);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      //const REAL liftedEdgeContribDiff = -std::min(0.0,repamPot[i_]) + std::min(0.0,msg);
      //repamPot.GetFactor()->LiftedEdgeContrib() += liftedEdgeContribDiff;
      //const REAL liftedEdgeForcedContribDiff = -std::max(0.0,repamPot[i_]) + std::max(0.0,msg);
      //repamPot.GetFactor()->LiftedEdgeForcedContrib() += liftedEdgeForcedContribDiff;
      repamPot.update_lifted_edge_contrib(i_, msg);
      return;

      repamPot.LiftedEdgeContrib() -= std::min(0.0,repamPot[i_]);
      repamPot.LiftedEdgeForcedContrib() -= std::max(0.0,repamPot[i_]);
      repamPot[i_] += msg; 
      repamPot.LiftedEdgeContrib() += std::min(0.0,repamPot[i_]);
      repamPot.LiftedEdgeForcedContrib() += std::max(0.0,repamPot[i_]);
   }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
    {
       assert(right[i_] == unknownState);
       //std::cout << "lifted edge index = " << i_ << ", primal value = " << INDEX(left[0]) << "\n";
       right[i_] = left[0];
    }

   template<typename EXTERNAL_SOLVER, typename EDGE_FACTOR, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, EDGE_FACTOR& l, VECTOR left_vars, LiftedMulticutCutFactor& r, VECTOR right_vars)
   {
	   s.make_equal(left_vars[0], right_vars[r.NoCutEdges() + i_]); 
   }

private:
   INDEX i_; // index of lifted edge in the multicut factor. Must be bigger than number of cut edges, smaller than size of cut factor
};

} // end namespace LP_MP

#endif // LP_MP_LIFTED_MULTICUT_FACTOR_MESSAGES_HXX

