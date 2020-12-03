#include "cut_base/lifted_factor.h"
#include <vector>
#include <cassert>

namespace LPMP {

LiftedMulticutCutFactor::LiftedMulticutCutFactor(const std::size_t noCutEdges) 
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

// the feasible set is: When the ordinary edges are all one (cut), then the lifted edges must be one as well
// if at least one ordinary edge is zero (not cut), then the lifted edges may be arbitrary
double LiftedMulticutCutFactor::LowerBound() const 
{
   assert(noCutEdges_ > 0);
   assert(noLiftedEdges_ > 0);
   return cutEdgeContrib_ + liftedEdgeContrib_ + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeVal_));


   //std::cout << "lifted factor repam: ";
   //for(std::size_t i=0; i<noCutEdges_; ++i) {
   //   std::cout << repam[i] << ", ";
   //}
   //std::cout << ";;; ";
   //for(std::size_t i=0; i<noLiftedEdges_; ++i) {
   //   std::cout << repam[i+noCutEdges_] << ", ";
   //}
   //std::cout << "\n";
   // do zrobienia: recompute all statistics here from scratch, once they are introduced

   // do zrobienia: start at index 1, initialize with zero value
   double maxCutEdgeVal = -std::numeric_limits<double>::max();
   double cutEdgeContrib = 0.0;
   for(std::size_t i=0; i<noCutEdges_; ++i) {
      maxCutEdgeVal = std::max(maxCutEdgeVal, (*this)[i]);
      cutEdgeContrib += std::min(0.0,(*this)[i]);
   }

   // do zrobienia: start at index 1, initialize with zero value
   double liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
   double liftedEdgeForcedContrib = 0.0;
   for(std::size_t i=0; i<noLiftedEdges_; ++i) {
      liftedEdgeContrib += std::min(0.0,(*this)[i + noCutEdges_]);
      liftedEdgeForcedContrib += std::max(0.0,(*this)[i + noCutEdges_]);
   }

   assert(std::abs(liftedEdgeContrib - liftedEdgeContrib_) < 1e-8);
   assert(std::abs(liftedEdgeForcedContrib - liftedEdgeForcedContrib_) < 1e-8);
   assert(std::abs(maxCutEdgeVal - maxCutEdgeVal_) < 1e-8);
   assert(std::abs(cutEdgeContrib - cutEdgeContrib_) < 1e-8);

   // all one <=> maxCutEdgeVal <= 0
   return cutEdgeContrib + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));

}

double LiftedMulticutCutFactor::EvaluatePrimal() const
{
   return std::numeric_limits<double>::infinity();
   std::size_t noCutEdgesOne = 0;
   double x = 0.0;
   for(std::size_t i=0; i<noCutEdges_; ++i) {
      noCutEdgesOne += primal_[i];
      x += primal_[i]*(*this)[i];
   }
   std::size_t noLiftedEdgesOne = 0;
   for(std::size_t i=0; i<noLiftedEdges_; ++i) {
      noLiftedEdgesOne += primal_[i + noCutEdges_];
      x += primal_[i+noCutEdges_]*(*this)[i+noCutEdges_];
   }
   if(noCutEdgesOne < noCutEdges_ || noLiftedEdgesOne == noLiftedEdges_) {
      return x;
   } else {
      assert(false);
      //std::cout << "solution infeasible: #cut edges = 1: " << noCutEdgesOne << ", #cut edges = " << noCutEdges_ << ", #lifted edges = 1: " << noLiftedEdgesOne << ", #lifted edges = " << noLiftedEdges_ << "\n\n";
      return std::numeric_limits<double>::infinity();
   }
}

// do zrobienia: use own reparametrization storage and make it bigger here as well.
void LiftedMulticutCutFactor::IncreaseLifted()
{
   ++noLiftedEdges_;
   this->push_back(0.0); 
   primal_.push_back(true);
}

// do zrobienia: make this more efficient by keeping track of needed informations
// compute by how much we can change the cut edge's cost
double LiftedMulticutCutFactor::CutEdgeBreakpoint(const std::size_t c) const 
{ 
   assert(c < noCutEdges_);

   const std::size_t edgeIndex = c;
   if((*this)[edgeIndex] < maxCutEdgeVal_ || liftedEdgeForcedContrib_ <= -maxCutEdgeVal_) { // do zrobienia: take into account std::max(0.0,...)
      return (*this)[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeVal_));
   } else { // we must recompute maxCutEdgeVal_ without the active repam value
      double maxCutEdgeValExcl = -std::numeric_limits<double>::max();
      for(std::size_t i=0; i<noCutEdges_; ++i) {
         if(i != c) {
            maxCutEdgeValExcl = std::max(maxCutEdgeValExcl, (*this)[i]);
         }
      }
      return (*this)[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeValExcl));
   }
   assert(false);

   //const std::size_t edgeIndex = c;
   //const double zeroAssignment = cutEdgeContrib - std::min(0.0,repam[edgeIndex]) + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
   //const double oneAssignment = cutEdgeContrib - std::min(0.0,repam[edgeIndex]) + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeValExcl_[edgeIndex]));
   //return oneAssignment - zeroAssignment;


   // do zrobienia: start at index 1, initialize with zero value
   /*
      double maxCutEdgeVal = -std::numeric_limits<double>::max();
      double cutEdgeContrib = 0.0;
      for(std::size_t i=0; i<noCutEdges_; ++i) {
      if(i != c) {
      maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
      cutEdgeContrib += std::min(0.0,repam[i]);
      }
      }

   // do zrobienia: start at index 1, initialize with 0 value
   double liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
   double liftedEdgeForcedContrib = 0.0;
   for(std::size_t i=0; i<noLiftedEdges_; ++i) {
   liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
   liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
   }


   // approach: compute the cost for assignment =1 and assignment -0 of current variable c and let change be difference of those two values
   const double zeroAssignment = cutEdgeContrib + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
   const double oneAssignment = cutEdgeContrib + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
   const double diff = repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
   assert(std::abs(oneAssignment - zeroAssignment - diff) < eps);
   return diff;
    */ 

   assert(false);
}

double LiftedMulticutCutFactor::LiftedEdgeBreakpoint(const std::size_t c) const
{
   assert(c < noCutEdges_+noLiftedEdges_);
   assert(c >= noCutEdges_);
   const std::size_t edgeIndex = c;

   return (*this)[edgeIndex] + std::min(0.0,maxCutEdgeVal_) + std::max(0.0,std::min(liftedEdgeForcedContrib_ - std::max(0.0,(*this)[edgeIndex]),-maxCutEdgeVal_));

   /*
      double maxCutEdgeVal = -std::numeric_limits<double>::max();
      double cutEdgeContrib = 0.0;
      for(std::size_t i=0; i<noCutEdges_; ++i) {
      maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
      cutEdgeContrib += std::min(0.0,repam[i]);
      }

      double liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      double liftedEdgeForcedContrib = 0.0;
      for(std::size_t i=0; i<noLiftedEdges_; ++i) {
      if(i + noCutEdges_ != c) {
      liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
      liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
      }
      }


      const double zeroAssignment = cutEdgeContrib + liftedEdgeContrib - std::min(0.0,maxCutEdgeVal); // this assignment allows at most noCutEdges_-1 cut edges to be active
      const double oneAssignment = cutEdgeContrib + liftedEdgeContrib + repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      const double diff = repam[edgeIndex] + std::min(0.0,maxCutEdgeVal) + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      assert(std::abs(oneAssignment - zeroAssignment - diff) < eps);
      return oneAssignment - zeroAssignment;
    */
}

// possibly update summary values implicitly when writing to them by going through a proxy object. Too complicated?
void LiftedMulticutCutFactor::update_cut_edge_contrib(const std::size_t c, const double msg)
{
   assert(c < NoCutEdges());

   CutEdgeContrib() -= std::min(0.0,(*this)[c]);
   const double prevRepam = (*this)[c];
   (*this)[c] += msg; 
   CutEdgeContrib() += std::min(0.0,(*this)[c]);

   // update maxCutEdgeVal_, if it needs to be updated.
   if((*this)[c] > MaxCutEdgeVal()) {
      MaxCutEdgeVal() = (*this)[c];
   } else if(prevRepam == MaxCutEdgeVal() && msg < 0.0) {
      // search through all indices to find possibly new maxCutEdgeVal_
      double maxCutEdgeVal = -std::numeric_limits<double>::max();
      for(std::size_t i=0; i<NoCutEdges(); ++i) {
         maxCutEdgeVal = std::max((*this)[i],maxCutEdgeVal);
      }
      MaxCutEdgeVal() = maxCutEdgeVal;
   }

}

void LiftedMulticutCutFactor::update_lifted_edge_contrib(const std::size_t c, const double msg)
{
   assert(c < NoLiftedEdges());

   LiftedEdgeContrib() -= std::min(0.0,(*this)[c]);
   LiftedEdgeForcedContrib() -= std::max(0.0,(*this)[c]);
   (*this)[c] += msg; 
   LiftedEdgeContrib() += std::min(0.0,(*this)[c]);
   LiftedEdgeForcedContrib() += std::max(0.0,(*this)[c]);
}

} // end namespace LP_MP
