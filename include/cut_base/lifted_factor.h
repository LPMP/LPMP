#ifndef LP_MP_LIFTED_MULTICUT_FACTOR_H
#define LP_MP_LIFTED_MULTICUT_FACTOR_H

#include <vector>
#include <limits>
#include <tuple>
#include <algorithm>

namespace LPMP {

class LiftedMulticutCutFactor : public std::vector<double> {
public:
   using repam_storage = std::vector<double>; // investigate using std::list with pool allocator. Possibly, the number of cut edges can be stored in a vector, and the lifted edges can be added via a list
   LiftedMulticutCutFactor(const std::size_t noCutEdges);

   double LowerBound() const;
   double EvaluatePrimal() const;

   void IncreaseLifted();

   double CutEdgeBreakpoint(const std::size_t c) const;
   double LiftedEdgeBreakpoint(const std::size_t c) const;
   void update_cut_edge_contrib(const std::size_t c, const double msg);
   void update_lifted_edge_contrib(const std::size_t c, const double msg);

   // statistics with which evaluation is fast.
   double& CutEdgeContrib() { return cutEdgeContrib_; }
   double& LiftedEdgeContrib() { return liftedEdgeContrib_; }
   double& LiftedEdgeForcedContrib() { return liftedEdgeForcedContrib_; }
   double& MaxCutEdgeVal() { return maxCutEdgeVal_; }

   std::size_t NoLiftedEdges() const { return noLiftedEdges_; }
   std::size_t NoCutEdges() const { return noCutEdges_; }

   void init_primal() {}
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<repam_storage*>(this) ); ar( maxCutEdgeVal_, cutEdgeContrib_, liftedEdgeContrib_, liftedEdgeForcedContrib_ ); } 
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

   auto export_variables() { return std::tie( *static_cast<repam_storage*>(this) ); }

   template<typename EXTERNAL_SOLVER, typename VECTOR>
   void construct_constraints(EXTERNAL_SOLVER& s, VECTOR vars);

   template<typename EXTERNAL_SOLVER, typename VECTOR>
   void convert_primal(EXTERNAL_SOLVER& s, VECTOR vars);

private:
   std::vector<char> primal_;
   // edges are arranged as follows: first come the lifted edges, then the original ones.
   // do zrobienia: make const
   std::size_t noLiftedEdges_; // number of lifted edges that have endpoints in different components
   std::size_t noCutEdges_; // number of cut edges in the original graph. Possibly can be automatically inferred from this->size() - noLiftedEdges_

   // hold these quantities explicitly to avoid having to recompute them after every update -> constant time operations
   double maxCutEdgeVal_;
   double cutEdgeContrib_;
   double liftedEdgeContrib_;
   double liftedEdgeForcedContrib_;
};

template<typename EXTERNAL_SOLVER, typename VECTOR>
void LiftedMulticutCutFactor::construct_constraints(EXTERNAL_SOLVER& s, VECTOR vars)
{
   auto cut_active = s.max(vars.begin(), vars.begin() + noCutEdges_);
   for(std::size_t i=0; i<noLiftedEdges_; ++i) {
      // if lifted edge is cut, at least one cut edge must be on as well
      s.add_implication( vars[noCutEdges_ + i], cut_active );
   }
}

template<typename EXTERNAL_SOLVER, typename VECTOR>
void LiftedMulticutCutFactor::convert_primal(EXTERNAL_SOLVER& s, VECTOR vars)
{
   for(std::size_t i=0; i<vars.size(); ++i) {
      primal_[i] = s.solution(vars[i]);
   } 
}

} // end namespace LP_MP

#endif // LP_MP_LIFTED_MULTICUT_FACTOR_H
