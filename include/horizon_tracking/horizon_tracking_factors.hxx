#ifndef LPMP_HORIZON_TRACKING_FACTORS_HXX
#define LPMP_HORIZON_TRACKING_FACTORS_HXX

#include "horizon_tracking_util.hxx"
#include "vector.hxx"
#include <queue>
#include <limits>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"
#include <unordered_set>
#include <cmath>
#include "omp.h"

namespace LPMP {
class max_potential_on_nodes {
struct MaxPotentialInUnary { REAL value; INDEX nodeIndex; INDEX labelIndex; };
public:
    max_potential_on_nodes() {}
    max_potential_on_nodes(const std::vector<Marginals>& marginals)
    {
        for(INDEX currentNodeIndex = 0; currentNodeIndex < marginals.size(); ++currentNodeIndex) {
            for(INDEX currentLabel = 0; currentLabel < marginals[currentNodeIndex].Get().size(); ++currentLabel ) {
                MaxPotentials.push_back( { marginals[currentNodeIndex].Get(currentLabel).MaxCost, currentNodeIndex, currentLabel } );
            }
        }
        SortingOrder = GetMaxPotsSortingOrder(MaxPotentials);
    }

    std::vector<INDEX> ComputeBestLabels(const std::vector<Marginals>& marginals) const {
        INDEX numCovered = 0;
        INDEX numNodes = marginals.size();
        std::vector<bool> coveredNodes(numNodes, false);
        REAL s = 0;
        std::vector<REAL> l(numNodes, std::numeric_limits<REAL>::max());
        double bestObjective = INFINITY;
        std::vector<INDEX> bestLabelsForNodes(numNodes);
        std::vector<INDEX> lablesForNodes(numNodes);
        for(const auto& currentElementToInsert : SortingOrder) {
            const INDEX currentNodeIndex = MaxPotentials[currentElementToInsert].nodeIndex;
            const INDEX currentLabelIndex = MaxPotentials[currentElementToInsert].labelIndex;
            const REAL currentLinearCost = marginals[currentNodeIndex].Get(currentLabelIndex).LinearCost;
            if (currentLinearCost == std::numeric_limits<REAL>::max())
                continue; // infeasible label.
            const REAL currentMaxCost =  MaxPotentials[currentElementToInsert].value;
            assert(currentMaxCost == marginals[currentNodeIndex].Get(currentLabelIndex).MaxCost);
            // If the edge is not yet covered:
            if (!coveredNodes[currentNodeIndex])  {
                coveredNodes[currentNodeIndex] = true;
                numCovered++;
                s += currentLinearCost;
                l[currentNodeIndex] = currentLinearCost;
                lablesForNodes[currentNodeIndex] = currentLabelIndex;
            }
            // If edge has been added, but current label has lower linear cost. We have two scenarios:
            // 1. Graph has not been covered completely in which case we want to increase our max pot. threshold anyway. Thus, if we are gaining 
            //      an improvement in linear cost take it.
            // 2. Graph is covered completely, but we want to see adding current label can possibly give us any benefit, which will be 
            //    checked by the 3rd if condition.
            if (currentLinearCost < l[currentNodeIndex]) {
                s = s - l[currentNodeIndex] + currentLinearCost;
                l[currentNodeIndex] = currentLinearCost;
                lablesForNodes[currentNodeIndex] = currentLabelIndex;
            }

            if (numCovered == numNodes && bestObjective > s + currentMaxCost) {
                // Found another solution which is better than the previous one, in which case mark current solution as the best so far.
                bestObjective = s + currentMaxCost;
                bestLabelsForNodes = lablesForNodes;
            }
        }
        assert(std::abs(bestObjective - CostOfLabelling(bestLabelsForNodes, marginals)) <= eps);
        BestCost = bestObjective;
        return bestLabelsForNodes;
    }

    REAL CostOfLabelling(const std::vector<INDEX>& labels, const std::vector<Marginals>& marginals) const {
        REAL linearCost = 0;
        REAL maxCost = std::numeric_limits<REAL>::lowest();
        for(INDEX currentNodeIndex = 0; currentNodeIndex < marginals.size(); ++currentNodeIndex) {
            maxCost = std::max(maxCost, marginals[currentNodeIndex].Get(labels[currentNodeIndex]).MaxCost);
            linearCost += marginals[currentNodeIndex].Get(labels[currentNodeIndex]).LinearCost;
        }
        return maxCost + linearCost;
    }

    REAL GetBestCost() const {return BestCost;}

private:
    std::vector<MaxPotentialInUnary> MaxPotentials;
    std::vector<INDEX> SortingOrder;
    mutable REAL BestCost;
  
    std::vector<INDEX> GetMaxPotsSortingOrder(const std::vector<MaxPotentialInUnary>& pots) const {
        std::vector<INDEX> idx(pots.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(),
            [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});
        return idx;
    }
};

class max_potential_on_multiple_chains {
struct MaxPotentialInChain { INDEX Edge; INDEX L1; INDEX L2; REAL Value; };
private:
    mutable std::vector<three_dimensional_variable_array<REAL>> LinearPotentials;
    std::vector<three_dimensional_variable_array<REAL>> MaxPotentials;
    std::vector<std::vector<INDEX>> NumLabels;
    INDEX NumChains;
    std::vector<INDEX> NumNodes;
    mutable max_potential_on_nodes UnarySolver;
    mutable bool UnarySolverInitialized = false;

    mutable std::vector<bool> MarginalsValid;
    mutable bool SolutionValid;
    mutable two_dim_variable_array<INDEX> Solution;
    mutable std::vector<Marginals> MarginalsChains;
    two_dim_variable_array<MaxPotentialInChain> MaxPotentialsOfChains;
    two_dim_variable_array<INDEX> MaxPotentialsOfChainsOrder;
    mutable std::vector<INDEX> BestChainMarginalIndices;
public:
    max_potential_on_multiple_chains(std::vector<three_dimensional_variable_array<REAL>>& linearPotentials,
                                     std::vector<three_dimensional_variable_array<REAL>>& maxPotentials, 
                                     std::vector<std::vector<INDEX>>& numLabels) {
        LinearPotentials = linearPotentials;
        MaxPotentials = maxPotentials;
        NumLabels = numLabels;
        assert(LinearPotentials.size() == MaxPotentials.size());
        assert(MaxPotentials.size() == NumLabels.size());
        NumChains = NumLabels.size();
        for (INDEX c = 0; c < NumChains; c++) { NumNodes.push_back(NumLabels[c].size()); }
        Solution.resize(NumNodes.begin(), NumNodes.end(), std::numeric_limits<INDEX>::max());
        MarginalsChains.resize(NumChains);
        init_primal();
        std::vector<std::vector<MaxPotentialInChain>> maxPotentialsChains(NumChains); 
        std::vector<std::vector<INDEX>> maxPotentialsChainsOrder(NumChains); 
        for (INDEX c = 0; c < NumChains; c++) {
            for (INDEX e = 0; e < MaxPotentials[c].size(); e++) {
                for (INDEX l1 = 0; l1 < MaxPotentials[c].dim2(e); l1++) {
                    for (INDEX l2 = 0; l2 < MaxPotentials[c].dim3(e); l2++) {
                        maxPotentialsChains[c].push_back({e, l1, l2, MaxPotentials[c](e, l1, l2)});
                    }
                }
            }
            maxPotentialsChainsOrder[c] = GetMaxPotentialSortingOrder(maxPotentialsChains[c]);
        }
        MaxPotentialsOfChains = two_dim_variable_array<MaxPotentialInChain>(maxPotentialsChains);
        MaxPotentialsOfChainsOrder = two_dim_variable_array<INDEX>(maxPotentialsChainsOrder);
        BestChainMarginalIndices.resize(NumChains);
        MarginalsValid.resize(NumChains);
        init_primal();
    }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }
    auto export_variables() { return std::tie( ); }
    template<typename EXTERNAL_SOLVER> void construct_constraints(EXTERNAL_SOLVER& s) const { assert(false); }
    template<typename EXTERNAL_SOLVER> void convert_primal(EXTERNAL_SOLVER& s) { assert(false); } 

    void init_primal() {
        std::fill(MarginalsValid.begin(), MarginalsValid.end(), false);
        for (INDEX c = 0; c < NumChains; c++) {
            std::fill(Solution[c].begin(), Solution[c].end(), std::numeric_limits<INDEX>::max());
        }
        std::fill(BestChainMarginalIndices.begin(), BestChainMarginalIndices.end(), std::numeric_limits<INDEX>::max());
        SolutionValid = false;
    }

    REAL LowerBound() const {
        if (*std::max_element(BestChainMarginalIndices.begin(), BestChainMarginalIndices.end()) >= std::numeric_limits<INDEX>::max())
            BestChainMarginalIndices = Solve(); 

        return UnarySolver.CostOfLabelling(BestChainMarginalIndices, MarginalsChains);
    }

    REAL EvaluatePrimal() const {
        max_linear_costs obj = ComputeChainLabellingObjective(Solution);
        assert(std::max(obj.MaxCost, obj.LinearCost) < std::numeric_limits<REAL>::max());
        return obj.TotalCost();
    }

    void MaximizePotentialAndComputePrimal() { 
        bool bestIndicesValid = *std::max_element(BestChainMarginalIndices.begin(), BestChainMarginalIndices.end()) < std::numeric_limits<INDEX>::max();
        if (bestIndicesValid && SolutionValid) return;
        if (bestIndicesValid) {
            Solution = ComputeLabelling(BestChainMarginalIndices);
            SolutionValid = true;
            return; 
        }
        if (!bestIndicesValid && !SolutionValid) { 
            BestChainMarginalIndices = Solve(); 
            Solution = ComputeLabelling(BestChainMarginalIndices);
            assert(std::abs(UnarySolver.CostOfLabelling(BestChainMarginalIndices, MarginalsChains) - ComputeChainLabellingObjective(Solution).TotalCost()) <= eps);
            SolutionValid = true;
            return;
        }
        if (!bestIndicesValid && !SolutionValid) {
            assert(false);
        }
    }

    REAL GetLinearPotential(INDEX chainIndex, chain_edge e) const {
        assert(chainIndex < NumChains); 
        return LinearPotentials[chainIndex](e.n1, e.l1, e.l2); 
    }

    void SetLinearPotential(INDEX chainIndex, chain_edge e, REAL val) const {
        assert(chainIndex < NumChains); 
        if (LinearPotentials[chainIndex](e.n1, e.l1, e.l2) == val) return;
        LinearPotentials[chainIndex](e.n1, e.l1, e.l2) = val;
        MarginalsValid[chainIndex] = false;
        SolutionValid = false;
        std::fill(BestChainMarginalIndices.begin(), BestChainMarginalIndices.end(), std::numeric_limits<INDEX>::max());
    }

    INDEX NumNodeLabels(INDEX chainIndex, INDEX nodeIndex) const { return NumLabels[chainIndex][nodeIndex]; }

    INDEX GetSolution(INDEX chainIndex, INDEX nodeIndex) const { return Solution[chainIndex][nodeIndex]; }
    void SetSolution(INDEX chainIndex, INDEX nodeIndex, INDEX val) const {
        if (Solution[chainIndex][nodeIndex] == val)
            return;
                
        SolutionValid = false; // invalidate the solution computed through the factor itself, assuming this call is made from outside
        Solution[chainIndex][nodeIndex] = val; 
    }

    std::vector<REAL> ComputeMessageForEdge(INDEX chain, INDEX e, REAL OMEGA = 0.9) const {
        REAL lb = LowerBound(); // Get Lower Bound and Populate Marginals of all Chains (if invalid).
        std::vector<REAL> message(LinearPotentials[chain].dim2(e) * LinearPotentials[chain].dim3(e));
        std::vector<Marginals> leftNodeLeftMarginals(LinearPotentials[chain].dim2(e));
        std::vector<Marginals> rightNodeRightMarginals(LinearPotentials[chain].dim3(e));
        for (INDEX l1 = 0; l1 < leftNodeLeftMarginals.size(); l1++) 
            ComputeChainMarginals<false, true>(leftNodeLeftMarginals[l1], chain, e, l1);

        for (INDEX l2 = 0; l2 < rightNodeRightMarginals.size(); l2++) 
            ComputeChainMarginals<true, true>(rightNodeRightMarginals[l2], chain, e + 1, l2);

        MarginalsValid[chain] = false;

        INDEX i = 0;
        for (INDEX l1 = 0; l1 < LinearPotentials[chain].dim2(e); l1++) {
            for (INDEX l2 = 0; l2 < LinearPotentials[chain].dim3(e); l2++) {
                // Copy the linear costs of current edge into marginals of current chain, as the original 
                // ones are not of any use now and will need to be updated anyway.
                MarginalsChains[chain] = MergeNodeMarginals(leftNodeLeftMarginals[l1], rightNodeRightMarginals[l2],
                                                             {MaxPotentials[chain](e, l1, l2), LinearPotentials[chain](e, l1, l2)});
                UnarySolver.ComputeBestLabels(MarginalsChains); // TO DO: Do not need to store labels.
                REAL minMarginal = UnarySolver.GetBestCost();

                message[i] = OMEGA * (minMarginal - lb);
                if (message[i] < 0 && message[i] > -eps)
                    message[i] = 0; // get rid of numerical errors. 
                assert(message[i] >= 0);
                i++;
            }
        }
        assert(*std::min_element(message.begin(), message.end()) <= eps); // one edge should be a part of lower bound solution.
        return message;
    }
 
private:
    std::vector<INDEX> Solve() const {
#pragma omp parallel for schedule(dynamic,3)
        for (INDEX c = 0; c < NumChains; c++) {
            if (MarginalsValid[c]) continue;
            ComputeChainMarginals(MarginalsChains[c], c);
            MarginalsValid[c] = true;
        }
        if (!UnarySolverInitialized) {
            UnarySolver = max_potential_on_nodes(MarginalsChains);
            UnarySolverInitialized = true;
        }
        return UnarySolver.ComputeBestLabels(MarginalsChains);
    }

    two_dim_variable_array<INDEX> ComputeLabelling(const std::vector<INDEX>& marginalIndices) const { 
        two_dim_variable_array<INDEX> solution(NumNodes.begin(), NumNodes.end(), std::numeric_limits<INDEX>::max());
        for (INDEX c = 0; c < NumChains; c++) {
            solution[c] = ComputeLabellingForOneChain(c, MarginalsChains[c].Get(marginalIndices[c]).MaxCost);
        }
        return solution;
    }

    std::vector<INDEX> ComputeLabellingForOneChain(INDEX chainIndex, REAL maxPotentialThresh) const {
        shortest_distance_calculator<true> distCalc(LinearPotentials[chainIndex], MaxPotentials[chainIndex], NumLabels[chainIndex]);
        distCalc.CalculateDistances(maxPotentialThresh);
        return distCalc.ShortestPath(maxPotentialThresh);
    }

    max_linear_costs ComputeChainLabellingObjective(const two_dim_variable_array<INDEX>& chainsLabelling) const {
        REAL maxPotValue = std::numeric_limits<REAL>::lowest();
        REAL linearCost = 0;
        for (INDEX c = 0; c < NumChains; c++) {
            for (INDEX n1 = 0; n1 < LinearPotentials[c].dim1(); n1++)
            {
                if (chainsLabelling[c][n1] == std::numeric_limits<INDEX>::max() || 
                    chainsLabelling[c][n1 + 1] == std::numeric_limits<INDEX>::max())
                    return {std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max()}; // Solution not ready yet

                maxPotValue = std::max(maxPotValue, MaxPotentials[c](n1, chainsLabelling[c][n1], chainsLabelling[c][n1 + 1]));
                linearCost += LinearPotentials[c](n1, chainsLabelling[c][n1], chainsLabelling[c][n1 + 1]);
            }
        }
        return {maxPotValue, linearCost};
    }

    template <bool doForward = true, bool directionalNodeMarginal = false>
    void ComputeChainMarginals(Marginals& marginals, const INDEX chainIndex, INDEX n = 0, INDEX l = 0) const
    {
        shortest_distance_calculator<doForward, directionalNodeMarginal> distCalc
                                    (LinearPotentials[chainIndex], MaxPotentials[chainIndex], NumLabels[chainIndex], n, l);

        for (const auto& e : MaxPotentialsOfChainsOrder[chainIndex]) {
            const auto n1 = MaxPotentialsOfChains[chainIndex][e].Edge;
            const auto l1 = MaxPotentialsOfChains[chainIndex][e].L1;
            const auto l2 = MaxPotentialsOfChains[chainIndex][e].L2;
            const auto bottleneckCost = MaxPotentialsOfChains[chainIndex][e].Value;
               
            distCalc.AddEdgeWithUpdate(n1, l1, l2, bottleneckCost);

            REAL l = distCalc.ShortestDistance();
            marginals.insert({bottleneckCost, l}); //storing infinities as well, to ensure consistency with min-marginal computation.
        }
        marginals.Populated();
    }

    std::vector<INDEX> GetMaxPotentialSortingOrder(const std::vector<MaxPotentialInChain>& pots) const {
        std::vector<INDEX> idx(pots.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&pots](INDEX i1, INDEX i2) {return pots[i1].Value < pots[i2].Value;});
        return idx;
    }

    Marginals MergeNodeMarginals(const Marginals& leftNode,const  Marginals& rightNode, max_linear_costs edgeCost) const {
        Marginals merged;
        assert(leftNode.size() == rightNode.size());
        for (INDEX i = 0; i < leftNode.size(); i++) {
            const REAL leftMaxCost = leftNode.Get(i).MaxCost;
            const REAL rightMaxCost = rightNode.Get(i).MaxCost;
            assert(leftMaxCost == rightMaxCost);
            const REAL leftLinearCost = leftNode.Get(i).LinearCost;
            const REAL rightLinearCost = rightNode.Get(i).LinearCost;
            if (std::max(leftLinearCost, rightLinearCost) >= std::numeric_limits<REAL>::max()  // no path exists at this bottleneck threshold
                || leftMaxCost < edgeCost.MaxCost) // below the bottleneck value of edge.
                merged.insert({leftMaxCost, std::numeric_limits<REAL>::max()});
            else
                merged.insert({leftMaxCost, leftLinearCost + rightLinearCost + edgeCost.LinearCost});
        }
        merged.Populated();
        return merged;
    }
};

class pairwise_max_potential_on_multiple_chains_message {
    public:
        pairwise_max_potential_on_multiple_chains_message(const INDEX chain_index, const INDEX pairwise_index, // TODO: simplify constructor, remove unary_1, unary_2
                                                          const INDEX unary_1, const INDEX unary_2) :
        ChainIndex(chain_index), EdgeIndex(pairwise_index), N1(unary_1), N2(unary_2) {}

        template<typename FACTOR, typename MSG>
        void RepamRight(FACTOR& r, const MSG& msgs) const
        {
            INDEX i = 0;
            for(INDEX l1=0; l1<r.NumNodeLabels(ChainIndex, N1); ++l1) {
                for(INDEX l2=0; l2<r.NumNodeLabels(ChainIndex, N2); ++l2, ++i) {
                    chain_edge e = {EdgeIndex, l1, l2};
                    REAL currentPot = r.GetLinearPotential(ChainIndex, e);
                    r.SetLinearPotential(ChainIndex, e, currentPot + msgs[i]);
                }
            }
        }

        template<typename FACTOR, typename MSG>
        void RepamLeft(FACTOR& l, const MSG& msgs) const
        {
            INDEX c=0;
            for(INDEX i=0; i<l.dim1(); ++i) {
                for(INDEX j=0; j<l.dim2(); ++j) {
                    l.cost(i,j) += msgs[c++];
                }
            } 
        }

        template<typename LEFT_FACTOR, typename MSG>
        void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0) const
        {
            vector<REAL> m(l.size());
            INDEX c=0;
            for(INDEX i=0; i<l.dim1(); ++i) {
                for(INDEX j=0; j<l.dim2(); ++j) {
                    m[c++] = l(i,j);
                }
            }
            msg -= omega*m; 
        }

        template<typename RIGHT_FACTOR, typename MSG>
        void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0) const
        {
            std::vector<REAL> message = r.ComputeMessageForEdge(ChainIndex, EdgeIndex);
            vector<REAL> m(message.size());
            for (INDEX i = 0; i < m.size(); i++)
                m[i] = message[i];

            msg -= omega*m; 
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
        {
            bool changed_N1 = false;
            if(r.GetSolution(ChainIndex, N1) < l.dim1()) {
                changed_N1 = (r.GetSolution(ChainIndex, N1) != l.primal()[0]);
                l.primal()[0] = r.GetSolution(ChainIndex, N1);
            }

            bool changed_N2 = false;
            if(r.GetSolution(ChainIndex, N2) < l.dim2()) {
                changed_N2 = (r.GetSolution(ChainIndex, N2) != l.primal()[1]);
                l.primal()[1] = r.GetSolution(ChainIndex, N2);
            }

            return changed_N1 || changed_N2;
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
        {
            bool changed_N1 = false;
            if(l.primal()[0] < l.dim1()) {
                changed_N1 = (r.GetSolution(ChainIndex, N1) != l.primal()[0]);
                r.SetSolution(ChainIndex, N1, l.primal()[0]);
            }

            bool changed_N2 = false;
            if(l.primal()[1] < l.dim2()) {
                changed_N2 = (r.GetSolution(ChainIndex, N2) != l.primal()[1]);
                r.SetSolution(ChainIndex, N2, l.primal()[1]);
            }

            return changed_N1 || changed_N2;
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
        {
            return r.GetSolution(ChainIndex, N1) == l.primal()[0] && r.GetSolution(ChainIndex, N2) == l.primal()[1];
        } 

        template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
        void construct_constraints(SOLVER& s, LEFT_FACTOR& l, 
        typename SOLVER::vector l_left_msg_variables, typename SOLVER::vector l_right_msg_variables, 
        typename SOLVER::matrix l_pairwise_variables, RIGHT_FACTOR& r) {}

    private:
        const INDEX ChainIndex;
        const INDEX EdgeIndex;
        const INDEX N1;
        const INDEX N2;
};
}

#endif // LPMP_HORIZON_TRACKING_FACTORS_NEW_HXX
