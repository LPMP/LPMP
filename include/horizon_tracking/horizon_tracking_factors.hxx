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
#include <chrono>

namespace LPMP {
class max_potential_on_nodes {
struct MaxPotentialInUnary { REAL value; INDEX nodeIndex; INDEX labelIndex; };
public:

    max_potential_on_nodes() : SolveIndependently(false) {}

    max_potential_on_nodes(const std::vector<Marginals>& marginals, const bool solveIndependently = false) : SolveIndependently(solveIndependently)
    {
        for(INDEX currentNodeIndex = 0; currentNodeIndex < marginals.size(); ++currentNodeIndex) {
            for(INDEX currentLabel = 0; currentLabel < marginals[currentNodeIndex].Get().size(); ++currentLabel ) {
                MaxPotentials.push_back( { marginals[currentNodeIndex].Get(currentLabel).MaxCost, currentNodeIndex, currentLabel } );
            }
        }
        SortingOrder = GetMaxPotsSortingOrder(MaxPotentials);
    }

    template <bool computeMarginals = false>
    std::vector<INDEX> ComputeBestLabels(const std::vector<Marginals>& marginals) const {
        if (SolveIndependently)
            return ComputeBestLabelsIndependent(marginals);

        INDEX numCovered = 0;
        INDEX numNodes = marginals.size();
        if (computeMarginals) {
            INDEX size = 0;
            for (INDEX i = 0; i < marginals.size(); i++)
                size += marginals[i].size();
            AllMarginals.Reserve(size);
        }
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

            if (numCovered == numNodes) {
                if (computeMarginals) {
                    AllMarginals.insert<true>({currentMaxCost, s});
                }
                if (bestObjective > s + currentMaxCost) {
                    // Found another solution which is better than the previous one, in which case mark current solution as the best so far.
                    bestObjective = s + currentMaxCost;
                    bestLabelsForNodes = lablesForNodes;
                }
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
            if (currentNodeIndex == 0)
                maxCost = marginals[currentNodeIndex].Get(labels[currentNodeIndex]).MaxCost;
            else {
                if (!SolveIndependently)
                    maxCost = std::max(maxCost, marginals[currentNodeIndex].Get(labels[currentNodeIndex]).MaxCost);
                else 
                    maxCost += marginals[currentNodeIndex].Get(labels[currentNodeIndex]).MaxCost;
            }
            linearCost += marginals[currentNodeIndex].Get(labels[currentNodeIndex]).LinearCost;
        }
        return maxCost + linearCost;
    }

    REAL GetBestCost() const {return BestCost;}
    Marginals& GetAllMarginals() const {return AllMarginals;}

private:
    std::vector<MaxPotentialInUnary> MaxPotentials;
    std::vector<INDEX> SortingOrder;
    mutable REAL BestCost;
    mutable Marginals AllMarginals;
    bool SolveIndependently;
    std::vector<INDEX> GetMaxPotsSortingOrder(const std::vector<MaxPotentialInUnary>& pots) const {
        std::vector<INDEX> idx(pots.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(),
            [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});
        return idx;
    }

    std::vector<INDEX> ComputeBestLabelsIndependent(const std::vector<Marginals>& marginals) const {
        INDEX numNodes = marginals.size();
        std::vector<REAL> bestCostOfNodes(numNodes, std::numeric_limits<REAL>::max());
        std::vector<INDEX> bestLabelsForNodes(numNodes);
        for(const auto& currentElementToInsert : SortingOrder) {
            const INDEX currentNodeIndex = MaxPotentials[currentElementToInsert].nodeIndex;
            const INDEX currentLabelIndex = MaxPotentials[currentElementToInsert].labelIndex;
            const REAL currentLinearCost = marginals[currentNodeIndex].Get(currentLabelIndex).LinearCost;
            if (currentLinearCost == std::numeric_limits<REAL>::max())
                continue; // infeasible label.
            const REAL currentMaxCost =  MaxPotentials[currentElementToInsert].value;
            assert(currentMaxCost == marginals[currentNodeIndex].Get(currentLabelIndex).MaxCost);
            if (bestCostOfNodes[currentNodeIndex] > currentMaxCost + currentLinearCost) {
                bestCostOfNodes[currentNodeIndex] = currentMaxCost + currentLinearCost;
                bestLabelsForNodes[currentNodeIndex] = currentLabelIndex;
            }
        }
        return bestLabelsForNodes;
    }
};

class max_potential_on_multiple_chains {
struct MaxPotentialInChain { INDEX Edge; INDEX L1; INDEX L2; REAL Value; };
private:
    mutable std::vector<three_dimensional_variable_array<REAL>> LinearPotentials;
    mutable std::vector<three_dimensional_variable_array<REAL>> MaxPotentials;
    const std::vector<std::vector<INDEX>> NumLabels;
    const bool SolveChainsIndependently;
    const INDEX NumChains;
    const two_dim_variable_array<INDEX> ChainNodeToOriginalNode;
    std::vector<INDEX> NumNodes; // TODO: make const
    mutable max_potential_on_nodes UnarySolver;
    mutable bool UnarySolverInitialized = false;

    mutable std::vector<bool> MarginalsValid;
    mutable bool SolutionValid;
    mutable two_dim_variable_array<INDEX> Solution;
    mutable std::vector<Marginals> MarginalsChains;
    mutable std::vector<std::vector<MaxPotentialInChain>> MaxPotentialsOfChains;
    mutable std::vector<std::vector<INDEX>> MaxPotentialsOfChainsOrder;
    mutable std::vector<INDEX> BestChainMarginalIndices;
    mutable INDEX numEdges;
    mutable INDEX messageNormalizer;
    mutable three_dimensional_variable_array<REAL> BackwardMessageFromChain;
    mutable INDEX BackwardMessageFromChainIndex;
    mutable INDEX NumChainsReparametrizedBackwards;

public:
    // TODO: use move semantics
    max_potential_on_multiple_chains(std::vector<three_dimensional_variable_array<REAL>>& linearPotentials,
                                     std::vector<three_dimensional_variable_array<REAL>>& maxPotentials, 
                                     std::vector<std::vector<INDEX>>& numLabels,
                                     two_dim_variable_array<INDEX>& chainNodeToOriginalNode,
                                     bool solveChainsIndependently = false) 
       : LinearPotentials(linearPotentials), 
        MaxPotentials(maxPotentials),
        NumLabels(numLabels),
        NumChains(NumLabels.size()),
        ChainNodeToOriginalNode(chainNodeToOriginalNode),
        SolveChainsIndependently(solveChainsIndependently)
   {
        assert(LinearPotentials.size() == MaxPotentials.size());
        assert(MaxPotentials.size() == NumLabels.size());
        for (INDEX c = 0; c < NumChains; c++) { NumNodes.push_back(NumLabels[c].size()); }
        Solution.resize(NumNodes.begin(), NumNodes.end(), std::numeric_limits<INDEX>::max());
        MarginalsChains.resize(NumChains);
        init_primal();
        numEdges = 0;
        MaxPotentialsOfChains.resize(NumChains);
        MaxPotentialsOfChainsOrder.resize(NumChains);
        for (INDEX c = 0; c < NumChains; c++) {
            numEdges += MaxPotentials[c].size();
            for (INDEX e = 0; e < MaxPotentials[c].size(); e++) {
                for (INDEX l1 = 0; l1 < MaxPotentials[c].dim2(e); l1++) {
                    for (INDEX l2 = 0; l2 < MaxPotentials[c].dim3(e); l2++) {
                        if (std::isinf(MaxPotentials[c](e, l1, l2))) continue;
                        MaxPotentialsOfChains[c].push_back({e, l1, l2, MaxPotentials[c](e, l1, l2)});
                    }
                }
            }
            MaxPotentialsOfChainsOrder[c] = GetMaxPotentialSortingOrder(MaxPotentialsOfChains[c]);
        }
        messageNormalizer = numEdges;
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
#pragma omp parallel for
        for (INDEX c = 0; c < NumChains; c++) {
            std::fill(Solution[c].begin(), Solution[c].end(), std::numeric_limits<INDEX>::max());
        }
        std::fill(BestChainMarginalIndices.begin(), BestChainMarginalIndices.end(), std::numeric_limits<INDEX>::max());
        SolutionValid = false;
        NumChainsReparametrizedBackwards = 0;
        BackwardMessageFromChainIndex = std::numeric_limits<INDEX>::max();
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

    std::vector<REAL> GetMessageForEdge(INDEX chain, INDEX e) const {
        if (chain != BackwardMessageFromChainIndex) {
            BackwardMessageFromChain = ComputeBackwardMessagesOnChain(chain);
            BackwardMessageFromChainIndex = chain;
            NumChainsReparametrizedBackwards++;
            if (NumChainsReparametrizedBackwards == NumChains) { // reset after all chains have been iterated over.
                NumChainsReparametrizedBackwards = 0;
                BackwardMessageFromChainIndex = std::numeric_limits<INDEX>::max();
            }
        }
        std::vector<REAL> message(BackwardMessageFromChain.dim2(e) * BackwardMessageFromChain.dim3(e));
        INDEX i = 0;
        for (INDEX l1 = 0; l1 < BackwardMessageFromChain.dim2(e); l1++) {
            for (INDEX l2 = 0; l2 < BackwardMessageFromChain.dim3(e); l2++, i++) {
                message[i] = BackwardMessageFromChain(e, l1, l2);
            }
        }
        return message;
    }

    void ComputeAndSetPrimal() const {
        Solution = ComputePrimal();
    }

private:
    std::vector<INDEX> Solve() const {
#pragma omp parallel for schedule(dynamic,3)
        for (INDEX c = 0; c < NumChains; c++) {
            if (MarginalsValid[c]) continue;
            ComputeChainMarginals(MarginalsChains[c], MaxPotentials[c], LinearPotentials[c], 
            MaxPotentialsOfChains[c], MaxPotentialsOfChainsOrder[c], c);
            MarginalsValid[c] = true;
        }
        if (!UnarySolverInitialized) {
            UnarySolver = max_potential_on_nodes(MarginalsChains, SolveChainsIndependently);
            UnarySolverInitialized = true;
        }
        return UnarySolver.ComputeBestLabels(MarginalsChains);
    }

    two_dim_variable_array<INDEX> ComputeLabelling(const std::vector<INDEX>& marginalIndices) const { 
        two_dim_variable_array<INDEX> solution(NumNodes.begin(), NumNodes.end(), std::numeric_limits<INDEX>::max());
        for (INDEX c = 0; c < NumChains; c++) {
            solution[c] = ComputeLabellingForOneChain(c, MarginalsChains[c].Get(marginalIndices[c]).MaxCost, MaxPotentials[c], LinearPotentials[c]);
        }
        return solution;
    }

    template <bool useFixedNode = false>
    std::vector<INDEX> ComputeLabellingForOneChain(const INDEX chainIndex, const REAL maxPotentialThresh,
                                                    const three_dimensional_variable_array<REAL>& maxPots, 
                                                    const three_dimensional_variable_array<REAL>& linearPots,
                                                    const INDEX fixedNode = 0, const INDEX fixedNodeLabel = 0) const {
        shortest_distance_calculator<true, false, useFixedNode> distCalc(linearPots, maxPots, NumLabels[chainIndex], 0, fixedNode, fixedNodeLabel);
        distCalc.CalculateDistances(maxPotentialThresh);
        return distCalc.ShortestPath(maxPotentialThresh);
    }

    max_linear_costs ComputeChainLabellingObjective(const two_dim_variable_array<INDEX>& chainsLabelling) const {
        REAL maxPotValue = std::numeric_limits<REAL>::lowest();
        REAL linearCost = 0;
        for (INDEX c = 0; c < NumChains; c++) {
            REAL currentChainMax = std::numeric_limits<REAL>::lowest();
            for (INDEX n1 = 0; n1 < LinearPotentials[c].dim1(); n1++)
            {
                if (chainsLabelling[c][n1] == std::numeric_limits<INDEX>::max() || 
                    chainsLabelling[c][n1 + 1] == std::numeric_limits<INDEX>::max())
                    return {std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max()}; // Solution not ready yet

                currentChainMax = std::max(currentChainMax, MaxPotentials[c](n1, chainsLabelling[c][n1], chainsLabelling[c][n1 + 1]));
                linearCost += LinearPotentials[c](n1, chainsLabelling[c][n1], chainsLabelling[c][n1 + 1]);
            }
            if (c == 0) {
                maxPotValue = currentChainMax;
            } else {
                if (SolveChainsIndependently)
                    maxPotValue += currentChainMax;
                else
                    maxPotValue = std::max(maxPotValue, currentChainMax);
            }
        }
        return {maxPotValue, linearCost};
    }

    template <bool doForward = true, bool onlyComputeB = false, bool useFixedNode = false>
    REAL ComputeChainMarginals(Marginals& marginals, const three_dimensional_variable_array<REAL>& maxPotentials, 
                                const three_dimensional_variable_array<REAL>& linearPotentials, 
                                const std::vector<MaxPotentialInChain>& maxPotentials1D,
                                const std::vector<INDEX>& maxPotentials1DOrder,
                                const INDEX chainIndex, const INDEX fixedNode = 0, const INDEX fixedNodeLabel = 0) const
    {
        assert(fixedNodeLabel < NumLabels[chainIndex][fixedNode]);
        REAL optimalB = std::numeric_limits<REAL>::max();
        REAL optimalCost = std::numeric_limits<REAL>::max();
        shortest_distance_calculator<doForward, false, useFixedNode> distCalc
                                    (linearPotentials, maxPotentials, NumLabels[chainIndex], 0, fixedNode, fixedNodeLabel);

        for (const auto& e : maxPotentials1DOrder) {
            const auto n1 = maxPotentials1D[e].Edge;
            const auto l1 = maxPotentials1D[e].L1;
            const auto l2 = maxPotentials1D[e].L2;
            const auto bottleneckCost = maxPotentials1D[e].Value;
               
            distCalc.AddEdgeWithUpdate(n1, l1, l2, bottleneckCost);

            REAL l = distCalc.ShortestDistance();
            if (optimalCost > bottleneckCost + l) {
                optimalCost = bottleneckCost + l;
                optimalB = bottleneckCost;
            }
            if (!onlyComputeB)
                marginals.insert({bottleneckCost, l}); //storing infinities as well, to ensure consistency with min-marginal computation.
        }
        if (!onlyComputeB)
            marginals.Populated();

        return optimalB;
    }

    // Assumes the solution in chain 'c' should be propagated
    void PropagateChainSolutionToOtherChains(const INDEX c, std::vector<INDEX>& gridSolution, std::vector<std::vector<INDEX>>& allChainsSolution) const {
        // Set the solution in grid first:
        for (INDEX n = 0; n < NumNodes[c]; n++) {
            // Either the solution is not set, or if it is set then it must be the same as incoming chain Solution:
            assert(gridSolution[ChainNodeToOriginalNode[c][n]] == std::numeric_limits<INDEX>::max() || 
                    gridSolution[ChainNodeToOriginalNode[c][n]] == allChainsSolution[c][n]);
            if (gridSolution[ChainNodeToOriginalNode[c][n]] == allChainsSolution[c][n])
                continue;

            gridSolution[ChainNodeToOriginalNode[c][n]] = allChainsSolution[c][n];
        }

        // Now propagate the grid solution to chains:
        for (INDEX otherC = 0; otherC < NumChains; otherC++) {
            for (INDEX n = 0; n < NumNodes[otherC]; n++) {
                if (gridSolution[ChainNodeToOriginalNode[otherC][n]] == std::numeric_limits<INDEX>::max()) continue;
                
                assert(allChainsSolution[otherC][n] == std::numeric_limits<INDEX>::max() ||
                allChainsSolution[otherC][n] == gridSolution[ChainNodeToOriginalNode[otherC][n]]);

                allChainsSolution[otherC][n] = gridSolution[ChainNodeToOriginalNode[otherC][n]];
            }
        }
    }

    std::vector<std::vector<INDEX>> ComputePrimal() const {
        std::vector<std::vector<INDEX>> allChainsLabels(NumChains);
        for (INDEX c = 0; c < NumChains; c++) {
            allChainsLabels[c].resize(NumNodes[c], std::numeric_limits<INDEX>::max());
        }
        INDEX totalNumNodes = std::accumulate(NumNodes.begin(), NumNodes.end(), 0);
        std::vector<INDEX> gridSolution(totalNumNodes, std::numeric_limits<INDEX>::max());
        ChainsInfo chainInfo(ChainNodeToOriginalNode);

        // Compute solution on - and | chain intersecting at node of minimum cardinality
        std::pair<INDEX, INDEX> h_v_chains = chainInfo.GetSeedChains(NumLabels);
        INDEX hStartingChain = h_v_chains.first;
        INDEX vStartingChain = h_v_chains.second;
        INDEX fixedNodeIndexInHorizontalChains = vStartingChain - chainInfo.NumHorizontal();

        // Solve on horizontal:
        REAL optimalBH = ComputeChainMarginals<true, true>(MarginalsChains[hStartingChain], MaxPotentials[hStartingChain], LinearPotentials[hStartingChain],
                                                            MaxPotentialsOfChains[hStartingChain], MaxPotentialsOfChainsOrder[hStartingChain], hStartingChain);

        allChainsLabels[hStartingChain] = ComputeLabellingForOneChain(hStartingChain, optimalBH, MaxPotentials[hStartingChain], LinearPotentials[hStartingChain]);   
        if (NumChains == 1)
            return allChainsLabels;
            
        PropagateChainSolutionToOtherChains(hStartingChain, gridSolution, allChainsLabels);

        // Solve on vertical:
        REAL optimalBV = ComputeChainMarginals<true, true>(MarginalsChains[vStartingChain], MaxPotentials[vStartingChain], LinearPotentials[vStartingChain],
                                                            MaxPotentialsOfChains[vStartingChain], MaxPotentialsOfChainsOrder[vStartingChain], vStartingChain);

        allChainsLabels[vStartingChain] = ComputeLabellingForOneChain(vStartingChain, optimalBV, MaxPotentials[vStartingChain], LinearPotentials[vStartingChain]);   
        PropagateChainSolutionToOtherChains(vStartingChain, gridSolution, allChainsLabels);  
        INDEX seedGridLoc = ChainNodeToOriginalNode[hStartingChain][fixedNodeIndexInHorizontalChains];

        // Solve Downward:         
        for (long int currentYOffset = hStartingChain - 1, fixedGridLoc = seedGridLoc - chainInfo.HorizontalSize(); 
            currentYOffset >= 0; currentYOffset--, fixedGridLoc -= chainInfo.HorizontalSize()) 
        {
            INDEX cDest = chainInfo.GetHorizontalChainAtOffset(currentYOffset);
            INDEX cSource = chainInfo.GetHorizontalChainAtOffset(currentYOffset + 1);
            auto [modifiedMaxP, modifiedLinearP] = ReparametrizeHorizontalChainFromFixedLabels(currentYOffset + 1, currentYOffset, chainInfo, allChainsLabels[cSource]);
            auto modifiedMaxP1D = GetMaxPotentials1D(modifiedMaxP);
            auto modifiedMaxP1DSortingOrder = GetMaxPotentialSortingOrder(modifiedMaxP1D);
            REAL optimalB = ComputeChainMarginals<true, true, true>(MarginalsChains[cDest], modifiedMaxP, modifiedLinearP, modifiedMaxP1D, modifiedMaxP1DSortingOrder,
                                                                    cDest, fixedNodeIndexInHorizontalChains, gridSolution[fixedGridLoc]);
            allChainsLabels[cDest] = ComputeLabellingForOneChain<true>(cDest, optimalB, modifiedMaxP, modifiedLinearP, 
                                                                        fixedNodeIndexInHorizontalChains, gridSolution[fixedGridLoc]);   

            PropagateChainSolutionToOtherChains(cDest, gridSolution, allChainsLabels);  
        }

        // Solve Upward:
        for (long int currentYOffset = hStartingChain + 1, fixedGridLoc = seedGridLoc + chainInfo.HorizontalSize(); 
            currentYOffset < chainInfo.NumHorizontal(); currentYOffset++, fixedGridLoc += chainInfo.HorizontalSize()) 
        {
            INDEX cDest = chainInfo.GetHorizontalChainAtOffset(currentYOffset);
            INDEX cSource = chainInfo.GetHorizontalChainAtOffset(currentYOffset - 1);
            auto [modifiedMaxP, modifiedLinearP] = ReparametrizeHorizontalChainFromFixedLabels(currentYOffset - 1, currentYOffset, chainInfo, allChainsLabels[cSource]);
            auto modifiedMaxP1D = GetMaxPotentials1D(modifiedMaxP);
            auto modifiedMaxP1DSortingOrder = GetMaxPotentialSortingOrder(modifiedMaxP1D);
            REAL optimalB = ComputeChainMarginals<true, true, true>(MarginalsChains[cDest], modifiedMaxP, modifiedLinearP, modifiedMaxP1D, modifiedMaxP1DSortingOrder,
                                                                    cDest, fixedNodeIndexInHorizontalChains, gridSolution[fixedGridLoc]);
            allChainsLabels[cDest] = ComputeLabellingForOneChain<true>(cDest, optimalB, modifiedMaxP, modifiedLinearP, 
                                                                        fixedNodeIndexInHorizontalChains, gridSolution[fixedGridLoc]);   
                                                                        
            PropagateChainSolutionToOtherChains(cDest, gridSolution, allChainsLabels);  
        }
        return allChainsLabels;
    }

    std::pair<three_dimensional_variable_array<REAL>, three_dimensional_variable_array<REAL>>
    ReparametrizeHorizontalChainFromFixedLabels(const INDEX vOffsetSource, const INDEX vOffsetDest, const ChainsInfo& chainInfo, const std::vector<INDEX>& sourceLabels) const 
    {
        assert(std::abs((long)vOffsetSource - (long)vOffsetDest) == 1);
        INDEX sourceC = chainInfo.GetHorizontalChainAtOffset(vOffsetSource);
        INDEX destC = chainInfo.GetHorizontalChainAtOffset(vOffsetDest);
    
        three_dimensional_variable_array<REAL> modifiedLinearPots(LinearPotentials[destC]);
        three_dimensional_variable_array<REAL> modifiedMaxPots(MaxPotentials[destC]); 

        assert(NumNodes[sourceC] == NumNodes[destC]);
        for (INDEX n = 0; n < NumNodes[destC]; n++) {
            std::vector<max_linear_costs> unaryPot(NumLabels[destC][n]);
            INDEX gridIndex = ChainNodeToOriginalNode[destC][n]; //Grid Index of node to reparameterize
            INDEX vertC = chainInfo.GetVerticalChainIndexAtGridLoc(gridIndex);
            INDEX destNodeInVert = vOffsetDest;
            if (vOffsetDest < vOffsetSource) { // Reparametrize downwards
                assert(NumLabels[vertC][destNodeInVert] == unaryPot.size());
                for (INDEX ld = 0; ld < NumLabels[vertC][destNodeInVert]; ld++) {
                    unaryPot[ld] = { MaxPotentials[vertC](destNodeInVert, ld, sourceLabels[n]), LinearPotentials[vertC](destNodeInVert, ld, sourceLabels[n])};
                }
            }
            else { // Reparametrize upwards
                assert(NumLabels[vertC][destNodeInVert] == unaryPot.size());
                for (INDEX ld = 0; ld < NumLabels[vertC][destNodeInVert]; ld++) {
                    unaryPot[ld] = { MaxPotentials[vertC](destNodeInVert - 1, sourceLabels[n], ld), LinearPotentials[vertC](destNodeInVert - 1, sourceLabels[n], ld)};
                }
            }
            ReparameterizePairwiseFromUnary(modifiedMaxPots, modifiedLinearPots, unaryPot, destC, n);
        }
        return std::make_pair(modifiedMaxPots, modifiedLinearPots);
    }

    void ReparameterizePairwiseFromUnary(three_dimensional_variable_array<REAL>& maxPots, three_dimensional_variable_array<REAL>& linearPots,
                                         const std::vector<max_linear_costs>& unary, INDEX c, INDEX n) const {
        REAL normalizer = 2.0;
        if (n == 0 || n == NumNodes[c] - 1) 
            normalizer = 1;

        if (n < NumNodes[c] - 1) {
            // reparam right:
            assert(NumLabels[c][n] == unary.size());
            for (INDEX l1 = 0; l1 < NumLabels[c][n]; l1++) {
                for (INDEX l2 = 0; l2 < NumLabels[c][n+1]; l2++) {
                    maxPots(n, l1, l2) = std::max(maxPots(n, l1, l2), unary[l1].MaxCost);
                    linearPots(n, l1, l2) = linearPots(n, l1, l2) + (unary[l1].LinearCost / normalizer);
                }
            }
        } 
        if (n > 0) { // reparam left:
            assert(NumLabels[c][n] == unary.size());
            for (INDEX l1 = 0; l1 < NumLabels[c][n-1]; l1++) {
                for (INDEX l2 = 0; l2 < NumLabels[c][n]; l2++) {
                    maxPots(n-1, l1, l2) = std::max(maxPots(n-1, l1, l2), unary[l2].MaxCost);
                    linearPots(n-1, l1, l2) = linearPots(n-1, l1, l2) + (unary[l2].LinearCost / normalizer);
                }
            }
        }
    }

    std::vector<MaxPotentialInChain> GetMaxPotentials1D(const three_dimensional_variable_array<REAL>& maxPotentials) const {
        std::vector<MaxPotentialInChain> maxPotentials1D;
        for (INDEX e = 0; e < maxPotentials.size(); e++) {
            for (INDEX l1 = 0; l1 < maxPotentials.dim2(e); l1++) {
                for (INDEX l2 = 0; l2 < maxPotentials.dim3(e); l2++) {
                    if (std::isinf(maxPotentials(e, l1, l2))) continue;
                    maxPotentials1D.push_back({e, l1, l2, maxPotentials(e, l1, l2)});
                }
            }
        }
        return maxPotentials1D;
    }

    std::vector<INDEX> ComputeGridSolution(const two_dim_variable_array<INDEX>& chainsSolution) const {
        INDEX numNodesGrid = 0;
        for (INDEX c = 0; c < NumChains; c++) {
            numNodesGrid = std::max(numNodesGrid, *std::max_element(ChainNodeToOriginalNode[c].begin(), ChainNodeToOriginalNode[c].end()));
        }
        std::vector<INDEX> gridSolution(numNodesGrid, std::numeric_limits<INDEX>::max());
        for (INDEX c = 0; c < NumChains; c++) {
            for (INDEX n = 0; n < NumNodes[c]; n++) {
                INDEX& currentGridLabel = gridSolution[ChainNodeToOriginalNode[c][n]];
                if (currentGridLabel == std::numeric_limits<INDEX>::max())      // solution not set
                    currentGridLabel = chainsSolution[c][n];
                else if (currentGridLabel != chainsSolution[c][n])
                    currentGridLabel = std::numeric_limits<INDEX>::infinity(); // Store disagreements as infinities   
            }
        }
        return gridSolution;
    }

    template <bool leftToRight>
    std::vector<Marginals> ComputeNodeMarginals(const INDEX chainIndex, INDEX node) const
    {
        INDEX numL = NumLabels[chainIndex][node];
        std::vector<Marginals> nodeMarginals(numL);
        shortest_distance_calculator<leftToRight, true> distCalc(LinearPotentials[chainIndex], MaxPotentials[chainIndex], NumLabels[chainIndex], node);

        for (const auto& e : MaxPotentialsOfChainsOrder[chainIndex]) {
            const auto n1 = MaxPotentialsOfChains[chainIndex][e].Edge;
            const auto l1 = MaxPotentialsOfChains[chainIndex][e].L1;
            const auto l2 = MaxPotentialsOfChains[chainIndex][e].L2;
            const auto bottleneckCost = MaxPotentialsOfChains[chainIndex][e].Value;
               
            distCalc.AddEdgeWithUpdate(n1, l1, l2, bottleneckCost);
            for (INDEX currentL = 0; currentL < numL; currentL++) {
                REAL l = distCalc.GetDistance(node, currentL);
                nodeMarginals[currentL].insert({bottleneckCost, l});
                //storing infinities as well, to ensure consistency with min-marginal computation.
            }
        }
        return nodeMarginals;
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

    //TODO: Maybe dealing with Infs explicitly can help.
    three_dimensional_variable_array<REAL> ComputeBackwardMessagesOnChain(INDEX c) const {
        //1. Initialize messages by zero.
        three_dimensional_variable_array<REAL> messages(LinearPotentials[c].size_begin(), LinearPotentials[c].size_end());
        three_dimensional_variable_array<INDEX> bestMIndices(LinearPotentials[c].size_begin(), LinearPotentials[c].size_end());
        for (INDEX e = 0; e < messages.dim1(); e++) {
            for (INDEX n1 = 0; n1 < messages.dim2(e); n1++) {
                for (INDEX n2 = 0; n2 < messages.dim3(e); n2++) {
                    messages(e, n1, n2) = std::numeric_limits<REAL>::max();
                    bestMIndices(e, n1, n2) = std::numeric_limits<INDEX>::max();
                }
            }
        }
        REAL normalizer = messages.dim1() * (NumChains - NumChainsReparametrizedBackwards);
        //2. Populate chain solution costs.
        REAL lb = LowerBound();
        //3. Merge the chain solution costs, by the node solver.
        Marginals otherChainsMergedMarginals;
        if (NumChains > 1) {
            std::vector<Marginals> otherChainMarginals;
            for (INDEX oc = 0; oc < NumChains; oc++) {
                if (oc == c) continue;
                otherChainMarginals.push_back(MarginalsChains[oc]);
            }
            max_potential_on_nodes otherChainsSolver = max_potential_on_nodes(otherChainMarginals);
            otherChainsSolver.ComputeBestLabels<true>(otherChainMarginals);
            otherChainsMergedMarginals = otherChainsSolver.GetAllMarginals();
        }
        // 4. Initialize two node solver, where node1 containing the merged marginals of all chains except of c:
        max_potential_on_two_nodes twoNodeSolver(otherChainsMergedMarginals);
        // 5. Create left and right shortest path calculators:
        shortest_distance_calculator<true> leftDistCalc(LinearPotentials[c], MaxPotentials[c], NumLabels[c]);
        shortest_distance_calculator<false> rightDistCalc(LinearPotentials[c], MaxPotentials[c], NumLabels[c]);
        for (const auto& e : MaxPotentialsOfChainsOrder[c]) {
            const auto n1 = MaxPotentialsOfChains[c][e].Edge;
            const auto l1 = MaxPotentialsOfChains[c][e].L1;
            const auto l2 = MaxPotentialsOfChains[c][e].L2;
            const auto bottleneckCost = MaxPotentialsOfChains[c][e].Value;

            // Calculate the min-marginal for given edge:
            const REAL n1LeftDistance = leftDistCalc.GetDistance(n1, l1);
            const REAL n2RightDistance = rightDistCalc.GetDistance(n1 + 1, l2);
            REAL currentEdgeMinMarginal = ComputeMinMarginal(n1LeftDistance, n2RightDistance, twoNodeSolver, bestMIndices, c, n1, l1, l2, bottleneckCost);
            assert(currentEdgeMinMarginal >= lb - eps);
            messages(n1, l1, l2) = std::min(messages(n1, l1, l2), (currentEdgeMinMarginal - lb) / normalizer);

            // Each edge addition can change the marginals of all edges on its right in the forward(left) distance calculator,
            // and all the edges on the left of backward(right) distance calculator.
            // 6. Update the edges on the right of left updated nodes:
            const std::vector<std::array<INDEX, 2>> leftUpdatedNodeLabels = leftDistCalc.AddEdgeWithUpdate<true>(n1, l1, l2, bottleneckCost);
            for (INDEX i = 0; i < leftUpdatedNodeLabels.size(); i++) {
                const auto [n_n1, n_l1] = leftUpdatedNodeLabels[i];
                if (n_n1 >= messages.dim1()) continue;
                for (INDEX n_l2 = 0; n_l2 < NumLabels[c][n_n1 + 1]; n_l2++) {
                    if (MaxPotentials[c](n_n1, n_l1, n_l2) > bottleneckCost) continue; // this edge is going to come, so can be computed right then.
                    const REAL n1LeftDistance = leftDistCalc.GetDistance(n_n1, n_l1);
                    const REAL n2RightDistance = rightDistCalc.GetDistance(n_n1 + 1, n_l2);
                    REAL currentEdgeMinMarginal = ComputeMinMarginal(n1LeftDistance, n2RightDistance, twoNodeSolver, bestMIndices, c, n_n1, n_l1, n_l2, bottleneckCost);
                    assert(currentEdgeMinMarginal >= lb - eps);
                    messages(n_n1, n_l1, n_l2) = std::min(messages(n_n1, n_l1, n_l2), (currentEdgeMinMarginal - lb) / normalizer);
                }
            }
            // 7. Update the edges on the left of right updated nodes:
            const std::vector<std::array<INDEX, 2>> rightUpdatedNodeLabels = rightDistCalc.AddEdgeWithUpdate<true>(n1, l1, l2, bottleneckCost);
            for (INDEX i = 0; i < rightUpdatedNodeLabels.size(); i++) {
                const auto [p_n2, p_l2] = rightUpdatedNodeLabels[i];
                if (p_n2 == 0) continue;
                for (INDEX p_l1 = 0; p_l1 < NumLabels[c][p_n2 - 1]; p_l1++) {
                    if (MaxPotentials[c](p_n2 - 1, p_l1, p_l2) > bottleneckCost) continue; // this edge is going to come, so can be computed right then.
                    const REAL n1LeftDistance = leftDistCalc.GetDistance(p_n2 - 1, p_l1);
                    const REAL n2RightDistance = rightDistCalc.GetDistance(p_n2, p_l2);
                    REAL currentEdgeMinMarginal = ComputeMinMarginal(n1LeftDistance, n2RightDistance, twoNodeSolver, bestMIndices, c, p_n2 - 1, p_l1, p_l2, bottleneckCost);
                    assert(currentEdgeMinMarginal >= lb - eps);
                    messages(p_n2 - 1, p_l1, p_l2) = std::min(messages(p_n2 - 1, p_l1, p_l2), (currentEdgeMinMarginal - lb) / normalizer);
                }
            }
        }
        return messages;
    }

    REAL ComputeMinMarginal(const REAL n1LeftDistance, const REAL n2RightDistance, 
                        const max_potential_on_two_nodes& twoNodeSolver, three_dimensional_variable_array<INDEX>& bestMIndices,
                        const INDEX c, const INDEX n1, const INDEX l1, const INDEX l2, REAL bottleneckCost) const {
        const max_linear_costs currentEdgeCost = {bottleneckCost, n1LeftDistance + n2RightDistance + LinearPotentials[c](n1, l1, l2)};
        const auto& pair = twoNodeSolver.ComputeBestIndexAndCost(currentEdgeCost, bestMIndices(n1, l1, l2)); 
        bestMIndices(n1, l1, l2) = pair.first;
        #ifndef NDEBUG
            const auto& pairD = twoNodeSolver.ComputeBestIndexAndCost(currentEdgeCost, 0); 
            assert(pairD.first == pair.first);
            assert(pairD.second.TotalCost() == pair.second.TotalCost());
        #endif
        return pair.second.TotalCost();
    }
};

class pairwise_max_potential_on_multiple_chains_message {
    public:
        pairwise_max_potential_on_multiple_chains_message(const INDEX chain_index, const INDEX pairwise_index, // TODO: simplify constructor, remove unary_1, unary_2
                                                          const INDEX unary_1, const INDEX unary_2) :
        ChainIndex(chain_index), EdgeIndex(pairwise_index), N1(unary_1), N2(unary_2) {}

        void TurnOffUpdateFlags() {
            ToUpdate = false;
        }

        template<typename FACTOR, typename MSG>
        void RepamRight(FACTOR& r, const MSG& msgs) const
        {
            INDEX i = 0;
            for(INDEX l1=0; l1<r.NumNodeLabels(ChainIndex, N1); ++l1) {
                for(INDEX l2=0; l2<r.NumNodeLabels(ChainIndex, N2); ++l2, ++i) {
                    chain_edge e = {EdgeIndex, l1, l2};
                    REAL currentPot = r.GetLinearPotential(ChainIndex, e);
                    if (std::isinf(currentPot))
                        continue;
                    else if (std::isinf(msgs[i]))
                        r.SetLinearPotential(ChainIndex, e, msgs[i]);
                    else                         
                        r.SetLinearPotential(ChainIndex, e, currentPot + msgs[i]);
                }
            }
        }

        template<typename FACTOR, typename MSG>
        void RepamLeft(FACTOR& l, const MSG& msgs) const
        {
            INDEX c=0;
            for(INDEX i=0; i<l.dim1(); ++i) {
                for(INDEX j=0; j<l.dim2(); ++j, ++c) {
                    if (std::isinf(l.cost(i,j)))
                        continue;
                    else if (std::isinf(msgs[c]))
                        l.cost(i,j) = msgs[c];
                    else
                        l.cost(i,j) += msgs[c];
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
           // TODO; avoid extra copy, write with a matrix
            std::vector<REAL> message = r.GetMessageForEdge(ChainIndex, EdgeIndex);
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

            return (changed_N1 || changed_N2) && ToUpdate;
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

            return (changed_N1 || changed_N2) && ToUpdate;
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

        INDEX GetChainIndex() const { return ChainIndex; }

    private:
        const INDEX ChainIndex;
        const INDEX EdgeIndex;
        const INDEX N1;
        const INDEX N2;
        bool ToUpdate = true;
};
}

#endif // LPMP_HORIZON_TRACKING_FACTORS_HXX
