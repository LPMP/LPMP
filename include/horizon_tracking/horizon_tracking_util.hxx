#ifndef LPMP_HORIZON_TRACKING_UTIL_HXX
#define LPMP_HORIZON_TRACKING_UTIL_HXX

#include <queue>
#include <limits>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"
#include <unordered_set>
#include <cmath>

namespace LPMP {
struct chain_edge { INDEX n1, l1, l2; };
struct max_linear_costs { 
    REAL MaxCost; REAL LinearCost; 
    REAL TotalCost() const { return MaxCost + LinearCost; }
};
class Marginals {
private:
    std::vector<max_linear_costs> M;
    INDEX i = 0;
    bool populated = false;

public:
    void insert(max_linear_costs m) {
        if (!populated) {
    #ifndef NDEBUG 
            if(M.size() > 0) assert(M.back().MaxCost <= m.MaxCost);
    #endif
            if (M.size() == 0) M.push_back(m);
            else if (M.back().MaxCost == m.MaxCost) M.back().LinearCost = std::min(M.back().LinearCost, m.LinearCost);
            else M.push_back(m);
        }
        else {
            if (i > 0 && M[i - 1].MaxCost == m.MaxCost) {
                M[i - 1].LinearCost = std::min(M[i - 1].LinearCost, m.LinearCost);
            }
            else { 
                assert(M[i].MaxCost == m.MaxCost);
                M[i++].LinearCost = m.LinearCost;
            }
        }
    }

    const std::vector<max_linear_costs>& Get() const {return M;}
    void Populated() {
        i = 0;
        populated = true;
    }

    const max_linear_costs Get(const INDEX j) const { assert(j < M.size()); return M[j]; }
    INDEX size() const { return M.size();}
    // const max_linear_costs operator[](const INDEX j) const { Get(j); }

    void Clear() { M.clear(); }
};

template <bool DoForward, bool useStartingNode = false>
class shortest_distance_calculator {
private: 
    const three_dimensional_variable_array<REAL>& LinearPairwisePotentials;
    const three_dimensional_variable_array<REAL>& MaxPairwisePotentials;
    const std::vector<INDEX>& NumLabels;
    two_dim_variable_array<REAL> distance;
    REAL shortestPathDistance;
    struct edge { INDEX n1, l1, l2; }; 
    const INDEX EndingNodeIndex;

public:
    shortest_distance_calculator(const three_dimensional_variable_array<REAL>& linearPairwisePotentials, 
                                const three_dimensional_variable_array<REAL>& maxPairwisePotentials,
                                const std::vector<INDEX>& numLabels, INDEX endingNode = 0):
                                LinearPairwisePotentials(linearPairwisePotentials),
                                MaxPairwisePotentials(maxPairwisePotentials), NumLabels(numLabels),
                                EndingNodeIndex(endingNode) {
        distance.resize(numLabels.begin(), numLabels.end(), std::numeric_limits<REAL>::max());
        init();
    }

    void init() {
        shortestPathDistance = std::numeric_limits<REAL>::max();
        if(DoForward) { std::fill(distance[0].begin(), distance[0].end(), 0); }
        else { std::fill(distance[distance.size() - 1].begin(), distance[distance.size() - 1].end(), 0); }
    }

    REAL GetDistance(INDEX n, INDEX l) const { return distance[n][l]; }

    //TODO: Only does forward distance calculation!
    void AddEdge(INDEX n1, INDEX l1, INDEX l2) {
        if (DoForward) {
            distance[n1 + 1][l2] = std::min(distance[n1 + 1][l2], distance[n1][l1] + LinearPairwisePotentials(n1, l1, l2));
            if (n1 + 1 == distance.size() - 1) // reached last node.
                shortestPathDistance = *std::min_element(distance[n1 + 1].begin(), distance[n1 + 1].end());
        }
        else {
            assert(false); // not required so far.
        }
    }
    
    void AddEdgeWithUpdate(INDEX n1, INDEX l1, INDEX l2, REAL bottleneckThreshold) {
        assert(MaxPairwisePotentials(n1, l1, l2) <= bottleneckThreshold);
        if (!ToAddEdge(n1, l1, l2))
            return;

        std::queue<edge> queue;
        queue.push({n1, l1, l2});

        while(!queue.empty()) {
            edge e = queue.front(); queue.pop();
            assert(e.n1 < LinearPairwisePotentials.dim1());
            REAL currentLinearPot = LinearPairwisePotentials(e.n1, e.l1, e.l2);
            INDEX currentNode = DoForward ? e.n1 : e.n1 + 1;
            INDEX nextNode = DoForward ? e.n1 + 1 : e.n1;
            INDEX currentLabel = DoForward ? e.l1 : e.l2;
            INDEX nextLabel = DoForward ? e.l2 : e.l1;
            REAL offeredDistance = distance[currentNode][currentLabel] + currentLinearPot;
            if (distance[nextNode][nextLabel] <= offeredDistance) continue; // TODO: do not add to queue if that condition holds

            distance[nextNode][nextLabel] = offeredDistance;
            if ((DoForward && currentNode == LinearPairwisePotentials.dim1() - 1) || (!DoForward && nextNode == 0))  {
                shortestPathDistance = std::min(shortestPathDistance, offeredDistance);
                continue;
            }
            INDEX childNode = nextNode + (DoForward ?  + 1 : -1);
            for (INDEX childLabel = 0; childLabel < NumLabels[childNode]; ++childLabel) {
                if (DoForward && MaxPairwisePotentials(nextNode, nextLabel, childLabel) <= bottleneckThreshold
                    && ToAddEdge(nextNode, nextLabel, childLabel)) {
                    queue.push({nextNode, nextLabel, childLabel});
                }
                else if (!DoForward && MaxPairwisePotentials(childNode, childLabel, nextLabel) <= bottleneckThreshold
                        && ToAddEdge(childNode, childLabel, nextLabel)) {
                    queue.push({childNode, childLabel, nextLabel});
                }
            }
        }
    }

    bool ToAddEdge(INDEX n1, INDEX l1, INDEX l2) const {
        if (!useStartingNode) return true;
        else if (DoForward && (n1 >= EndingNodeIndex))  {
            return false;
        }
        else if (!DoForward && (n1 < EndingNodeIndex)) {
            return false;
        }
        return true;
    }

    REAL ShortestDistance() const { return shortestPathDistance; }
    
    template <bool ToSpecificNode = false>
    std::vector<INDEX> ShortestPath(REAL bottleneckThreshold, INDEX endingNode = 0, INDEX endingLabel = 0) const {
        std::vector<INDEX> path;
        if (DoForward) {
            if(!ToSpecificNode)
                endingNode = NumLabels.size() - 1;

            INDEX numNodes = endingNode + 1;
            path.resize(numNodes, std::numeric_limits<INDEX>::max());
            INDEX n2 = endingNode;
            INDEX n1 = n2 - 1;
            path[n2] = endingLabel; 
            if (!ToSpecificNode)
                path[n2] = std::min_element(distance[n2].begin(),distance[n2].end()) - distance[n2].begin();
            while (n2 > 0) {
                for (INDEX l1 = 0; l1 < NumLabels[n1]; ++l1) {
                    if (MaxPairwisePotentials(n1, l1, path[n2]) > bottleneckThreshold) continue;
                    REAL currentLinearPot = LinearPairwisePotentials(n1, l1, path[n2]);
                    if (distance[n1][l1] + currentLinearPot == distance[n2][path[n2]]) {
                        path[n1] = l1;
                        break;
                    }
                }
                n1--;
                n2--;
            }
        }
        else {
            if(!ToSpecificNode)
                endingNode = 0;

            INDEX numNodes = NumLabels.size() - endingNode;
            path.resize(numNodes, std::numeric_limits<INDEX>::max());
            INDEX n2 = endingNode + 1;
            INDEX n1 = endingNode;
            path[0] = endingLabel;
            if (!ToSpecificNode)
                path[0] = std::min_element(distance[n1].begin(),distance[n1].end()) - distance[n1].begin();

            INDEX pathNode = 0;
            while (n2 <= NumLabels.size() - 1) {
                for (INDEX l2 = 0; l2 < NumLabels[n2]; ++l2) {
                    if (MaxPairwisePotentials(n1, path[pathNode], l2) > bottleneckThreshold) continue;
                    REAL currentLinearPot = LinearPairwisePotentials(n1, path[pathNode], l2);
                    if (distance[n1][path[pathNode]] == distance[n2][l2] + currentLinearPot) {
                        path[pathNode + 1] = l2;
                        break;
                    }
                }
                assert(path[pathNode + 1] < std::numeric_limits<INDEX>::max());
                pathNode++;
                n1++;
                n2++;
            }
        }
        assert(*std::max_element(path.begin(), path.end()) < std::numeric_limits<INDEX>::max());
        return path;
    }

    template <bool specificNode = false>
    void CalculateDistances(REAL bottleneckThreshold, INDEX startNode = 0, INDEX startLabel = 0) {
        init();
        INDEX numNodes = NumLabels.size();

        if (DoForward) {
            for (INDEX n1 = startNode, n2 = startNode + 1; n2 < numNodes; n1++, n2++) {
                for (INDEX l1 = startLabel; l1 < NumLabels[n1]; l1++) {
                    for (INDEX l2 = 0; l2 < NumLabels[n2]; l2++) {
                        if(MaxPairwisePotentials(n1, l1, l2) > bottleneckThreshold) continue;
                        if (distance(n2,l2) > distance(n1,l1) + LinearPairwisePotentials(n1, l1, l2))
                            distance(n2,l2) = distance(n1,l1) + LinearPairwisePotentials(n1, l1, l2);
                    }
                    if (specificNode && n1 == startNode)
                        break;
                }
            }
            shortestPathDistance = *std::min_element(distance[numNodes - 1].begin(), distance[numNodes - 1].end());
        }
        else {
            if (!specificNode) startNode = numNodes - 1;  
            for (long int n2 = startNode , n1 = startNode - 1; n1 >= 0; n2--, n1--) {
            // TODO: exchanging iteration order slows down code, even tough it should result in more contiguous access
               for (INDEX l1 = 0; l1 < NumLabels[n2 - 1]; l1++) {
                    for (INDEX l2 = startLabel; l2 < NumLabels[n2]; l2++) {
                        if(MaxPairwisePotentials(n1, l1, l2) > bottleneckThreshold) continue;
                        if (distance(n1,l1) > distance(n2,l2) + LinearPairwisePotentials(n1, l1, l2))
                            distance(n1,l1) = distance(n2,l2) + LinearPairwisePotentials(n1, l1, l2);
                    }
                    if (specificNode && n2 == startNode)
                        break;
                }
            }
            shortestPathDistance = *std::min_element(distance[0].begin(), distance[0].end());
        }
    }
};

class ChainsInfo {
public:
    ChainsInfo(const two_dim_variable_array<INDEX>& chainNodeToOriginalNode) : 
    ChainNodeToOriginalNode(chainNodeToOriginalNode),
    NumChains(chainNodeToOriginalNode.size()) {
        IsHorizontalChain.resize(NumChains);
        NumHorizontalChains = 0;
        NumVerticalChains = 0;
        for (INDEX c = 0; c < NumChains; c++) {
            if (ChainNodeToOriginalNode[c][1] - ChainNodeToOriginalNode[c][0] == 1) {
                IsHorizontalChain[c] = true;
                NumHorizontalChains++;
            } else {
                IsHorizontalChain[c] = false;
                NumVerticalChains++;
            }
        }
        HChainIndices.resize(NumHorizontalChains);
        VChainIndices.resize(NumVerticalChains);
        for (INDEX c = 0; c < NumChains; c++) {
            if (IsHorizontalChain[c]) {
                INDEX yOffset = ChainNodeToOriginalNode[c][0] / NumVerticalChains;
                HChainIndices[yOffset] = c;
            } else {
                INDEX xOffset = ChainNodeToOriginalNode[c][0] % NumVerticalChains;
                VChainIndices[xOffset] = c;
            }
        }
    }

    INDEX GetVerticalChainIndexAtGridLoc(const INDEX gridLoc) const {
        return VChainIndices[gridLoc % NumVerticalChains];
    }

    INDEX GetHorizontalChainIndexAtGridLoc(const INDEX gridLoc) const {
        return HChainIndices[gridLoc / NumVerticalChains];
    }

    // Returns the y-offset for a horizontal chain having a node with minimum number of labels, for 
    // starting primal computation.
    INDEX GetPrimalStartChainOffset(const std::vector<std::vector<INDEX>>& numLabels) const {
        INDEX minIndex;
        INDEX minValue = std::numeric_limits<INDEX>::max();
        for (INDEX c = 0; c < NumChains; c++) {
            if (IsVertical(c)) continue;
            INDEX currentChainMinLabels = *std::min_element(numLabels[c].begin(), numLabels[c].end());
            if (minValue > currentChainMinLabels) {
                minValue = currentChainMinLabels;
                minIndex = c;
            }
        }
        return GetVerticalOffset(ChainNodeToOriginalNode[minIndex][0]);
    }

    INDEX GetHorizontalChainAtOffset(const INDEX hOffset) const { return HChainIndices[hOffset]; }
    INDEX GetVerticalChainAtOffset(const INDEX vOffset) const { return VChainIndices[vOffset]; }

    INDEX GetHorizontalOffset(const INDEX gridLoc) const { return gridLoc % NumVerticalChains; }
    INDEX GetVerticalOffset(const INDEX gridLoc) const { return gridLoc / NumVerticalChains; }
    bool IsHorizontal(const INDEX c) const { return IsHorizontalChain[c]; }
    bool IsVertical(const INDEX c) const { return !IsHorizontal(c); }
    INDEX NumHorizontal() const { return NumHorizontalChains; }
    INDEX NumVertical() const { return NumVerticalChains; }

private:
    const two_dim_variable_array<INDEX> ChainNodeToOriginalNode;
    const INDEX NumChains;
    std::vector<bool> IsHorizontalChain;   // 1 -> horizontal, 0 -> vertical
    std::vector<INDEX> HChainIndices; // gives horizontal chain index at the given y-offset
    std::vector<INDEX> VChainIndices; // gives vertical chain index at the given x-offset
    INDEX NumHorizontalChains;
    INDEX NumVerticalChains;

};
}
#endif //HORIZON_TRACKING_UTIL_HXX
