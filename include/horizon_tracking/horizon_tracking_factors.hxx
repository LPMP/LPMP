#ifndef LPMP_HORIZON_TRACKING_FACTORS_HXX
#define LPMP_HORIZON_TRACKING_FACTORS_HXX

#include "vector.hxx"
#include <queue>
#include <limits>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"
#include <unordered_set>
#include <cmath>

namespace LPMP {
struct max_linear_costs { REAL MaxCost; REAL LinearCost; };
struct max_linear_rep_costs { REAL MaxCost; REAL LinearCost; REAL ReparamCost; };

class max_potential_on_graph {
struct MaxPotentialElement { REAL value; INDEX tableIndex; INDEX labelIndex; };
public:
    // indices correspond to: chain index, entry in chain, max pot value, linear cost.
    max_potential_on_graph(const std::vector<std::vector<max_linear_costs>>& marginals_collection)
    : marginals_collection_(marginals_collection) {
        for(INDEX currentTableIndex = 0; currentTableIndex < marginals_collection_.size(); ++currentTableIndex) {
            for(INDEX currentLabel = 0; currentLabel < marginals_collection_[currentTableIndex].size(); ++currentLabel ) {
                MaxPotentials.push_back( { marginals_collection_[currentTableIndex][currentLabel].MaxCost, currentTableIndex, currentLabel } );
            }
        }
        SortingOrder = GetMaxPotsSortingOrder(MaxPotentials);
        max_potential_index_.resize(marginals_collection_.size());
        init_primal();
    }

    void init_primal() {
        std::fill(max_potential_index_.begin(), max_potential_index_.end(), std::numeric_limits<INDEX>::max());
        chainMessageNormalizer = marginals_collection_.size();                
    }

    REAL LowerBound() const { return Solve(); }

    REAL EvaluatePrimal() const {
        const bool primal_computed = *std::max_element(max_potential_index_.begin(), max_potential_index_.end()) < std::numeric_limits<INDEX>::max();
        assert(primal_computed);
        // return cost of current solution //TODO: Replace by CostOfLabelling
        REAL linearCost = 0;
        REAL maxCost = std::numeric_limits<REAL>::lowest();
        for(INDEX currentTableIndex = 0; currentTableIndex < marginals_collection_.size(); ++currentTableIndex)
        {
            assert(max_potential_index_[currentTableIndex] < marginals_collection_[currentTableIndex].size());
            maxCost = std::max(maxCost, marginals_collection_[currentTableIndex][max_potential_index_[currentTableIndex]].MaxCost);
            linearCost += marginals_collection_[currentTableIndex][max_potential_index_[currentTableIndex]].LinearCost;
        }
        return maxCost + linearCost;
    }

    void MaximizePotentialAndComputePrimal() { Solve(); }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

    auto export_variables() { return std::tie( ); }
                
    template<typename EXTERNAL_SOLVER>
    void construct_constraints(EXTERNAL_SOLVER& s) const { assert(false); }
    template<typename EXTERNAL_SOLVER>
    void convert_primal(EXTERNAL_SOLVER& s) { assert(false); } 

    void set_max_potential_index(const INDEX tableIndex, const INDEX newValue) {
        assert(tableIndex < max_potential_index_.size()); 
        max_potential_index_[tableIndex] = newValue;
    }

    INDEX max_potential_index(const INDEX tableIndex) const {
        assert(tableIndex < max_potential_index_.size());
        return max_potential_index_[tableIndex]; 
    }

    two_dim_variable_array<max_linear_costs>& marginals_collection() { return marginals_collection_; }

    const two_dim_variable_array<max_linear_costs>& marginals_collection() const { return marginals_collection_; }

    std::vector<REAL> messages_to_chain(INDEX chainIndex) const { return ComputeMessagesToChains(chainIndex); }

private:
    two_dim_variable_array<max_linear_costs> marginals_collection_;
    std::vector<MaxPotentialElement> MaxPotentials;
    mutable std::vector<INDEX> max_potential_index_;
    std::vector<INDEX> SortingOrder;
    mutable INDEX chainMessageNormalizer;
    template <bool fixLabels = false> 
    REAL Solve(INDEX fixedTableIndex = 0, INDEX fixedLabelIndex = 0) const {
        INDEX numCovered = 0;
        INDEX numTables = marginals_collection_.size();
        std::vector<bool> coveredTables(numTables, false);
        REAL s = 0;
        std::vector<REAL> l(numTables, INFINITY);
        double bestObjective = INFINITY;
        std::vector<INDEX> bestLabelsForTables(numTables);
        std::vector<INDEX> lablesForTables(numTables);
        for(const auto& currentElementToInsert : SortingOrder) {
            INDEX currentTableIndex = MaxPotentials[currentElementToInsert].tableIndex;
            INDEX currentLabelIndex = MaxPotentials[currentElementToInsert].labelIndex;
            if (fixLabels && currentTableIndex == fixedTableIndex && currentLabelIndex != fixedLabelIndex) { continue; }

            REAL currentLinearCost = marginals_collection_[currentTableIndex][currentLabelIndex].LinearCost;
            REAL currentMaxCost =  MaxPotentials[currentElementToInsert].value;
            assert(currentMaxCost == marginals_collection_[currentTableIndex][currentLabelIndex].MaxCost);
            // If the edge is not yet covered:
            if (!coveredTables[currentTableIndex])  {
                coveredTables[currentTableIndex] = true;
                numCovered++;
                s += currentLinearCost;
                l[currentTableIndex] = currentLinearCost;
                lablesForTables[currentTableIndex] = currentLabelIndex;
            }
            // If edge has been added, but current label has lower linear cost. We have two scenarios:
            // 1. Graph has not been covered completely in which case we want to increase our max pot. threshold anyway. Thus, if we are gaining 
            //      an improvement in linear cost take it.
            // 2. Graph is covered completely, but we want to see adding current label can possibly give us any benefit, which will be 
            //    checked by the 3rd if condition.
            if (currentLinearCost < l[currentTableIndex]) {
                s = s - l[currentTableIndex] + currentLinearCost;
                l[currentTableIndex] = currentLinearCost;
                lablesForTables[currentTableIndex] = currentLabelIndex;
            }

            if (numCovered == numTables && bestObjective > s + currentMaxCost) {
                // Found another solution which is better than the previous one, in which case mark current solution as the best so far.
                    bestObjective = s + currentMaxCost;
                    bestLabelsForTables = lablesForTables;
            }
        }
        if (!fixLabels) {
            max_potential_index_ = bestLabelsForTables;
#ifndef NDEBUG 
            assert(std::abs(bestObjective - EvaluatePrimal()) <= eps);
#endif          
        }
        return bestObjective;
    }

    std::vector<INDEX> ComputeSolution(REAL maxPotValue) const {
        INDEX numTables = marginals_collection_.size();
        std::vector<INDEX> solution(numTables);
        for (INDEX tableIndex = 0; tableIndex < numTables; tableIndex++) {
            REAL currentTableCost = std::numeric_limits<REAL>::max();
            auto numLabels = marginals_collection_[tableIndex].size();
            for (INDEX labelIndex = 0; labelIndex < numLabels; labelIndex++) {
                REAL currentMaxPotValue = marginals_collection_[tableIndex][labelIndex].MaxCost;
                if (currentMaxPotValue > maxPotValue) // as marginals are sorted so we can break.
                    break;
                
                REAL currentCost = marginals_collection_[tableIndex][labelIndex].LinearCost;
                if (currentCost < currentTableCost) {
                    currentTableCost = currentCost;
                    solution[tableIndex] = labelIndex;
                }
            }
        }
        return solution;
    }
    
    std::vector<INDEX> GetMaxPotsSortingOrder(const std::vector<MaxPotentialElement>& pots) const {
        std::vector<INDEX> idx(pots.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::sort(idx.begin(), idx.end(),
            [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});
        return idx;
    }

    std::vector<REAL> ComputeMessagesToChains(INDEX chainIndex, REAL retentionFactor = 1e-2) const {
        std::vector<REAL> messagesToChain(marginals_collection_[chainIndex].size());
        assert(chainMessageNormalizer > 0);
        for(INDEX currentLabelIndex = 0; currentLabelIndex < messagesToChain.size(); ++currentLabelIndex) {
            messagesToChain[currentLabelIndex] = Solve<true>(chainIndex, currentLabelIndex) / (chainMessageNormalizer + retentionFactor);
        }
        chainMessageNormalizer--;
        return messagesToChain;
    }

    REAL CostOfLabelling(const std::vector<INDEX>& labels) const {
        REAL linearCost = 0;
        REAL maxCost = std::numeric_limits<REAL>::lowest();
        for(INDEX currentTableIndex = 0; currentTableIndex < marginals_collection_.size(); ++currentTableIndex) {
            assert(labels[currentTableIndex] < marginals_collection_[currentTableIndex].size());
            maxCost = std::max(maxCost, marginals_collection_[currentTableIndex][max_potential_index_[currentTableIndex]].MaxCost);
            linearCost += marginals_collection_[currentTableIndex][labels[currentTableIndex]].LinearCost;
        }
        return maxCost + linearCost;
    }
};

template <bool DoForward>
class shortest_distance_calculator {
private: 
    three_dimensional_variable_array<REAL> LinearPairwisePotentials;
    three_dimensional_variable_array<REAL> MaxPairwisePotentials;
    std::vector<INDEX> NumLabels;
    two_dim_variable_array<REAL> distance;
    REAL shortestPathDistance;
    struct edge { INDEX n1, l1, l2; }; 

public:
    shortest_distance_calculator(const three_dimensional_variable_array<REAL>& linearPairwisePotentials, 
                                const three_dimensional_variable_array<REAL>& maxPairwisePotentials,
                                const std::vector<INDEX>& numLabels):
                                LinearPairwisePotentials(linearPairwisePotentials),
                                MaxPairwisePotentials(maxPairwisePotentials),
                                NumLabels(numLabels) {
        distance.resize(numLabels.begin(), numLabels.end(), std::numeric_limits<REAL>::max());
        init();
    }

    void init() {
        if(DoForward) { std::fill(distance[0].begin(), distance[0].end(), 0); }
        else { std::fill(distance[distance.size() - 1].begin(), distance[distance.size() - 1].end(), 0); }
        shortestPathDistance = std::numeric_limits<REAL>::max();
    }

    REAL GetDistance(INDEX n, INDEX l) const { return distance[n][l]; }

    template <bool checkCollision = false>
    void AddEdgeWithUpdate(INDEX n1, INDEX l1, INDEX l2, REAL bottleneckThreshold, INDEX cn = 0, INDEX cl1 = 0, INDEX cl2 = 0) {
        assert(MaxPairwisePotentials(n1, l1, l2) <= bottleneckThreshold);
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
            if (distance[nextNode][nextLabel] <= offeredDistance) continue;

            distance[nextNode][nextLabel] = offeredDistance;
            if ((DoForward && currentNode == LinearPairwisePotentials.dim1() - 1) || (!DoForward && nextNode == 0))  {
                shortestPathDistance = std::min(shortestPathDistance, offeredDistance);
                continue;
            }
            INDEX childNode = nextNode + (DoForward ?  + 1 : -1);
            for (INDEX childLabel = 0; childLabel < NumLabels[childNode]; ++childLabel) {
                if (DoForward) {
                    if (MaxPairwisePotentials(nextNode, nextLabel, childLabel) > bottleneckThreshold || 
                        (checkCollision && nextNode == cn && (cl1 != nextLabel || cl2 != childLabel))) 
                        continue;
                    
                    queue.push({nextNode, nextLabel, childLabel});
                }
                else {
                    if (MaxPairwisePotentials(childNode, childLabel, nextLabel) > bottleneckThreshold || 
                        (checkCollision && childNode == cn && (cl1 != childLabel || cl2 != nextLabel)))
                        continue;

                    queue.push({childNode, childLabel, nextLabel});
                }
            }
        }
    }

    void AddEdge(INDEX n1, INDEX l1, INDEX l2) {
        distance[n1 + 1][l2] = std::min(distance[n1 + 1][l2], distance[n1][l1] + LinearPairwisePotentials(n1, l1, l2));
        if (n1 + 1 == LinearPairwisePotentials.dim1() - 1)
            shortestPathDistance = *std::min_element(distance[n1 + 1].begin(), distance[n1 + 1].end());
    }

    REAL shortestDistance() const { return shortestPathDistance; }

    std::vector<INDEX> shortestPath(INDEX endingNode, INDEX endingLabel, REAL bottleneckThreshold) const {
        std::vector<INDEX> path;

        if (DoForward) {
            INDEX numNodes = endingNode + 1;
            path.resize(numNodes, std::numeric_limits<INDEX>::max());
            INDEX n2 = endingNode;
            INDEX n1 = n2 - 1;
            path[n2] = endingLabel; //std::min_element(distance[n2].begin(),distance[n2].end()) - distance[n2].begin();
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
            INDEX numNodes = NumLabels.size() - endingNode;
            path.resize(numNodes, std::numeric_limits<INDEX>::max());
            INDEX n2 = endingNode + 1;
            INDEX n1 = endingNode;
            path[0] = endingLabel; // std::min_element(distance[n1].begin(),distance[n1].end()) - distance[n1].begin();
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

    void CalculateDistances(REAL bottleneckThreshold, bool specificNode = false, INDEX startNode = 0, INDEX startLabel = 0) {
        init();
        INDEX numNodes = NumLabels.size();

        if (DoForward) {
            for (INDEX n1 = startNode, n2 = startNode + 1; n2 < numNodes; n1++, n2++) {
                for (INDEX l1 = startLabel; l1 < NumLabels[n1]; l1++) {
                    for (INDEX l2 = 0; l2 < NumLabels[n2]; l2++) {
                        if(MaxPairwisePotentials(n1, l1, l2) > bottleneckThreshold) continue;
                        if (distance[n2][l2] > distance[n1][l1] + LinearPairwisePotentials(n1, l1, l2))
                            distance[n2][l2] = distance[n1][l1] + LinearPairwisePotentials(n1, l1, l2);
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
                for (INDEX l2 = startLabel; l2 < NumLabels[n2]; l2++) {
                    for (INDEX l1 = 0; l1 < NumLabels[n2 - 1]; l1++) {
                        if(MaxPairwisePotentials(n1, l1, l2) > bottleneckThreshold) continue;
                        if (distance[n1][l1] > distance[n2][l2] + LinearPairwisePotentials(n1, l1, l2))
                            distance[n1][l1] = distance[n2][l2] + LinearPairwisePotentials(n1, l1, l2);
                    }
                    if (specificNode && n2 == startNode)
                        break;
                }
            }
            shortestPathDistance = *std::min_element(distance[0].begin(), distance[0].end());
        }
    }
};
class max_potential_on_chain {       
public:     
    struct MaxPairwisePotentialEntry {
        REAL value;
        INDEX edgeIndex;
        INDEX n1;
        INDEX n2; 
        INDEX l1;
        INDEX l2; 
    };
    struct edge { INDEX n1, l1, l2; }; 
    struct EdgePriority { REAL value; INDEX index;};
    struct EdgePriorityComparison { bool operator() (const EdgePriority& lhs, const EdgePriority& rhs) { lhs.value < rhs.value; }};

    mutable three_dimensional_variable_array<REAL> LinearPairwisePotentials;  

    max_potential_on_chain() {}
    max_potential_on_chain(const three_dimensional_variable_array<REAL>& maxPairwisePotentials, const three_dimensional_variable_array<REAL>& linearPairwisePotentials, const std::vector<INDEX>& numLabels, INDEX chainIndex, bool useEdgeDeletion = false)
    :
        LinearPairwisePotentials(linearPairwisePotentials), MaxPairwisePotentials(maxPairwisePotentials),
        NumNodes(numLabels.size()), NumLabels(numLabels), ChainIndex_(chainIndex)
    {
        assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
        assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());
        for(INDEX n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
            assert(maxPairwisePotentials.dim2(n1) == linearPairwisePotentials.dim2(n1) && maxPairwisePotentials.dim3(n1) == linearPairwisePotentials.dim3(n1));
            for(INDEX i=0; i<maxPairwisePotentials.dim2(n1); ++i) {
                for(INDEX j=0; j<maxPairwisePotentials.dim3(n1); ++j) {
                    MaxPotentials1D.push_back( {maxPairwisePotentials(n1,i,j), n1, n1, n1+1, i, j} );
                }
            }
        }
        MaxPotsSortingOrder = GetPairwisePotsSortingOrder(MaxPotentials1D);
        solution_.assign(NumNodes, 0);
        init_primal();
    }

    INDEX ChainIndex() const {return ChainIndex_; }

    REAL LowerBound() const {
        if (!max_potential_marginals_valid_) Solve();
            
        INDEX bestIndex = GetBestMarginal();
        return max_potential_marginals_[bestIndex].LinearCost + max_potential_marginals_[bestIndex].ReparamCost;
    }

    REAL EvaluatePrimal() const {
        assert(max_potential_index_ != std::numeric_limits<INDEX>::max());
        auto objFromPath = PathSolutionObjective(solution_, LinearPairwisePotentials, MaxPairwisePotentials);
        INDEX indexFromSolution = GetMarginalIndexFromSolution(solution_);
        assert(objFromPath.MaxCost ==  max_potential_marginals_[max_potential_index_].MaxCost);
        return objFromPath.LinearCost + max_potential_marginals_[indexFromSolution].ReparamCost;
    }

    void MaximizePotentialAndComputePrimal() {
        const bool labels_computed = *std::max_element(solution_.begin(), solution_.end()) < std::numeric_limits<INDEX>::max();
        const bool max_potential_index_computed = max_potential_index_ != std::numeric_limits<INDEX>::max();
        
        if (max_potential_index_computed && labels_computed) return;
        if (max_potential_index_computed) {
            solution_ = ComputeSolution(max_potential_index_);
            assert(GetMarginalIndexFromSolution(solution_) == max_potential_index_);
            return; 
        }
        if (!max_potential_index_computed && !labels_computed) { 
            if (!max_potential_marginals_valid_) {
                Solve();
            }
            set_max_potential_index(GetBestMarginal());
            solution_ = ComputeSolution(max_potential_index_);
            // assert(GetMarginalIndexFromSolution(solution_) == max_potential_index_);
        }
        if (labels_computed && !max_potential_index_computed) {
            assert(false);
        }
    }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( LinearPairwisePotentials.data() ); }

    auto export_variables() { return std::tie( ); }
    void init_primal() {
        max_potential_index_ = std::numeric_limits<INDEX>::max();
        std::fill(solution_.begin(), solution_.end(), std::numeric_limits<INDEX>::max());
        pairwiseMessagesComputed_ = false;
        pairwiseMessageNormalizer = LinearPairwisePotentials.dim1();
    } 

    void invalidate_marginals() { max_potential_marginals_valid_ = false; }

    template<typename ARRAY>
    void apply(ARRAY& a) const { 
        INDEX offset = 0;
        for(INDEX n=0; n<solution_.size()-1; ++n) {
            a[offset + solution_[n]* LinearPairwisePotentials.dim3(n) + solution_[n+1]];
            offset += LinearPairwisePotentials.dim2(n) * LinearPairwisePotentials.dim3(n); 
        }
    }
    template<typename EXTERNAL_SOLVER> void construct_constraints(EXTERNAL_SOLVER& s) const { assert(false); }
    template<typename EXTERNAL_SOLVER> void convert_primal(EXTERNAL_SOLVER& s) { assert(false); } 

    INDEX& solution(const INDEX i) { assert(i < solution_.size()); return solution_[i]; }

    INDEX solution(const INDEX i) const { assert(i < solution_.size()); return solution_[i]; }

    //TO ADDRESS: Move solution_ computation to messages.
    void set_max_potential_index(const INDEX index) const {
        assert(max_potential_marginals_valid_);
        if (max_potential_index_ == index)
            return;

        max_potential_index_ = index;
        solution_ = ComputeSolution(max_potential_index_); //TODO: DO OR DO NOT ?
    }

    INDEX max_potential_index() const { return max_potential_index_; }
    void SetMaxPotentialIndexFromSolution() const {
        if (*std::max_element(solution_.begin(), solution_.end()) >= std::numeric_limits<INDEX>::max())
            return;
            
        max_potential_index_ = GetMarginalIndexFromSolution(solution_); 
    }

    max_linear_rep_costs& max_potential_marginal(const INDEX i) {
        assert(i < max_potential_marginals_.size()); 
        assert(max_potential_marginals_valid_);
        return max_potential_marginals_[i]; 
    }

    const max_linear_rep_costs& max_potential_marginal(const INDEX i) const {
        assert(i < max_potential_marginals_.size());
        assert(max_potential_marginals_valid_);
        return max_potential_marginals_[i];
    }

    void set_marginal_slack(std::vector<REAL> slack) const {
        marginalSlack_ = slack;
    }
    
    const INDEX max_potential_marginals_size() const { return max_potential_marginals_.size(); }

    // // Computes the messages to pairwise MRF and tries to distributes uniformly, this one 
    // // can be used for parallel update, down-side is that some slack can still remain inside.
    // three_dimensional_variable_array<REAL> edge_marginals() const {
    //     if (!pairwiseMessagesComputed_) {
    //         edgePairwiseMessages_ = ComputeMessagesToPairwise();
    //         pairwiseMessagesComputed_ = true;
    //     }
        
    //     return edgePairwiseMessages_;
    // }

    // // Computes the messages to pairwise MRF only to the specific edgeIndex, this one has a variable
    // // which tells the amount of edges already reparameterized and therefore adjusts the normalization
    // // accordingly. This one can move all the slack back to pairwise MRF when pairwiseMessageNormalizer is 1
    // // i.e. the last edge.
    std::vector<REAL> ComputeMessagesToPairwiseEdge(INDEX edgeIndex, REAL retentionFactor = 1e-2) const {
        std::vector<REAL> edge_messages(LinearPairwisePotentials.dim2(edgeIndex) * LinearPairwisePotentials.dim3(edgeIndex));
        assert(pairwiseMessageNormalizer > 0);
        INDEX i = 0;
        for (INDEX l1 = 0; l1 < LinearPairwisePotentials.dim2(edgeIndex); l1++) {
            for (INDEX l2 = 0; l2 < LinearPairwisePotentials.dim3(edgeIndex); l2++) {
                auto m = ComputeMinMarginalForEdge({edgeIndex, l1, l2});
                edge_messages[i] = m / (retentionFactor + pairwiseMessageNormalizer);  
                i++;
            }
        }
        pairwiseMessageNormalizer--;
        return edge_messages;
    }

    void ConvertMarginalSlackToPairwiseSlack() const
    {
        assert(marginalSlack_.size() == max_potential_marginals_.size());
        std::vector<INDEX> slackOrder(marginalSlack_.size());
        std::iota(slackOrder.begin(), slackOrder.end(), 0);
        std::sort(slackOrder.begin(), slackOrder.end(), [&slack = marginalSlack_](INDEX i1, INDEX i2) {return slack[i1] < slack[i2];});
        three_dimensional_variable_array<uint8_t> locked_edges(LinearPairwisePotentials.size_begin(), LinearPairwisePotentials.size_end(), 0);
        three_dimensional_variable_array<REAL> delta(LinearPairwisePotentials.size_begin(), LinearPairwisePotentials.size_end(), 0);
        for (const auto& index : slackOrder) {
            //TODO: Can be optimized:
            std::vector<INDEX> path = ComputeSolution(index);
            INDEX highestMaxPotentialN1;
            REAL highestMaxPotentialValue = std::numeric_limits<REAL>::lowest();
            for (INDEX n1 = 0; n1 < locked_edges.dim1(); n1++) {
                if (MaxPairwisePotentials(n1, path[n1], path[n1 + 1]) > highestMaxPotentialValue && locked_edges(n1, path[n1], path[n1 + 1]) == 0) {
                    highestMaxPotentialValue = MaxPairwisePotentials(n1, path[n1], path[n1 + 1]);
                    highestMaxPotentialN1 = n1;
                }
            }
            if(highestMaxPotentialValue == std::numeric_limits<REAL>::lowest()) continue; // All edges locked!
            REAL currentSlack = marginalSlack_[index];
            assert(currentSlack >= 0);
            INDEX numLocked = 0;
            for (INDEX n1 = 0; n1 < path.size() - 1; n1++) {
                currentSlack -= delta(n1, path[n1], path[n1 + 1]); // Subtract the slack which is already covered, because of some edges getting reparameterized previously.
                numLocked += locked_edges(n1, path[n1], path[n1 + 1]);
            }
            if (currentSlack < 0) // Too much slack pushed into current path, remove equally.
            {
                REAL totalSlack = currentSlack;
                for (INDEX n1 = 0; n1 < path.size() - 1; n1++) {
                    if (locked_edges(n1, path[n1], path[n1 + 1]) == 0) continue;
                    delta(n1, path[n1], path[n1 + 1]) += totalSlack / numLocked; 
                    currentSlack -= totalSlack / numLocked;
                }
                assert(std::abs(currentSlack) <= eps);
                marginalSlack_[index] = 0;
                continue;
            }
            // Move slack into the pairwise potential with highest bottleneck potential value. However lock all the edges of the path to ensure that the edges of 
            // the path are not reparametrized again.
            delta(highestMaxPotentialN1, path[highestMaxPotentialN1], path[highestMaxPotentialN1 + 1]) += currentSlack;
            for (INDEX n1 = 0; n1 < locked_edges.dim1(); n1++) { locked_edges(n1, path[n1], path[n1 + 1]) = 1; }
            auto pathSolutionSlack = PathSolutionObjective(path, delta, MaxPairwisePotentials); 
            assert(std::abs(pathSolutionSlack.LinearCost - marginalSlack_[index]) <= eps);
            marginalSlack_[index] = 0;
            //TODO: Move delta to LinearPots.
        }
        for (INDEX n1 = 0; n1 < LinearPairwisePotentials.dim1(); n1++) {
            for (INDEX l1 = 0; l1 < NumLabels[n1]; l1++) {
                for (INDEX l2 = 0; l2 < NumLabels[n1 + 1]; l2++) {
                    LinearPairwisePotentials(n1, l1, l2) += delta(n1, l1, l2);
                    delta(n1, l1, l2) = 0;
                }
            }
        }
    }

protected:
    std::vector<INDEX> NumLabels;
    // primal solution
    mutable std::vector<INDEX> solution_;
    mutable INDEX max_potential_index_;
    three_dimensional_variable_array<REAL> MaxPairwisePotentials;
    mutable std::vector<MaxPairwisePotentialEntry> MaxPotentials1D;
    std::vector<INDEX> MaxPotsSortingOrder;
    mutable std::vector<max_linear_rep_costs> max_potential_marginals_; 
    mutable std::vector<edge> edge_in_max_potential_marginals_;
    mutable bool max_potential_marginals_valid_ = false;
    mutable bool MaxPotMarginalsInitialized = false;
    INDEX ChainIndex_;

    INDEX NumNodes;
    mutable three_dimensional_variable_array<REAL> edgePairwiseMessages_;
    mutable std::vector<REAL> marginalSlack_;
    mutable bool pairwiseMessagesComputed_ = false;
    mutable INDEX pairwiseMessageNormalizer;

    INDEX GetBestMarginal() const {
        assert(max_potential_marginals_valid_);
        REAL bestCost = std::numeric_limits<REAL>::max();
        INDEX bestIndex;
        for (INDEX currentMaxPotIndex = 0; currentMaxPotIndex < max_potential_marginals_.size(); ++currentMaxPotIndex) {
            auto currentMarginal = max_potential_marginals_[currentMaxPotIndex];
            auto currentCost = currentMarginal.MaxCost + currentMarginal.LinearCost + currentMarginal.ReparamCost;
            if (currentCost < bestCost) {
                bestIndex = currentMaxPotIndex;
                bestCost = currentCost;
            }
        }
        assert(bestCost < std::numeric_limits<REAL>::max());
        assert(max_potential_marginals_.size() > 0);
        return bestIndex;
    }

    template <bool insertEnd>
    bool InsertMarginal(std::vector<max_linear_rep_costs>& marginals, REAL maxPotValue, INDEX insertionIndex, REAL currentLinearCost, bool marginals_populated) const {
        if (!marginals_populated) {
            if(insertEnd) {
                if(marginals.size() > 0 && marginals.back().MaxCost == maxPotValue) {
                    marginals.back().LinearCost = std::min(marginals.back().LinearCost, currentLinearCost);
                    return false;
                } else { 
                    marginals.push_back({maxPotValue, currentLinearCost, 0});
                    return true;
                }
            }
            else {
                if(marginals.size() > 0 && marginals.front().MaxCost == maxPotValue) {
                    marginals.front().LinearCost = std::min(marginals.front().LinearCost, currentLinearCost);
                    return false;
                }
                else {
                    marginals.insert(marginals.begin(), {maxPotValue, currentLinearCost, 0});
                    return true;
                }
            }
        }
        else {                       
            // Check if the current max potential value is also present in the marginals (can only be present at adjacent index as they are sorted),
            // if yes just take the minimum of the linear costs.
            INDEX adjacentIndex = insertionIndex + (insertEnd ? -1:1);
            if (adjacentIndex >= 0 && adjacentIndex <= marginals.size() - 1 &&
                marginals[adjacentIndex].MaxCost == maxPotValue) {
                    marginals[adjacentIndex].LinearCost = std::min(marginals[adjacentIndex].LinearCost, currentLinearCost);
                    return false;
            }
            
            assert(maxPotValue == marginals[insertionIndex].MaxCost);                   
            marginals[insertionIndex].LinearCost = currentLinearCost;

            return true;
        }
    }

    std::vector<INDEX> GetPairwisePotsSortingOrder(const std::vector<MaxPairwisePotentialEntry>& pots) const {
        std::vector<INDEX> idx(pots.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&pots](INDEX i1, INDEX i2) {return pots[i1].value < pots[i2].value;});
        return idx;
    }

private:
    virtual void Solve() const {
        SolveByEdgeAddition();
        MaxPotMarginalsInitialized = true;
        max_potential_marginals_valid_ = true;
    }

    // Computes the labelling given a max pot index.
    virtual std::vector<INDEX> ComputeSolution(INDEX maxPotIndex) const {
        assert(max_potential_marginals_valid_);

        shortest_distance_calculator<true> forwardCalc(LinearPairwisePotentials, MaxPairwisePotentials, NumLabels);
        //TODO: Distance calculation can be optimized:
        forwardCalc.CalculateDistances(max_potential_marginals_[maxPotIndex].MaxCost); //, true, edge_in_max_potential_marginals_[maxPotIndex].n1,
                                                                                       //     edge_in_max_potential_marginals_[maxPotIndex].l1);

        std::vector<INDEX> solF = forwardCalc.shortestPath(edge_in_max_potential_marginals_[maxPotIndex].n1,
                                                           edge_in_max_potential_marginals_[maxPotIndex].l1,
                                                            max_potential_marginals_[maxPotIndex].MaxCost); 

        shortest_distance_calculator<false> backwardCalc(LinearPairwisePotentials, MaxPairwisePotentials, NumLabels);
        backwardCalc.CalculateDistances(max_potential_marginals_[maxPotIndex].MaxCost); // true, edge_in_max_potential_marginals_[maxPotIndex].n1 + 1,
                                                                                       //     edge_in_max_potential_marginals_[maxPotIndex].l2);

        std::vector<INDEX> solB = backwardCalc.shortestPath(edge_in_max_potential_marginals_[maxPotIndex].n1 + 1,
                                                           edge_in_max_potential_marginals_[maxPotIndex].l2,
                                                            max_potential_marginals_[maxPotIndex].MaxCost);

        solF.insert(solF.end(), solB.begin(), solB.end());
        #ifndef NDEBUG 
            auto pathSolObj = PathSolutionObjective(solF, LinearPairwisePotentials, MaxPairwisePotentials);
            assert(pathSolObj.MaxCost <= max_potential_marginals_[maxPotIndex].MaxCost);
        #endif 
        return solF;
    }

    virtual max_linear_costs PathSolutionObjective(std::vector<INDEX> sol, three_dimensional_variable_array<REAL> linearPots, three_dimensional_variable_array<REAL> maxPots) const
    {
        REAL maxPotValue = std::numeric_limits<REAL>::lowest();
        REAL linearCost = 0;
        for (int n1 = 0; n1 < linearPots.dim1(); n1++)
        {
            if (sol[n1] == std::numeric_limits<INDEX>::max() || 
                sol[n1 + 1] == std::numeric_limits<INDEX>::max())
                return {std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max()}; // Solution not ready yet

            maxPotValue = std::max(maxPotValue, maxPots(n1, sol[n1], sol[n1 + 1]));
            linearCost += linearPots(n1, sol[n1], sol[n1 + 1]);
        }
        return {maxPotValue, linearCost};
    }

    virtual INDEX GetMarginalIndexFromSolution(std::vector<INDEX> sol) const 
    {
        REAL maxPotValue = std::numeric_limits<REAL>::lowest();
        REAL linearCost = 0;
        for (INDEX n1 = 0; n1 < LinearPairwisePotentials.dim1(); n1++) {
            const auto& l1 = sol[n1];
            const auto& l2 = sol[n1 + 1];
            maxPotValue = std::max(maxPotValue, MaxPairwisePotentials(n1, l1, l2));
            linearCost += LinearPairwisePotentials(n1, l1, l2);
        }
        INDEX m = std::numeric_limits<INDEX>::max();
        INDEX m2;
        REAL linearCostDiff = std::numeric_limits<REAL>::max();
        for (INDEX i = 0; i < max_potential_marginals_.size(); i++) {
            bool fullMatch =  sol[edge_in_max_potential_marginals_[i].n1] == edge_in_max_potential_marginals_[i].l1 &&  
                 sol[edge_in_max_potential_marginals_[i].n1 + 1] == edge_in_max_potential_marginals_[i].l2 && 
                 std::abs(max_potential_marginals_[i].LinearCost - linearCost) < linearCostDiff &&
                 max_potential_marginals_[i].MaxCost == maxPotValue;

            if (fullMatch) {
                m = i;
                linearCostDiff = std::abs(max_potential_marginals_[i].LinearCost - linearCost);
            }
            // Its possible that current solution does not have 
            // any correspondence with marginals, then find the index corresponding to bottleneck value.
            if (max_potential_marginals_[i].MaxCost == maxPotValue) { m2 = i; }
        }
        return m < std::numeric_limits<INDEX>::max() ? m : m2;
    }

    virtual max_linear_costs SolveByEdgeAddition(bool computeMarginals = true) const { 
        max_linear_costs bestSolutionCost = {std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max()};
        shortest_distance_calculator<true> distanceCalcForward(LinearPairwisePotentials, MaxPairwisePotentials, NumLabels);
        shortest_distance_calculator<false> distanceCalcBackward(LinearPairwisePotentials, MaxPairwisePotentials, NumLabels);
        INDEX currentMarginalIndex = 0;
        for(const auto& currentEdgeToInsert : MaxPotsSortingOrder) {
            const auto& n1 = MaxPotentials1D[currentEdgeToInsert].n1;
            const auto& n2 = MaxPotentials1D[currentEdgeToInsert].n2;
            const auto& l1 = MaxPotentials1D[currentEdgeToInsert].l1;
            const auto& l2 = MaxPotentials1D[currentEdgeToInsert].l2;
            const auto& bottleneckCost = MaxPotentials1D[currentEdgeToInsert].value;
            distanceCalcForward.AddEdgeWithUpdate(n1, l1, l2, bottleneckCost);
            distanceCalcBackward.AddEdgeWithUpdate(n1, l1, l2, bottleneckCost);

            REAL forwardDistance = distanceCalcForward.GetDistance(n1, l1);
            REAL backwardDistance = distanceCalcBackward.GetDistance(n2, l2);
            if (std::max(forwardDistance, backwardDistance) >= std::numeric_limits<REAL>::max()) continue;

            REAL linearCost = forwardDistance + backwardDistance + LinearPairwisePotentials(n1, l1, l2);
            if (!MaxPotMarginalsInitialized) {
                max_potential_marginals_.push_back({bottleneckCost, linearCost, 0});
                edge_in_max_potential_marginals_.push_back({n1, l1, l2});
            }
            else {
                max_potential_marginals_[currentMarginalIndex].LinearCost = linearCost; 
                currentMarginalIndex++;
            }
            
            if (bottleneckCost + linearCost < bestSolutionCost.MaxCost + bestSolutionCost.LinearCost) {
                bestSolutionCost.MaxCost = bottleneckCost;
                bestSolutionCost.LinearCost = linearCost;
            }
        }
        return bestSolutionCost;
    }

    // // For primal propagation back to pairwise potentials.
    // three_dimensional_variable_array<REAL> ComputeMessagesToPairwise() const
    // {
    //     three_dimensional_variable_array<REAL> edge_messages(LinearPairwisePotentials.size_begin(), LinearPairwisePotentials.size_end(), 0);
    //     auto lb = SolveByEdgeAddition(false);
    //     for (INDEX n1 = 0; n1 < edge_messages.dim1(); ++n1) {
    //         REAL minMinMarginal = std::numeric_limits<REAL>::max();
    //         for (INDEX l1 = 0; l1 < edge_messages.dim2(n1); l1++) {
    //             for (INDEX l2 = 0; l2 < edge_messages.dim3(n1); l2++) {
    //                 edge_messages(n1, l1, l2) = ComputeMinMarginalForEdge({n1, l1, l2});  
    //                 minMinMarginal = std::min(minMinMarginal, edge_messages(n1, l1, l2));
    //             }
    //         }
    //         assert(std::abs(lb.MaxCost + lb.LinearCost - minMinMarginal) <= eps);   
    //     }

    //     for (INDEX n1 = 0; n1 < edge_messages.dim1(); ++n1) {
    //         for (INDEX l1 = 0; l1 < edge_messages.dim2(n1); l1++) {
    //             for (INDEX l2 = 0; l2 < edge_messages.dim3(n1); l2++) {
    //                 edge_messages(n1, l1, l2) = (edge_messages(n1, l1, l2) - lb.MaxCost - lb.LinearCost) / edge_messages.dim1(); 
    //                 // Equal distribution
    //             }
    //         }
    //     }

    //     return edge_messages;
    // }


    REAL ComputeMinMarginalForEdge(edge e) const 
    {
        REAL maxPotV = MaxPairwisePotentials(e.n1, e.l1, e.l2);
        REAL minMarginal = std::numeric_limits<REAL>::max();
        shortest_distance_calculator<true> distCalc(LinearPairwisePotentials, MaxPairwisePotentials, NumLabels);

        for (INDEX n1 = 0; n1 < LinearPairwisePotentials.dim1(); ++n1) {         
            for (INDEX l1 = 0; l1 < NumLabels[n1]; l1++) {
                for (INDEX l2 = 0; l2 < NumLabels[n1+1]; l2++) {
                    if (n1 == e.n1 && (l1 != e.l1 || l2 != e.l2)) continue;
                    // Directly add all the edges which have max pot value lower than e:   
                    if (MaxPairwisePotentials(n1, l1, l2) <= maxPotV) 
                        distCalc.AddEdge(n1, l1, l2);
                }
            }
        }

        if (distCalc.shortestDistance() < std::numeric_limits<REAL>::max())
            minMarginal = maxPotV + distCalc.shortestDistance();

        // Iterative shortest path for higher edges:
        for(const auto& currentEdgeToInsert : MaxPotsSortingOrder) {
            const auto& n1 = MaxPotentials1D[currentEdgeToInsert].n1; 
            const auto& l1 = MaxPotentials1D[currentEdgeToInsert].l1; 
            const auto& l2 = MaxPotentials1D[currentEdgeToInsert].l2; 
            const auto& currentMaxCost = MaxPotentials1D[currentEdgeToInsert].value;

            if (n1 == e.n1 && (l1 != e.l1 || l2 != e.l2) || currentMaxCost <= maxPotV)
                continue; // do not consider colliding edges, lesser max pots edges already inserted

            distCalc.AddEdgeWithUpdate<true>(n1, l1, l2,  currentMaxCost, e.n1, e.l1, e.l2);
            if (minMarginal > distCalc.shortestDistance() + currentMaxCost)
                minMarginal = distCalc.shortestDistance() + currentMaxCost;
        }

        return minMarginal;
    }
};


class max_potential_on_tree_iterative : public max_potential_on_chain {
public:
// Assuming that the pairwise potentials are stored in the same order as messagePassingSchedule(mps), where edge \in mps points from tail to head
// Assuming that l1 = three_dimensional_variable_array.dim2(edge) corresponds to labels of node min(edge[0], edge[1]) and l2 to node labels of other node.
    max_potential_on_tree_iterative(const three_dimensional_variable_array<REAL>& maxPairwisePotentials, const three_dimensional_variable_array<REAL>& linearPairwisePotentials,  
        const std::vector<INDEX>& numLabels, const std::vector<std::array<INDEX, 2>>& messagePassingSchedule)
    {
        LinearPairwisePotentials = linearPairwisePotentials;
        MaxPairwisePotentials = maxPairwisePotentials;
        NumLabels = numLabels;
        NumNodes = numLabels.size();
        MessagePassingSchedule = messagePassingSchedule;
        assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
        assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());
        solution_.resize(NumNodes);
        INDEX edgeIndex = 0;
        for (const auto& currentEdge : messagePassingSchedule){
            assert(maxPairwisePotentials.dim2(edgeIndex) == linearPairwisePotentials.dim2(edgeIndex) &&
                   maxPairwisePotentials.dim3(edgeIndex) == linearPairwisePotentials.dim3(edgeIndex));

            for(INDEX l1=0; l1<maxPairwisePotentials.dim2(edgeIndex); ++l1) {
                for(INDEX l2=0; l2<maxPairwisePotentials.dim3(edgeIndex); ++l2) {
                    if (!isReverse(edgeIndex))
                        MaxPotentials1D.push_back( {maxPairwisePotentials(edgeIndex,l1,l2), edgeIndex, currentEdge[0], currentEdge[1], l1, l2} );
                    else
                        MaxPotentials1D.push_back( {maxPairwisePotentials(edgeIndex,l1,l2), edgeIndex, currentEdge[0], currentEdge[1], l2, l1} );
                }
            }
            edgeIndex++;
        }
        MaxPotsSortingOrder = GetPairwisePotsSortingOrder(MaxPotentials1D);
        BuildAdjacency();
        init_primal();
    }

protected:
    std::vector<std::array<INDEX, 2>> MessagePassingSchedule;
    struct NodeAdjacency {
        INDEX OutgoingEdge = std::numeric_limits<INDEX>::max(); //single outgoing edge index, (cannot be more than one outgoing edge)
        std::vector<INDEX> IncomingEdges; // its incoming edge indices
    };

    std::vector<NodeAdjacency> NodeIncidentEdges;

    // Given message passing order, populates outgoing edge index -1 corresponds to root node
    void BuildAdjacency() {
        NodeIncidentEdges.resize(NumNodes);
        INDEX edgeIndex = 0;
        for (const auto& currentEdge : MessagePassingSchedule) {
            NodeIncidentEdges[currentEdge[0]].OutgoingEdge = edgeIndex;
            NodeIncidentEdges[currentEdge[1]].IncomingEdges.push_back(edgeIndex);
            edgeIndex++;
        }

        // Validation:
        bool rootFound = false;
        for (INDEX n = 0; n < NodeIncidentEdges.size(); n++)
        {
            if (NodeIncidentEdges[n].OutgoingEdge < std::numeric_limits<INDEX>::max())
                continue;
            
            assert(!rootFound);
            rootFound = true;
        }
        assert(rootFound);
    }

    void Solve() const {
        REAL bestSolutionCost = INFINITY;
        std::vector<std::vector<std::vector<REAL>>> distanceFromSource(NumNodes); // node index, label index, incoming edge index
        REAL initialValue = 0; // For leaf nodes
        for (INDEX i = 0 ; i < NumNodes ; i++) {
            if (NodeIncidentEdges[i].IncomingEdges.size() > 0) // not leaf
                initialValue = std::numeric_limits<REAL>::max();
            else
                initialValue = 0;

            distanceFromSource[i].resize(NumLabels[i]);
            auto numIncoming = NodeIncidentEdges[i].IncomingEdges.size(); 
            if (numIncoming < 1) numIncoming = 1; // For leaves

            for (INDEX l = 0; l < distanceFromSource[i].size(); l++)
                distanceFromSource[i][l].resize(numIncoming, initialValue);
        }

        INDEX currentMaxPotIndex = 0;
        for(const auto& currentEdgeToInsert : MaxPotsSortingOrder)
        {
            REAL rootCost = UpdateDistances(MaxPotentials1D[currentEdgeToInsert], distanceFromSource, MaxPotentials1D[currentEdgeToInsert].value);

            if (rootCost == std::numeric_limits<REAL>::max())
                continue; 

            // Insert the marginal, and do not increment the index if the max pot was already present
            // at previous index in which case the marginal was not inserted and we only took min:
            if (InsertMarginal<true>(max_potential_marginals_, MaxPotentials1D[currentEdgeToInsert].value, currentMaxPotIndex, rootCost, MaxPotMarginalsInitialized))
                currentMaxPotIndex++;
        }
        MaxPotMarginalsInitialized = true;
        max_potential_marginals_valid_ = true;
    }

    REAL UpdateDistances(MaxPairwisePotentialEntry newElement, std::vector<std::vector<std::vector<REAL>>>& distanceFromSource, REAL maxPotThresh) const
    {
        REAL bestRootCost = std::numeric_limits<REAL>::max();
        std::queue<MaxPairwisePotentialEntry> queue;
        queue.push(newElement);

        while(!queue.empty())
        {
            auto currentElement = queue.front();
            queue.pop();
            const auto& edgeIndex = currentElement.edgeIndex;
            const auto& n1 = currentElement.n1;
            const auto& n2 = currentElement.n2;
            const auto& l1 = currentElement.l1;
            const auto& l2 = currentElement.l2;
            REAL worstDistanceToN1L1 = *std::max_element(distanceFromSource[n1][l1].begin(), distanceFromSource[n1][l1].end());
            if (worstDistanceToN1L1 >= std::numeric_limits<REAL>::max())
                continue; // There is no incoming path from atleast one branch

            REAL currentLinearPot = isReverse(edgeIndex) ? LinearPairwisePotentials(edgeIndex, l2, l1) : LinearPairwisePotentials(edgeIndex, l1, l2);
            
            // Sum up the cost of reaching n2, l2 from all of the incoming branches of n1, l1:
            REAL offeredDistanceToN2L2 = std::accumulate(distanceFromSource[n1][l1].begin(), distanceFromSource[n1][l1].end(), currentLinearPot);
            auto incomingEdgeIndexN2L2It = std::find(NodeIncidentEdges[n2].IncomingEdges.begin(), NodeIncidentEdges[n2].IncomingEdges.end(), edgeIndex);
#ifndef NDEBUG 
            assert(incomingEdgeIndexN2L2It != NodeIncidentEdges[n2].IncomingEdges.end()); // Should be able to find
#endif
            INDEX incomingEdgeIndexN2L2 = std::distance(NodeIncidentEdges[n2].IncomingEdges.begin(), incomingEdgeIndexN2L2It);
            auto currentDistanceToN2L2 = distanceFromSource[n2][l2][incomingEdgeIndexN2L2];
            if (offeredDistanceToN2L2 >= currentDistanceToN2L2)
                continue;

            distanceFromSource[n2][l2][incomingEdgeIndexN2L2] = offeredDistanceToN2L2;
            INDEX outgoingEdgeFromN2 = NodeIncidentEdges[n2].OutgoingEdge;
            
            if (outgoingEdgeFromN2 >= std::numeric_limits<INDEX>::max()) // Reached Root
            {
                REAL worstDistanceToRoot = *std::max_element(distanceFromSource[n2][l2].begin(), distanceFromSource[n2][l2].end());
                if (worstDistanceToRoot >= std::numeric_limits<REAL>::max())
                    continue;

                REAL currentRootCost = std::accumulate(distanceFromSource[n2][l2].begin(), distanceFromSource[n2][l2].end(), 0.0);
                bestRootCost = std::min(bestRootCost, currentRootCost);
                continue;
            }

            INDEX n3 = MessagePassingSchedule[outgoingEdgeFromN2][1];
            
            for (INDEX l3 = 0; l3 < NumLabels[n3]; l3++) {
                REAL currentMaxPot = isReverse(outgoingEdgeFromN2) ? MaxPairwisePotentials(outgoingEdgeFromN2, l3, l2) : MaxPairwisePotentials(outgoingEdgeFromN2, l2, l3);

                if (currentMaxPot > maxPotThresh)
                    continue;
                
                queue.push({currentMaxPot, outgoingEdgeFromN2, n2, n3, l2, l3});
            }
        }
        return bestRootCost;
    }

    std::vector<INDEX> ComputeSolution(INDEX maxPotIndex) const 
    {
        std::vector<INDEX> solution(NumNodes, std::numeric_limits<INDEX>::max());
        REAL maxPotThresh = max_potential_marginals_[maxPotIndex].MaxCost;
        std::vector<std::vector<REAL>> messages(NumNodes);
        std::vector<bool> messageSent(NumNodes, false);
        std::vector<bool> messageReceived(NumNodes, false);
        for (INDEX i = 0; i < NumNodes; i++) { messages[i].resize(NumLabels[i], 0); }

        INDEX edgeIndex = 0;
        for (const auto & currentEdge : MessagePassingSchedule)
        {
            INDEX tail = currentEdge[0];
            INDEX head = currentEdge[1];
            if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
            {
                assert(!messageSent[tail]);
            }
            for (INDEX lh = 0; lh < messages[head].size(); lh++)
            {
                REAL minMessage = ComputeMessageValue(messages[tail], edgeIndex, lh, false, maxPotThresh);
                messages[head][lh] += minMessage;
            }
            messageSent[tail] = true;
            messageReceived[head] = true;
            edgeIndex++;
        }

        for (auto it = MessagePassingSchedule.rbegin(); it != MessagePassingSchedule.rend(); ++it)           
        {
            edgeIndex--;
            const auto & currentEdge = *it;
            INDEX tail = currentEdge[1];
            INDEX head = currentEdge[0];
            REAL tailNodeBestMessage = INFINITY;
            solution[tail] = std::distance(messages[tail].begin(), std::min_element(messages[tail].begin(), messages[tail].end()));
            
            for (INDEX lh = 0; lh < messages[head].size(); lh++)
            {
                REAL minMessage = ComputeMessageValue(messages[tail], edgeIndex, lh, true, maxPotThresh, solution[tail]);
                messages[head][lh] += minMessage;
            }
        }

        // Set solution for leaf nodes:
        for (INDEX n = 0; n < NumNodes; n++) {
            if (solution[n] == std::numeric_limits<INDEX>::max())
                solution[n] = std::distance(messages[n].begin(), std::min_element(messages[n].begin(), messages[n].end()));
        }
#ifndef NDEBUG 
        auto solutionObj = PathSolutionObjective(solution);
        assert(solutionObj.MaxCost <= maxPotThresh);
        assert(std::abs(solutionObj.MaxCost + solutionObj.LinearCost  - 
        (max_potential_marginals_[maxPotIndex].MaxCost + max_potential_marginals_[maxPotIndex].LinearCost)) <= eps);
#endif 
        return solution;
    }

    REAL ComputeMessageValue(const std::vector<REAL>& tailMessages, INDEX edgeIndex, 
    INDEX lh, bool oppositeDirection, REAL maxPotThresh, INDEX labelTail = std::numeric_limits<INDEX>::max()) const
    {
        REAL bestMessageValue = INFINITY;
        INDEX bestlt;
        INDEX ltStart = 0;
        INDEX ltEnd = tailMessages.size() - 1;
        
        if (labelTail < std::numeric_limits<INDEX>::max())
        {
            ltStart = labelTail;
            ltEnd = labelTail;
        }
        bool takeReverse = !oppositeDirection ? isReverse(edgeIndex) : !isReverse(edgeIndex);
        for (INDEX lt = ltStart; lt <= ltEnd; lt++)
        {
            REAL currentEdgeMaxPot = !takeReverse ? MaxPairwisePotentials(edgeIndex, lt, lh) :  MaxPairwisePotentials(edgeIndex, lh, lt);
            if (currentEdgeMaxPot > maxPotThresh)
                continue;

            REAL currentEdgePot = !takeReverse ? LinearPairwisePotentials(edgeIndex, lt, lh) :  LinearPairwisePotentials(edgeIndex, lh, lt);
            REAL currentIncomingMsg = tailMessages[lt];

            REAL msg = currentIncomingMsg + currentEdgePot;
            if (msg < bestMessageValue)
            {
                bestMessageValue = msg;
                bestlt = lt;
            }
        }

        return bestMessageValue;
    }

    max_linear_costs PathSolutionObjective(std::vector<INDEX> sol) const
    {
        REAL maxPotValue = std::numeric_limits<REAL>::lowest();
        REAL linearCost = 0;
        INDEX edgeIndex = 0;
        for (const auto& currentEdge : MessagePassingSchedule) {
            INDEX n1 = currentEdge[0]; 
            INDEX n2 = currentEdge[1];
            if (isReverse(edgeIndex)) {
                n1 = currentEdge[1];
                n2 = currentEdge[0];
            }
            if (sol[n1] == std::numeric_limits<INDEX>::max() || 
                sol[n2] == std::numeric_limits<INDEX>::max())
                return {std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max()}; // Solution not ready yet

            maxPotValue = std::max(maxPotValue, MaxPairwisePotentials(edgeIndex, sol[n1], sol[n2]));
            linearCost += LinearPairwisePotentials(edgeIndex, sol[n1], sol[n2]);
            edgeIndex++;
        }
        return {maxPotValue, linearCost};
    }

    bool isReverse(INDEX edgeIndex) const { return MessagePassingSchedule[edgeIndex][1] < MessagePassingSchedule[edgeIndex][0]; }
};

class pairwise_max_factor_tree_message {
    public:
        pairwise_max_factor_tree_message(const INDEX _pairwise_entry, INDEX _unary_1, INDEX _unary_2) :
            pairwise_entry(_pairwise_entry), 
            unary_1(_unary_1),
            unary_2(_unary_2)
            {}

        template<typename FACTOR, typename MSG>
        void RepamRight(FACTOR& r, const MSG& msgs)
        {
            INDEX i = 0;
            for(INDEX l1=0; l1<r.LinearPairwisePotentials.dim2(pairwise_entry); ++l1) {
                for(INDEX l2=0; l2<r.LinearPairwisePotentials.dim3(pairwise_entry); ++l2, ++i) {
                    r.LinearPairwisePotentials(pairwise_entry, l1, l2) += msgs[i];
                }
            }
            r.invalidate_marginals();
        }

        template<typename FACTOR, typename MSG>
        void RepamLeft(FACTOR& l, const MSG& msgs)
        {
            INDEX c=0;
            for(INDEX i=0; i<l.dim1(); ++i) {
                for(INDEX j=0; j<l.dim2(); ++j) {
                    l.cost(i,j) += msgs[c++];
                }
            } 
        }

        template<typename LEFT_FACTOR, typename MSG>
        void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
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
        void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
        {
            r.ConvertMarginalSlackToPairwiseSlack();
            std::vector<REAL> m = r.ComputeMessagesToPairwiseEdge(pairwise_entry);
            const auto min = *std::min_element(m.begin(), m.end());
            vector<REAL> mm(m.size());
            for (INDEX i = 0; i < m.size(); i++) { mm[i] = m[i] - min; }
            msg -= omega * mm;
        }

        // template<typename RIGHT_FACTOR, typename MSG_ARRAY>
        // void SendMessagesToLeft(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
        // {
        //     assert(false);
        // }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
        {
            bool changed_unary_1 = false;
            if(r.solution(unary_1) < l.dim1()) {
                changed_unary_1 = (r.solution(unary_1) != l.primal()[0]);
                l.primal()[0] = r.solution(unary_1);
            }

            bool changed_unary_2 = false;
            if(r.solution(unary_2) < l.dim2()) {
                changed_unary_2 = (r.solution(unary_2) != l.primal()[1]);
                l.primal()[1] = r.solution(unary_2);
            }

            return changed_unary_1 || changed_unary_2;
        }

        // TODO: Should we update the max pot indices of chain based on the solution we just set 
        // from MRF factors, so that when we propagate from chain to graph the indices correspond to
        // the actual solution we got from unaries!
        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
        {
            bool changed_unary_1 = false;
            if(l.primal()[0] < l.dim1()) {
                changed_unary_1 = (r.solution(unary_1) != l.primal()[0]);
                r.solution(unary_1) = l.primal()[0];
                r.SetMaxPotentialIndexFromSolution();
            }

            bool changed_unary_2 = false;
            if(l.primal()[1] < l.dim2()) {
                changed_unary_2 = (r.solution(unary_2) != l.primal()[1]);
                r.solution(unary_2) = l.primal()[1];
                r.SetMaxPotentialIndexFromSolution();
            }

            return changed_unary_1 || changed_unary_2;
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
        {
            return r.solution(unary_1) == l.primal()[0] && r.solution(unary_2) == l.primal()[1];
        } 

        template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
        void construct_constraints(SOLVER& s, 
                LEFT_FACTOR& l, typename SOLVER::vector l_left_msg_variables, typename SOLVER::vector l_right_msg_variables, typename SOLVER::matrix l_pairwise_variables, 
                RIGHT_FACTOR& r)
        {
        }

    private:
        const INDEX pairwise_entry;
        const INDEX unary_1, unary_2;
};

// left factor is chain/tree
// right factor is max_potential_on_graph
class max_factor_tree_graph_message {
    public:

        max_factor_tree_graph_message(const INDEX _chain_index) : chain_index(_chain_index) {}

        template<typename FACTOR, typename MSG>
        void RepamRight(FACTOR& r, const MSG& msg)
        {
            assert(r.marginals_collection()[chain_index].size() == msg.size());
            for(INDEX i=0; i<r.marginals_collection()[chain_index].size(); ++i) {
                r.marginals_collection()[chain_index][i].LinearCost += msg[i];
            }
        }

        template<typename FACTOR, typename MSG>
        void RepamLeft(FACTOR& l, const MSG& msg)
        {
            assert(msg.size() == l.max_potential_marginals_size());
            std::vector<REAL> m(msg.size());
            for(INDEX i=0; i<msg.size(); ++i) {
                l.max_potential_marginal(i).ReparamCost += msg[i]; 
                m[i] = msg[i];
            }
            l.set_marginal_slack(m);
        }

        template<typename LEFT_FACTOR, typename MSG>
        void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
        {
            std::vector<REAL> m(l.max_potential_marginals_size());
            REAL min = std::numeric_limits<REAL>::max();
            for(INDEX i=0; i<m.size(); ++i) {
                m[i] = l.max_potential_marginal(i).LinearCost + l.max_potential_marginal(i).ReparamCost;
                min = std::min(min, m[i]);
            }
            vector<REAL> mm(m.size());           
            for(INDEX i=0; i<m.size(); ++i) { mm[i] = m[i] - min; }
            msg -= omega*mm; 
        }

        template<typename RIGHT_FACTOR, typename MSG>
        void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
        {
            std::vector<REAL> m = r.messages_to_chain(chain_index); 
            REAL min = *std::min_element(m.begin(), m.end());
            vector<REAL> mm(m.size());
            for(INDEX i=0; i<m.size(); ++i) {
                mm[i] = m[i] - min;
            }
            msg -= omega*mm;
        }

        // template<typename RIGHT_FACTOR, typename MSG_ARRAY>
        // void SendMessagesToLeft(const RIGHT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const REAL omega)
        // {
        //     assert(false);
        // }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
        {
            if(r.max_potential_index(chain_index) != std::numeric_limits<INDEX>::max()) {
                const bool changed = (l.max_potential_index() != r.max_potential_index(chain_index));
                l.set_max_potential_index(r.max_potential_index(chain_index));
                return changed;
            } else {
                return false;
            }
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
        {
            if(l.max_potential_index() != std::numeric_limits<INDEX>::max()) {
                const bool changed = (r.max_potential_index(chain_index) != l.max_potential_index());
                r.set_max_potential_index(chain_index, l.max_potential_index());
                return changed;
            } else {
                return false;
            }
        }

        template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
        bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
        {
            return r.max_potential_index(chain_index) == l.max_potential_index();
        } 

        template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
        void construct_constraints(SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r)
        {
        }

    private:
    const INDEX chain_index;
};

class LabelStateSpace {
    private:
        REAL MaxPotentialUpperBound;
        REAL MaxPotentialLowerBound;
        bool lowerBoundAdded = false;
        bool cleared = false;
        std::map<REAL, REAL> potentials; 

    public:
        LabelStateSpace()
        {}

        LabelStateSpace(REAL maxPotLowerBound, REAL maxPotUpperBound) :
        MaxPotentialLowerBound(maxPotLowerBound), MaxPotentialUpperBound(maxPotUpperBound)
        {
            assert(MaxPotentialLowerBound <= MaxPotentialUpperBound);
        }

        const std::map<REAL, REAL>& GetPotentials() const
        {
            assert(!cleared);
            return potentials;
        }

        void CheckAndInsert(REAL maxPot, REAL linearPot)
        {
            assert(!cleared);
            if (maxPot > MaxPotentialUpperBound)
                return; // This edge has value greater than ub, and we hope to find something <= ub, thus discarding this.

            else if (maxPot < MaxPotentialLowerBound)
            {
                if (!lowerBoundAdded)
                {
                    potentials.insert(std::pair<REAL, REAL>(MaxPotentialLowerBound, linearPot));
                    lowerBoundAdded = true;
                }
                else
                {
                    auto itr = potentials.find(MaxPotentialLowerBound);
                    assert(itr != potentials.end());
                    itr->second = std::min(linearPot, itr->second);
                }
            }

            else
            {
                auto lbKey = potentials.lower_bound(maxPot);
                if (lbKey != potentials.end() || potentials.size() > 0)
                {
                    if (lbKey->first == maxPot) // Current key already exists in the map
                    {
                        lbKey->second = std::min(linearPot, lbKey->second);
                        return;
                    }

                    //Need to insert ALL max pots between lb and ub even if for now they dont seem optimal, they CAN become optimal after reparametrization!
                    // else if (lbKey != potentials.begin()) // There is a key before the current element which needs to be checked with current max pot.
                    // {
                    //     assert(maxPot > std::prev(lbKey)->first);
                    //     if (!(linearPot < std::prev(lbKey)->second)) // Inserting has no benefit as linearPot is not strictly less than previous linear pot.
                    //         return;
                    // }
                }
                potentials.insert(lbKey, std::pair<REAL, REAL>(maxPot, linearPot)); // Insert with hint 
                if (maxPot == MaxPotentialLowerBound)
                    lowerBoundAdded = true;
            }
        }

        REAL GetMaxPotLowerBound() const
        {
            return MaxPotentialLowerBound;
        }

        REAL GetMaxPotUpperBound() const
        {
            return MaxPotentialUpperBound;
        }

        void TakeUnion(LabelStateSpace newStateSpace)
        {
            assert(!isCleared());
            assert(!newStateSpace.isCleared());
            const auto& newPotentials = newStateSpace.GetPotentials();
            for (auto const& currentPot : newPotentials)
            {
                CheckAndInsert(currentPot.first, currentPot.second);
            }
        }

        static LabelStateSpace MergeStateSpaceFromDifferentNeighbours(LabelStateSpace s1, LabelStateSpace s2)
        {
            assert(!s1.isCleared());
            assert(!s2.isCleared());
            const auto& s1Potentials = s1.GetPotentials();
            const auto& s2Potentials = s2.GetPotentials();

            if (s1Potentials.size() == 0)
                return s2;

            if (s2Potentials.size() == 0)
                return s1;

            LabelStateSpace merged = LabelStateSpace(s1.GetMaxPotLowerBound(), s1.GetMaxPotUpperBound());
            assert(s1.GetMaxPotLowerBound() == s2.GetMaxPotLowerBound());
            assert(s1.GetMaxPotUpperBound() == s2.GetMaxPotUpperBound());
            for (auto const& currentS1Pot : s1Potentials)
            {
                for (auto const& currentS2Pot : s2Potentials)
                {
                    merged.CheckAndInsert(std::max(currentS1Pot.first, currentS2Pot.first), currentS1Pot.second + currentS2Pot.second);
                    //Probably can be made more efficient, also validate the logic.
                }
            }
            return merged;
        }

        void ClearPotentials()
        {
            potentials.clear();
            cleared = true;
        }

        bool isCleared() const
        {
            return cleared;
        }

};

class max_potential_on_tree_dynamic_prog {
    public:       
        max_potential_on_tree_dynamic_prog(const three_dimensional_variable_array<REAL>& maxPairwisePotentials, const three_dimensional_variable_array<REAL>& linearPairwisePotentials,
         const std::vector<INDEX>& numLabels, const std::vector<std::array<INDEX, 2>>& messagePassingSchedule, const std::vector<INDEX>& numEdgesForNode)
        :
            LinearPairwisePotentials(linearPairwisePotentials),
            MaxPairwisePotentials(maxPairwisePotentials),
            NumLabels(numLabels),
            NumNodes(numLabels.size()),
            NumEdgesForNode(numEdgesForNode),
            MessagePassingSchedule(messagePassingSchedule)
        {
            assert(maxPairwisePotentials.dim1() + 1 == numLabels.size()); 
            assert(maxPairwisePotentials.dim1() == linearPairwisePotentials.dim1());
            ComputeMaxPotLowerBound();
            solution_.assign(NumNodes, 0);
        }

        // Call this function whenever linear/max potentials get changed so the bounds need to be recomputed.
        void LinearPotsChanged()
        {
            boundsDirty = true;
        }

        REAL GetMaxPotLowerBound()
        {
            return MaxPotentialLowerBoundForAllTrees;
        }

        REAL GetMaxPotUpperBound()
        {
            if (boundsDirty)
            {
                ComputeMaxPotUpperBound();
                boundsDirty = false;
            }
            return MaxPotentialUpperBound;
        }

        REAL LowerBound() const
        {
            Solve();
            return solutionObjective;
            // compute optimal solution and return its cost
        }

        REAL EvaluatePrimal() const
        {
            return solutionObjective;
            // return cost of current solution
        }

        void MaximizePotentialAndComputePrimal() 
        {
            Solve();
            ComputeSolution();
            // compute optimal solution and store it
        }

        INDEX& solution(const INDEX i) { assert(i < solution_.size()); return solution_[i]; }
        INDEX solution(const INDEX i) const { assert(i < solution_.size()); return solution_[i]; }

        void set_max_potential_index(const INDEX index)
        {
            if (max_potential_index_ == index)
                return;

            max_potential_index_ = index;
            ComputeSolution();
        }

        INDEX max_potential_index() const { return max_potential_index_; }
                        
        std::array<REAL,3>& max_potential_marginal(const INDEX i) { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }
        std::array<REAL,3> max_potential_marginal(const INDEX i) const { assert(i < max_potential_marginals_.size()); return max_potential_marginals_[i]; }
        std::vector<std::array<REAL,3>> max_potential_marginals() const { return max_potential_marginals_; }

    private:
        mutable std::vector<INDEX> solution_;
        three_dimensional_variable_array<REAL> MaxPairwisePotentials;
        mutable three_dimensional_variable_array<REAL> LinearPairwisePotentials;
        INDEX NumNodes;
        std::vector<INDEX> NumLabels;
        std::vector<INDEX> NumEdgesForNode;
        
        std::vector<std::array<INDEX, 2>> MessagePassingSchedule;
        static REAL MaxPotentialLowerBoundForAllTrees;    // Computed by max potential message passing.
        mutable REAL MaxPotentialUpperBound;
        mutable REAL LinearPotentialLowerBound; // Computed by conventional message passing.
        // mutable REAL LinearPotentialUpperBound; // Can be computed from max potential message passing and by used to prune paths longer than this bound.
        mutable INDEX max_potential_index_;
        mutable std::vector<std::array<REAL,3>> max_potential_marginals_; // (i) max potential, (ii) minimum linear potential, (iii) cost of configuration 
        mutable bool MaxPotMarginalsInitialized = false;

        bool boundsDirty = true;
        mutable REAL solutionObjective = std::numeric_limits<REAL>::max();

        void ComputeMaxPotLowerBound() const
        {
            std::array<REAL, 2> bounds = MessagePassingForOnePotential(MaxPairwisePotentials, LinearPairwisePotentials, false);
            MaxPotentialLowerBoundForAllTrees = std::max(bounds[0], MaxPotentialLowerBoundForAllTrees); // To take max over all trees.
            // LinearPotentialUpperBound = bounds[1]; //If this is not useful remove it and compute max potential lb only once and store it, as it will not change.
            // Seems like LinearPotentialUpperBound wont help in anything, so we dont need to call this function again and again after updated linear pots.
        }

        void ComputeMaxPotUpperBound() const
        {
            std::array<REAL, 2> bounds = MessagePassingForOnePotential(LinearPairwisePotentials, MaxPairwisePotentials, true);
            LinearPotentialLowerBound = bounds[0];
            MaxPotentialUpperBound = bounds[1];
        }

        // Computes and stores the marginals of the root node.
        void Solve() const
        {
            std::vector<std::vector<LabelStateSpace>> messages(NumNodes);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                messages[i].resize(NumLabels[i]);
                for (INDEX l = 0; l < NumLabels[i]; l++) 
                    messages[i][l] = LabelStateSpace(MaxPotentialLowerBoundForAllTrees, MaxPotentialUpperBound);
            }

            std::vector<INDEX> totalSentAndReceivedMessages(NumNodes, 0);
            std::vector<bool> messageReceived(NumNodes, false);
            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule)
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                bool isTailLeaf = false;
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    isTailLeaf = true;
                    assert(totalSentAndReceivedMessages[tail] == 0);
                }
                for (INDEX lh = 0; lh < NumLabels[head]; lh++)
                {
                    LabelStateSpace lhStateSpace = GetHeadLabelStateSpaceFromCurrentTail(messages[tail], edgeIndex, lh, head < tail, isTailLeaf);
                    LabelStateSpace lhPrevStateSpace = messages[head][lh];
                    messages[head][lh] = LabelStateSpace::MergeStateSpaceFromDifferentNeighbours(lhStateSpace, lhPrevStateSpace);
                }

                totalSentAndReceivedMessages[head]++;
                totalSentAndReceivedMessages[tail]++;

                // Clear the potentials of all nodes except root node to conserve memory as they have 'messaged' their beliefs already.
                if (currentEdge != MessagePassingSchedule.back())
                {
                    if (totalSentAndReceivedMessages[tail] == NumEdgesForNode[tail])
                    {   
                        for (INDEX lt = 0; lt < NumLabels[tail]; lt++)
                            messages[tail][lt].ClearPotentials();
                    }
                    if (totalSentAndReceivedMessages[head] == NumEdgesForNode[head])
                    {
                        for (INDEX lh = 0; lh < NumLabels[head]; lh++)
                            messages[head][lh].ClearPotentials();
                    }
                }
                edgeIndex++;
            }

            // Merge the marginals of all the labels of root node.
            const auto& lastEdge = MessagePassingSchedule.back();
            INDEX rootNode = lastEdge[1];
            LabelStateSpace mergedRootNodeStateSpace;
            for (INDEX l = 0; l < NumLabels[rootNode]; l++)
            {
                mergedRootNodeStateSpace.TakeUnion(messages[rootNode][l]);
            }
            auto rootPots = mergedRootNodeStateSpace.GetPotentials();
            INDEX maxPotIndex = 0;
            for (const auto &currentPair : rootPots)
            {
                if (!MaxPotMarginalsInitialized)
                    max_potential_marginals_.push_back({currentPair.first, currentPair.second, 0});
                else
                {
                    // Due to std::map the rootPots should be already sorted w.r.t increasing max pots:
                    assert(currentPair.first == max_potential_marginals_[maxPotIndex][0]); 
                    max_potential_marginals_[maxPotIndex][1] = currentPair.second;
                }
                maxPotIndex++;
            }
        }

        LabelStateSpace GetHeadLabelStateSpaceFromCurrentTail(const std::vector<LabelStateSpace>& tailMessages, INDEX edgeIndex, INDEX lh, bool isReverse, bool isTailLeaf = false) const 
        {
            LabelStateSpace lhStateSpace(MaxPotentialLowerBoundForAllTrees, MaxPotentialUpperBound);
            for (INDEX lt = 0; lt < tailMessages.size(); lt++)
            {
                // Assuming that the 2D matrix of potentials for each edge is stored always in the order in which tail always counts as l2 and head as l1.
                assert(!isReverse && MaxPairwisePotentials.dim2(edgeIndex) == tailMessages.size() || 
                        isReverse && MaxPairwisePotentials.dim3(edgeIndex) == tailMessages.size());
                assert(!isReverse && LinearPairwisePotentials.dim2(edgeIndex) == tailMessages.size() ||
                        isReverse && LinearPairwisePotentials.dim3(edgeIndex) == tailMessages.size()); 

                REAL edgeMaxPot = !isReverse ? MaxPairwisePotentials(edgeIndex, lt, lh) : MaxPairwisePotentials(edgeIndex, lh, lt);
                REAL edgeLinearPot = !isReverse ? LinearPairwisePotentials(edgeIndex, lt, lh) : LinearPairwisePotentials(edgeIndex, lh, lt);
                if (!isTailLeaf)
                {
                    for (const auto& currentMessageTolt: tailMessages[lt].GetPotentials()) // Iterator over all messages incoming to l1.
                    {
                        lhStateSpace.CheckAndInsert(std::max(currentMessageTolt.second, edgeMaxPot), currentMessageTolt.first + edgeLinearPot);
                    }
                }
                else
                {
                    lhStateSpace.CheckAndInsert(edgeMaxPot, edgeLinearPot);                    
                }
            }
            return lhStateSpace;
        }

        void ComputeSolution() const
        {
            std::vector<std::vector<REAL>> messages(NumNodes);
            std::vector<bool> messageSent(NumNodes, false);
            std::vector<bool> messageReceived(NumNodes, false);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                messages[i].resize(NumLabels[i], 0);
            }

            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule)
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    assert(!messageSent[tail]);
                }
                for (INDEX lh = 0; lh < messages[head].size(); lh++)
                {
                    REAL minMessage = ComputeMessageValue(messages[tail], edgeIndex, lh, head < tail);
                    messages[head][lh] += minMessage;
                }
                messageSent[tail] = true;
                messageReceived[head] = true;
                edgeIndex++;
            }

            bool solutionObjectiveSet = false;

            for (auto it = MessagePassingSchedule.rbegin(); it != MessagePassingSchedule.rend(); ++it)           
            {
                edgeIndex--;
                const auto & currentEdge = *it;
                INDEX tail = currentEdge[1];
                INDEX head = currentEdge[0];
                REAL tailNodeBestMessage = INFINITY;
                for (INDEX lt = 0; lt < messages[tail].size(); lt++)
                {
                    if (messages[tail][lt] < tailNodeBestMessage)
                    {
                        solution_[tail] = lt;
                        tailNodeBestMessage = messages[tail][lt];
                    }
                }
                
                if (!solutionObjectiveSet) // Works for root node.
                {
                    solutionObjective = tailNodeBestMessage + max_potential_marginals_[max_potential_index_][2];
                    solutionObjectiveSet = true;
                }

                for (INDEX lh = 0; lh < messages[head].size(); lh++)
                {
                    REAL minMessage = ComputeMessageValue(messages[tail], edgeIndex, lh, head < tail, solution_[tail]);
                    messages[head][lh] += minMessage;
                }
            }
        }

        REAL ComputeMessageValue(const std::vector<REAL>& tailMessages, INDEX edgeIndex, 
        INDEX lh, bool isReverse, bool labelTail = -1) const
        {
            REAL bestMessageValue = INFINITY;
            INDEX bestlt;
            INDEX ltStart = 0;
            INDEX ltEnd = tailMessages.size() - 1;
            
            if (labelTail >= 0)
            {
                ltStart = labelTail;
                ltEnd = labelTail;
            }

            for (INDEX lt = ltStart; lt <= ltEnd; lt++)
            {
                REAL currentEdgeMaxPot = !isReverse ? MaxPairwisePotentials(edgeIndex, lt, lh) :  MaxPairwisePotentials(edgeIndex, lh, lt);
                if (currentEdgeMaxPot > max_potential_marginals_[max_potential_index_][2])
                    continue;

                REAL currentEdgePot = !isReverse ? LinearPairwisePotentials(edgeIndex, lt, lh) :  LinearPairwisePotentials(edgeIndex, lh, lt);
                REAL currentIncomingMsg = tailMessages[lt];

                REAL msg = currentIncomingMsg + currentEdgePot;
                if (msg < bestMessageValue)
                {
                    bestMessageValue = msg;
                    bestlt = lt;
                }
            }

            return bestMessageValue;
        }

        std::array<REAL, 2> MessagePassingForOnePotential(const three_dimensional_variable_array<REAL>& mainPairwisePots, const three_dimensional_variable_array<REAL>& otherPairwisePots, bool doConventional) const
        {
            std::vector<std::vector<REAL>> mainMessages(NumNodes);
            std::vector<std::vector<REAL>> otherMessages(NumNodes);
            std::vector<bool> messageSent(NumNodes, false);
            std::vector<bool> messageReceived(NumNodes, false);
            for (INDEX i = 0; i < NumNodes; i++)
            {
                mainMessages[i].resize(NumLabels[i], 0);
                otherMessages[i].resize(NumLabels[i], 0);
            }

            INDEX edgeIndex = 0;
            for (const auto & currentEdge : MessagePassingSchedule)
            {
                INDEX tail = currentEdge[0];
                INDEX head = currentEdge[1];
                if (!messageReceived[tail]) // Only leaf nodes can send a message (only one) and only before receiving any.
                {
                    assert(!messageSent[tail]);
                }
                for (INDEX lh = 0; lh < mainMessages[head].size(); lh++)
                {
                    std::array<REAL, 2> minMessage = ComputeMessageValuePair(mainMessages[tail], otherMessages[tail], edgeIndex, lh, mainPairwisePots, otherPairwisePots, doConventional, head < tail);
                    if (doConventional)
                    {
                        mainMessages[head][lh] += minMessage[0];
                        otherMessages[head][lh] = std::max(otherMessages[head][lh], minMessage[1]);
                    }
                    else
                    {
                        mainMessages[head][lh] = std::max(mainMessages[head][lh], minMessage[0]);                                            
                        otherMessages[head][lh] += minMessage[1];
                    }
                }
                messageSent[tail] = true;
                messageReceived[head] = true;
                edgeIndex++;
            }

            const auto& lastEdge = MessagePassingSchedule.back();
            INDEX rootNode = lastEdge[1];
            assert(messageReceived[rootNode]);
            assert(!messageSent[rootNode]);
            REAL mainBound = INFINITY;
            REAL otherBound = INFINITY;

            for (int l = 0; l < NumLabels[rootNode]; l++)
            {
                if (mainMessages[rootNode][l] < mainBound)
                {
                    mainBound = mainMessages[rootNode][l];
                    otherBound = otherMessages[rootNode][l];
                }
            }

            return std::array<REAL, 2>({mainBound, otherBound});
        }

        std::array<REAL, 2> ComputeMessageValuePair(const std::vector<REAL>& tailMainMessages, const std::vector<REAL>& tailOtherMessages, INDEX edgeIndex, INDEX lh, 
        const three_dimensional_variable_array<REAL>& mainPairwisePots, const three_dimensional_variable_array<REAL>& otherPairwisePots, bool doConventional, bool isReverse) const
        {
            REAL bestMainMessage = INFINITY;
            REAL bestOtherMessage = INFINITY;
            INDEX bestlt;
            for (INDEX lt = 0; lt < tailMainMessages.size(); lt++)
            {
                REAL currentEdgePot = !isReverse ? mainPairwisePots(edgeIndex, lt, lh) :  mainPairwisePots(edgeIndex, lh, lt);
                REAL currentIncomingMsg = tailMainMessages[lt];

                if (doConventional)
                {
                    REAL msg = currentIncomingMsg + currentEdgePot;
                    if (msg < bestMainMessage)
                    {
                        bestMainMessage = msg;
                        bestlt = lt;
                    }
                }

                else
                {
                    REAL msg = std::max(currentIncomingMsg, currentEdgePot);
                    if (msg < bestMainMessage)
                    {
                        bestMainMessage = msg;
                        bestlt = lt;
                    }
                }
            }

            if (doConventional)
                bestOtherMessage = std::max(!isReverse ? otherPairwisePots(edgeIndex, bestlt, lh) :  otherPairwisePots(edgeIndex, lh, bestlt), tailOtherMessages[bestlt]);

            else
                bestOtherMessage = (!isReverse ? otherPairwisePots(edgeIndex, bestlt, lh) : otherPairwisePots(edgeIndex, lh, bestlt)) + tailOtherMessages[bestlt];

            return std::array<REAL, 2>({bestMainMessage, bestOtherMessage});
        }

};
REAL max_potential_on_tree_dynamic_prog::MaxPotentialLowerBoundForAllTrees = std::numeric_limits<REAL>::lowest(); //Will need to be re-initialized for a new instance.

// class unary_max_potential_on_chain_message {
//     public:
//         // UNUSED CLASS:
//         unary_max_potential_on_chain_message(const INDEX nodeIndex) : variable(nodeIndex) {}
//         template<typename FACTOR, typename MSG>
//         void RepamRight(FACTOR& r, const MSG& msgs)
//         {
//             if(variable < r.LinearPairwisePotentials.dim1()) {
//                 for(INDEX i=0; i<r.LinearPairwisePotentials.dim2(variable); ++i) {
//                     for(INDEX j=0; j<r.LinearPairwisePotentials.dim3(variable); ++j) {
//                         r.LinearPairwisePotentials(variable,i,j) += msgs[i];
//                     }
//                 }
//             } else {
//                 for(INDEX i=0; i<r.LinearPairwisePotentials.dim2(variable-1); ++i) {
//                     for(INDEX j=0; j<r.LinearPairwisePotentials.dim3(variable-1); ++j) {
//                         r.LinearPairwisePotentials(variable,i,j) += msgs[j];
//                     }
//                 } 
//             }
//             r.invalidate_marginals();
//         }

//         template<typename FACTOR, typename MSG>
//         void RepamLeft(FACTOR& l, const MSG& msgs)
//         {
//             for(INDEX i=0; i<l.size(); ++i) {
//                 l[i] += msgs[i];
//             } 
//         }

//         template<typename LEFT_FACTOR, typename MSG>
//         void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0)
//         {
//             msg -= omega*l;
//         }

//         template<typename RIGHT_FACTOR, typename MSG>
//         void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0)
//         {
//             assert(false);
//         }

//         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
//         bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
//         {
//             if(l.primal() < l.size()) {
//                 const bool changed = (l.primal() != r.solution[variable]);
//                 l.primal() = r.solution[variable];
//                 return changed;
//             } else {
//                 return false;
//             }
//         }

//         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
//         bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
//         {
//             if(r.primal()[variable] < l.size()) {
//                 const bool changed = (l.primal() != r.solution[variable]);
//                 r.solution[variable] = l.primal();
//                 return changed;
//             } else {
//                 return false;
//             }
//         }

//         template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
//         bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
//         {
//             return l.primal() == r.solution[variable];
//         } 

//     private:
//         const INDEX variable;
// };


    
}

#endif // LPMP_HORIZON_TRACKING_FACTORS_HXX
