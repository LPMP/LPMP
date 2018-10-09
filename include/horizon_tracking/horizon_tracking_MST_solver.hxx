#ifndef LPMP_HORIZON_TRACKING_MST_SOLVER_HXX
#define LPMP_HORIZON_TRACKING_MST_SOLVER_HXX

#include "vector.hxx"
#include <queue>
#include <limits>
#include "two_dimensional_variable_array.hxx"
#include "three_dimensional_variable_array.hxx"
#include <unordered_set>
#include <cmath>
#include <sstream>
#include <fstream>

namespace LPMP {
class horizon_tracking_MST_solver {
    struct NodeWithPriority {
        INDEX N;
        INDEX Label;
        REAL Cost; // higher means less priority
    };

private:
    const mrf_input Input;
    std::vector<std::vector<INDEX>> NodeToEdges; // Given a unary index, provides the edges connected to it.
    std::vector<INDEX> Solution;

public:
    horizon_tracking_MST_solver(const mrf_input& input) : Input(input) {
        NodeToEdges.resize(Input.no_variables());
        Solution.resize(Input.no_variables(), std::numeric_limits<INDEX>::max());
        for(INDEX p=0; p<input.no_pairwise_factors(); ++p) {
            auto [i, j] = Input.get_pairwise_variables(p);
            NodeToEdges[i].push_back(p);
            NodeToEdges[j].push_back(p);
        }
    }

    void ComputeSolution() {
        auto cmp = [](const NodeWithPriority& left, const NodeWithPriority& right) { return left.Cost > right.Cost; };
        std::priority_queue<NodeWithPriority, std::vector<NodeWithPriority>, decltype(cmp)> pQ(cmp);
        for (INDEX i = 0; i < Input.no_variables(); i++) {
            if (Input.cardinality(i) > 1) continue;
            pQ.push({i, 0, Input.get_unary(i)[0]});
        }
        if (pQ.empty()) {
            std::runtime_error("No seed present");
            std::exit(0);
        }
        while (!pQ.empty()) {
            const auto bestN = pQ.top();
            pQ.pop();
            if (Solution[bestN.N] < std::numeric_limits<INDEX>::max()) {
                continue; // already set.
            }
            Solution[bestN.N] = bestN.Label;
            for (const INDEX& p : NodeToEdges[bestN.N]) {
                auto [i, j] = Input.get_pairwise_variables(p);
                assert(i!=j);
                assert(i == bestN.N || j == bestN.N);
                INDEX neighbourNode = i == bestN.N ? j : i;
                if (Solution[neighbourNode] < std::numeric_limits<INDEX>::max()) continue;
                // find the best label for current node and enqueue:
                const auto& bestLabelAndCost = ComputeNeighbourBestLabel(bestN.N, p);
                pQ.push({neighbourNode, bestLabelAndCost.first, bestLabelAndCost.second});
            }
        }
        assert(*std::min_element(Solution.begin(), Solution.end()) < std::numeric_limits<INDEX>::max());
    }

    void PrintPrimal() const {
        assert(*std::min_element(Solution.begin(), Solution.end()) < std::numeric_limits<INDEX>::max());
        REAL linearCost = 0;
        REAL bottleneckCost = 0;
        for (INDEX i = 0; i < Input.no_variables(); i++) {
            linearCost += Input.get_unary(i)[Solution[i]];
        }
        for (INDEX p = 0; p < Input.no_pairwise_factors(); p++) {
            const auto [i, j] = Input.get_pairwise_variables(p);
            assert(i < j);
            REAL currentPairwisePot = Input.get_pairwise_potential(p)(Solution[i], Solution[j]);
            linearCost += currentPairwisePot;
            bottleneckCost = std::max(bottleneckCost, currentPairwisePot);
        }

        std::cout<<"\nLinear Cost: "<<linearCost<<std::endl;
        std::cout<<"Bottleneck Cost: "<<bottleneckCost<<std::endl;
    }

    std::pair<INDEX, REAL> ComputeNeighbourBestLabel(const INDEX& root, const INDEX& p) const {
        auto [i, j] = Input.get_pairwise_variables(p);
        assert(i < j);
        assert(i == root || j == root);
        const INDEX neighbour = i == root ? j : i;
        const auto& pairwisePotentials = Input.get_pairwise_potential(p);
        INDEX bestLabel;
        REAL bestCost = std::numeric_limits<REAL>::max();
        for (INDEX ln = 0; ln < Input.cardinality(neighbour); ln++) {
            REAL unaryPot = Input.unaries[neighbour][ln];
            REAL pairwisePot;
            if (i == root)
                pairwisePot = pairwisePotentials(Solution[root], ln);
            else
                pairwisePot = pairwisePotentials(ln, Solution[root]);
            if (bestCost > unaryPot + pairwisePot) {
                bestCost = unaryPot + pairwisePot;
                bestLabel = ln;
            }
        }
        return std::make_pair(bestLabel, bestCost);
    }
    
    void WritePrimal(std::string fileName) const {
        std::stringstream s;
        for(INDEX i=0; i<Solution.size()-1; ++i) {
            s << Solution[i] << ", ";
        }
        s << Solution[Solution.size()-1] << "\n";

        std::string sol = s.str();
        std::ofstream output_file;
        output_file.open(fileName, std::ofstream::out);
        if(!output_file.is_open()) {
            throw std::runtime_error("could not open file " + fileName);
        }
        output_file << sol; 
    }
};
}

#endif // LPMP_HORIZON_TRACKING_MST_SOLVER_HXX