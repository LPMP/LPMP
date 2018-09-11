#ifndef LPMP_SPANNING_TREE_HXX
#define LPMP_SPANNING_TREE_HXX

#include "arboricity.h"
#include "mrf/mrf_input.h"

namespace LPMP {
class spanning_trees {
    
private:
mrf_input input_;
std::vector<std::vector<std::array<INDEX, 2>>> message_passing_schedules_;
std::map<std::pair<INDEX, INDEX>, INDEX> edge_to_pairwise_entry_;

void Compute()
{
    UndirectedGraph g = UndirectedGraph(input_.no_variables(), input_.no_pairwise_factors());

    INDEX idx = 0;
    for(auto e : input_.pairwise_indices) {
        const auto i = e[0];
        const auto j = e[1];
        g.AddEdge(i,j,1);
    }

    const INDEX forest_num = g.Solve();
    std::vector<std::vector<std::array<INDEX,2>>> spTreesPairwiseIndices(forest_num);

    std::vector<int> forestEdges(input_.no_variables());

    for (INDEX k=0; k < forest_num; k++) {     
        g.GetForestEdges(k, forestEdges.data());
        union_find uf(input_.no_variables());

        for (const auto& currentEdge: forestEdges) {
            if (currentEdge < 0)
                continue;
            
            spTreesPairwiseIndices[k].push_back(input_.pairwise_indices[currentEdge]);
            uf.merge(input_.pairwise_indices[currentEdge][0], input_.pairwise_indices[currentEdge][1]);
        }

    // Convert the forest to spanning tree by iterating over the edges not in the forest and add which does not create a cycle:
        for(auto e : input_.pairwise_indices) {
            const auto i = e[0];
            const auto j = e[1];
            if (uf.find(i) != uf.find(j))
            {
                spTreesPairwiseIndices[k].push_back(e);
                uf.merge(i, j);
            }
        }
    }

    for (auto const &currentTree : spTreesPairwiseIndices)
    {
        //Assure that the trees have the correct size to be a spanning tree.
        std::map<INDEX, std::vector<INDEX>> currentAdjacency; // Node-Node Adjacency relations
        assert(currentTree.size() == input_.no_variables() - 1);
        for (const auto& currentEdge : currentTree)
        {
            INDEX n1 = currentEdge[0];
            INDEX n2 = currentEdge[1];
            if (currentAdjacency.count(n1) == 0) {
                std::vector<INDEX> neighborVector{n2};
                currentAdjacency.insert(std::pair<INDEX, std::vector<INDEX>>(n1, neighborVector));
            }
            else currentAdjacency[n1].push_back(n2);

            if (currentAdjacency.count(n2) == 0) {
                std::vector<INDEX> neighborVector{n1};
                currentAdjacency.insert(std::pair<INDEX, std::vector<INDEX>>(n2, neighborVector));
            }
            else currentAdjacency[n2].push_back(n1);
        }
        message_passing_schedules_.push_back(compute_leaf_to_root_order(currentAdjacency));   
    }
}

std::vector<std::array<INDEX, 2>> compute_leaf_to_root_order(const std::map<INDEX, std::vector<INDEX>>& adjacency)
{
    std::queue<INDEX> nodeQ;
    std::vector<INDEX> nodeStatus(adjacency.size(), 0); // 0 -> Not processed, 1-> In queue, 2-> Popped from queue
    for (const auto& currentNode : adjacency) {
        if (currentNode.second.size() == 1) { // Find leaf nodes and add into queue
            nodeQ.push(currentNode.first);
            nodeStatus[currentNode.first] = 1;
        }
    }
    std::vector<std::array<INDEX, 2>> order;
    while (nodeQ.size() > 0) {
        INDEX currentN = nodeQ.front();
        nodeQ.pop();
        nodeStatus[currentN] = 2;
        std::vector<INDEX> neighbors = adjacency.at(currentN);
        for (const auto& neighbourNode : neighbors) {
            if (nodeStatus[neighbourNode] == 2) // Node already processed and should have already added the edge=(neighbour->currentNode)
                continue;

            order.push_back({currentN, neighbourNode});
            if (nodeStatus[neighbourNode] == 0) {
                nodeQ.push(neighbourNode);
                nodeStatus[neighbourNode] = 1;
            }
        }
    }
    return order;
}

//   struct pair_hash {
//     template <class T1, class T2>
//     std::size_t operator () (const std::pair<T1,T2> &p) const {
//         auto h1 = std::hash<T1>{}(p.first);
//         auto h2 = std::hash<T2>{}(p.second);

//         // Mainly for demonstration purposes, i.e. works but is overly simple
//         // In the real world, use sth. like boost.hash_combine
//         return h1 ^ h2;  
//     }
//   };
//   std::vector<std::array<INDEX, 2>> compute_leaf_to_root_order(const std::map<INDEX, std::vector<INDEX>>& adjacency, std::map<INDEX, INDEX>& nodeDegrees)
//   {
//     std::unordered_set<INDEX> leaves;
//     for (const auto& currentNode : adjacency)
//     {
//         nodeDegrees.emplace(currentNode.first, currentNode.second.size());
//         if (currentNode.second.size() == 1)
//             leaves.insert(currentNode.first);
//     }
//     std::queue<INDEX> nodeQ;
//     // std::unordered_set<INDEX> coveredNodes(leaves);

//     for (const auto& currentL : leaves) {
//         nodeQ.push(currentL);
//     }
//     std::unordered_set<std::pair<INDEX, INDEX>,pair_hash> coveredEdges;
//     std::vector<std::array<INDEX, 2>> order;
//     while (nodeQ.size() > 0) {
//         INDEX currentN = nodeQ.front();
//         nodeQ.pop();
//         std::vector<INDEX> neighbors = adjacency.at(currentN);
//         for (const auto& neighbourNode : neighbors) {
//             if (coveredEdges.count({currentN, neighbourNode}) > 0 || 
//                 coveredEdges.count({neighbourNode, currentN}) > 0 )
//                 continue;

//             order.push_back({currentN, neighbourNode});
//             coveredEdges.insert({currentN, neighbourNode});
//             // if (coveredNodes.count(neighbourNode) > 0) // Neighbouring node has already added its incident edges which will already include the edge {currentN, neighbourNode}.
//             //     continue;

//             nodeQ.push(neighbourNode);
//             // coveredNodes.insert(neighbourNode);
//         }
//     }
//     return order;
//   }
public:
spanning_trees(mrf_input input) : input_(input) {
    for (INDEX p = 0; p < input.no_pairwise_factors(); p++) {
        std::array<std::size_t,2> current_pairwise = input.get_pairwise_variables(p);
        edge_to_pairwise_entry_.insert(std::make_pair(std::make_pair(current_pairwise[0], current_pairwise[1]), p));
    }
    Compute();
}

INDEX number_trees() const {
    return message_passing_schedules_.size();
}

std::vector<std::array<INDEX, 2>> get_schedule(INDEX tree_index) const {
    assert(tree_index < message_passing_schedules_.size());
    return message_passing_schedules_[tree_index];
}

INDEX get_pairwise_index(INDEX n1, INDEX n2) const {
    if (edge_to_pairwise_entry_.count(std::make_pair(n1, n2))) {
        return edge_to_pairwise_entry_.at(std::make_pair(n1, n2));
    }

    else if (edge_to_pairwise_entry_.count(std::make_pair(n2, n1))) {
        return edge_to_pairwise_entry_.at(std::make_pair(n2, n1));
    }
    assert(false); // Edge should be present somewhere
}
};
}

#endif // LPMP_SPANNING_TREE_HXX
