#ifndef LPMP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
#define LPMP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX

#include "horizon_tracking_uai_input.h"
#include "three_dimensional_variable_array.hxx"
#include "grid.hxx"

namespace LPMP {

template<typename MRF_CONSTRUCTOR, typename MAX_MULTIPLE_CHAINS_FACTOR, typename PAIRWISE_MULTIPLE_CHAINS_MESSAGE>
class max_multiple_chains_constructor : public MRF_CONSTRUCTOR {

public:
using FMC = typename MRF_CONSTRUCTOR::FMC;
using mrf_constructor = MRF_CONSTRUCTOR;
using max_multiple_chains_factor_container = MAX_MULTIPLE_CHAINS_FACTOR;
using pairwise_multiple_chains_message_container = PAIRWISE_MULTIPLE_CHAINS_MESSAGE;
using mrf_constructor::mrf_constructor;

std::vector<max_multiple_chains_factor_container*> max_multiple_chains_factors() const { return max_multiple_chains_factors_; }
std::vector<pairwise_multiple_chains_message_container*> pairwise_to_multiple_chain_messages() const { return pairwise_multiple_chains_messages_; }

max_multiple_chains_factor_container* add_multiple_chains_factor(const std::vector<std::vector<INDEX>>& numLabels,
                                    const std::vector<three_dimensional_variable_array<REAL>>& maxPairwisePotentials,
                                    const two_dim_variable_array<INDEX>& chainNodeToOriginalNode,
                                    factor_tree<FMC>* t = nullptr)
{
    INDEX numChains = numLabels.size();
    assert(numChains == maxPairwisePotentials.size());
    std::vector<three_dimensional_variable_array<REAL>> linearPairwisePotentials;
    for (std::size_t c=0; c<maxPairwisePotentials.size();c++) {
        three_dimensional_variable_array<REAL> currentChainLinearPotentials(maxPairwisePotentials[c]);
        for(std::size_t n1=0; n1<maxPairwisePotentials[c].dim1(); ++n1) {
            for(std::size_t l1=0; l1<maxPairwisePotentials[c].dim2(n1); ++l1) {
                for(std::size_t l2=0; l2<maxPairwisePotentials[c].dim3(n1); ++l2) {
                    currentChainLinearPotentials(n1, l1, l2) = 0.0;
                }
            }
        }
        linearPairwisePotentials.push_back(currentChainLinearPotentials);
    }
    auto* multiple_chains_factor = this->lp_->template add_factor<max_multiple_chains_factor_container>(linearPairwisePotentials, maxPairwisePotentials, numLabels);

    assert(numChains == chainNodeToOriginalNode.size());
    for (std::size_t chain_index=0; chain_index<numChains;chain_index++) {
        INDEX pairwise_index=0;
        INDEX node1_index = 0;
        INDEX node2_index = 1;
        for (INDEX i = 0; i < chainNodeToOriginalNode[chain_index].size() - 1; i++, pairwise_index++, node1_index++, node2_index++) {
            const INDEX original_node_1 = chainNodeToOriginalNode[chain_index][i];
            const INDEX original_node_2 = chainNodeToOriginalNode[chain_index][i + 1];
            auto* pairwise_factor = this->get_pairwise_factor(original_node_1, original_node_2);
            auto* msg = this->lp_->template add_message<pairwise_multiple_chains_message_container>
                        (pairwise_factor, multiple_chains_factor, chain_index, pairwise_index, node1_index, node2_index);
            
            pairwise_multiple_chains_messages_.push_back(msg);
            if(t != nullptr) {
                t->add_message(msg, Chirality::right); 
            }
        }
    }
    max_multiple_chains_factors_.push_back(multiple_chains_factor);
    return multiple_chains_factor;
}

void order_factors() const
{
   class cardinality_traversal {
      public:
         cardinality_traversal(std::size_t num_nodes) {
            nodes.resize(num_nodes);
            nodes_covered.resize(num_nodes);
         }

         struct adjacency_list {
            std::size_t num_labels;
            std::vector<std::size_t> neighbors;
         };

         struct node {
            std::size_t node_index;
            std::size_t num_labels;
            bool operator<(const node& o) const { return num_labels > o.num_labels; }
         };

         void add_node(std::size_t index, std::size_t num_labels) {
            nodes[index].num_labels = num_labels;
            nodes_covered[index] = false;
            if (num_labels < min_num_labels)
            {
               min_num_labels = num_labels;
               min_label_node_index = index;
            }
         }

         void add_edge(std::size_t i, std::size_t j) {
            nodes[i].neighbors.push_back(j);
            nodes[j].neighbors.push_back(i);
        }

        std::vector<std::size_t> get_traversal_order() const {
           std::vector<std::size_t> traversal_order;
           std::priority_queue<node, std::vector<node>> nodeIndexAndNodeLabels;

           nodeIndexAndNodeLabels.push({min_label_node_index, min_num_labels});
           nodes_covered[min_label_node_index] = true;

           while (!nodeIndexAndNodeLabels.empty()) {
              auto currentBestNode = nodeIndexAndNodeLabels.top();
              nodeIndexAndNodeLabels.pop();
              traversal_order.push_back(currentBestNode.node_index);

              for (const auto& currentNeighbourIndex : nodes[currentBestNode.node_index].neighbors) {
                 if (nodes_covered[currentNeighbourIndex])
                    continue;

                 nodeIndexAndNodeLabels.push({currentNeighbourIndex, nodes[currentNeighbourIndex].num_labels});
                 nodes_covered[currentNeighbourIndex] = true;
              }
           }
           return traversal_order;
        }

      private:

        std::vector<adjacency_list> nodes;
        mutable std::vector<bool> nodes_covered;
        std::size_t min_label_node_index;
        std::size_t min_num_labels = std::numeric_limits<std::size_t>::max();
   };

    cardinality_traversal traversal(this->get_number_of_variables());
    for(std::size_t i=0; i<this->get_number_of_variables(); ++i) {
        const auto no_labels = this->get_number_of_labels(i);
        traversal.add_node(i, no_labels); 
    }

    for(std::size_t p=0; p<this->get_number_of_pairwise_factors(); ++p) {
        const auto [i,j] = this->get_pairwise_variables(p);
        traversal.add_edge(i,j);
    }

    std::vector<std::size_t> order = traversal.get_traversal_order();
    std::vector<std::size_t> inverse_order(order.size());
    for(std::size_t i=0; i<order.size(); ++i) {
        inverse_order[order[i]] = i;
    }

    for (std::size_t i=0; i<order.size()-1; ++i ) {
        auto* f1 = this->get_unary_factor(order[i]);
        auto* f2 = this->get_unary_factor(order[i+1]);
        this->get_lp()->add_factor_relation(f1,f2);
    }

    for(std::size_t p=0; p<this->get_number_of_pairwise_factors(); ++p) {
        const auto [i,j] = this->get_pairwise_variables(p);
        auto* f_p = this->get_pairwise_factor(p);
        auto f_i = this->get_unary_factor(i);
        auto f_j = this->get_unary_factor(j);
        if (inverse_order[i] < inverse_order[j]) {
            this->get_lp()->add_factor_relation(f_i, f_p);
            this->get_lp()->add_factor_relation(f_p, f_j);
        } else {
            this->get_lp()->add_factor_relation(f_j, f_p);
            this->get_lp()->add_factor_relation(f_p, f_i);
        }
    }
}

private: 
    std::vector<max_multiple_chains_factor_container*> max_multiple_chains_factors_;
    std::vector<pairwise_multiple_chains_message_container*> pairwise_multiple_chains_messages_;
};

template<typename SOLVER, typename HORIZON_TRACKING_CONSTRUCTOR>
void construct_horizon_tracking_problem_on_grid_to_chains(const horizon_tracking_input& input, SOLVER& solver, HORIZON_TRACKING_CONSTRUCTOR& multiple_chain_constructor)
{
    // construct mrf part
    multiple_chain_constructor.construct(input.mrf); 
    auto trees = multiple_chain_constructor.compute_forest_cover();
    for(auto& tree : trees) {
        solver.GetLP().add_tree(tree);
    }

    for(const auto& bottleneck_potential : input.bottleneck_potentials) {
        // we assume that bottleneck potential and mrf potential have same number of variables and in same order. TODO: Check for this!

        // check whether bottleneck potentials come from grid graph and if so, construct horizontal and vertical chains
        auto grid = recognize_grid(bottleneck_potential.pairwise_indices);

        // allocate space for max potentials on chains
        std::vector<three_dimensional_variable_array<REAL>> max_pairwise_potentials_on_chains(grid.number_of_chains());
        std::vector<std::vector<INDEX>> nodeIndices;
        std::vector<std::vector<INDEX>> numLabels(grid.number_of_chains());
        for(std::size_t i=0; i<grid.number_of_chains(); ++i) {
            const auto nodes = grid.chain(i);
            nodeIndices.push_back(nodes);
            assert(nodes.size() > 0);
            std::vector<std::array<std::size_t,2>> function_table_size;
            function_table_size.reserve(nodes.size()-1);
            for(auto node_it = nodes.begin(); node_it!=std::prev(nodes.end()); ++node_it) {
                const std::size_t l1Size = bottleneck_potential.cardinality(*node_it);
                const std::size_t l2Size = bottleneck_potential.cardinality(*std::next(node_it, 1));
                function_table_size.push_back({l1Size, l2Size});
            }
            for(auto node_it = nodes.begin(); node_it!=nodes.end(); ++node_it) {
                numLabels[i].push_back(bottleneck_potential.cardinality(*node_it));
            }
            max_pairwise_potentials_on_chains[i].resize(function_table_size.begin(), function_table_size.end());
        }

        // Populate max potentials
        for(std::size_t i=0; i<bottleneck_potential.no_pairwise_factors(); ++i) {
            const auto pairwise_variables = bottleneck_potential.get_pairwise_variables(i);
            const auto [chain_number, chain_position] = grid.edge_to_chain(pairwise_variables[0], pairwise_variables[1]);
            assert(chain_number < max_pairwise_potentials_on_chains.size());

            auto pairwise_potential = bottleneck_potential.get_pairwise_potential(i);
            for(std::size_t l1=0; l1<bottleneck_potential.cardinality(pairwise_variables[0]); ++l1) {
                for(std::size_t l2=0; l2<bottleneck_potential.cardinality(pairwise_variables[1]); ++l2) {
                    max_pairwise_potentials_on_chains[chain_number](chain_position, l1, l2) = pairwise_potential(l1,l2);
                }
            }
        }

        using FMC = typename SOLVER::FMC;
        factor_tree<FMC> tree;
        two_dim_variable_array<INDEX> nodeIndices2D(nodeIndices);
        multiple_chain_constructor.add_multiple_chains_factor(numLabels, max_pairwise_potentials_on_chains, nodeIndices2D, &tree);
        solver.GetLP().add_tree(tree);
    }
}


} // namespace LPMP 

#endif // LPMP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
