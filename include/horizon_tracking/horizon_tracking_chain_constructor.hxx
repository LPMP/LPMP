#ifndef LPMP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
#define LPMP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX

#include "horizon_tracking_uai_input.h"
#include "three_dimensional_variable_array.hxx"
#include "grid.hxx"

namespace LPMP {

template<typename MRF_CONSTRUCTOR, typename MAX_CHAIN_FACTOR, typename MAX_POTENTIAL_FACTOR, typename PAIRWISE_MAX_CHAIN_MESSAGE, typename MAX_CHAIN_MAX_POTENTIAL_MESSAGE>
class max_chain_constructor : public MRF_CONSTRUCTOR {

public:
using FMC = typename MRF_CONSTRUCTOR::FMC;
using mrf_constructor = MRF_CONSTRUCTOR;
using max_chain_factor_container = MAX_CHAIN_FACTOR;
using max_potential_factor_container = MAX_POTENTIAL_FACTOR;
using pairwise_max_factor_message_container = PAIRWISE_MAX_CHAIN_MESSAGE;
using max_chain_max_potential_message_container = MAX_CHAIN_MAX_POTENTIAL_MESSAGE;
using mrf_constructor::mrf_constructor;

std::vector<max_chain_factor_container*> max_chain_factors() const {
    return max_chain_factors_;
}

std::vector<max_potential_factor_container*> max_graph_factors() const {
    return max_graph_factors_;
}

template <typename ITERATOR>
max_chain_factor_container* add_max_chain(ITERATOR node_var_begin, ITERATOR node_var_end,
                                    const three_dimensional_variable_array<REAL>& maxPairwisePotentials,
                                    const std::size_t chainIndex,
                                    factor_tree<FMC>* t = nullptr)
{
    std::vector<INDEX> num_labels;
    for(auto it = node_var_begin; it!=node_var_end; ++it) {
        const INDEX i = (*it);
        num_labels.push_back( this->get_number_of_labels(i) );
    }
    three_dimensional_variable_array<REAL> linearPairwisePotentials(maxPairwisePotentials);
    for(std::size_t n1=0; n1<maxPairwisePotentials.dim1(); ++n1) {
        for(std::size_t l1=0; l1<maxPairwisePotentials.dim2(n1); ++l1) {
            for(std::size_t l2=0; l2<maxPairwisePotentials.dim3(n1); ++l2) {
                linearPairwisePotentials(n1, l1, l2) = 0.0;
            }
        }
    }
    auto* chain_factor = this->lp_->template add_factor<max_chain_factor_container>(maxPairwisePotentials, linearPairwisePotentials, num_labels, chainIndex);

    INDEX pairwise_index=0;
    INDEX node1_index = 0;
    INDEX node2_index = 1;
    for(auto it = node_var_begin; std::next(it, 1)!=node_var_end; ++it, ++pairwise_index, ++node1_index, ++node2_index) {
        const INDEX i = (*it);
        const INDEX j = *std::next(it, 1);
        auto* pairwise_factor = this->get_pairwise_factor(i,j);
        auto* msg = this->lp_->template add_message<pairwise_max_factor_message_container>(pairwise_factor, chain_factor, pairwise_index, node1_index, node2_index);

        if(t != nullptr) {
            t->add_message(msg, Chirality::right); 
        }
    }
    max_chain_factors_.push_back(chain_factor);
    return chain_factor;
}

template<typename ITERATOR>
max_potential_factor_container* add_max_potential(ITERATOR max_chain_begin, ITERATOR max_chain_end, factor_tree<FMC>* t = nullptr)
{
   // ugly: we first build a std::vector<std::vector<..>> and then convert it to two_dim_variable_array
    std::vector<std::vector<max_linear_costs>> all_marginals;
    for(auto max_chain_it = max_chain_begin; max_chain_it!=max_chain_end; ++max_chain_it) {
        auto* f = (*max_chain_it)->get_factor();
        f->MaximizePotentialAndComputePrimal();
        std::vector<max_linear_costs> current_chain_marginals_max;
        for (INDEX i = 0; i < f->max_potential_marginals_size(); i++) {
            current_chain_marginals_max.push_back({f->max_potential_marginal(i).MaxCost, f->max_potential_marginal(i).LinearCost});
        }
        all_marginals.push_back(current_chain_marginals_max);
    }

    auto* max_factor = this->lp_->template add_factor<max_potential_factor_container>(all_marginals);
    for(auto max_chain_it = max_chain_begin; max_chain_it!=max_chain_end; ++max_chain_it) {
        const auto chain_index = std::distance(max_chain_begin, max_chain_it);
        auto* current_chain = *max_chain_it;

        auto* msg = this->lp_->template add_message<max_chain_max_potential_message_container>(current_chain, max_factor, chain_index);

        if(t != nullptr) {
            t->add_message(msg, Chirality::right);
        } 
    }
    max_graph_factors_.push_back(max_factor);
    return max_factor;
}

private: 
    std::vector<max_chain_factor_container*> max_chain_factors_;
    std::vector<max_potential_factor_container*> max_graph_factors_;
};

template<typename SOLVER, typename HORIZON_TRACKING_CONSTRUCTOR>
void construct_horizon_tracking_problem_on_grid_to_chains(const horizon_tracking_input& input, SOLVER& solver, HORIZON_TRACKING_CONSTRUCTOR& chain_constructor)
{
    // construct mrf part
    chain_constructor.construct(input.mrf); 
    auto trees = chain_constructor.compute_forest_cover();
    for(auto& tree : trees) {
        solver.GetLP().add_tree(tree);
    }

    for(const auto& bottleneck_potential : input.bottleneck_potentials) {
        // we assume that bottleneck potential and mrf potential have same number of variables and in same order. TODO: Check for this!

        // check whether bottleneck potentials come from grid graph and if so, construct horizontal and vertical chains
        auto grid = recognize_grid(bottleneck_potential.pairwise_indices);

        // allocate space for max potentials on chains
        std::vector<three_dimensional_variable_array<REAL>> max_pairwise_potentials_on_chains(grid.number_of_chains());
        for(std::size_t i=0; i<grid.number_of_chains(); ++i) {
            const auto nodes = grid.chain(i);
            assert(nodes.size() > 0);
            std::vector<std::array<std::size_t,2>> function_table_size;
            function_table_size.reserve(nodes.size()-1);
            for(auto node_it = nodes.begin(); node_it!=std::prev(nodes.end()); ++node_it) {
                const std::size_t l1Size = bottleneck_potential.cardinality(*node_it);
                const std::size_t l2Size = bottleneck_potential.cardinality(*std::next(node_it, 1));
                function_table_size.push_back({l1Size, l2Size});
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

        // add chain potentials
        using FMC = typename SOLVER::FMC;
        factor_tree<FMC> tree;
        std::vector<typename std::remove_reference_t<decltype(chain_constructor)>::max_chain_factor_container*> max_chain_potentials;
        max_chain_potentials.reserve(grid.number_of_chains());

        for(std::size_t i=0; i<grid.number_of_chains(); ++i) {
            const auto nodes = grid.chain(i);
            auto* f = chain_constructor.add_max_chain(nodes.begin(), nodes.end(), max_pairwise_potentials_on_chains[i], i, &tree);
            max_chain_potentials.push_back(f);
        }

        auto* f = chain_constructor.add_max_potential(max_chain_potentials.begin(), max_chain_potentials.end(), &tree);
        solver.GetLP().add_tree(tree);
    }
}

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

template<typename MRF_CONSTRUCTOR>
void order_nodes_by_label_space_cadinality(MRF_CONSTRUCTOR& mrf)
{
    cardinality_traversal traversal(mrf.get_number_of_variables());
    for(std::size_t i=0; i<mrf.get_number_of_variables(); ++i) {
        const auto no_labels = mrf.get_number_of_labels(i);
        traversal.add_node(i, no_labels); 
    }

    for(std::size_t p=0; p<mrf.get_number_of_pairwise_factors(); ++p) {
        const auto [i,j] = mrf.get_pairwise_variables(p);
        traversal.add_edge(i,j);
    }

    std::vector<std::size_t> order = traversal.get_traversal_order();
    std::vector<std::size_t> inverse_order(order.size());
    for(std::size_t i=0; i<order.size(); ++i) {
        inverse_order[order[i]] = i;
    }

    for (std::size_t i=0; i<order.size()-1; ++i ) {
        auto* f1 = mrf.get_unary_factor(order[i]);
        auto* f2 = mrf.get_unary_factor(order[i+1]);
        mrf.get_lp()->add_factor_relation(f1,f2);
    }

    for(std::size_t p=0; p<mrf.get_number_of_pairwise_factors(); ++p) {
        const auto [i,j] = mrf.get_pairwise_variables(p);
        auto* f_p = mrf.get_pairwise_factor(p);
        auto f_i = mrf.get_unary_factor(i);
        auto f_j = mrf.get_unary_factor(j);
        if (inverse_order[i] < inverse_order[j]) {
            mrf.get_lp()->add_factor_relation(f_i, f_p);
            mrf.get_lp()->add_factor_relation(f_p, f_j);
        } else {
            mrf.get_lp()->add_factor_relation(f_j, f_p);
            mrf.get_lp()->add_factor_relation(f_p, f_i);
        }
    }
}
} // namespace LPMP 

#endif // LPMP_HORIZON_TRACKING_CHAIN_CONSTRUCTOR_HXX
