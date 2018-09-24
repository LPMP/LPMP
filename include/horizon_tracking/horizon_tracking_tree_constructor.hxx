// #ifndef LPMP_HORIZON_TRACKING_TREE_CONSTRUCTOR_HXX
// #define LPMP_HORIZON_TRACKING_TREE_CONSTRUCTOR_HXX

// #include "horizon_tracking_uai_input.h"
// #include "three_dimensional_variable_array.hxx"
// #include "spanning_trees.hxx"
// #include "horizon_tracking_util.hxx"

// namespace LPMP {

// template<typename MRF_CONSTRUCTOR, typename MAX_TREE_FACTOR, typename MAX_POTENTIAL_FACTOR, typename PAIRWISE_MAX_TREE_MESSAGE, typename MAX_TREE_MAX_POTENTIAL_MESSAGE>
// class max_tree_constructor : public MRF_CONSTRUCTOR {
// public:
// using FMC = typename MRF_CONSTRUCTOR::FMC;
// using mrf_constructor = MRF_CONSTRUCTOR;
// using max_tree_factor_container = MAX_TREE_FACTOR;
// using max_potential_factor_container = MAX_POTENTIAL_FACTOR;
// using pairwise_max_factor_message_container = PAIRWISE_MAX_TREE_MESSAGE;
// using max_tree_max_potential_message_container = MAX_TREE_MAX_POTENTIAL_MESSAGE;
// using mrf_constructor::mrf_constructor;

// std::vector<max_tree_factor_container*> max_tree_factors() const { return max_tree_factors_; }
// std::vector<max_potential_factor_container*> max_graph_factors() const { return max_graph_factors_; }

// max_tree_factor_container* add_max_tree(const std::vector<std::array<std::size_t,2>>& messagePassingSchedule,
//                                     const three_dimensional_variable_array<double>& maxPairwisePotentials,
//                                     std::vector<std::size_t>& numLabels, factor_tree<FMC>* t = nullptr)
// {
//     three_dimensional_variable_array<double> linearPairwisePotentials(maxPairwisePotentials);
//     for(std::size_t e=0; e<linearPairwisePotentials.dim1(); ++e) {
//         for(std::size_t l1=0; l1<linearPairwisePotentials.dim2(e); ++l1) {
//             for(std::size_t l2=0; l2<linearPairwisePotentials.dim3(e); ++l2) {
//                 linearPairwisePotentials(e, l1, l2) = 0.0;
//             }
//         }
//     }

//     auto* tree_factor = this->lp_->template add_factor<max_tree_factor_container>(maxPairwisePotentials, linearPairwisePotentials, numLabels, messagePassingSchedule);

//     std::size_t pairwise_index=0;
//     for(const auto& currentEdge : messagePassingSchedule) {
//         std::size_t i = fmin(currentEdge[0], currentEdge[1]);
//         std::size_t j = fmax(currentEdge[0], currentEdge[1]);
//         auto* pairwise_factor = this->get_pairwise_factor(i, j);
//         auto* msg = this->lp_->template add_message<pairwise_max_factor_message_container>(pairwise_factor, tree_factor, pairwise_index, i, j);

//         if(t != nullptr) {
//             t->add_message(msg, Chirality::right); 
//         }
//         pairwise_index++;
//     }
//     max_tree_factors_.push_back(tree_factor);
//     return tree_factor;
// }

// template<typename ITERATOR>
// max_potential_factor_container* add_max_potential(ITERATOR max_tree_begin, ITERATOR max_tree_end, factor_tree<FMC>* t = nullptr)
// {
//     std::vector<std::vector<max_linear_costs>> all_marginals;
//     for(auto max_tree_it = max_tree_begin; max_tree_it!=max_tree_end; ++max_tree_it) {
//         auto* f = (*max_tree_it)->get_factor();
//         f->MaximizePotentialAndComputePrimal();
//         std::vector<max_linear_rep_costs> current_tree_marginals = f->max_potential_marginals();
//         std::vector<max_linear_costs> current_tree_marginals_max;
//         for (auto current_marginal_item : current_tree_marginals) {
//             current_tree_marginals_max.push_back({current_marginal_item.MaxCost, current_marginal_item.LinearCost});   // Ignoring the third column in the first iteration. 
//         }
//         all_marginals.push_back(current_tree_marginals_max);
//     }

//     auto* max_factor = this->lp_->template add_factor<max_potential_factor_container>(all_marginals);
//     for(auto max_tree_it = max_tree_begin; max_tree_it!=max_tree_end; ++max_tree_it) {
//         const auto tree_index = std::distance(max_tree_begin, max_tree_it);
//         auto* current_tree = *max_tree_it;

//         auto* msg = this->lp_->template add_message<max_tree_max_potential_message_container>(current_tree, max_factor, tree_index);

//         if(t != nullptr) {
//             t->add_message(msg, Chirality::right);
//         } 
//     }
//     max_graph_factors_.push_back(max_factor);
//     return max_factor;
// }

// private: 
//     std::vector<max_tree_factor_container*> max_tree_factors_;
//     std::vector<max_potential_factor_container*> max_graph_factors_;
// };

// template<typename SOLVER, typename HORIZON_TRACKING_CONSTRUCTOR>
// void construct_horizon_tracking_problem_on_grid_to_trees(const horizon_tracking_input& input, SOLVER& solver, HORIZON_TRACKING_CONSTRUCTOR& tree_constructor)
// {
//     // construct mrf part
//     tree_constructor.construct(input.mrf);

//     auto trees = tree_constructor.compute_forest_cover();
//     for(auto& tree : trees) {
//         solver.GetLP().add_tree(tree);
//     }
//     std::size_t numNodes = input.mrf.no_variables();
//     std::vector<std::size_t> numLabels(numNodes);
//     // we assume that bottleneck potential and mrf potential have same number of variables and in same order. 
//     for (std::size_t n = 0; n < numNodes; n++) { numLabels[n] = input.mrf.cardinality(n); }

//     for(const auto& bottleneck_potential : input.bottleneck_potentials) {
//         spanning_trees sp_trees(input.mrf);

//         // allocate space for max potentials on trees
//         std::vector<three_dimensional_variable_array<double>> max_pairwise_potentials_on_trees(sp_trees.number_trees());

//         for(std::size_t t=0; t<sp_trees.number_trees(); ++t) {
//             std::vector<std::array<std::size_t,2>> function_table_size;

//             auto current_tree_schedule = sp_trees.get_schedule(t);
//             for (const auto& current_pairwise : current_tree_schedule) {
//                 std::size_t pairwise_index = sp_trees.get_pairwise_index(current_pairwise[0], current_pairwise[1]);
//                 auto l1Size = bottleneck_potential.cardinality(std::min(current_pairwise[0], current_pairwise[1]));
//                 auto l2Size = bottleneck_potential.cardinality(std::max(current_pairwise[0], current_pairwise[1]));
//                 function_table_size.push_back({l1Size, l2Size});
//             }
//             max_pairwise_potentials_on_trees[t].resize(function_table_size.begin(), function_table_size.end());
//         }

//         for(std::size_t t=0; t<sp_trees.number_trees(); ++t) {
//             std::size_t p = 0;
//             auto current_tree_schedule = sp_trees.get_schedule(t);
//             for (const auto& current_pairwise : current_tree_schedule) {
//                 std::size_t pairwise_index = sp_trees.get_pairwise_index(current_pairwise[0], current_pairwise[1]);
//                 auto l1Size = bottleneck_potential.cardinality(std::min(current_pairwise[0], current_pairwise[1]));
//                 auto l2Size = bottleneck_potential.cardinality(std::max(current_pairwise[0], current_pairwise[1]));
//                 auto pairwise_potential = bottleneck_potential.get_pairwise_potential(pairwise_index);

//                 for(std::size_t l1=0; l1<l1Size; ++l1) {
//                     for(std::size_t l2=0; l2<l2Size; ++l2) {
//                         max_pairwise_potentials_on_trees[t](p, l1, l2) = pairwise_potential(l1,l2);
//                     }
//                 }
//                 p++;
//             }
//         }
//        // add tree factors
//         using FMC = typename SOLVER::FMC;
//         factor_tree<FMC> tree;
//         std::vector<typename std::remove_reference_t<decltype(tree_constructor)>::max_tree_factor_container*> max_tree_factors;
//         max_tree_factors.reserve(sp_trees.number_trees());

//         for(std::size_t t=0; t<sp_trees.number_trees(); ++t) {
//             auto* f = tree_constructor.add_max_tree(sp_trees.get_schedule(t), max_pairwise_potentials_on_trees[t], numLabels, &tree);
//             max_tree_factors.push_back(f);
//         }

//         auto* f = tree_constructor.add_max_potential(max_tree_factors.begin(), max_tree_factors.end(), &tree);
//         solver.GetLP().add_tree(tree);
//     }
// }
// }


// #endif // LPMP_HORIZON_TRACKING_TREE_CONSTRUCTOR_HXX
