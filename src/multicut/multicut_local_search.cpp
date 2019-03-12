#include "multicut/multicut_local_search.h"
#include "multicut/multicut_greedy_additive_edge_contraction.h"
#include "union_find.hxx"
#include "sequence_compression.h"
#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

namespace LPMP {

    multicut_local_search::multicut_local_search(const multicut_instance& instance, const multicut_node_labeling& labeling)
        : instance_(instance)
    {
        construct(instance, labeling);
    }

    void multicut_local_search::construct(const multicut_instance& instance, const multicut_node_labeling& labeling)
    {
        g_.construct(instance.edges().begin(), instance.edges().end(), 
                [](const auto& e) { return e.cost[0]; } 
                );

        std::cout << "initial objective = " << instance.evaluate(labeling) << "\n";

        total_edge_sum_.clear();
        total_edge_sum_.resize(instance.no_nodes(), 0.0);
        same_component_edge_sum_.clear();
        same_component_edge_sum_.resize(instance.no_nodes(), 0.0);
        labeling_ = labeling; 
        assert(labeling_.size() == g_.no_nodes());

        // TODO: make more efficient by doing below computations in one pass
        for(std::size_t i=0; i<g_.no_nodes(); ++i)
            for(auto edge_it=g_.begin(i); edge_it!=g_.end(i); ++edge_it)
                total_edge_sum_[i] += edge_it->edge();

        for(std::size_t i=0; i<g_.no_nodes(); ++i)
            same_component_edge_sum_[i] = compute_cut_value(i, labeling_[i]);

        lower_bound_ = instance.evaluate(labeling_);
        empty_label_ = *std::max_element(labeling_.begin(), labeling_.end()) + 1;
    }

    // compute cost of moving i to partition k
    double multicut_local_search::compute_cut_value(const std::size_t i, const std::size_t k) const
    {
        double cut_cost = 0.0;
        for(auto edge_it_=g_.begin(i); edge_it_!=g_.end(i); ++edge_it_)
            if(labeling_[edge_it_->head()] == k)
                cut_cost += edge_it_->edge();
        return cut_cost; 
    }

    void multicut_local_search::move(const std::size_t i, const std::size_t new_label)
    {
        assert(i < g_.no_nodes());
        assert(new_label != labeling_[i]);
        assert(std::abs(same_component_edge_sum_[i] - compute_cut_value(i, labeling_[i])) < 1e-8);
        const std::size_t prev_label = labeling_[i];
        labeling_[i] = new_label;
        // TODO: make in one go
        same_component_edge_sum_[i] = compute_cut_value(i, labeling_[i]);; 
        for(auto edge_it=g_.begin(i); edge_it!=g_.end(i); ++edge_it) {
            const std::size_t j = edge_it->head();
            if(labeling_[j] == prev_label)
                same_component_edge_sum_[j] -= edge_it->edge();
            if(labeling_[j] == new_label)
                same_component_edge_sum_[j] += edge_it->edge();
        }

        empty_label_ = std::max(new_label+1, empty_label_);
    }

    double multicut_local_search::move_cost(const std::size_t i, const std::size_t new_label) const
    {
        assert(std::abs(same_component_edge_sum_[i] - compute_cut_value(i, labeling_[i])) < 1e-8);
        const double cut_value = compute_cut_value(i, new_label);
        return same_component_edge_sum_[i] - cut_value;
    }

    double multicut_local_search::perform_1_swaps()
    {
        const double prev_lower_bound_ = lower_bound_;
        std::size_t last_update = 0;
        for(std::size_t iter=0; ; ++iter) {
            if(iter - last_update >= g_.no_nodes())
                break;
            const std::size_t i = iter%g_.no_nodes();
            // check whether it is profitable to make node an isolated component. If so, do so
            if(same_component_edge_sum_[i] < 0.0) {
                lower_bound_ += same_component_edge_sum_[i];
                move(i, iter + g_.no_nodes());
                last_update = iter;
                assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
            }
            // check whether sum of all cut edges is greater than sum of uncut edges. Only if so, proceed
            // idea: comput all cut costs by going over edges at once
            const double cut_edge_sum = total_edge_sum_[i] - same_component_edge_sum_[i];
            if(true || cut_edge_sum > same_component_edge_sum_[i]) {
                for(auto edge_it=g_.begin(i); edge_it!=g_.end(i); ++edge_it) {
                    const std::size_t j = edge_it->head();
                    if(labeling_[i] == labeling_[j])
                        continue;
                    const double delta = move_cost(i, labeling_[j]); // TODO: reuse cut_cost by storing and checking whether next element belongs to same component and no move occurred
                    if(delta < 0.0) {
                        lower_bound_ += delta;
                        move(i, labeling_[j]);
                        last_update = iter; 
                        assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
                    } 
                }
            }
        } 
        std::cout << "improvement in 1 swaps = " << lower_bound_ - prev_lower_bound_ << "\n";
        assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
        return lower_bound_ - prev_lower_bound_;
    }

    double multicut_local_search::compute_2_swap_cost(const std::size_t i, const std::size_t j, const double edge_cost) const
    {
        assert(std::abs(same_component_edge_sum_[i] - compute_cut_value(i, labeling_[i])) < 1e-8);
        assert(std::abs(same_component_edge_sum_[j] - compute_cut_value(j, labeling_[j])) < 1e-8);
        if(labeling_[i] == labeling_[j]) { // compute cost to make new partition out of i and j
            return same_component_edge_sum_[i] + same_component_edge_sum_[j] - 2*edge_cost;
        } else { // compute cost to exchange node labels of j and j
            const double cut_i = compute_cut_value(i, labeling_[j]);
            const double cut_j = compute_cut_value(j, labeling_[i]); 
            return same_component_edge_sum_[i] - cut_i + same_component_edge_sum_[j] - cut_j + 2*edge_cost;
        } 
    }

    double multicut_local_search::perform_2_swaps()
    {
        std::size_t last_update = 0;
        std::size_t iter = 0;
        const double prev_lower_bound_ = lower_bound_;
        for(std::size_t i=0; i<g_.no_nodes(); ++i) {
            for(auto edge_it=g_.begin(i); edge_it!=g_.end(i); ++edge_it, ++iter) {
                if(iter - last_update >= 2*g_.no_edges())
                    break;
                const std::size_t j = edge_it->head();
                const double swap_cost = compute_2_swap_cost(i, j, edge_it->edge());
                if(swap_cost < 0.0) {
                    lower_bound_ += swap_cost;
                    // TODO: reuse computations
                    if(labeling_[i] == labeling_[j]) { // not nice!
                        move(i, iter + g_.no_nodes());
                        move(j, iter + g_.no_nodes());
                        assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
                    } else {
                        const std::size_t prev_label_i = labeling_[i];
                        move(i, labeling_[j]);
                        move(j, prev_label_i);
                        assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
                    }
                    last_update = iter;
                }
            }
        }
        std::cout << "improvement in 2 swaps = " << lower_bound_ - prev_lower_bound_ << "\n";
        assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
        return lower_bound_ - prev_lower_bound_;
    }

    std::pair<double, std::array<std::size_t,3>> multicut_local_search::compute_triangle_swap_cost(std::array<std::size_t,3> nodes, std::array<double,3> edge_costs) const
    {
        auto evaluate_labeling = [&](const auto& c, const auto& labeling) {
            double current_cost = 0.0;
            constexpr std::size_t n = std::tuple_size<std::remove_reference_t<decltype(labeling)>>::value;
            for(std::size_t i=0; i<n; ++i)
                for(std::size_t j=i+1; j<n; ++j)
                    if(labeling[i] != labeling[j]) {
                        assert(std::isfinite(c[i][j]));
                        current_cost += c[i][j];
                    }
            return current_cost;
        }; 

        if(labeling_[nodes[0]] == labeling_[nodes[1]] && labeling_[nodes[1]] == labeling_[nodes[2]]) {
            // check wether putting all nodes in new components is advantageous. Putting one or two in new component should be detected already
            const double delta = same_component_edge_sum_[nodes[0]] + same_component_edge_sum_[nodes[1]] + same_component_edge_sum_[nodes[2]] - (edge_costs[0] + edge_costs[1] + edge_costs[2]);
            return {delta, {empty_label(), empty_label(), empty_label()}};
        } else if(labeling_[nodes[0]] == labeling_[nodes[1]] || labeling_[nodes[0]] == labeling_[nodes[2]] || labeling_[nodes[1]] == labeling_[nodes[2]]) { // two different labels present
            const std::size_t l1 = labeling_[nodes[0]];
            const std::size_t l2 = labeling_[nodes[0]] == labeling_[nodes[1]] ? labeling_[nodes[2]] : labeling_[nodes[1]];
            assert(l1 != l2);
            const std::size_t no_l1_nodes = 1 + (labeling_[nodes[1]] == l1) + (labeling_[nodes[2]] == l1);
            const std::size_t no_l2_nodes = (labeling_[nodes[1]] == l2) + (labeling_[nodes[2]] == l2);
            assert(no_l1_nodes >= 1 && no_l1_nodes <= 2);
            assert(no_l2_nodes >= 1 && no_l2_nodes <= 2);
            assert(no_l1_nodes + no_l2_nodes == 3);

            std::array<std::size_t, 5> current_labeling = {labeling_[nodes[0]], labeling_[nodes[1]], labeling_[nodes[2]], l1, l2};
            std::array<std::array<double,5>,5> c;

            c[0][1] = edge_costs[0];
            c[0][2] = edge_costs[1];
            c[1][2] = edge_costs[2];

            if(labeling_[nodes[0]] == labeling_[nodes[1]]) {
                c[0][3] = same_component_edge_sum_[nodes[0]] - edge_costs[0];
                c[1][3] = same_component_edge_sum_[nodes[1]] - edge_costs[0];
                c[0][4] = compute_cut_value(nodes[0], l2) - edge_costs[1];
                c[1][4] = compute_cut_value(nodes[1], l2) - edge_costs[2];
                c[2][3] = compute_cut_value(nodes[2], l1) - edge_costs[1] - edge_costs[2];
                c[2][4] = same_component_edge_sum_[nodes[2]];
            } else if(labeling_[nodes[0]] == labeling_[nodes[2]]) {
                assert(labeling_[nodes[0]] == l1 && labeling_[nodes[2]] == l1 && labeling_[nodes[1]] == l2);
                c[0][3] = same_component_edge_sum_[nodes[0]] - edge_costs[1];
                c[2][3] = same_component_edge_sum_[nodes[2]] - edge_costs[1];
                c[0][4] = compute_cut_value(nodes[0], l2) - edge_costs[0];
                c[2][4] = compute_cut_value(nodes[2], l2) - edge_costs[2];
                c[1][3] = compute_cut_value(nodes[1], l1) - edge_costs[0] - edge_costs[2];
                c[1][4] = same_component_edge_sum_[nodes[1]];
            } else {
                assert(labeling_[nodes[1]] == l2 && labeling_[nodes[1]] == l2);
                c[1][3] = compute_cut_value(nodes[1], l1) - edge_costs[0];
                c[2][3] = compute_cut_value(nodes[2], l1) - edge_costs[1];
                c[1][4] = same_component_edge_sum_[nodes[1]] - edge_costs[2];
                c[2][4] = same_component_edge_sum_[nodes[2]] - edge_costs[2];
                c[0][3] = same_component_edge_sum_[nodes[0]];
                c[0][4] = compute_cut_value(nodes[0], l2) - edge_costs[0] - edge_costs[1];
            }

            c[3][4] = 0.0;

            const double current_cost = evaluate_labeling(c, current_labeling);
            // solve the multicut problem with given costs by enumeration
            double best_cost = current_cost;
            std::array<std::size_t, 5> best_labeling = current_labeling;

            std::array<std::size_t,5> labeling = current_labeling;
            std::array<std::size_t,3> possible_labels = {l1,l2,empty_label()};
            for(const std::size_t l0 : possible_labels) {
                for(const std::size_t l1 : possible_labels) {
                    for(const std::size_t l2 : possible_labels) {
                        labeling[0] = l0;
                        labeling[1] = l1;
                        labeling[2] = l2;
                        double cost = evaluate_labeling(c,labeling);
                        if(cost < best_cost) {
                            best_cost = cost;
                            best_labeling = labeling;
                        }
                    }
                }
            }
            assert(best_cost <= current_cost);
            return {best_cost - current_cost, {best_labeling[0], best_labeling[1], best_labeling[2]}};

        } else if(labeling_[nodes[0]] != labeling_[nodes[1]] && labeling_[nodes[0]] != labeling_[nodes[2]] && labeling_[nodes[1]] != labeling_[nodes[2]]) {
            // first three nodes are from triangles
            // last three nodes are the partitions that are connected to triangle
            std::array<std::size_t, 6> current_labeling = {labeling_[nodes[0]], labeling_[nodes[1]], labeling_[nodes[2]], labeling_[nodes[0]], labeling_[nodes[1]], labeling_[nodes[2]]};

            std::array<std::array<double,6>,6> c;
            c[0][1] = edge_costs[0];
            c[0][2] = edge_costs[1];
            c[1][2] = edge_costs[2];

            c[0][3] = same_component_edge_sum_[nodes[0]];
            c[1][4] = same_component_edge_sum_[nodes[1]];
            c[2][5] = same_component_edge_sum_[nodes[2]];

            // TODO: eliminate double counting for same labeling value!
            c[0][4] = compute_cut_value(nodes[0], labeling_[nodes[1]]) - edge_costs[0];
            c[0][5] = compute_cut_value(nodes[0], labeling_[nodes[2]]) - edge_costs[1];

            c[1][3] = compute_cut_value(nodes[1], labeling_[nodes[0]]) - edge_costs[0];
            c[1][5] = compute_cut_value(nodes[1], labeling_[nodes[2]]) - edge_costs[2];

            c[2][3] = compute_cut_value(nodes[2], labeling_[nodes[0]]) - edge_costs[1];
            c[2][4] = compute_cut_value(nodes[2], labeling_[nodes[1]]) - edge_costs[2];

            if(labeling_[nodes[0]] == labeling_[nodes[1]]) {
                c[3][4] = std::numeric_limits<double>::infinity();
            } else {
                c[3][4] = 0.0;
            }

            if(labeling_[nodes[0]] == labeling_[nodes[2]]) {
                c[3][5] = std::numeric_limits<double>::infinity();
            } else {
                c[3][5] = 0.0;
            }

            if(labeling_[nodes[1]] == labeling_[nodes[2]]) {
                c[4][5] = std::numeric_limits<double>::infinity();
            } else {
                c[4][5] = 0.0;
            }

            const double current_cost = evaluate_labeling(c,current_labeling);

            // solve the multicut problem with given costs by enumeration
            double best_cost = current_cost;
            std::array<std::size_t, 6> best_labeling = current_labeling;

            std::array<std::size_t,4> possible_labels = {current_labeling[0], current_labeling[1], current_labeling[2], empty_label()};
            std::array<std::size_t,6> labeling = current_labeling;
            for(const std::size_t l0 : possible_labels) {
                for(const std::size_t l1 : possible_labels) {
                    for(const std::size_t l2 : possible_labels) {
                        labeling[0] = l0;
                        labeling[1] = l1;
                        labeling[2] = l2;
                        double cost = evaluate_labeling(c,labeling);
                        if(cost < best_cost) {
                            best_cost = cost;
                            best_labeling = labeling;
                        }
                    }
                }
            }
            assert(best_cost <= current_cost);
            return {best_cost - current_cost, {best_labeling[0], best_labeling[1], best_labeling[2]}};
        }
        assert(false);
        return {0.0, {labeling_[nodes[0]], labeling_[nodes[1]], labeling_[nodes[2]]}};
    }

    double multicut_local_search::perform_3_swaps()
    {
        double prev_lower_bound_ = lower_bound_;
        std::size_t no_triangles = 0;
        g_.for_each_triangle( [&](const std::size_t i, const std::size_t j, const std::size_t k, const double c01, const double c02, const double c12) {
                const auto [delta, improved_labeling] = compute_triangle_swap_cost({i, j, k}, {c01, c02, c12});
                if(delta < 0.0) {
                if(labeling_[i] != improved_labeling[0])
                move(i, improved_labeling[0]);
                if(labeling_[j] != improved_labeling[1])
                move(j, improved_labeling[1]);
                if(labeling_[k] != improved_labeling[2])
                move(k, improved_labeling[2]);

                lower_bound_ += delta;
                assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
                }
                ++no_triangles; 
                });
        std::cout << "improvement in 3 swaps = " << lower_bound_ - prev_lower_bound_ << "\n";
        assert(std::abs(instance_.evaluate(labeling_) - lower_bound_) < 1e-8);
        return lower_bound_ - prev_lower_bound_;
    }

    void multicut_local_search::perform_swaps()
    {
        double prev_lower_bound = lower_bound_;
        std::bitset<4> perform_swap;
        perform_swap.set();
        while(perform_swap.count() > 0) {
            if(perform_swap[0]) {
                if(perform_1_swaps() < -1e-8)
                    perform_swap.set();
                perform_swap[0] = false; 
            } else if(perform_swap[1]) {
                if(perform_2_swaps() < -1e-8)
                    perform_swap.set();
                perform_swap[1] = false; 
            } else if(perform_swap[2]) {
                if(perform_3_swaps() < -1e-8)
                    perform_swap.set();
                perform_swap[2] = false; 
            } else if(perform_swap[3]) {
                if(perform_joins() < -1e-8)
                    perform_swap.set();
                perform_swap[3] = false;
            } 
        }
        std::cout << "total improvement = " <<  lower_bound_ - prev_lower_bound << "\n";
        perform_joins();
    }

    double multicut_local_search::perform_joins()
    {
        // possibly put reduction directly into gaec code
        double prev_lower_bound = lower_bound_;
        // build instance of connected components and solve it with greedy additive edge contraction
        union_find uf(g_.no_nodes());
        g_.for_each_edge([&](const std::size_t i, const std::size_t j, const auto& edge_info) {
                if(labeling_[i] == labeling_[j])
                uf.merge(i,j);
                });

        auto edge = [](const std::size_t i, const std::size_t j) -> std::array<std::size_t,2> {
            return {std::min(i,j), std::max(i,j)};
        };

        tsl::robin_map<std::array<std::size_t,2>, double> connectivity_graph;
        g_.for_each_edge([&](const std::size_t i, const std::size_t j, const double edge_cost) {
                if(labeling_[i] != labeling_[j]) {
                assert(!uf.connected(i,j));
                const auto compressed_edge = edge(uf.find(i), uf.find(j));
                auto edge_it = connectivity_graph.find(compressed_edge);
                if(edge_it == connectivity_graph.end())
                connectivity_graph.insert(std::make_pair(compressed_edge, edge_cost));
                else
                edge_it.value() += edge_cost;
                }
                });

        sequence_compression connectivity_graph_node_compression(g_.no_nodes()); // TODO: reuse
        for(auto edge_it : connectivity_graph) {
            const std::size_t i = edge_it.first[0];
            connectivity_graph_node_compression.add_index(i);
            const std::size_t j = edge_it.first[1];
            connectivity_graph_node_compression.add_index(j);
        } 
        multicut_instance reduced_instance;
        for(auto edge_it : connectivity_graph) {
            const std::size_t i = connectivity_graph_node_compression.orig_to_compressed_index( edge_it.first[0] );
            const std::size_t j = connectivity_graph_node_compression.orig_to_compressed_index( edge_it.first[1] );
            reduced_instance.add_edge(i, j, edge_it.second);
        }

        multicut_edge_labeling improved_sol_e = greedy_additive_edge_contraction(reduced_instance);
        multicut_node_labeling improved_sol = improved_sol_e.transform_to_node_labeling(reduced_instance); 

        // transform from reduced graph solution to original one
        multicut_node_labeling improved_sol_full;
        improved_sol_full.reserve(labeling_.size());
        for(std::size_t i=0; i<g_.no_nodes(); ++i)
            improved_sol_full.push_back( improved_sol[connectivity_graph_node_compression.orig_to_compressed_index(uf.find(i))] );

        lower_bound_ = instance_.evaluate(labeling_);
        std::cout << "improvement in joins = " <<  lower_bound_ - prev_lower_bound << "\n";
        return lower_bound_ - prev_lower_bound;

    }

    multicut_node_labeling multicut_local_search::get_labeling() const 
    { 
        return labeling_; 
    };

    std::size_t multicut_local_search::empty_label() const
    {
        return empty_label_;
    }

} // namespace LPMP
