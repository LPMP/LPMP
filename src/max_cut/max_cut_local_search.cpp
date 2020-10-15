#include "max_cut/max_cut_local_search.h"
#include <cassert>
#include <bitset>

namespace LPMP {

    max_cut_local_search::max_cut_local_search(const max_cut_instance& instance, const max_cut_node_labeling& labeling)
        : instance_(instance),
        lower_bound(instance_.evaluate(labeling))
    {
        g.construct(instance.edges().begin(), instance.edges().end(), 
                [](const auto& e) { return e.cost[0]; } 
                );

        assert(labeling.size() == instance.no_nodes());
        label.reserve(instance.no_nodes()); 
        for(const auto l : labeling ) 
            label.push_back(l == 1 ? 1 : 0);

        cut_values.reserve(instance.no_nodes());
        for(std::size_t i=0; i<g.no_nodes(); ++i)
            cut_values.push_back(swap_1_cost(i));

    }

    double max_cut_local_search::swap_1_cost(const std::size_t i) const
    {
        double cost = 0.0;
        for(auto edge_it=g.begin(i); edge_it!=g.end(i); ++edge_it) {
            const std::size_t j = edge_it->head();
            const double sign = label[i] == label[j] ? +1.0 : -1.0;
            cost += sign * edge_it->edge();
        }
        return cost;
    }

    void max_cut_local_search::swap(const std::size_t i)
    {
        assert(label[i] == 0 || label[i] == 1);
        const char orig_i_label = label[i];
        label[i] = 1 - label[i];
        for(auto edge_it=g.begin(i); edge_it!=g.end(i); ++edge_it) {
            const std::size_t j = edge_it->head();
            const double sign = orig_i_label == label[j] ? +2.0 : -2.0;
            const double cost = edge_it->edge();
            cut_values[i] -= sign*cost;
            cut_values[j] -= sign*cost;
            assert(std::abs(swap_1_cost(j) - cut_values[j]) <= 1e-8);
        }
        assert(std::abs(swap_1_cost(i) - cut_values[i]) <= 1e-8);
    }

    double max_cut_local_search::swap_2_cost(const std::size_t i, const std::size_t j, const double edge_cost) const
    {
        if(label[i] == label [j])
            return cut_values[i] + cut_values[j] - 2*edge_cost;
        else
            return cut_values[i] + cut_values[j] + 2*edge_cost;
    }

    double max_cut_local_search::perform_1_swaps()
    {
        double prev_lower_bound = lower_bound;
        std::size_t last_swap = 0;
        for(std::size_t iter=0;; ++iter) {
            if(iter - last_swap >= g.no_nodes())
                break;
            std::size_t i=iter%g.no_nodes();
            const double delta = cut_values[i];
            if(delta < 0.0) {
                lower_bound += delta;
                last_swap = iter;
                swap(i);
            }
        }
        std::cout << "1 swaps improvement = " << prev_lower_bound - lower_bound << "\n";
        return lower_bound - prev_lower_bound;
    }

    double max_cut_local_search::perform_2_swaps()
    {
        double prev_lower_bound = lower_bound;
        std::size_t last_swap = 0;
        std::size_t iter = 0;
        for(std::size_t i=0; i<g.no_nodes(); ++i) {
            for(auto edge_it=g.begin(i); edge_it!=g.end(i); ++edge_it, ++iter) {
                if(iter - last_swap >= 2*g.no_edges())
                    break;

                const double delta = swap_2_cost(i, edge_it->head(), edge_it->edge());
                if(delta < 0.0) {
                    lower_bound += delta;
                    swap(i); 
                    swap(edge_it->head());
                }
            }
        }
        std::cout << "2 swaps improvement = " << prev_lower_bound - lower_bound << "\n";
        return lower_bound - prev_lower_bound;
    }

    std::tuple<double, std::array<std::size_t,3>> max_cut_local_search::best_3_swap(const std::array<std::size_t,3>& nodes, const std::array<double,3>& edge_costs) const
    {
        assert(std::abs(cut_values[nodes[0]] - swap_1_cost(nodes[0])) < 1e-8);
        assert(std::abs(cut_values[nodes[1]] - swap_1_cost(nodes[1])) < 1e-8);
        assert(std::abs(cut_values[nodes[2]] - swap_1_cost(nodes[2])) < 1e-8);

        std::array<std::size_t,4> labeling{label[nodes[0]], label[nodes[1]], label[nodes[2]], 0};

        std::array<std::array<double,4>,4> c;

        c[0][1] = edge_costs[0];
        c[0][2] = edge_costs[1];
        c[1][2] = edge_costs[2];

        if(labeling[0] == 0 && labeling[1] == 0 && labeling[2] == 0) {
            c[0][3] = +cut_values[nodes[0]] - edge_costs[0] - edge_costs[1];
            c[1][3] = +cut_values[nodes[1]] - edge_costs[0] - edge_costs[2];
            c[2][3] = +cut_values[nodes[2]] - edge_costs[1] - edge_costs[2];
        } else if(labeling[0] == 0 && labeling[1] == 0 && labeling[2] == 1) {
            c[0][3] = +cut_values[nodes[0]] - edge_costs[0] + edge_costs[1];
            c[1][3] = +cut_values[nodes[1]] - edge_costs[0] + edge_costs[2];
            c[2][3] = -(cut_values[nodes[2]] + edge_costs[1] + edge_costs[2]);
        } else if(labeling[0] == 0 && labeling[1] == 1 && labeling[2] == 0) {
            c[0][3] = +cut_values[nodes[0]] + edge_costs[0] - edge_costs[1];
            c[1][3] = -(cut_values[nodes[1]] + edge_costs[0] + edge_costs[2]);
            c[2][3] = +cut_values[nodes[2]] - edge_costs[1] + edge_costs[2];
        } else if(labeling[0] == 0 && labeling[1] == 1 && labeling[2] == 1) {
            c[0][3] = +cut_values[nodes[0]] + edge_costs[0] + edge_costs[1];
            c[1][3] = -(cut_values[nodes[1]] + edge_costs[0] - edge_costs[2]);
            c[2][3] = -(cut_values[nodes[2]] + edge_costs[1] - edge_costs[2]);
        } else if(labeling[0] == 1 && labeling[1] == 0 && labeling[2] == 0) {
            c[0][3] = -(cut_values[nodes[0]] + edge_costs[0] + edge_costs[1]);
            c[1][3] = +cut_values[nodes[1]] + edge_costs[0] - edge_costs[2];
            c[2][3] = +cut_values[nodes[2]] + edge_costs[1] - edge_costs[2];
        } else if(labeling[0] == 1 && labeling[1] == 0 && labeling[2] == 1) {
            c[0][3] = -(cut_values[nodes[0]] + edge_costs[0] - edge_costs[1]);
            c[1][3] = +cut_values[nodes[1]] + edge_costs[0] + edge_costs[2];
            c[2][3] = -(cut_values[nodes[2]] - edge_costs[1] + edge_costs[2]);
        } else if(labeling[0] == 1 && labeling[1] == 1 && labeling[2] == 0) {
            c[0][3] = -(cut_values[nodes[0]] - edge_costs[0] + edge_costs[1]);
            c[1][3] = -(cut_values[nodes[1]] - edge_costs[0] + edge_costs[2]);
            c[2][3] = +cut_values[nodes[2]] + edge_costs[1] + edge_costs[2];
        } else {
            assert(labeling[0] == 1 && labeling[1] == 1 && labeling[2] == 1);
            c[0][3] = -cut_values[nodes[0]] + edge_costs[0] + edge_costs[1];
            c[1][3] = -cut_values[nodes[1]] + edge_costs[0] + edge_costs[2];
            c[2][3] = -cut_values[nodes[2]] + edge_costs[1] + edge_costs[2];
        }

        // TODO: make more efficient by only holding upper triangular part of c
        auto evaluate_labeling = [&](const auto& cost, const auto labeling) {
            double c = 0.0;
            for(std::size_t i=0; i<4; ++i) {
                for(std::size_t j=i+1; j<4; ++j) {
                    if(labeling[i] != labeling[j]) {
                        c += cost[i][j];
                    }
                }
            }
            return c;
        };

        const double current_cost = evaluate_labeling(c,labeling);

        double best_cost = current_cost;
        std::array<std::size_t,3> best_labeling{labeling[0], labeling[1], labeling[2]};

        const std::array<std::size_t,2> labelset{0,1};
        for(const auto l0 : labelset) {
            for(const auto l1 : labelset) {
                for(const auto l2 : labelset) {
                    const double cost = evaluate_labeling(c, std::array<std::size_t,4>{l0,l1,l2,0});
                    if(cost < best_cost) {
                        best_cost = cost;
                        best_labeling = {l0,l1,l2};
                    }
                }
            }
        }

        assert(best_cost <= current_cost);
        return {best_cost - current_cost, best_labeling};
    }

    double max_cut_local_search::perform_3_swaps()
    {
        double prev_lower_bound = lower_bound;
        g.for_each_triangle([&](const std::size_t i, const std::size_t j, const std::size_t k, const double e01, const double e02, const double e12) {
                const auto [delta, l] = best_3_swap({i,j,k}, {e01, e02, e12});
                if(delta < 0.0) {
                lower_bound += delta;
                if(label[i] != l[0])
                swap(i);
                if(label[j] != l[1])
                swap(j);
                if(label[k] != l[2])
                swap(k);
                }
                assert(std::abs(instance_.evaluate(label) - lower_bound) < 1e-8);
                });

        std::cout << "3 swaps improvement = " << prev_lower_bound - lower_bound << "\n";
        return lower_bound - prev_lower_bound;
    }

    double max_cut_local_search::perform_swaps()
    {
        double prev_lower_bound = lower_bound;
        std::bitset<3> actions;
        actions.set();

        while(actions.count() > 0) {
            if(actions[2]) {
                const double improvement = perform_3_swaps();
                if(improvement < -1e-8)
                    actions.set();
                actions[2] = false;
            } else if(actions[1]) {
                const double improvement = perform_2_swaps();
                if(improvement < -1e-8)
                    actions.set();
                actions[1] = false;
            } else if(actions[0]) {
                const double improvement = perform_1_swaps();
                if(improvement < -1e-8)
                    actions.set();
                actions[0] = false;
            }
        }

        return lower_bound - prev_lower_bound;
    }

    max_cut_node_labeling max_cut_local_search::get_labeling() const
    {
        max_cut_node_labeling output;
        output.reserve(g.no_nodes());
        for(const auto l : label)
            output.push_back(l == 1 ? 1 : 0);
        return output;
    }

} // namespace LPMP
