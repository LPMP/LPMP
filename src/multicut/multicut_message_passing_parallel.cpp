#include "cut_base/cut_base_instance.hxx"
#include "cut_base/cut_base_apply_packing.hxx"
#include "multicut/multicut_factors_messages.h"
#include "multicut/multicut_message_passing_parallel.h"
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "graph.hxx"
#include "union_find.hxx"
#include <array>
#include <chrono>
#include <atomic>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <taskflow/taskflow.hpp>
#include "multicut/multicut_cycle_packing_parallel.h"
#include "atomic_helper.h"
#include "dynamic_graph_thread_safe.hxx"
#define ITERATION 5

namespace LPMP {

    void send_weights_to_triplets(edge_item& e, std::vector<triangle_item>& triangle_to_edge) {
        assert(e.nodes[0] < e.nodes[1]);
        auto tmp_cost = e.cost; 
        e.cost.store(0.0);
        auto nr_triangles = e.triangle_indices.size();
        assert(nr_triangles > 0);
        for (auto t: e.triangle_indices){
            std::array<std::pair<size_t, int>,3> sorted_nodes = {std::make_pair(e.nodes[0],0),std::make_pair(e.nodes[1],1),std::make_pair(t[0],2)};
            std::sort(sorted_nodes.begin(), sorted_nodes.end(), [](std::pair<size_t, double> n1, std::pair<size_t, double> n2) {return n1.first < n2.first;});
            std::size_t i=sorted_nodes[0].first, j=sorted_nodes[1].first, k=sorted_nodes[2].first;
            if (sorted_nodes[0].second == 2){
                atomic_addition(triangle_to_edge[t[1]].weights[1], tmp_cost/nr_triangles);
            } else if (sorted_nodes[1].second == 2){
                atomic_addition(triangle_to_edge[t[1]].weights[2], tmp_cost/nr_triangles);
            } else {
                atomic_addition(triangle_to_edge[t[1]].weights[0], tmp_cost/nr_triangles);
            }
        }
    }

    void marginalize(CopyableAtomic<double>& edge_cost, std::array<CopyableAtomic<double>,3>& cost, const int option, const double omega){
        auto marginal = std::min({cost[option]+cost[(option+1)%3], cost[option]+cost[(option+2)%3], cost[0]+cost[1]+cost[2]}) 
                      - std::min(0.0, cost[(option+1)%3]+cost[(option+2)%3]);
        atomic_addition(edge_cost, omega*marginal);
        cost[option] = cost[option] - omega*marginal;
    }

    void send_triplets_to_edge(triangle_item& t, std::vector<edge_item>& edge_to_triangle){
        // ij: 0 jk: 1 ik:2
        marginalize(edge_to_triangle[t.edge_indices[0]].cost, t.weights, 0, 1.0/3.0);
        marginalize(edge_to_triangle[t.edge_indices[2]].cost, t.weights, 2, 1.0/2.0);
        marginalize(edge_to_triangle[t.edge_indices[1]].cost, t.weights, 1, 1.0/1.0);
        marginalize(edge_to_triangle[t.edge_indices[0]].cost, t.weights, 0, 1.0/2.0);
        marginalize(edge_to_triangle[t.edge_indices[2]].cost, t.weights, 2, 1.0/1.0);
        marginalize(edge_to_triangle[t.edge_indices[0]].cost, t.weights, 0, 1.0/1.0);
        auto marginal = std::min({t.weights[0]+t.weights[1], t.weights[0]+t.weights[2], t.weights[0]+t.weights[1]+t.weights[2]}) 
                      - std::min(0.0, t.weights[1]+t.weights[2]);
        if (marginal > 1e-8) std::cout << "Incorrect marginal.\n";       
    }

    double compute_lower_bound(std::vector<edge_t>& other_edges, std::vector<edge_item>& edge_to_triangle, 
                               std::vector<triangle_item>& triangle_to_edge){
        double lb = 0.0;
        for(auto e: other_edges)
            lb += std::min(0.0, e.cost);
        for (auto e: edge_to_triangle)
            lb += std::min(0.0, e.cost.load());
        for (auto&t:  triangle_to_edge)
            lb += std::min({0.0, t.weights[0]+t.weights[1], t.weights[0]+t.weights[2], t.weights[1]+t.weights[2],
                            t.weights[0]+t.weights[1]+t.weights[2]});
        return lb;
    }

    void send_weights_to_triplets_parallel(tf::Taskflow& taskflow, std::vector<edge_item>& edge_to_triangle, std::vector<triangle_item>& triangle_to_edge, 
        const int nr_threads){
        auto send = taskflow.for_each_index(0, nr_threads, 1, [&](const std::size_t thread_no){
            const std::size_t batch_size = edge_to_triangle.size()/nr_threads + 1;
            int first_edge = thread_no*batch_size;
            int last_edge = std::min((thread_no+1)*batch_size, edge_to_triangle.size());
            for(auto i = first_edge; i != last_edge; ++i){
                send_weights_to_triplets(edge_to_triangle[i], triangle_to_edge);
            }
        });
    }

    void send_triplets_to_edge_parallel(tf::Taskflow& taskflow, std::vector<edge_item>& edge_to_triangle, std::vector<triangle_item>& triangle_to_edge, 
        const int nr_threads){
         auto send = taskflow.for_each_index(0, nr_threads, 1, [&](const std::size_t thread_no){
            const std::size_t batch_size = triangle_to_edge.size()/nr_threads + 1;
            int first_triangle = thread_no*batch_size;
            int last_triangle = std::min((thread_no+1)*batch_size, triangle_to_edge.size());
            for(auto i = first_triangle; i != last_triangle; ++i){
                send_triplets_to_edge(triangle_to_edge[i], edge_to_triangle);
            }
        });
    }
    

    std::pair<multicut_instance, double> multicut_message_passing_parallel(const multicut_instance& input, const bool record_cycles, const int nr_threads)
    {
        const auto begin_time = std::chrono::steady_clock::now();
        tf::Executor executor(nr_threads);
        tf::Taskflow taskflow;
        std::size_t no_nodes = input.no_nodes();
        std::cout << "Message Passing\n";

        // Cycle Packing
        cycle_packing cp;
        std::vector<triangle_item> triangle_to_edge;
        std::vector<edge_item> edge_to_triangle;
        std::vector<edge_t> other_edges;

        cp = cycle_packing_triangulation_parallel(input, nr_threads, triangle_to_edge, edge_to_triangle, other_edges);
        const auto MP_begin_time = std::chrono::steady_clock::now();
        double lower_bound;

        std::cout << "#Edges in triangle: " << edge_to_triangle.size() << std::endl;

        // Message Passing
        for (int i=0; i < ITERATION; ++i){
            taskflow.clear();

            send_weights_to_triplets_parallel(taskflow, edge_to_triangle, triangle_to_edge, nr_threads);
            executor.run(taskflow);
            executor.wait_for_all();

            taskflow.clear();
            lower_bound = compute_lower_bound(other_edges, edge_to_triangle, triangle_to_edge);
            std::cout << "Lower bound after MP step 1: " << lower_bound << std::endl;

            send_triplets_to_edge_parallel(taskflow, edge_to_triangle, triangle_to_edge, nr_threads);
            executor.run(taskflow);
            executor.wait_for_all();
            lower_bound = compute_lower_bound(other_edges, edge_to_triangle, triangle_to_edge);
            std::cout << "Lower bound after MP step 2: " << lower_bound << std::endl;
        }

        const auto MP_end_time = std::chrono::steady_clock::now();
        std::cout << "MP took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(MP_end_time - MP_begin_time).count() << " milliseconds\n";

        multicut_instance final_instance;
        for(auto e: other_edges)
            final_instance.add_edge(e[0], e[1], e.cost);
        for (auto e: edge_to_triangle)
            final_instance.add_edge(e.nodes[0], e.nodes[1], e.cost);

        return std::make_pair(final_instance, lower_bound);
    }

}
