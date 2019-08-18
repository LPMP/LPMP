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
#include "copyable_atomic.h"
#include "dynamic_graph_thread_safe.hxx"

namespace LPMP {

    void send_weights_to_triplets(std::size_t node1, std::size_t node2, CopyableAtomic<double> cost, edge_to_triangle_map& M, multicut_triangle_factor& T){
        assert(node1 < node2);
        if (M.find({node1,node2}) == M.end()){
   //         std::cout << "send_weights_to_triplets: no triangle for edge " << e.first[0] << " " << e.first[1] << std::endl;
        } else {
            auto nr_triangles = M[{node1,node2}].size();
            for (auto third: M[{node1,node2}]){
                std::array<std::pair<size_t, int>,3> nodes = {std::make_pair(node1,0),std::make_pair(node2,1),std::make_pair(third,2)};
                std::sort(nodes.begin(), nodes.end(), [](std::pair<size_t, double> n1, std::pair<size_t, double> n2) {return n1.first < n2.first;});
                std::size_t i=nodes[0].first, j=nodes[1].first, k=nodes[2].first;
                if (T.find({i,j,k}) == T.end()) {
                    std::cout << "send_weights_to_triplets: triangle" << i << " "<< j << " " << k << " missing. \n";
                } else {
                    if (nodes[0].second == 2){
                        T[{i,j,k}][1].store(T[{i,j,k}][1]+cost/nr_triangles);
                    } else if (nodes[1].second == 2){
                        T[{i,j,k}][2].store(T[{i,j,k}][2]+cost/nr_triangles);
                    } else {
                        T[{i,j,k}][0].store(T[{i,j,k}][0]+cost/nr_triangles);
                    }
                }
            }
        }
        
    }

    void marginalize(CopyableAtomic<double>& edge_cost, std::array<CopyableAtomic<double>,3>& cost, int option, double omega){
        auto marginal = std::min({cost[option]+cost[(option+1)%3], cost[option]+cost[(option+2)%3], cost[0]+cost[1]+cost[2]}) 
                      - std::min(0.0, cost[(option+1)%3]+cost[(option+2)%3]);
        edge_cost = edge_cost + omega*marginal;
        cost[option] = cost[option] - omega*marginal;
    }

    void send_triplets_to_edge(atomic_edge_container& edges, std::size_t i, std::size_t j, std::size_t k, std::array<CopyableAtomic<double>,3>& triangle_cost){
        // ij: 0 jk: 1 ik:2
        marginalize(edges[{i, j}], triangle_cost, 0, 1.0/3.0);
        marginalize(edges[{i, k}], triangle_cost, 2, 1.0/2.0);
        marginalize(edges[{j, k}], triangle_cost, 1, 1.0/1.0);
        marginalize(edges[{i, j}], triangle_cost, 0, 1.0/2.0);
        marginalize(edges[{i, k}], triangle_cost, 2, 1.0/1.0);
        marginalize(edges[{i, j}], triangle_cost, 0, 1.0/1.0);
        auto marginal = std::min({triangle_cost[0]+triangle_cost[1], triangle_cost[0]+triangle_cost[2], triangle_cost[0]+triangle_cost[1]+triangle_cost[2]}) 
                      - std::min(0.0, triangle_cost[2]+triangle_cost[3]);
        // assert(marginal <= tolerance);
    }

    double compute_lower_bound(atomic_edge_container& edges, multicut_triangle_factor& T){
        double lb = 0.0;
        for (auto& e: edges)
            lb += std::min(0.0, e.second.load());
        for (auto&t:  T)
            lb += std::min({0.0, t.second[0]+t.second[1], t.second[0]+t.second[2], t.second[1]+t.second[1],
                            t.second[0]+t.second[1]+t.second[2]});
        return lb;
    }
    

    multicut_edge_labeling multicut_message_passing_parallel(const multicut_instance& input, const bool record_cycles, const int nr_threads)
    {
        const auto begin_time = std::chrono::steady_clock::now();
        tf::Executor executor(nr_threads);
        tf::Taskflow taskflow;
        std::size_t no_nodes = input.no_nodes();
        std::cout << "Message Passing\n";

        // Cycle Packing
        cycle_packing cp;
        multicut_triangle_factor T;
        edge_to_triangle_map M;
        atomic_edge_container all_edges;

        cp = cycle_packing_triangulation_parallel(input, nr_threads, T, M, all_edges);

        // Message Passing
        const auto MP_begin_time = std::chrono::steady_clock::now();
        
        auto [send_weights_start, send_weights_end] = taskflow.parallel_for(0, nr_threads, 1, [&](const std::size_t thread_no){
            const std::size_t batch_size = all_edges.size()/nr_threads + 1;
            int first_edge= thread_no*batch_size;
            int last_edge = std::min((thread_no+1)*batch_size, all_edges.size());
            atomic_edge_container::iterator iter_begin = std::next(all_edges.begin(), first_edge);
            atomic_edge_container::iterator iter_end = std::next(iter_begin, last_edge);
            for(auto it = iter_begin; it != iter_end; ++it){
                send_weights_to_triplets(std::ref(it->first[0]), std::ref(it->first[1]),
                    std::ref(it->second), std::ref(M), std::ref(T));
            }
        });

        auto [send_triplets_start, send_triplets_end] = taskflow.parallel_for(0, nr_threads, 1, [&](const std::size_t thread_no){
            const std::size_t batch_size = T.size()/nr_threads + 1;
            int first_triangle = thread_no*batch_size;
            int last_triangle = std::min((thread_no+1)*batch_size, T.size());
            multicut_triangle_factor::iterator iter_begin = std::next(T.begin(), first_triangle);
            multicut_triangle_factor::iterator iter_end = std::next(iter_begin, last_triangle);
            for(auto it = iter_begin; it != iter_end; ++it){
                send_triplets_to_edge(std::ref(all_edges), std::ref(it->first[0]), std::ref(it->first[1]), std::ref(it->first[2]), 
                    std::ref(it->second));
            }
        });

        send_triplets_start.gather(send_weights_end);
        
        executor.run(taskflow);
        executor.wait_for_all();

        const auto MP_end_time = std::chrono::steady_clock::now();
        std::cout << "lower_bound after MP: " << compute_lower_bound(all_edges, T) << std::endl;
        std::cout << "MP took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(MP_end_time - MP_begin_time).count() << " milliseconds\n";
  
        return gaec_parallel_non_blocking(input, nr_threads, all_edges);
    }

}
