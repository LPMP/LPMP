#include <vector>
#include <queue>
#include <cassert>
#include <functional>
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "union_find.hxx"
#include "dynamic_graph_thread_safe.hxx"
#include <thread>
#include <mutex>
#include <chrono>
#define SHUFFLE 0

namespace LPMP {

    struct edge_type {
        double cost;
        std::size_t stamp;
    };

    struct edge_type_q : public std::array<std::size_t,2> {
        double cost;
        std::size_t stamp;
        std::size_t no_neighbors;
    };

    // std::mutex partition_mutex;
    auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost/e1.no_neighbors < e2.cost/e2.no_neighbors;  };
    using pq_t = std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)>;

    tf::Task extract_positive_edges(tf::Taskflow& taskflow, std::vector<edge_type_q>& positive_edges, const multicut_instance& instance, const bool& sort, const dynamic_graph_thread_safe<edge_type>& g, tf::Task prev) {
        auto extract = [&]() {
            for (auto& edge : instance.edges()){
                if (edge.cost > 0.0)
                    positive_edges.push_back(edge_type_q{edge[0], edge[1], edge.cost, 0, g.edges(edge[0]).size() + g.edges(edge[1]).size()});
            }
        };
        auto Q = taskflow.emplace(extract);
        Q.gather(prev);
        if (sort) {
            auto sort_edge = [&](){std::sort(positive_edges.begin(), positive_edges.end(), pq_cmp);};
            auto T = taskflow.emplace(sort_edge);
            T.gather(Q);
            return T;
        } else {
            return Q;
        }
    }

    void distribute_edges_by_endnodes(tf::Taskflow& taskflow, const dynamic_graph_thread_safe<edge_type>& g, const multicut_instance& instance, std::pair<tf::Task, tf::Task>& Q, tf::Task prev, std::vector<pq_t>& queues, const int& nr_threads) {
        Q = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            const std::size_t nodes_batch_size = instance.no_nodes()/nr_threads + 1;
            const std::size_t first_node = thread_no*nodes_batch_size;
            const std::size_t last_node = std::min((thread_no+1)*nodes_batch_size, instance.no_nodes());
            for(std::size_t n=first_node; n<last_node; ++n){
                for (auto& neighbor_edge : g.edges(n)){
                    const std::size_t neighbor_node = neighbor_edge.first;
                    if (neighbor_node > n) {
                        auto edge_cost = g.edge(n, neighbor_node).cost;
                        if (edge_cost > 0.0)
                            queues[thread_no].push(edge_type_q{n, neighbor_node, edge_cost, 0, g.edges(n).size() + g.edges(neighbor_node).size()});
                    }
                }
            }
        });
        Q.first.gather(prev);
    }

    tf::Task distribute_edges_no_conflicts(tf::Taskflow& taskflow, const dynamic_graph_thread_safe<edge_type>& g, const multicut_instance& instance, tf::Task prev, std::vector<pq_t>& queues, const int& nr_threads, dynamic_graph_thread_safe<edge_type>& partial_graph, std::vector<std::tuple<std::size_t, std::size_t, double>>& remaining_edges) {
        auto extract = [&](){
            const std::size_t nodes_batch_size = instance.no_nodes()/nr_threads + 1;
            for(std::size_t n=0; n < instance.no_nodes(); ++n){
                bool within = true;
                for (auto& neighbor_edge : g.edges(n)){
                    const std::size_t neighbor_node = neighbor_edge.first;
                    if (neighbor_node > n && !((int)(n/nodes_batch_size) == (int)(neighbor_node/nodes_batch_size))){
                        within = false;
                        break;
                    }
                }
                for (auto& neighbor_edge : g.edges(n)){
                    const std::size_t neighbor_node = neighbor_edge.first;
                    if (neighbor_node > n) {
                        auto edge_cost = g.edge(n, neighbor_node).cost;
                        if (within) {
                            partial_graph.insert_edge(n, neighbor_node, {edge_cost,0});
                            if (edge_cost > 0.0) queues[n/nodes_batch_size].push(edge_type_q{n, neighbor_node, edge_cost, 0, g.edges(n).size() + g.edges(neighbor_node).size()});
                        } else {
                            remaining_edges.push_back(std::make_tuple(n, neighbor_node, edge_cost));
                            if (edge_cost > 0.0) queues[nr_threads].push(edge_type_q{n, neighbor_node, edge_cost, 0, g.edges(n).size() + g.edges(neighbor_node).size()});
                        }
                    }
                }
            }
            std::cout << "Edges in final queue: " << queues[nr_threads].size() << std::endl;
        };
        auto Q = taskflow.emplace(extract);
        Q.gather(prev);
        return Q;
    }


    void distribute_edges_round_robin(tf::Taskflow& taskflow, std::vector<edge_type_q>& positive_edges, std::pair<tf::Task, tf::Task>& Q, tf::Task prev, std::vector<pq_t>& queues, const int& nr_threads) {
        if (SHUFFLE) std::random_shuffle(positive_edges.begin(), positive_edges.end());
        Q = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::size_t e = thread_no;
            while(e < positive_edges.size()) {
                queues[thread_no].push(positive_edges[e]);
                e += nr_threads;
            }
        });
        Q.first.gather(prev);
    }

    void distribute_edges_in_chunks(tf::Taskflow& taskflow, std::vector<edge_type_q>& positive_edges, std::pair<tf::Task, tf::Task>& Q, tf::Task prev, std::vector<pq_t>& queues, const int& nr_threads) {
        if (SHUFFLE) std::random_shuffle(positive_edges.begin(), positive_edges.end());
        Q = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::size_t positive_edges_batch_size = positive_edges.size()/nr_threads + 1;
            const std::size_t first_edge = thread_no*positive_edges_batch_size;
            const std::size_t last_edge = std::min((thread_no+1)*positive_edges_batch_size, positive_edges.size());
            for(std::size_t e=first_edge; e<last_edge; ++e)
                queues[thread_no].push(positive_edges[e]);
        });
        Q.first.gather(prev);
    }

    template<bool SYNCHRONIZE=true>
    void gaec_single(dynamic_graph_thread_safe<edge_type>& g, pq_t& Q, union_find& partition, std::vector<std::atomic_flag>& mask)
    {
        const auto begin_time = std::chrono::steady_clock::now();
	    std::vector<edge_type_q> conflicted;
        std::vector<std::pair<std::array<std::size_t,2>, edge_type>> insert_candidates;
        std::vector<std::size_t> marked_nodes;
        std::vector<std::size_t> neighbor;

main_loop:
        while(!Q.empty() || !conflicted.empty()) {
            if (Q.empty()) {
                if constexpr (SYNCHRONIZE) {
                    std::cout << std::this_thread::get_id() << ": number of edges in conflicted set: " << conflicted.size() << std::endl;
                    for (auto & c: conflicted) Q.push(c);
                    conflicted.clear();
                    return;
                }
            }

            const edge_type_q e_q = Q.top();
            Q.pop();
            const std::size_t i = e_q[0];
            const std::size_t j = e_q[1];

            if constexpr(SYNCHRONIZE) {
                if (mask[i].test_and_set(std::memory_order_relaxed)){
                    conflicted.push_back(e_q);
                    continue;
                }
                if (mask[i].test_and_set(std::memory_order_relaxed)){
                    mask[i].clear(std::memory_order_relaxed);
                    conflicted.push_back(e_q);
                    continue;
                }
            }

            if(!g.edge_present(i,j)){
                if constexpr(SYNCHRONIZE) {
                    mask[i].clear(std::memory_order_relaxed);
                    mask[j].clear(std::memory_order_relaxed);
                }
                continue;
            }

            const auto& e = g.edge(i,j);
            if(e_q.stamp < e.stamp){
                if constexpr(SYNCHRONIZE) {
                    mask[i].clear(std::memory_order_relaxed);
                    mask[j].clear(std::memory_order_relaxed);
                }
                continue;
            }
            if(e.cost <= 0.0){
                if constexpr(SYNCHRONIZE) {
                    mask[i].clear(std::memory_order_relaxed);
                    mask[j].clear(std::memory_order_relaxed);
                }
                break;
            }

            // Mark relevant nodes
            if constexpr(SYNCHRONIZE) {
                neighbor.clear();
                for (auto& e : g.edges(i))
                    neighbor.push_back(e.first);
                for (auto& e : g.edges(j))
                    neighbor.push_back(e.first);
                std::sort(neighbor.begin(), neighbor.end());
                neighbor.erase(std::unique(neighbor.begin(), neighbor.end()), neighbor.end());

                marked_nodes.clear();
                marked_nodes.push_back(i);
                marked_nodes.push_back(j);
                for (auto& n: neighbor){
                    if (n == i || n == j)
                        continue;
                    if (mask[n].test_and_set(std::memory_order_relaxed)){
                        conflicted.push_back(e_q);
                        for (auto& m: marked_nodes) mask[m].clear(std::memory_order_relaxed);
                        goto main_loop;
                    } else {
                        marked_nodes.push_back(n);
                    }
                }
                assert(marked_nodes.size() == neighbor.size());
            }

            partition.merge(i,j);

            const auto node_pair  = [&]() -> std::array<std::size_t,2> {
                if(g.no_edges(i) < g.no_edges(j))
                    return {j,i};
                else
                    return {i,j};
            }();
            auto stable_node = node_pair[0];
            auto merge_node = node_pair[1];

            for (auto& e: g.edges(merge_node)){
                std::size_t head = e.first;
                if(head == stable_node)
                    continue;
                auto& p = g.edge(merge_node,head);

                if(g.edge_present(stable_node, head)) {
                    auto& pp = g.edge(stable_node, head);
                    pp.cost += p.cost;
                    pp.stamp++;
                    auto& pp_reverse = g.edge(head, stable_node);
                    pp_reverse.cost += p.cost;
                    pp_reverse.stamp++;
                    if(pp.cost >= 0.0)
                        Q.push(edge_type_q{stable_node, head, pp.cost, pp.stamp, g.edges(stable_node).size() + g.edges(head).size()});
                } else {
                    if(p.cost >= 0.0)
                        Q.push(edge_type_q{stable_node, head, p.cost, 0, g.edges(stable_node).size() + g.edges(head).size()});
                    insert_candidates.push_back({{stable_node, head}, {p.cost, 0}});
                }
            }

            g.remove_node(merge_node);
            for(const auto& e : insert_candidates)
                g.insert_edge(e.first[0], e.first[1], e.second);
            insert_candidates.clear();

            if constexpr(SYNCHRONIZE)
                for (auto& n: neighbor)
                    mask[n].clear(std::memory_order_relaxed);;
        }
        const auto end_time = std::chrono::steady_clock::now();
        //std::cout << "Parallel gaec time for thread " << std::this_thread::get_id() << " : "<<
        //std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }


    multicut_edge_labeling greedy_additive_edge_contraction_parallel(const multicut_instance& instance, const int nr_threads, const std::string option)
    {
        const auto begin_time = std::chrono::steady_clock::now();

        union_find partition(instance.no_nodes());

        tf::Executor executor;
        tf::Taskflow taskflow;

        std::vector<std::atomic_flag>  mask(instance.no_nodes());
        auto fill_mask = [&]() {
            for (auto& m: mask)
                m.clear(std::memory_order_relaxed);
        };
        auto M = taskflow.emplace(fill_mask);

        dynamic_graph_thread_safe<edge_type> g(instance.no_nodes()), partial_graph(instance.no_nodes());
        auto construct_graph = [&]() { g.construct(instance.edges().begin(), instance.edges().end(), [](const auto& e) -> edge_type { return {e.cost, 0}; }); };
        auto G = taskflow.emplace(construct_graph);

	    std::vector<pq_t> queues(nr_threads+1, pq_t(pq_cmp));
        std::vector<edge_type_q> positive_edges;
        std::vector<std::tuple<std::size_t, std::size_t, double>> remaining_edges;
        //std::unordered_map<std::size_t, std::pair<std::size_t, double>> remaining_edges;
        std::pair<tf::Task, tf::Task> Q;
        tf::Task E, prev;

        if(option == "round_robin_not_sorted"){
            E = extract_positive_edges(taskflow, positive_edges, instance, false, g, G);
            distribute_edges_round_robin(taskflow, positive_edges, Q, E, queues, nr_threads);
            prev = Q.second;
        } else if (option == "endnodes"){
            distribute_edges_by_endnodes(taskflow, g, instance, Q, G, queues, nr_threads);
            prev = Q.second;
        } else if (option == "chunk_not_sorted"){
            E = extract_positive_edges(taskflow, positive_edges, instance, false, g, G);
            distribute_edges_in_chunks(taskflow, positive_edges, Q, E, queues, nr_threads);
            prev = Q.second;
        } else if (option == "round_robin_sorted"){
            E = extract_positive_edges(taskflow, positive_edges, instance, true, g, G);
            distribute_edges_round_robin(taskflow, positive_edges, Q, E, queues, nr_threads);
            prev = Q.second;
        } else if (option == "chunk_sorted"){
            E = extract_positive_edges(taskflow, positive_edges, instance, true, g, G);
            distribute_edges_in_chunks(taskflow, positive_edges, Q, E, queues, nr_threads);
            prev = Q.second;
        } else { // non-blocking
            E = distribute_edges_no_conflicts(taskflow, g, instance, G, queues, nr_threads, partial_graph, remaining_edges);
            prev = E;
        }

        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Initialization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        auto [GAEC_begin,GAEC_end] = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::cout << "Edges in queue " << thread_no << ": " << queues[thread_no].size() << std::endl;
            if (option == "non-blocking")
                gaec_single<false>(std::ref(partial_graph), std::ref(queues[thread_no]), std::ref(partition), std::ref(mask));
            else
                gaec_single<true>(std::ref(g), std::ref(queues[thread_no]), std::ref(partition), std::ref(mask));
        });
        GAEC_begin.gather(prev);
        executor.run(taskflow);
        executor.wait_for_all();

        if (option == "non-blocking" && !queues[nr_threads].empty()) {
            std::cout << "Edges in the final thread " << nr_threads << ": " << queues[nr_threads].size() << std::endl;
            // for (auto& e : remaining_edges) {
            //     partial_graph.insert_edge(std::get<0>(e),std::get<1>(e),{std::get<2>(e),0});
            // }
            //gaec_single<false>(std::ref(g), std::ref(queues[nr_threads]), std::ref(partition), std::ref(mask));
        }
	    return multicut_edge_labeling(instance, partition);
   }

} // namespace LPMP
