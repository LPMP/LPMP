#include <vector>
#include <queue>
#include <cassert>
#include <functional>
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "union_find.hxx"
#include "dynamic_graph_thread_safe.hxx"
#include <boost/dynamic_bitset.hpp>
#include <thread>
#include <mutex>
#include <chrono>

namespace LPMP {

    struct edge_type {
        double cost;
        std::size_t stamp;
    };

    struct edge_type_q : public std::array<std::size_t,2> {
        double cost;
        std::size_t stamp;
        //std::size_t no_neighbors;
    };

    //std::mutex partition_mutex;
    //auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost/e1.no_neighbors < e2.cost/e2.no_neighbors;  };
    auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost < e2.cost;  };
    using pq_t = std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)>;

    //char bool_clear = false, bool_mark = true;
    void unmark(std::atomic<char>& node) {
        assert(node == 1);
        node.store(0); 
        //char bool_mark = true;
        //while (!node.compare_exchange_weak(bool_mark, 0, std::memory_order_relaxed));
    }
    bool marked(std::atomic<char>& node){
        char bool_clear = false;
        return !(node.compare_exchange_strong(bool_clear, 1, std::memory_order_relaxed));
    }

    void gaec_single(dynamic_graph_thread_safe<edge_type>& g, pq_t& Q, union_find& partition, std::vector<std::atomic<char>>& mask)
    {
        const auto begin_time = std::chrono::steady_clock::now();
	    std::vector<edge_type_q> conflicted; 
        std::vector<std::pair<std::array<std::size_t,2>, edge_type>> insert_candidates; 
        std::vector<std::size_t> marked_nodes;
        std::vector<std::size_t> neighbor;

main_loop:
        while(!Q.empty() || !conflicted.empty()) {
            if (Q.empty()) {
                //std::cout << std::this_thread::get_id() << ": number of edges in conflicted set: " << conflicted.size() << std::endl;
                for (auto & c: conflicted) Q.push(c);
                conflicted.clear();
                return;
            }

            const edge_type_q e_q = Q.top();
            Q.pop();
            const std::size_t i = e_q[0];
            const std::size_t j = e_q[1];
//            std::cout << "\nProcessing edge " << i << " " << j << std::endl;

            if (marked(mask[i])){
                conflicted.push_back(e_q);
                continue;
            }
            if (marked(mask[j])){
                unmark(mask[i]);
                conflicted.push_back(e_q);
                continue;
            }
            if(!g.edge_present(i,j)){
                unmark(mask[i]);
                unmark(mask[j]);
                continue;
            }
            const auto& e = g.edge(i,j);
            if(e_q.stamp < e.stamp){
                unmark(mask[i]);
                unmark(mask[j]);
                continue;
            }
            if(e.cost <= 0.0){
                unmark(mask[i]);
                unmark(mask[j]);
                break;
            }

            // Mark relevant nodes
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
                if (marked(mask[n])){
                    conflicted.push_back(e_q);
                    for (auto& m: marked_nodes) unmark(mask[m]);
                    goto main_loop;
                } else {
                    marked_nodes.push_back(n);
                }
            }

            assert(marked_nodes.size() == neighbor.size());
            
            {
                //std::lock_guard<std::mutex> lck(partition_mutex);
                partition.merge(i,j);
            }

            std::size_t stable_node, merge_node;
            {
                const auto node_pair  = [&]() -> std::array<std::size_t,2> {
                    if(g.no_edges(i) < g.no_edges(j))
                        return {j,i};
                    else
                        return {i,j};
                }();
                stable_node = node_pair[0];
                merge_node = node_pair[1];
                    
                for (auto& e: g.edges(merge_node)){
                    std::size_t head = e.first;
                    if(head == stable_node)
                        continue;
                    auto& p = g.edge(merge_node,head);
                    
                    if(g.edge_present(stable_node, head)) {
                        auto& pp = g.edge(stable_node, head);
                        auto& pp_reverse = g.edge(head, stable_node);
                        pp.cost += p.cost;
                        pp.stamp++;
                        pp_reverse.cost += p.cost;
                        pp_reverse.stamp++; 
                        if(pp.cost >= 0.0)
                            Q.push(edge_type_q{stable_node, head, pp.cost/neighbor.size(), pp.stamp});

                    } else {
                        if(p.cost >= 0.0)
                            Q.push(edge_type_q{stable_node, head, p.cost/neighbor.size(), 0});
//                        std::cout << "Insert node " << stable_node << " " << head << std::endl;
                        insert_candidates.push_back({{stable_node, head}, {p.cost, 0}});
                    } 
                }
            }
            {
                g.remove_node(merge_node);
                for(const auto& e : insert_candidates)
                    g.insert_edge(e.first[0], e.first[1], e.second);
                insert_candidates.clear();
            }
            
            for (auto& n: neighbor) unmark(mask[n]);
        }
        const auto end_time = std::chrono::steady_clock::now();
        //std::cout << "Parallel gaec time for thread " << std::this_thread::get_id() << " : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }


    multicut_edge_labeling greedy_additive_edge_contraction_parallel(const multicut_instance& instance, const int nr_threads)
    {
        const auto begin_time = std::chrono::steady_clock::now();

        tf::Executor executor;
        tf::Taskflow taskflow;

        union_find partition;
        auto create_partition = [&]() { partition = union_find(instance.no_nodes()); };
        auto P = taskflow.emplace(create_partition);

        std::vector<std::atomic<char>>  mask(instance.no_nodes());
        auto fill_mask = [&]() { std::fill(mask.begin(), mask.end(), false); };
        auto M = taskflow.emplace(fill_mask); 

        dynamic_graph_thread_safe<edge_type> g;
        auto construct_graph = [&]() { g.construct(instance.edges().begin(), instance.edges().end(), [](const auto& e) -> edge_type { return {e.cost, 0}; }); };
        auto G = taskflow.emplace(construct_graph); 

	    std::vector<pq_t> queues(nr_threads, pq_t(pq_cmp));
        std::size_t batch_size = instance.edges().size()/nr_threads + 1;
        auto [Q_begin,Q_end] = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
                const std::size_t first_edge = thread_no*batch_size;
                const std::size_t last_edge = std::min((thread_no+1)*batch_size, instance.no_edges());
                for(std::size_t e=first_edge; e<last_edge; ++e) {
                auto& edge = instance.edges()[e];
                if(edge.cost > 0.0)
                queues[thread_no].push(edge_type_q{edge[0], edge[1], edge.cost/(g.edges(edge[0]).size() + g.edges(edge[1]).size()), 0 });
                }
                });

        Q_begin.gather(G);


/*
        std::vector<edge_type_q> total_edges(instance.no_nodes());
        for(const auto& e : instance.edges()){
            if(e.cost > 0.0)
                total_edges.push_back(edge_type_q{e[0], e[1], e.cost, 0, g.edges(e[0]).size() + g.edges(e[1]).size()});
    	}
        std::random_shuffle(total_edges.begin(), total_edges.end()); // faster for knott but not for gm_large
        std::size_t batch_size = total_edges.size()/nr_threads + 1;
        for (const auto&e: total_edges){    
            queues[edge_no/batch_size].push(e);
            edge_no++;
        }
*/
        //const auto end_time = std::chrono::steady_clock::now();
        //std::cout << "Initialization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        auto [GAEC_begin,GAEC_end] = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
                //std::cout << "Edges in thread " << thread_no << ": " << queues[thread_no].size() << std::endl;
                gaec_single(g, queues[thread_no], partition, mask);
                });

        GAEC_begin.gather(Q_end);
        GAEC_begin.gather(P);
        GAEC_begin.gather(M);

        executor.run(taskflow);
        executor.wait_for_all();

	    return multicut_edge_labeling(instance, partition);
   }

} // namespace LPMP 
