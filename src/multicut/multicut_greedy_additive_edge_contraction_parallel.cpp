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

    std::mutex partition_mutex;
    auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost/e1.no_neighbors < e2.cost/e2.no_neighbors;  };
    using pq_t = std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)>;

    int bool_clear = false, bool_mark = true;
    void unmark(std::atomic<int>& node) {
        while (!node.compare_exchange_weak(bool_mark, 0, std::memory_order_relaxed));
    }
    int marked(std::atomic<int>& node){
        return !(node.compare_exchange_strong(bool_clear, 1, std::memory_order_relaxed));
    }

    void gaec_single(dynamic_graph_thread_safe<edge_type>& g, pq_t& Q, union_find& partition, std::vector<std::atomic<int>>& mask)
    {
	    std::vector<edge_type_q> conflicted; 
        std::vector<std::pair<std::array<std::size_t,2>, edge_type>> insert_candidates; 
        std::vector<std::size_t> marked_nodes;
        std::vector<std::size_t> neighbor;

main_loop:
        while(!Q.empty() || !conflicted.empty()) {
            if (Q.empty()) {
                std::cout << std::this_thread::get_id() << ": number of edges in conflicted set: " << conflicted.size() << std::endl;
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
                std::lock_guard<std::mutex> lck(partition_mutex);
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
                            Q.push(edge_type_q{stable_node, head, pp.cost, pp.stamp, neighbor.size()});

                    } else {
                        if(p.cost >= 0.0)
                            Q.push(edge_type_q{stable_node, head, p.cost, 0, neighbor.size()});
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
    }


    multicut_edge_labeling greedy_additive_edge_contraction_parallel(const multicut_instance& instance, const int nr_threads)
    {
        dynamic_graph_thread_safe<edge_type> g(instance.edges().begin(), instance.edges().end(), [](const auto& e) -> edge_type { return {e.cost, 0}; });

        union_find partition(instance.no_nodes());
    
    	std::vector<std::atomic<int>>  mask(instance.no_nodes());
        std::fill(mask.begin(), mask.end(), false);

	    std::vector<pq_t> queues(nr_threads, pq_t(pq_cmp));
        std::vector<std::thread> T(nr_threads);
	
	    // distribute edges among priority queues
	    int edge_no = 0;
/*
        std::vector<edge_type_q> total_edges(instance.no_nodes());
        for(const auto& e : instance.edges()){
            if(e.cost > 0.0)
                total_edges.push_back(edge_type_q{e[0], e[1], e.cost, 0, g.edges(e[0]).size() + g.edges(e[1]).size()});
    	}
        std::random_shuffle(total_edges.begin(), total_edges.end()); // faster for knott but not for gm_large
        int batch_size = total_edges.size()/nr_threads + 1;
        for (const auto&e: total_edges){    
            queues[edge_no/batch_size].push(e);
            edge_no++;
        }
*/
        int batch_size = instance.edges().size()/nr_threads + 1;
        for(const auto& e : instance.edges()){
            if(e.cost > 0.0)
	            queues[edge_no/batch_size].push(edge_type_q{e[0], e[1], e.cost, 0, g.edges(e[0]).size() + g.edges(e[1]).size()});
            edge_no++;
    	}

	    for (int i = 0; i < nr_threads; ++i){
            std::cout << "Edges in thread " << i << ": " << queues[i].size() << std::endl;
	        T[i] = std::thread(gaec_single, std::ref(g), std::ref(queues[i]), std::ref(partition), std::ref(mask));
	    }   

	    for (auto &t : T) {
	        t.join();
	    }   
	    return multicut_edge_labeling(instance, partition);

   }

} // namespace LPMP 
