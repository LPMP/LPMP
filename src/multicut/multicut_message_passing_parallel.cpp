#include "cut_base/cut_base_instance.hxx"
#include "cut_base/cut_base_apply_packing.hxx"
#include "multicut/multicut_factors_messages.h"
#include "multicut/multicut_message_passing_parallel.h"
#include "graph.hxx"
#include "union_find.hxx"
#include <array>
#include <chrono>
#include <atomic>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <thread>
#include <mutex>
#include <taskflow/taskflow.hpp>
#include "dynamic_graph_thread_safe.hxx"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    struct weighted_edge : public std::array<std::size_t,2> { double cost; };

    struct edge_type_cost_stamp {
        double cost;
        std::size_t stamp;
    };

    struct edge_type_q : public std::array<std::size_t,2> {
        double cost;
        std::size_t stamp;
        std::size_t no_neighbors;
    };

    auto pq_cmp = [](const edge_type_q& e1, const edge_type_q& e2) { return e1.cost/e1.no_neighbors < e2.cost/e2.no_neighbors;  };
    using pq_t = std::priority_queue<edge_type_q, std::vector<edge_type_q>, decltype(pq_cmp)>;

    
    using atomic_edge_container = std::map<std::array<std::size_t,2>, double>;

    template<class T>
    class CopyableAtomic : public std::atomic<T>
    {
    public:
        CopyableAtomic() :
            std::atomic<T>(T{})
        {}

        constexpr CopyableAtomic(T desired): std::atomic<T>(desired){}

        constexpr CopyableAtomic(const CopyableAtomic<T>& other): CopyableAtomic(other.load(std::memory_order_relaxed))
        {}

        CopyableAtomic& operator=(const CopyableAtomic<T>& other) {
            this->store(other, std::memory_order_relaxed);
            return *this;
        }

        CopyableAtomic& operator=(const T& other) {
            this->store(other, std::memory_order_relaxed);
            return *this;
        }
    };

    void compute_connectivity(union_find& uf, graph<CopyableAtomic<double>>& pos_edges_graph) {
        uf.reset();
        pos_edges_graph.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost) {
            if(cost >= tolerance) { uf.merge(i,j); }
        });
    };

    bool connected(union_find& uf, const std::size_t i, const std::size_t j) { return uf.connected(i,j); };

    tf::Task distribute_edges_round_robin(tf::Taskflow& taskflow, const multicut_instance& input, 
        std::vector<std::vector<weighted_edge>>& positive_edge_vec, 
        std::vector<std::vector<weighted_edge>>& repulsive_edge_vec, 
        const int& nr_threads, std::vector<double>& lower) {
        auto Q = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::size_t e = thread_no;
            while(e < input.edges().size()) {
                auto edge = input.edges()[e];
                auto smaller_node = std::min(edge[0], edge[1]);
                auto larger_node =  std::max(edge[0], edge[1]); 
                if(edge.cost < 0.0)
                    repulsive_edge_vec[thread_no].push_back({smaller_node, larger_node, edge.cost});
                else if(edge.cost > 0.0)
                    positive_edge_vec[thread_no].push_back({smaller_node, larger_node, edge.cost});
                lower[thread_no] += std::min(0.0, edge.cost[0]);
                e += nr_threads;
            }
        });
        return Q.second;
    }

   tf::Task distribute_edges_no_conflicts(tf::Taskflow& taskflow, atomic_edge_container& instance, 
    std::vector<pq_t>& queues, const int& nr_threads, const std::size_t& no_nodes,
    dynamic_graph_thread_safe<edge_type_cost_stamp>& partial_graph, 
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, double>>>& remaining_edges) { 
        std::vector<size_t> degree(no_nodes);
        for (auto& edge: instance){
            degree[edge.first[0]]++;
            degree[edge.first[1]]++;
        }
        partial_graph.pre_allocate(degree.begin(), degree.end());
        auto [extract_begin, extract_end] = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            const std::size_t nodes_batch_size = instance.size()/nr_threads + 1;
            std::cout<< "nodes_batch_size "<< nodes_batch_size << std::endl;
            for (auto& edge: instance){
                int nth = (int)(edge.first[0]/nodes_batch_size);
                if (nth == thread_no){
                    if (nth == (int)(edge.first[1]/nodes_batch_size)) {
                        partial_graph.insert_edge(edge.first[0], edge.first[1], {edge.second,0});
                        if (edge.second > 0.0)
                            queues[thread_no].emplace(edge_type_q{edge.first[0], edge.first[1], edge.second, 0, 1});
                    } else {
                        remaining_edges[thread_no].push_back(std::make_tuple(edge.first[0], edge.first[1], edge.second));
                    }
                }
            }
        });
        return extract_end;
    }

    using multicut_triangle_factor = std::map<std::array<std::size_t,3>, std::vector<double>>;
    using edge_to_triangle_map = std::unordered_map<std::array<std::size_t,2>, std::set<std::size_t>>;

    void add_triangles(std::vector<std::size_t>& cycle, double cycle_cap, edge_to_triangle_map& M, multicut_triangle_factor& T){
        for(std::size_t c=1; c<cycle.size()-1; ++c) {
            std::vector<std::pair<size_t, int>> nodes = {{cycle[0],0}, {cycle[c],1}, {cycle[c+1],2}};
            std::sort(nodes.begin(), nodes.end(), [](std::pair<size_t, double> n1, std::pair<size_t, double> n2) {return n1.first < n2.first;});
            std::size_t i=nodes[0].first, j=nodes[1].first, k=nodes[2].first;
            double w_ij, w_jk, w_ik;
            if (nodes[0].second == 1){
                w_ij = cycle_cap; w_jk = -cycle_cap; w_ik = cycle_cap;
            } else if (nodes[1].second == 1){
                w_ij = cycle_cap; w_jk = cycle_cap; w_ik = -cycle_cap;
            } else {
                w_ij = -cycle_cap; w_jk = cycle_cap; w_ik = cycle_cap;
            }
            if (T.find({i,j,k}) != T.end()){
     //           std::cout << "Triangle already exists: " << std::endl;
                T[{i,j,k}][0] += w_ij; T[{i,j,k}][1] += w_jk; T[{i,j,k}][2] += w_ik;
            } else {
                std::vector<double> weights = {w_ij, w_jk, w_ik};
                T[{i,j,k}] = weights;
            }
            M[{i,j}].insert(k); 
            M[{i,k}].insert(j); 
            M[{j,k}].insert(i); 
        }
    }

    void cp_single(std::vector<weighted_edge>& repulsive_edges,
        graph<CopyableAtomic<double>>& pos_edges_graph, 
        const std::size_t& cycle_length, const std::size_t& no_nodes, const bool& record_cycles, std::atomic<double>& lower_bound, 
        cycle_packing& cp, edge_to_triangle_map& M, multicut_triangle_factor& T){

        bfs_data<graph<CopyableAtomic<double>>> bfs(pos_edges_graph);
        std::size_t progress = 0;
        union_find uf(no_nodes);

        for(std::size_t i=0; i<repulsive_edges.size(); ++i) {
            auto& re = repulsive_edges[i];

            // periodically update component labeling for speed-up
            progress++;
            if (progress > 0.05 * repulsive_edges.size()) {
                compute_connectivity(uf, pos_edges_graph);
                progress = 0;
            }

            // check if conflicted cycle exists and repulsive edge has positive weight
            if(-re.cost<= tolerance || !connected(uf, re[0], re[1]) || std::max(re[0], re[1]) >= pos_edges_graph.no_nodes()) {
                const std::size_t imax = repulsive_edges.size() - 1;
                repulsive_edges[i] = repulsive_edges[imax];
                repulsive_edges.resize(imax);
                i--;
                continue;
            }

            // balance short cycles as long as available and positive weight left
            auto mask_small_edges = [cycle_length](const std::size_t i, const std::size_t j, const double cost, const std::size_t distance) {
                if(cost <= tolerance) return false;
                if(distance >= cycle_length) return false;
                return true;
            };
            std::vector<std::size_t> cycle;

            while(-re.cost > tolerance) {
                double cycle_cap = std::numeric_limits<double>::infinity();
                auto cycle_capacity = [&cycle_cap](const std::size_t i, const std::size_t j, const double cost) {
                        cycle_cap = std::min(cycle_cap, cost); };
                cycle = bfs.find_path(re[0], re[1], mask_small_edges, cycle_capacity); // possibly do not reallocate for every search but give as argument to find_pathatomic_weighted_edge
                cycle_cap = std::min(cycle_cap, -re.cost);
                if(cycle.size() > 0 && cycle_cap >= tolerance) {
                    if(record_cycles)
                        cp.add_cycle(cycle.begin(), cycle.end(), cycle_cap);
                    
                    // subtract minimum weight and remove edges of negligible weight
                    assert(-re.cost>= cycle_cap && re.cost<= 0.0);
                    re.cost += cycle_cap;
                    assert(re.cost <= 0.0);

                    if (cycle_length>=3) add_triangles(cycle, cycle_cap, M, T);
                    for(std::size_t c=1; c<cycle.size(); ++c) {
                        pos_edges_graph.edge(cycle[c-1], cycle[c]).store(pos_edges_graph.edge(cycle[c-1], cycle[c]).load() - cycle_cap);
                        pos_edges_graph.edge(cycle[c], cycle[c-1]).store(pos_edges_graph.edge(cycle[c], cycle[c-1]).load() - cycle_cap);
                        assert(pos_edges_graph.edge(cycle[c-1], cycle[c]).load() == pos_edges_graph.edge(cycle[c], cycle[c-1]).load());
                    }

                    lower_bound.store(lower_bound.load() + cycle_cap);

                    // remove repulsive edge if weight is negligible
                    if (-re.cost<= tolerance) {
                    //    std::cout << "weight is negligible:";atomic_weighted_edge
                        auto imax = repulsive_edges.size() - 1;
                        repulsive_edges[i] = repulsive_edges[imax];
                        repulsive_edges.resize(imax);
                        i--;
                        break;
                    }
                } else {
                    break;
                }
            }
        }
    }

    void send_weights_to_triplets(std::size_t node1, std::size_t node2, double cost, edge_to_triangle_map& M, multicut_triangle_factor& T){
        assert(node1 < node2);
        if (M.find({node1,node2}) == M.end()){
   //         std::cout << "send_weights_to_triplets: no triangle for edge " << e.first[0] << " " << e.first[1] << std::endl;
        } else {
            auto nr_triangles = M[{node1,node2}].size();
            for (auto third: M[{node1,node2}]){
                std::vector<std::pair<size_t, int>> nodes = {{node1,0},{node2,1},{third,2}};
                std::sort(nodes.begin(), nodes.end(), [](std::pair<size_t, double> n1, std::pair<size_t, double> n2) {return n1.first < n2.first;});
                std::size_t i=nodes[0].first, j=nodes[1].first, k=nodes[2].first;
                if (T.find({i,j,k}) == T.end()) {
                    std::cout << "send_weights_to_triplets: triangle" << i << " "<< j << " " << k << " missing. \n";
                } else {
                    if (nodes[0].second == 2){
                        T[{i,j,k}][1] += cost/nr_triangles;
                    } else if (nodes[1].second == 2){
                        T[{i,j,k}][2] += cost/nr_triangles;
                    } else {
                        T[{i,j,k}][0] += cost/nr_triangles;
                    }
                }
            }
        }
        
    }

    void marginalize(double& edge_cost, std::vector<double>& cost, int option, double omega){
        auto marginal = std::min({cost[option]+cost[(option+1)%3], cost[option]+cost[(option+2)%3], cost[0]+cost[1]+cost[2]}) 
                      - std::min(0.0, cost[(option+1)%3]+cost[(option+2)%3]);
        edge_cost = edge_cost + omega*marginal;
        cost[option] -= omega*marginal;
    }

    void send_triplets_to_edge(atomic_edge_container& edges, std::size_t i, std::size_t j, std::size_t k, std::vector<double>& triangle_cost){
        // ij: 0 jk: 1 ik:2
        marginalize(edges[{i, j}], triangle_cost, 0, 1/3);
        marginalize(edges[{i, k}], triangle_cost, 2, 1/2);
        marginalize(edges[{j, k}], triangle_cost, 1, 1/1);
        marginalize(edges[{i, j}], triangle_cost, 0, 1/2);
        marginalize(edges[{i, k}], triangle_cost, 2, 1/1);
        marginalize(edges[{i, j}], triangle_cost, 0, 1/1);
    }

    static std::vector<std::size_t> empty_partition_node = {};

    void gaec_single(dynamic_graph_thread_safe<edge_type_cost_stamp>& g, pq_t& Q, union_find& partition, 
        std::vector<std::size_t>& partition_to_node=empty_partition_node) {
        std::vector<std::pair<std::array<std::size_t,2>, edge_type_cost_stamp>> insert_candidates;

        while(!Q.empty()) {
            const edge_type_q e_q = Q.top();
            Q.pop();
            const std::size_t i = e_q[0];
            const std::size_t j = e_q[1];

            if(!g.edge_present(i,j)) continue;
            const auto& e = g.edge(i,j);
            if(e_q.stamp < e.stamp) continue;
            if(e.cost <= 0.0) break;

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

            if (!partition_to_node.empty()) {
                const std::size_t t = partition.find(stable_node);
                assert(t == partition.find(merge_node));
                partition_to_node[t] = stable_node;
            }
            
            for(const auto& e : insert_candidates)
            g.insert_edge(e.first[0], e.first[1], e.second);
            insert_candidates.clear();

        }
    }
   

    multicut_edge_labeling multicut_message_passing_parallel(const multicut_instance& input, const bool record_cycles, const int nr_threads)
    {
        const auto begin_time = std::chrono::steady_clock::now();
        tf::Executor executor(nr_threads);
        tf::Taskflow taskflow;


        // cycle packing initialization
        std::vector<double> lower(nr_threads, 0.0);
        std::atomic<double> lower_bound = 0.0;
        std::vector<weighted_edge> positive_edges;
        std::vector<std::vector<weighted_edge>> repulsive_edge_vec(nr_threads);
        std::vector<std::vector<weighted_edge>> positive_edge_vec(nr_threads);
        const std::size_t no_nodes = input.no_nodes();
        auto init_edge = distribute_edges_round_robin(taskflow, input, positive_edge_vec, repulsive_edge_vec, nr_threads, lower);

        graph<CopyableAtomic<double>> pos_edges_graph;
        auto construct_pos_edges_graph = [&]() {
            for (auto l: lower) lower_bound = lower_bound + l;
            for (auto& pos: positive_edge_vec) std::move(pos.begin(), pos.end(), std::back_inserter(positive_edges));
            pos_edges_graph.construct(positive_edges.begin(), positive_edges.end(),
                [](const weighted_edge& e) { return e.cost; }); };
        auto CPE = taskflow.emplace(construct_pos_edges_graph);
        CPE.gather(init_edge);

        std::future<void> fu = executor.run(taskflow);
        fu.get();

        const auto initialization_end_time = std::chrono::steady_clock::now();
        std::cout << "cycle packing initialization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(initialization_end_time - begin_time).count() << " milliseconds\n";
        std::cout << "initial lower bound = " << lower_bound << "\n";
        std::cout << "#attractive edges = " << positive_edges.size() << "\n";
        std::cout << "pos_edges_graph.no_nodes:" << pos_edges_graph.no_nodes() << std::endl;
        std::cout << "pos_edges_graph.no_edges:" << pos_edges_graph.no_edges() << std::endl;

        // Cycle Packing
        cycle_packing cp;
        multicut_triangle_factor T;
        edge_to_triangle_map M;

        const std::array<std::size_t,11> cycle_lengths = {1,2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
        for(const std::size_t cycle_length : cycle_lengths) {

            taskflow.clear();
            std::cout << "find cycles of length " << cycle_length << " lower bound = " << lower_bound << "\n";

            //shuffling can give great speed-up if edges with similar indices are spatially close in the graph
            for (auto& re: repulsive_edge_vec) std::random_shuffle(re.begin(), re.end());

            auto CP = taskflow.parallel_for(0, nr_threads, 1, [&](const std::size_t thread_no){
            //    std::cout << "#repulsive edges :" << repulsive_edge_vec[thread_no].size() << " remaining in thread " << thread_no << std::endl;
                cp_single(std::ref(repulsive_edge_vec[thread_no]), std::ref(pos_edges_graph), std::ref(cycle_length),  
                    std::ref(no_nodes), std::ref(record_cycles), std::ref(lower_bound), std::ref(cp), std::ref(M), std::ref(T));
            });

            executor.run(taskflow);
            executor.wait_for_all();
        }
        const auto CP_end_time = std::chrono::steady_clock::now();
        std::cout << "CP took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(CP_end_time - begin_time).count() << " milliseconds\n";
        std::cout << "final lower bound = " << lower_bound << "\n";
        

        // Message Passing
        // collect edges
        taskflow.clear();
        atomic_edge_container all_edges;
        pos_edges_graph.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost){
            all_edges[{i, j}] = cost;
        });
        for (auto& re: repulsive_edge_vec) { 
            for (auto& e: re){
                all_edges[{e[0], e[1]}] = e.cost;
            }
        }
        std::cout << "#Edges = " << all_edges.size() << "\n";
        std::cout << "Number of triangles:" << M.size() << std::endl;
        
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
        std::cout << "MP took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(MP_end_time - CP_end_time).count() << " milliseconds\n";

        // GAEC
        taskflow.clear();
        union_find partition(no_nodes);
        dynamic_graph_thread_safe<edge_type_cost_stamp>  g(no_nodes);
        std::vector<std::size_t> partition_to_node(no_nodes);
        std::vector<pq_t> queues(nr_threads+1, pq_t(pq_cmp));
        std::vector<std::vector<std::tuple<std::size_t, std::size_t, double>>> remaining_edges(nr_threads+1);
        for(std::size_t i=0; i<no_nodes; ++i)
            partition_to_node[partition.find(i)] = i;

        auto E = distribute_edges_no_conflicts(taskflow, all_edges, queues, nr_threads, no_nodes, g, remaining_edges);

        auto [GAEC_begin,GAEC_end] = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::cout << "Edges in queue " << thread_no << ": " << queues[thread_no].size() << std::endl;
            gaec_single(std::ref(g), std::ref(queues[thread_no]), std::ref(partition), std::ref(partition_to_node));
            
        });
        GAEC_begin.gather(E);
        executor.run(taskflow);
        executor.wait_for_all();

        const auto remaining_edges_begin_time = std::chrono::steady_clock::now();
        for (auto& v: remaining_edges) {
            for (auto& e: v){
                const double cost = std::get<2>(e);
                const std::size_t pi = partition.find(std::get<0>(e));
                const std::size_t i = partition_to_node[pi];

                const std::size_t pj = partition.find(std::get<1>(e));
                const std::size_t j = partition_to_node[pj];

                if (g.edge_present(i,j)){
                    g.edge(i,j).cost += cost;
                } else
                    g.insert_edge(i,j, {cost,0});

                if (cost > 0.0) queues[nr_threads].push(edge_type_q{i,j,cost,0,1});
            }
        }
        std::cout << "GAEC: Edges in the final thread " << nr_threads << ": " << queues[nr_threads].size() << std::endl;
        gaec_single(std::ref(g), std::ref(queues[nr_threads]), std::ref(partition));

        const auto remaining_edges_end_time = std::chrono::steady_clock::now();
        std::cout << "handling remaining edges took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(remaining_edges_end_time - remaining_edges_begin_time).count() << " milliseconds\n";
    
        return multicut_edge_labeling(input, partition);
    }

}
