#include <array>
#include <chrono>
#include <atomic>
#include <algorithm>
#include <cassert>
#include <mutex>
#include <taskflow/taskflow.hpp>
#include "graph.hxx"
#include "union_find.hxx"
#include "cut_base/cut_base_instance.hxx"
#include "multicut/multicut_cycle_packing_parallel.h"

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    std::mutex T_mutex, M_mutex;

    using triangle_dict = std::unordered_map<std::array<std::size_t,3>, int>;
    using edge_dict = std::unordered_map<std::array<std::size_t,2>, int>;

    static std::vector<edge_item> empty_e = {};
    static std::vector<triangle_item> empty_t = {};
    static std::vector<edge_t> empty_o = {};
    static triangle_dict empty_t_dict = {};
    static edge_dict empty_e_dict = {};

    void compute_connectivity(union_find& uf, graph<CopyableAtomic<double>>& pos_edges_graph) {
        uf.reset();
        pos_edges_graph.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost) {
            if(cost >= tolerance) { uf.merge(i,j); }
        });
    };

    bool connected(union_find& uf, const std::size_t i, const std::size_t j) { return uf.connected(i,j); };
    
    tf::Task distribute_edges_round_robin(tf::Taskflow& taskflow, const multicut_instance& input, 
        std::vector<std::vector<edge_t>>& positive_edge_vec, std::vector<std::vector<edge_t>>& repulsive_edge_vec, 
        const int& nr_threads, std::vector<double>& lower) 
    {
        auto Q = taskflow.for_each_index(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::size_t e = thread_no;
            while(e < input.edges().size()) {
                auto edge = input.edges()[e];
                if(edge.cost < 0.0)
                    repulsive_edge_vec[thread_no].push_back({edge[0], edge[1], edge.cost});
                else if(edge.cost > 0.0)
                    positive_edge_vec[thread_no].push_back({edge[0], edge[1], edge.cost});
                lower[thread_no] += std::min(0.0, edge.cost[0]);
                e += nr_threads;
            }
        });
        return Q;
    }


    tf::Task distribute_edges_in_chunks(tf::Taskflow& taskflow, const multicut_instance& input, 
        std::vector<std::vector<edge_t>>& positive_edge_vec, std::vector<std::vector<edge_t>>& repulsive_edge_vec, 
        const int& nr_threads, std::vector<double>& lower) 
    {
        auto Q = taskflow.for_each_index(0,nr_threads,1, [&](const std::size_t thread_no) {
            std::size_t edges_batch_size = input.edges().size()/nr_threads + 1;
            const std::size_t first_edge = thread_no*edges_batch_size;
            const std::size_t last_edge = std::min((thread_no+1)*edges_batch_size, input.edges().size());
            for(std::size_t e=first_edge; e<last_edge; ++e){
                auto edge = input.edges()[e];
                if(edge.cost < 0.0)
                    repulsive_edge_vec[thread_no].push_back({edge[0], edge[1], edge.cost});
                else if(edge.cost > 0.0)
                    positive_edge_vec[thread_no].push_back({edge[0], edge[1], edge.cost});
                lower[thread_no] += std::min(0.0, edge.cost[0]);
            }
        });
        return Q;
    }


    void add_edge(edge_dict& M, const std::size_t i, const std::size_t j, const std::size_t k, 
        std::vector<edge_item>& edge_to_triangle, std::size_t& edge_index, const std::size_t triangle_index)
    {
        auto search = M.find({i,j});
        if(search != M.end()){
            edge_index = search->second;
            edge_to_triangle[edge_index].triangle_indices.push_back({k, triangle_index});
        } else {
            M[{i,j}] = edge_to_triangle.size();
            edge_index = M[{i,j}];
            edge_item e({i,j});
            e.triangle_indices.push_back({k, triangle_index});
            edge_to_triangle.push_back(e);
        }
    }

    void add_triangles(const std::vector<std::size_t> cycle, double cycle_cap, triangle_dict& T, edge_dict& M, 
        std::vector<triangle_item>& triangle_to_edge, std::vector<edge_item>& edge_to_triangle)
    {
        for(std::size_t c=1; c<cycle.size()-1; ++c) {
            // std::cout << "Processing triangle: " << cycle[0] << "," << cycle[c] << "," << cycle[c+1] << std::endl;
            std::array<std::pair<size_t, int>,3> sorted_nodes = {std::make_pair(cycle[0],0), std::make_pair(cycle[c],1), 
                                                                 std::make_pair(cycle[c+1],2)};
            std::sort(sorted_nodes.begin(), sorted_nodes.end(), 
                    [](std::pair<size_t, double> n1, std::pair<size_t, double> n2) {return n1.first < n2.first;});
            std::size_t i=sorted_nodes[0].first, j=sorted_nodes[1].first, k=sorted_nodes[2].first;
            double w_ij, w_jk, w_ik;
            if (sorted_nodes[0].second == 1){
                w_ij = cycle_cap; w_jk = -cycle_cap; w_ik = cycle_cap;
            } else if (sorted_nodes[1].second == 1){
                w_ij = cycle_cap; w_jk = cycle_cap; w_ik = -cycle_cap;
            } else {
                w_ij = -cycle_cap; w_jk = cycle_cap; w_ik = cycle_cap;
            }
            std::size_t triangle_index = 0;
            std::array<std::size_t, 3> edge_indices = {0, 0, 0};
            {
                std::lock_guard<std::mutex> guard(T_mutex);
                auto search = T.find({i,j,k});
                if (search != T.end()){
                    triangle_index = search->second;
                    assert(i == triangle_to_edge[triangle_index].nodes[0]);
                    assert(j == triangle_to_edge[triangle_index].nodes[1]);
                    assert(k == triangle_to_edge[triangle_index].nodes[2]);
                    atomic_addition(triangle_to_edge[triangle_index].weights[0], w_ij);
                    atomic_addition(triangle_to_edge[triangle_index].weights[1], w_jk);
                    atomic_addition(triangle_to_edge[triangle_index].weights[2], w_ik);
                } else {
                    T[{i,j,k}] = triangle_to_edge.size();
                    triangle_index = T[{i,j,k}];
                    triangle_item t({i,j,k});
                    t.set_weights(w_ij, w_jk, w_ik);
                    triangle_to_edge.push_back(t);
                }
            }
            {
                std::lock_guard<std::mutex> guard(M_mutex);
                add_edge(M, i, j, k, edge_to_triangle, edge_indices[0], triangle_index);
                add_edge(M, j, k, i, edge_to_triangle, edge_indices[1], triangle_index);
                add_edge(M, i, k, j, edge_to_triangle, edge_indices[2], triangle_index);
            }
            {
                std::lock_guard<std::mutex> guard(T_mutex);
                triangle_to_edge[triangle_index].set_edge_indices(edge_indices);
            }
        }
    }

    template<bool ADD_TRIANGLES=true>
    void cp_single(std::vector<edge_t>& repulsive_edges, graph<CopyableAtomic<double>>& pos_edges_graph,
        const std::size_t& cycle_length, const std::size_t& no_nodes, const bool& record_cycles, std::atomic<double>& lower_bound, 
        cycle_packing& cp, triangle_dict& T=empty_t_dict, edge_dict& E=empty_e_dict, std::vector<triangle_item>& triangle_to_edge=empty_t, 
        std::vector<edge_item>& edge_to_triangle=empty_e, std::vector<edge_t>& other_edges=empty_o)
    {
        bfs_data<graph<CopyableAtomic<double>>> bfs(pos_edges_graph);
        std::size_t progress = 0;
        union_find uf(no_nodes);

    outer:
        for(std::size_t i=0; i<repulsive_edges.size(); ++i) {
            auto& re = repulsive_edges[i];

            // periodically update component labeling for speed-up
            progress++;
            if (progress > 0.05 * repulsive_edges.size()) {
                compute_connectivity(uf, pos_edges_graph);
                progress = 0;
            }

            // check if conflicted cycle exists and repulsive edge has positive weight
            if (!connected(uf, re[0], re[1])){
                if(-re.cost <= tolerance || std::max(re[0], re[1]) >= pos_edges_graph.no_nodes()) {} 
                else{
                    other_edges.push_back(edge_t{re[0], re[1], re.cost});
                }
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
                cycle = bfs.find_path(re[0], re[1], mask_small_edges, cycle_capacity); // possibly do not reallocate for every search but give as argument to find_pathatomic_edge_t
                cycle_cap = std::min(cycle_cap, -re.cost);
                if(cycle.size() > 0 && cycle_cap >= tolerance) {
                    if(record_cycles)
                        cp.add_cycle(cycle.begin(), cycle.end(), cycle_cap);
                    
                    // subtract minimum weight and remove edges of negligible weight
                    assert(-re.cost>= cycle_cap && re.cost<= 0.0);
                    re.cost += cycle_cap;
                    assert(re.cost <= 0.0);

                    for(std::size_t c=1; c<cycle.size(); ++c) {
                        atomic_addition(pos_edges_graph.edge(cycle[c-1], cycle[c]), -cycle_cap);
                        atomic_addition(pos_edges_graph.edge(cycle[c], cycle[c-1]), -cycle_cap);
                        if(pos_edges_graph.edge(cycle[c-1], cycle[c]) < 0.0 || pos_edges_graph.edge(cycle[c], cycle[c-1]) < 0.0){
                            // std::cout << "Edge is negative after subtracted by cycle_cap.\n";
                            for(std::size_t c_p=1; c_p <=c; ++c_p){
                                atomic_addition(pos_edges_graph.edge(cycle[c_p-1], cycle[c_p]), cycle_cap);
                                atomic_addition(pos_edges_graph.edge(cycle[c_p], cycle[c_p-1]), cycle_cap);
                            }
                            goto outer;
                        } 
                    }

                    if constexpr(ADD_TRIANGLES)
                        add_triangles(cycle, cycle_cap, T, E, triangle_to_edge, edge_to_triangle);
                    atomic_addition(lower_bound, cycle_cap);

                    // remove repulsive edge if weight is negligible
                    if (-re.cost<= tolerance) {
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

    void collect_edges(const triangle_dict& T, edge_dict& M, std::vector<edge_item>& edge_to_triangle, std::vector<edge_t>& other_edges, 
        const std::vector<std::vector<edge_t>>& other_edges_container, const graph<CopyableAtomic<double>>& pos_edges_graph){
        std::cout << "Number of triangles:" << T.size() << std::endl;

        for (auto edges: other_edges_container){
            for (auto& e: edges){
                auto i = std::min(e[0], e[1]); 
                auto j = std::max(e[0], e[1]);
                auto search = M.find({i,j});
                if(search != M.end()){
                    if (search->second != -1)
                       edge_to_triangle[search->second].cost = e.cost;
                } else {
                    M[{i,j}] = -1;
                    other_edges.push_back(edge_t{i, j, e.cost});
                }
            }
        }
        pos_edges_graph.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost){
            if (i < j) {
                auto search = M.find({i,j});
                if(search != M.end()){
                    if (search->second != -1)
                        edge_to_triangle[search->second].cost = cost;
                }
                else other_edges.push_back(edge_t{i,j,cost});
            }
        });
        std::cout << "Collecting edges not in the triangles = " << other_edges.size() << "\n";
    }

    cycle_packing multicut_cycle_packing_parallel_impl(const multicut_instance& input, const bool record_cycles, const int nr_threads, 
        const int triangulation, std::vector<triangle_item>& triangle_to_edge=empty_t, std::vector<edge_item>& edge_to_triangle=empty_e, 
        std::vector<edge_t>& other_edges=empty_o)
    {
        const auto begin_time = std::chrono::steady_clock::now();
        tf::Executor executor(nr_threads);
        tf::Taskflow taskflow;
        cycle_packing cp;

        std::vector<double> lower(nr_threads, 0.0);
        std::atomic<double> lower_bound = 0.0;
        std::vector<edge_t> positive_edges;
        triangle_dict T;
        edge_dict M;
        std::vector<std::vector<edge_t>> other_edges_container(nr_threads);
        std::vector<std::vector<edge_t>> repulsive_edge_vec(nr_threads);
        std::vector<std::vector<edge_t>> positive_edge_vec(nr_threads);
        std::size_t no_nodes = input.no_nodes();

        auto init_edge = distribute_edges_round_robin(taskflow, input, positive_edge_vec, repulsive_edge_vec, nr_threads, lower);
        //auto init_edge = distribute_edges_in_chunks(taskflow, input, positive_edge_vec, repulsive_edge_vec, nr_threads, lower);

        graph<CopyableAtomic<double>> pos_edges_graph;
        auto construct_pos_edges_graph = [&]() {
            for (auto l: lower) lower_bound = lower_bound + l;
            for (auto& pos: positive_edge_vec) std::move(pos.begin(), pos.end(), std::back_inserter(positive_edges));
            pos_edges_graph.construct(positive_edges.begin(), positive_edges.end(),
                [](const edge_t& e) { return e.cost; }); };
        auto CPE = taskflow.emplace(construct_pos_edges_graph);

        init_edge.precede(CPE);

        std::future<void> fu = executor.run(taskflow);
        fu.get();

        const auto initialization_end_time = std::chrono::steady_clock::now();
        std::cout << "CP initialization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(initialization_end_time - begin_time).count() << " milliseconds\n";
        std::cout << "Cycle packing: initial lower bound = " << lower_bound << "\n";
        std::cout << "#attractive edges = " << positive_edges.size() << "\n";
        std::cout << "pos_edges_graph.no_nodes:" << pos_edges_graph.no_nodes() << std::endl;
        std::cout << "pos_edges_graph.no_edges:" << pos_edges_graph.no_edges() << std::endl;

        // iteratively pack cycles of given length
        // const std::array<std::size_t,11> cycle_lengths = {1,2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
        const std::array<std::size_t,10> cycle_lengths = {1,2,3,4,5,6,7,8,9,10};

        for(const std::size_t cycle_length : cycle_lengths) {

            taskflow.clear();
            std::cout << "find cycles of length " << cycle_length << " lower bound = " << lower_bound << "\n";

            //shuffling can give great speed-up if edges with similar indices are spatially close in the graph
            // for (auto& re: repulsive_edge_vec)
            //     std::random_shuffle(re.begin(), re.end());

            auto CP = taskflow.for_each_index(0, nr_threads, 1, [&](const std::size_t thread_no){
            //    std::cout << "#repulsive edges :" << repulsive_edge_vec[thread_no].size() << " remaining in thread " << thread_no << std::endl;
                if (triangulation)
                    cp_single<true>(repulsive_edge_vec[thread_no], pos_edges_graph, cycle_length,  
                        no_nodes, record_cycles, lower_bound, cp, T, M, triangle_to_edge, edge_to_triangle, other_edges_container[thread_no]);
                else 
                    cp_single<false>(repulsive_edge_vec[thread_no], pos_edges_graph, cycle_length,  
                        no_nodes, record_cycles, lower_bound, cp);
            });

            executor.run(taskflow);
            executor.wait_for_all();
        }

        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "CP parallel optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        size_t remain_repulsive_edges = 0;
        for (auto v: repulsive_edge_vec) { 
            for (auto e: v)
                other_edges.push_back(edge_t{e[0], e[1], e.cost});
            remain_repulsive_edges += v.size(); 
        }

        std::cout << "Cycle Packing final lower bound = " << lower_bound << "\n";
        std::cout << "#repulsive edges = " << remain_repulsive_edges << "\n";

        if (triangulation) collect_edges(T, M, edge_to_triangle, other_edges, other_edges_container, pos_edges_graph);

        return cp;
    }

    void multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads)
    {
        multicut_cycle_packing_parallel_impl(input, false, nr_threads, false);
    }

    cycle_packing compute_multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads)
    {
        return multicut_cycle_packing_parallel_impl(input, true, nr_threads, false);
    }

    cycle_packing cycle_packing_triangulation_parallel(const multicut_instance& input, const int nr_threads, 
        std::vector<triangle_item>& triangle_to_edge, std::vector<edge_item>& edge_to_triangle, std::vector<edge_t>& other_edges){
        return multicut_cycle_packing_parallel_impl(input, false, nr_threads, true, triangle_to_edge, edge_to_triangle, other_edges);
    }


}
