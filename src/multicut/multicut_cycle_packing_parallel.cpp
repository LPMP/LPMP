#include "cut_base/cut_base_instance.hxx"
#include "multicut/multicut_cycle_packing_parallel.h"
#include "multicut/multicut_factors_messages.h"
#include "cut_base/cut_base_apply_packing.hxx"
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

namespace LPMP {

    constexpr static double tolerance = 1e-8;

    struct weighted_edge : public std::array<std::size_t,2> { double cost; };

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

    void compute_connectivity(union_find& uf,graph<CopyableAtomic<double>>& pos_edges_graph) {
        uf.reset();
        pos_edges_graph.for_each_edge([&](const std::size_t i, const std::size_t j, const double cost) {
            if(cost >= tolerance) { uf.merge(i,j); }
        });
    };

    bool connected(union_find& uf, const std::size_t i, const std::size_t j) { return uf.connected(i,j); };

    tf::Task distribute_edges_round_robin(tf::Taskflow& taskflow, const multicut_instance& input, std::vector<std::vector<weighted_edge>>& positive_edge_vec, std::vector<std::vector<weighted_edge>>& repulsive_edge_vec, const int& nr_threads, std::vector<double>& lower) {
        auto Q = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
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
        return Q.second;
    }

    tf::Task distribute_edges_in_chunks(tf::Taskflow& taskflow, const multicut_instance& input, std::vector<std::vector<weighted_edge>>& positive_edge_vec, std::vector<std::vector<weighted_edge>>& repulsive_edge_vec, const int& nr_threads, std::vector<double>& lower) {
        auto Q = taskflow.parallel_for(0,nr_threads,1, [&](const std::size_t thread_no) {
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
        return Q.second;
    }

    void cp_single(std::vector<weighted_edge>& repulsive_edges,graph<CopyableAtomic<double>>& pos_edges_graph, const std::size_t& cycle_length, std::size_t& no_nodes, const bool& record_cycles, std::atomic<double>& lower_bound, cycle_packing& cp){

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
                    if(record_cycles) cp.add_cycle(cycle.begin(), cycle.end(), cycle_cap);

                    // subtract minimum weight and remove edges of negligible weight
                    assert(-re.cost>= cycle_cap && re.cost<= 0.0);
                    re.cost += cycle_cap;
                    assert(re.cost <= 0.0);
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

    cycle_packing multicut_cycle_packing_parallel_impl(const multicut_instance& input, const bool record_cycles, const int nr_threads)
    {
        const auto begin_time = std::chrono::steady_clock::now();
        tf::Executor executor(nr_threads);
        tf::Taskflow taskflow;
        cycle_packing cp;

        std::vector<double> lower(nr_threads, 0.0);
        std::atomic<double> lower_bound = 0.0;
        std::vector<weighted_edge> positive_edges;
        std::vector<std::vector<weighted_edge>> repulsive_edge_vec(nr_threads);
        std::vector<std::vector<weighted_edge>> positive_edge_vec(nr_threads);
        std::size_t no_nodes = input.no_nodes();

        //auto init_edge = distribute_edges_round_robin(taskflow, input, positive_edge_vec, repulsive_edge_vec, nr_threads, lower);
        auto init_edge = distribute_edges_in_chunks(taskflow, input, positive_edge_vec, repulsive_edge_vec, nr_threads, lower);

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
        std::cout << "initialization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(initialization_end_time - begin_time).count() << " milliseconds\n";

        std::cout << "cycle packing\n";
        std::cout << "initial lower bound = " << lower_bound << "\n";
        std::cout << "#attractive edges = " << positive_edges.size() << "\n";
        std::cout << "pos_edges_graph.no_nodes:" << pos_edges_graph.no_nodes() << std::endl;
        std::cout << "pos_edges_graph.no_edges:" << pos_edges_graph.no_edges() << std::endl;

        // iteratively pack cycles of given length
        const std::array<std::size_t,11> cycle_lengths = {1,2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
        for(const std::size_t cycle_length : cycle_lengths) {

            taskflow.clear();
            std::cout << "find cycles of length " << cycle_length << " lower bound = " << lower_bound << "\n";

            //shuffling can give great speed-up if edges with similar indices are spatially close in the graph
            for (auto& re: repulsive_edge_vec)
                std::random_shuffle(re.begin(), re.end());

            auto CP = taskflow.parallel_for(0, nr_threads, 1, [&](const std::size_t thread_no){
            //    std::cout << "#repulsive edges :" << repulsive_edge_vec[thread_no].size() << " remaining in thread " << thread_no << std::endl;
                cp_single(std::ref(repulsive_edge_vec[thread_no]), std::ref(pos_edges_graph), std::ref(cycle_length),  std::ref(no_nodes), std::ref(record_cycles), std::ref(lower_bound), std::ref(cp));
            });

            executor.run(taskflow);
            executor.wait_for_all();
            // terminate if no conflicted cycle remain
            size_t remain_repulsive_edges = 0;
            for (auto& re: repulsive_edge_vec) { remain_repulsive_edges += re.size(); }
            if (remain_repulsive_edges == 0) break;
        }

        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Parallel Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        size_t remain_repulsive_edges = 0;
        for (auto& re: repulsive_edge_vec) { remain_repulsive_edges += re.size(); }

        std::cout << "final lower bound = " << lower_bound << "\n";
        std::cout << "#repulsive edges = " << remain_repulsive_edges << "\n";

        return cp;
    }

    void multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads)
    {
        multicut_cycle_packing_parallel_impl(input, false, nr_threads);
    }

    cycle_packing compute_multicut_cycle_packing_parallel(const multicut_instance& input, const int& nr_threads)
    {
        return multicut_cycle_packing_parallel_impl(input, true, nr_threads);
    }

    // apply the cycle packing
    triplet_multicut_instance pack_multicut_instance(const multicut_instance& input, const cycle_packing& cp)
    {
        return pack_cut_base_instance<multicut_edge_triplet_message_0, multicut_edge_triplet_message_1, multicut_edge_triplet_message_2 , triplet_multicut_instance>(input, cp);
    }

}
