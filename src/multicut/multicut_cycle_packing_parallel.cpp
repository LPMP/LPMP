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
    //struct atomic_weighted_edge : public std::array<std::size_t,2> { std::atomic<double> cost; };
    struct atomic_weighted_edge {
        std::array<size_t,2> edges;
        std::atomic<double> cost;
        atomic_weighted_edge(){}
        ~atomic_weighted_edge(){}
        atomic_weighted_edge(std::array<size_t,2> e, double c): edges(e), cost(c) {}
        atomic_weighted_edge(size_t e1, size_t e2, double c): edges({e1, e2}), cost(c) {}
        atomic_weighted_edge(const atomic_weighted_edge& we){ edges = we.edges; cost = we.cost.load(); }
        atomic_weighted_edge& operator=(const atomic_weighted_edge& we){
            edges = we.edges;
            cost = we.cost.load();
            return *this;
        }
    };

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
                //std::cout<<"update component labeling\n";
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

        std::atomic<double> lower_bound = 0.0;
        std::vector<weighted_edge> repulsive_edges;
        std::vector<weighted_edge> positive_edges;
        std::vector<std::vector<weighted_edge>> repulsive_edge_vec(nr_threads);
        std::size_t no_nodes = input.no_nodes();
        std::vector<atomic_weighted_edge> positive_atomic_weighted_edge;

        auto classify_edges = [&] () {
            for (const auto& e : input.edges()){
                if(e.cost < 0.0)
                    repulsive_edges.push_back({e[0], e[1], e.cost});
                else if(e.cost > 0.0){
                    std::atomic<double> atomic_cost;
                    atomic_cost.store(e.cost);
                    positive_edges.push_back({e[0], e[1], e.cost});
                    positive_atomic_weighted_edge.push_back({e[0], e[1], atomic_cost});
                }
                lower_bound = lower_bound + std::min(0.0, e.cost[0]);
            }
        };
        auto CE = taskflow.emplace(classify_edges);

        graph<CopyableAtomic<double>> pos_edges_graph;
        auto construct_pos_edges_graph = [&]() { pos_edges_graph.construct(positive_atomic_weighted_edge.begin(), positive_atomic_weighted_edge.end(),
                [](const atomic_weighted_edge& e) -> std::array<size_t,2> { return {e.edges[0], e.edges[1]}; },
                [](const atomic_weighted_edge& e) { return e.cost.load(); }); };
        auto CPE = taskflow.emplace(construct_pos_edges_graph);

        CPE.gather(CE);

        auto RV = taskflow.parallel_for(0, nr_threads, 1, [&](const std::size_t thread_no) {
            std::size_t e = thread_no;
            while (e < repulsive_edges.size()) {
                repulsive_edge_vec[thread_no].push_back(repulsive_edges[e]);
                e += nr_threads;
            }
        });
        RV.first.gather(CE);

        std::future<void> fu = executor.run(taskflow);
        fu.get();

        std::cout << "cycle packing\n";
        std::cout << "initial lower bound = " << lower_bound << "\n";
        std::cout << "#repulsive edges = " << repulsive_edges.size() << "\n";
        std::cout << "#attractive edges = " << positive_edges.size() << "\n";
        std::cout << "pos_edges_graph.no_nodes:" << pos_edges_graph.no_nodes() << std::endl;
        std::cout << "pos_edges_graph.no_edges:" << pos_edges_graph.no_edges() << std::endl;

        // iteratively pack cycles of given length
        const std::array<std::size_t,11> cycle_lengths = {1,2,3,4,5,6,7,8,9,10,std::numeric_limits<std::size_t>::max()};
        for(const std::size_t cycle_length : cycle_lengths) {

            tf::Executor executor2(nr_threads);
            tf::Taskflow taskflow2;

            std::cout << "find cycles of length " << cycle_length << " lower bound = " << lower_bound << "\n";

            // shuffling can give great speed-up if edges with similar indices are spatially close in the graph
            for (auto& re: repulsive_edge_vec)
                std::random_shuffle(re.begin(), re.end());

            auto CP = taskflow2.parallel_for(0, nr_threads, 1, [&](const std::size_t thread_no){
            //    std::cout << "#repulsive edges :" << repulsive_edge_vec[thread_no].size() << " remaining in thread " << thread_no << std::endl;
                cp_single(std::ref(repulsive_edge_vec[thread_no]), std::ref(pos_edges_graph), std::ref(cycle_length),  std::ref(no_nodes), std::ref(record_cycles), std::ref(lower_bound), std::ref(cp));
            });

            executor2.run(taskflow2);
            executor2.wait_for_all();
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
