#pragma once

#include <array>
#include <vector>
#include "maxflow/maxflow.h"
#include <map> // TODO: change to unordered_map

namespace LPMP {

    template<class BASE_CONSTRUCTOR,
        class CUT_CONSTRUCTOR,
        typename LIFTED_CUT_FACTOR,
        typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG>
            class lifted_constructor : public CUT_CONSTRUCTOR {
                public:
                    using base_constructor = BASE_CONSTRUCTOR;
                    using class_type = lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG>;
                    using lifted_cut_factor_container = LIFTED_CUT_FACTOR;
                    using cut_edge_lifted_factor_message_container = CUT_EDGE_LIFTED_FACTOR_MSG;
                    using cut_edge_lifted_factor_message = typename cut_edge_lifted_factor_message_container::MessageType;
                    using lifted_edge_lifted_factor_message_container = LIFTED_EDGE_LIFTED_FACTOR_MSG;
                    using lifted_edge_lifted_factor_message = typename lifted_edge_lifted_factor_message_container::MessageType;

                    // do zrobienia: use this everywhere instead of std::array<std::size_t,2>
                    struct Edge : public std::array<std::size_t,2> {
                        Edge(const std::size_t i, const std::size_t j) : std::array<std::size_t,2>({std::min(i,j), std::max(i,j)}) {}
                    };
                    using CutId = std::vector<Edge>;

                    template<typename SOLVER>
                        lifted_constructor(SOLVER& pd) : CUT_CONSTRUCTOR(pd) {}

                    typename CUT_CONSTRUCTOR::edge_factor* add_edge_factor(const std::size_t i1, const std::size_t i2, const double cost);
                    typename CUT_CONSTRUCTOR::edge_factor* add_lifted_edge_factor(const std::size_t i1, const std::size_t i2, const double cost);
                    bool has_cut_factor(const CutId& cut);
                    bool has_lifted_edge_in_cut_factor(const CutId& cut, const std::size_t i1, const std::size_t i2);

                    // do zrobienia: provide add_cut_factor(const CutId& cut, const std::size_t i1, const std::size_t i2) as well
                    lifted_cut_factor_container* add_cut_factor(const CutId& cut);
                    lifted_cut_factor_container* get_cut_factor(const CutId& cut);
                    void AddLiftedEdge(const CutId& cut, const std::size_t i1, const std::size_t i2);
                    std::size_t Tighten(const std::size_t max_factors_to_add);
                    double FindViolatedCutsThreshold(const std::size_t max_triplets_to_add);
                    std::size_t FindViolatedCuts(const std::size_t minDualIncrease, const std::size_t noConstraints);

                    // check if all lifted edges are primally consistent by asserting that a path of zero values exists in the ground graph whenever lifted edge is zero
                    bool CheckPrimalConsistency() const;

                private:
                    //struct Edge {std::size_t i; std::size_t j; double w;}; // replace by WeightedEdge
                    struct weighted_edge {
                        std::size_t i; 
                        std::size_t j; 
                        typename CUT_CONSTRUCTOR::edge_factor* f;
                        double weight() const { return (*f->get_factor())[0]; }
                    };
                    bool addingTighteningEdges = false; // controls whether edges are added to baseEdges_
                    std::vector<weighted_edge> baseEdges_;
                    std::vector<weighted_edge> liftedEdges_;

                    std::vector<std::vector<std::size_t>> cutEdgesLiftedFactors_;
                    std::vector<std::vector<std::size_t>> liftedEdgesLiftedFactors_;

                    std::map<CutId,std::pair<lifted_cut_factor_container*,std::vector<Edge>>> liftedFactors_;
            }; 

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        typename CUT_CONSTRUCTOR::edge_factor* lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::add_edge_factor(const std::size_t i1, const std::size_t i2, const double cost)
        {
            assert(i1<i2);
            auto* f = CUT_CONSTRUCTOR::add_edge_factor(i1,i2,cost);
            if(!addingTighteningEdges) {
                baseEdges_.push_back({i1, i2, f});
            } else { // all auxiliary edges may be regarded as lifted ones
                add_lifted_edge_factor(i1, i2, cost);
            }
            return f;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        typename CUT_CONSTRUCTOR::edge_factor* lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::add_lifted_edge_factor(const std::size_t i1, const std::size_t i2, const double cost)
        {
            auto* f = CUT_CONSTRUCTOR::add_edge_factor(i1, i2, cost);
            liftedEdges_.push_back({i1,i2,f}); 
            return f;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        bool lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::has_cut_factor(const CutId& cut) 
        {
            assert(std::is_sorted(cut.begin(), cut.end()));
            return liftedFactors_.find(cut) != liftedFactors_.end();
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        bool lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::has_lifted_edge_in_cut_factor(const CutId& cut, const std::size_t i1, const std::size_t i2)
        {
            assert(i1<i2);
            assert(has_cut_factor(cut));
            const auto& edgeList = liftedFactors_[cut].second;
            // speed this up by holding edge list sorted
            return std::find(edgeList.begin(),edgeList.end(),Edge({i1,i2})) != edgeList.end();
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        LIFTED_CUT_FACTOR* lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::add_cut_factor(const CutId& cut)
        {
            assert(!has_cut_factor(cut));
            //std::cout << "Add cut with edges ";
            //for(auto i : cut) { std::cout << "(" << std::get<0>(i) << "," << std::get<1>(i) << ");"; } std::cout << "\n";
            auto* f = new lifted_cut_factor_container(cut.size());
            CUT_CONSTRUCTOR::lp_->add_factor(f);
            // connect the cut edges
            for(std::size_t e=0; e<cut.size(); ++e) {
                auto* unaryFactor = CUT_CONSTRUCTOR::get_edge_factor(cut[e][0],cut[e][1]);
                //new cut_edge_lifted_factor_message_container(cut_edge_lifted_factor_message(e),unaryFactor,f);
                auto* m = CUT_CONSTRUCTOR::lp_->template add_message<cut_edge_lifted_factor_message_container>(unaryFactor, f, e);
            }
            liftedFactors_.insert(std::make_pair(cut,std::make_pair(f,std::vector<Edge>())));
            return f;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        LIFTED_CUT_FACTOR* lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::get_cut_factor(const CutId& cut)
        {
            assert(has_cut_factor(cut));
            return liftedFactors_[cut].first;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        void lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::AddLiftedEdge(const CutId& cut, const std::size_t i1, const std::size_t i2)
        {
            assert(!has_lifted_edge_in_cut_factor(cut,i1,i2));
            assert(has_cut_factor(cut));
            assert(i1<i2);
            //std::cout << "Add lifted edge (" << i1 << "," << i2 << ") to cut\n";
            auto& c = liftedFactors_[cut];
            auto* f = c.first;
            auto& edgeList = c.second;
            auto* unaryFactor = CUT_CONSTRUCTOR::get_edge_factor(i1,i2);
            f->get_factor()->IncreaseLifted();
            //new lifted_edge_lifted_factor_message_container(lifted_edge_lifted_factor_message(edgeList.size() + cut.size()), unaryFactor, f);
            auto* m = CUT_CONSTRUCTOR::lp_->template add_message<lifted_edge_lifted_factor_message_container>(unaryFactor, f, edgeList.size() + cut.size());
            c.second.push_back(Edge({i1,i2}));
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        std::size_t lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::Tighten(const std::size_t max_factors_to_add)
        {
            const bool prevMode = addingTighteningEdges;
            addingTighteningEdges = true;
            assert(max_factors_to_add > 5); //otherwise the below arrangement makes little sense.
            const std::size_t noBaseConstraints = CUT_CONSTRUCTOR::Tighten(std::ceil(0.8*max_factors_to_add));
            //return noBaseConstraints;
            std::size_t noLiftingConstraints = 0;
            if(debug()) {
                std::cout << "number of cut constraints = " << liftedFactors_.size() << "\n";
            }
            if(noBaseConstraints < max_factors_to_add) {
                double th = FindViolatedCutsThreshold(max_factors_to_add - noBaseConstraints);
                if(th >= 0.0) {
                    noLiftingConstraints = FindViolatedCuts(th, max_factors_to_add - noBaseConstraints);
                    if(diagnostics()) {
                        std::cout << "added " << noLiftingConstraints << " lifted cut factors.\n";
                    }
                }
            }
            addingTighteningEdges = prevMode;
            return noBaseConstraints + noLiftingConstraints;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        double lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::FindViolatedCutsThreshold(const std::size_t max_triplets_to_add)
        {
            // make one function to reuse allocated datastructures.
            union_find uf(this->noNodes_);
            std::vector<std::tuple<std::size_t,std::size_t,double>> edges;
            edges.reserve(baseEdges_.size());
            for(const auto& e : baseEdges_) {
                const std::size_t i = e.i;
                const std::size_t j = e.j;
                const double weight = e.weight();
                if(weight > 0) {
                    uf.merge(i,j);
                } else {
                    edges.push_back(std::make_tuple(i,j,weight));
                }
            }
            std::sort(edges.begin(),edges.end(), [] (auto& a, auto& b)->bool { return std::get<2>(a) > std::get<2>(b); });
            std::vector<std::list<std::tuple<std::size_t,double>>> liftedEdges(this->noNodes_);
            for(auto& e : liftedEdges_) {
                const std::size_t i = e.i;
                const std::size_t j = e.j;
                const double weight = e.weight();
                if(weight > 0.0) {
                    liftedEdges[i].push_front(std::make_tuple(j,weight));
                    liftedEdges[j].push_front(std::make_tuple(i,weight));
                }
            }
            auto edge_sort = [] (const std::tuple<std::size_t,double>& a, const std::tuple<std::size_t,double>& b)->bool { return std::get<0>(a) < std::get<0>(b); };
            for(std::size_t i=0; i<liftedEdges.size(); ++i) {
                liftedEdges[i].sort(edge_sort);
            }
            double prevWeight = 0.0;
            double maxTh = -std::numeric_limits<double>::infinity();
            for(std::size_t e=0; e<edges.size(); ++e) {
                const std::size_t i = std::get<0>(edges[e]);
                const std::size_t j = std::get<1>(edges[e]);
                const double weight = std::get<2>(edges[e]);
                if(uf.find(i) != uf.find(j)) {
                    uf.merge(i,j);
                    const std::size_t c = uf.find(i);
                    // (i) merge edges emanating from same connected component
                    if(i != c) {
                        liftedEdges[c].merge(liftedEdges[i],edge_sort);
                    } 
                    if(j != c) {
                        liftedEdges[c].merge(liftedEdges[j],edge_sort);
                    }
                    std::transform(liftedEdges[c].begin(), liftedEdges[c].end(), liftedEdges[c].begin(), 
                            [&uf,&maxTh,weight] (std::tuple<std::size_t,double> a)->std::tuple<std::size_t,double> {
                            const std::size_t cc = uf.find(std::get<0>(a));
                            if(cc != std::get<0>(a)) {
                            const double th = std::min(-weight, std::get<1>(a));
                            maxTh = std::max(th,maxTh); // do zrobienia: correct?
                            }
                            std::get<0>(a) = cc;
                            return a; 
                            });
                    // do zrobienia: unnecessary sort here: has been sorted more or less after merge, only parallel edges need to be sorted
                    liftedEdges[c].sort([] (const std::tuple<std::size_t,double>& a, const std::tuple<std::size_t,double>& b)->bool { 
                            if(std::get<0>(a) != std::get<0>(b)) {
                            return std::get<0>(a) < std::get<0>(b); 
                            }
                            return std::get<1>(a) > std::get<1>(b); // this ensures that remove removes copies with smaller weight. Thus, the largest one only remains.
                            });
                    // (ii) remove lifted edges that have been collapsed. do zrobienia: record weight of those
                    liftedEdges[c].remove_if([&uf,c,maxTh] (const auto& a)->bool { 
                            if(std::get<1>(a) < maxTh) return true; // we need not take this edge into account anymore, as it would lead to smaller minimum dual improvement.
                            // first check whether lifted edge belonged to different connected components before the last merge operation. If yes, update threshold
                            const std::size_t cc = uf.find(std::get<0>(a));
                            return uf.find(std::get<0>(a)) == c; 
                            });
                    // (iii) take maximal representative of edges that are parallel now
                    liftedEdges[c].unique([] (const auto& a, const auto& b)->bool { return std::get<0>(a) == std::get<0>(b); }); 
                }

            }

            //if(maxTh != -std::numeric_limits<double>::infinity()) {
            //   std::cout << "\nnon-trivial cut factor with weight = " << maxTh << "\n\n";
            //}
            return 0.1*maxTh;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        std::size_t lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::FindViolatedCuts(const std::size_t minDualIncrease, const std::size_t noConstraints)
        {
            // update weight of base edges
            //for(auto& e : baseEdges_) {
            //   e.w = CUT_CONSTRUCTOR::get_edge_factor(e.i,e.j)->operator[](0);
            //}
            //std::sort(baseEdges_.begin(), baseEdges_.end(), [](const Edge& e1, const Edge& e2) { return e1.weight() < e2.weight(); });
            union_find uf(CUT_CONSTRUCTOR::noNodes_);
            for(const auto& e : baseEdges_) {
                if(e.weight() >= -minDualIncrease) {
                    uf.merge(e.i,e.j);
                }
            }

            // build reduced graph with connected components as nodes and edge weights as number of edges with weight < -minDualIncrease

            // union find indices are not contiguous. Make them so, to use them as identifiers for connected components
            std::map<std::size_t,std::size_t> ufIndexToContiguous;
            for(std::size_t i=0; i<CUT_CONSTRUCTOR::noNodes_; ++i) {
                const std::size_t ufIndex = uf.find(i);
                if(ufIndexToContiguous.find(ufIndex) == ufIndexToContiguous.end()) {
                    ufIndexToContiguous.insert(std::make_pair(ufIndex,ufIndexToContiguous.size()));
                }
            }
            const std::size_t ccNodes = ufIndexToContiguous.size();

            std::map<std::size_t,std::size_t> origToCompressedNode; // union find index to compressed node indices // do zrobienia: use hash map
            for(std::size_t i=0; i<CUT_CONSTRUCTOR::noNodes_; ++i) {
                const std::size_t ufIndex = uf.find(i);
                const std::size_t collapsedIndex = ufIndexToContiguous[ufIndex];
                origToCompressedNode.insert(std::make_pair(i,collapsedIndex));
            }

            std::size_t ccEdges = 0;
            std::map<Edge,std::vector<Edge>> ccToBaseEdges;
            for(const auto& e : baseEdges_) {
                if(!uf.connected(e.i,e.j)) {
                    const std::size_t i = ufIndexToContiguous[uf.find(e.i)];
                    const std::size_t j = ufIndexToContiguous[uf.find(e.j)];
                    //if(ccEdgeCap.find({i,j}) == cc.EdgeCap.end()) {
                    //   ccEdgeCap.insert(std::make_pair(std::array<std::size_t,2>(i,j),1));
                    //}
                    //ccEdgeCap[std::array<std::size_t,2>(i,j)]++;
                    if(ccToBaseEdges.find(Edge(i,j)) == ccToBaseEdges.end()) {
                        ++ccEdges;
                        ccToBaseEdges.insert(std::make_pair(Edge(i,j),std::vector<Edge>()));
                    }
                    ccToBaseEdges[Edge(i,j)].push_back(Edge(e.i,e.j));
                }
            }
            maxflow::Graph<int,int,int> maxFlow(ccNodes, ccEdges); 
            maxFlow.add_node(ccNodes);
            for(auto& e : ccToBaseEdges) {
                const std::size_t i = e.first.operator[](0);
                const std::size_t j = e.first.operator[](1);
                const std::size_t cap = e.second.size();
                maxFlow.add_edge(i,j,cap,cap);
            }

            // note: this can possibly be made faster by caching the weight
            std::sort(liftedEdges_.begin(), liftedEdges_.end(), [](const weighted_edge& e1, const weighted_edge& e2) { return e1.weight() > e2.weight(); });
            std::size_t factorsAdded = 0;

            struct empty {};
            graph<empty> base_graph(ccToBaseEdges.begin(), ccToBaseEdges.end(), ccNodes);
            bfs_data<graph<empty>> bfs(base_graph);

            // note that currently possibly multiple max flow computations are performed, when lifted edges come from the same connected components. This is superfluous and searches could be remembered and reused.
            const int capacityMax = baseEdges_.size()+1;
            for(const auto& liftedEdge : liftedEdges_) {
                if(factorsAdded >= noConstraints) { 
                    if(debug()) {
                        std::cout << "maximal number of constraints to add reached";
                    }
                    break; 
                }
                if(liftedEdge.weight() > minDualIncrease) {
                    const std::size_t i = origToCompressedNode[liftedEdge.i];
                    const std::size_t j = origToCompressedNode[liftedEdge.j];
                    if(!uf.connected(i,j)) {
                        // find minimum cut in unweighted graph containing only base edges with weight < eps
                        maxFlow.add_tweights(i,capacityMax,0);
                        maxFlow.add_tweights(j,0,capacityMax);
                        const std::size_t noCutEdges = maxFlow.maxflow();
                        assert(noCutEdges > 0 && noCutEdges < baseEdges_.size()); // otherwise there is no path from i to j or all paths were collapsed
                        std::vector<Edge> minCut;
                        minCut.reserve(noCutEdges);
                        // now do a dfs from i on those vertices which are in the same segment as i. The edges from those to vertices in segment of j form a minimum cut.
                        auto cut_edge_op = [&](const std::size_t i, const std::size_t j, auto edge_info) -> bool {
                            const bool cut = maxFlow.what_segment(i) != maxFlow.what_segment(j);
                            if(cut)
                                minCut.push_back(Edge(i,j));
                            return cut;
                        };

                        bfs.traverse_current_component(i, cut_edge_op);
                        /*
                           std::stack<std::size_t> q;
                           std::vector<bool> visited(ccNodes,false);
                           q.push(i);
                           while(!q.empty()) {
                           const std::size_t v = q.top();
                           q.pop();
                           if(visited[v]) {
                           continue;
                           }
                           visited[v] = true;
                           auto* a = maxFlow.get_first_arc(v);
                           while(a != nullptr) {
                           int v_test, w; // do zrobienia: use proper type
                           maxFlow.get_arc_ends(a, v_test, w);
                           assert(v != w);
                           assert(v_test == v);
                           if(maxFlow.what_segment(std::size_t(v)) != maxFlow.what_segment(std::size_t(w))) {
                        // expand all edges that were collapsed into (v,w) in the original graph
                        //spdlog::get("logger")->info() << "edge in mincut: " << v << "," << w;
                        for(const auto& e : ccToBaseEdges[Edge(std::size_t(v),std::size_t(w))]) {
                        //spdlog::get("logger")->info() << " expanded edge : " << e[0] << "," << e[1];
                        minCut.push_back(Edge(e[0],e[1]));
                        }
                        } else if(!visited[std::size_t(w)]) {
                        q.push(std::size_t(w));
                        }
                        a = a->next;
                        }
                        }
                        */
                        assert(minCut.size() == noCutEdges);
                        std::sort(minCut.begin(),minCut.end()); // unique form of cut

                        // check if minimum cut is already present
                        if(!has_cut_factor(minCut)) {
                            auto* f = add_cut_factor(minCut);
                            AddLiftedEdge(minCut,liftedEdge.i,liftedEdge.j);
                            ++factorsAdded;
                        } else {
                            auto* MinCutFactor = get_cut_factor(minCut);
                            if(!has_lifted_edge_in_cut_factor(minCut,liftedEdge.i,liftedEdge.j)) {
                                AddLiftedEdge(minCut,liftedEdge.i,liftedEdge.j);
                                ++factorsAdded;
                            }
                        }

                        // restore original terminal weights
                        maxFlow.add_tweights(i,-capacityMax,0);
                        maxFlow.add_tweights(j,0,-capacityMax);
                    }
                }
            }

            return factorsAdded;
        }

    template< class BASE_CONSTRUCTOR, class CUT_CONSTRUCTOR, typename LIFTED_CUT_FACTOR, typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG >
        bool lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG >::CheckPrimalConsistency() const
        {
            const bool multicutConsistent = CUT_CONSTRUCTOR::CheckPrimalConsistency();
            if(!multicutConsistent) {
                return false;
            }

            //collect connectivity information with union find w.r.t. base edges
            union_find uf(CUT_CONSTRUCTOR::noNodes_);
            for(const auto& e : baseEdges_) {
                if(e.f->get_factor()->primal()[0] == false) {
                    uf.merge(e.i,e.j);
                }
            }
            for(const auto& e : liftedEdges_) {
                if(e.f->get_factor()->primal()[0] == false) {
                    if(!uf.connected(e.i,e.j)) {
                        return false;
                    }
                }
            }
            return true;
        }

} // namespace LPMP
