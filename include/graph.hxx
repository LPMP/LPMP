#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include <queue>
#include "two_dimensional_variable_array.hxx"
#include "union_find.hxx"
#include "help_functions.hxx"

namespace LPMP {

    // possibly templatize for support for sister pointers and masking edges out
    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER=true, bool SUPPORT_MASKING=false>
        class graph {
            public:

                using type = graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>;
                using edge_information = EDGE_INFORMATION;

                class edge_type {
                    public:
                        const std::size_t tail() const
                        {
                            assert(sister_->sister_ == this);
                            return sister_->head();
                        }
                        const std::size_t head() const { return head_; }
                        const edge_type& sister() const { return *sister_; }
                        edge_type& sister() { return *sister_; }

                        EDGE_INFORMATION& edge() { return edge_; }
                        const EDGE_INFORMATION& edge() const { return edge_; }

                        bool operator<(const edge_type& o) const { return head() < o.head(); }
                        template<typename T>
                        edge_type& operator=(const T& o) { edge() = o; }

                        std::size_t head_;
                        edge_type* sister_;
                        EDGE_INFORMATION edge_;
                };

                // compute sorted adjacency list representation of graph given by edges.
                template<typename EDGE_ITERATOR, typename EDGE_ENDPOINTS_LAMBDA, typename EDGE_INFORMATION_LAMBDA>
                    void construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_ENDPOINTS_LAMBDA e, EDGE_INFORMATION_LAMBDA f);

                template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
                    void construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA f);

                graph() {}

                template<typename EDGE_ITERATOR, typename EDGE_ENDPOINTS_LAMBDA, typename EDGE_INFORMATION_LAMBDA>
                    graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_ENDPOINTS_LAMBDA e, EDGE_INFORMATION_LAMBDA f);

                constexpr static auto return_edge_op = [](const auto& edge) -> EDGE_INFORMATION { return EDGE_INFORMATION{}; };
                constexpr static auto return_endpoints_op = [](const auto& edge) -> std::array<std::size_t,2>{ return {edge[0], edge[1]}; };

                template<typename EDGE_ITERATOR>
                    graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end)
                    : graph(edge_begin, edge_end, return_endpoints_op, return_edge_op)
                    {}

                template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
                    graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA f)
                    : graph(edge_begin, edge_end, return_endpoints_op, f)
                    {}

                void set_sister_pointers();

                // TODO: implement move operator
                graph& operator=(const graph& o);

                template<typename LAMBDA>
                    void for_each_edge(LAMBDA f) const;

                struct edge_comparator {
                    bool operator()(const std::size_t i, const edge_type& e) const { return i < e.head(); }
                    bool operator()(const edge_type& e, const std::size_t i) const { return e.head() < i; }
                };

                bool edge_present(const std::size_t i, const std::size_t j) const;
                EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j);
                const EDGE_INFORMATION& edge(const std::size_t i, const std::size_t j) const;

                std::size_t no_nodes() const { return edges_.size(); }
                std::size_t no_edges() const { return edges_.no_elements(); }
                std::size_t no_edges(const std::size_t i) const { return edges_[i].size(); }
                std::size_t edge_no(const edge_type* e) const;
                std::size_t edge_no(const edge_type& e) const;

                auto begin(const std::size_t i) { return edges_[i].begin(); }
                auto end(const std::size_t i) { return edges_[i].end(); }
                auto begin(const std::size_t i) const { return edges_[i].begin(); }
                auto end(const std::size_t i) const { return edges_[i].end(); }

                // enumerate all triangles and call f on each. f expects the three node indices (i<j<k) of triangles (sorted) and references to edges in lexicographical order (ij,ik,jk)
                template<typename LAMBDA>
                    void for_each_triangle(LAMBDA f) const;

                // we enumerate quadrangles with the method C4 from "Arboricity and subgraph listin algorithms" by Norishige Chiba and Takao Nishizeki
                template<typename LAMBDA>
                    void for_each_quadrangle(LAMBDA f);

                // return contracted graph and mapping from original nodes to contracted nodes
                template<typename EDGE_ITERATOR, typename MERGE_FUNC>
                    std::tuple<graph, std::vector<std::size_t>> contract(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, MERGE_FUNC merge_func) const;

                template<typename EDGE_ITERATOR>
                    std::tuple<graph, std::vector<std::size_t>> contract(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end) const;

                void check_graph() const;

            private:
                two_dim_variable_array<edge_type> edges_;
        };

    // bfs search
    template<typename GRAPH>
        struct bfs_data {
            using edge_type = typename GRAPH::edge_type;
            using edge_information = typename GRAPH::edge_information;

            struct item {
                edge_information e;
                std::size_t parent;
                std::size_t flag;
            };

            void update_graph();
            bfs_data(const GRAPH& _g);
            void reset() ;
            item& operator[](const std::size_t i) { return d[i]; }

            void label1(const std::size_t i) { d[i].flag = flag1; }
            void label2(const std::size_t i) { d[i].flag = flag2; }
            bool labelled(const std::size_t i) const { return labelled1(i) || labelled2(i); }
            bool labelled1(const std::size_t i) const { return d[i].flag == flag1; }
            bool labelled2(const std::size_t i) const { return d[i].flag == flag2; }

            std::size_t& parent(const std::size_t i) { return d[i].parent; }
            std::size_t parent(const std::size_t i) const { return d[i].parent; }

            template<typename EDGE_OP>
                std::vector<std::size_t> trace_path(const std::size_t i1, const std::size_t i2, EDGE_OP edge_op, const edge_information& edge_info) const;

            // do bfs with thresholded costs and iteratively lower threshold until enough cycles are found
            // only consider edges that have cost equal or larger than th
            constexpr static auto no_mask_op = [](const auto i, const auto j, const auto& edge_info, const std::size_t distance) -> bool { return true; };
            constexpr static auto no_edge_op = [](const auto i, const auto j, const edge_information& edge_info) {};
            auto find_path(const std::size_t start_node, const std::size_t end_node);

            template<typename MASK_OP, typename EDGE_OP>
                std::vector<std::size_t> find_path(const std::size_t start_node, const std::size_t end_node, MASK_OP mask_op, EDGE_OP edge_op);

            // traverse all edges that have at least one endpoint in current component
            template<typename CUT_EDGE_OP>
                void traverse_component(const std::size_t start_node, CUT_EDGE_OP cut_edge_op);

            private:
            std::vector<item> d;
            std::deque<std::array<std::size_t,2>> visit; // node number, distance from start or end
            std::size_t flag1, flag2;
            const GRAPH& g;
        };


    // implementation

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename EDGE_ITERATOR, typename EDGE_ENDPOINTS_LAMBDA, typename EDGE_INFORMATION_LAMBDA>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_ENDPOINTS_LAMBDA e, EDGE_INFORMATION_LAMBDA f)
        {
            std::vector<std::size_t> adjacency_list_count;
            // first determine size for adjacency_list
            for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
                const auto [i,j] = e(*edge_it);
                adjacency_list_count.resize(std::max({i+1,j+1,adjacency_list_count.size()}));
                adjacency_list_count[i]++;
                adjacency_list_count[j]++;
            }

            edges_.resize(adjacency_list_count.begin(), adjacency_list_count.end());
            std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);

            for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
                //const auto i = (*edge_it)[0];
                //const auto j = (*edge_it)[1];
                const auto [i,j] = e(*edge_it);

                edges_[i][adjacency_list_count[i]].head_ = j;
                edges_[i][adjacency_list_count[i]].edge() = f(*edge_it);
                adjacency_list_count[i]++;

                edges_[j][adjacency_list_count[j]].head_ = i;
                edges_[j][adjacency_list_count[j]].edge() = f(*edge_it);
                adjacency_list_count[j]++;
            }

            //#pragma omp parallel for schedule(guided)
            for(std::size_t i=0; i<adjacency_list_count.size(); ++i) {
                std::sort(edges_[i].begin(), edges_[i].end());
            }

            set_sister_pointers();
            check_graph();
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename EDGE_ITERATOR, typename EDGE_INFORMATION_LAMBDA>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::construct(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_INFORMATION_LAMBDA f)
        {
            return construct(edge_begin, edge_end, return_endpoints_op, f);
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename EDGE_ITERATOR, typename EDGE_ENDPOINTS_LAMBDA, typename EDGE_INFORMATION_LAMBDA>
        graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::graph(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, EDGE_ENDPOINTS_LAMBDA e, EDGE_INFORMATION_LAMBDA f)
        {
            construct(edge_begin, edge_end, e, f);
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::set_sister_pointers()
        {
            for(std::size_t i=0; i<edges_.size(); ++i) { assert(std::is_sorted(edges_[i].begin(), edges_[i].end())); }
            std::vector<std::size_t> adjacency_list_count(edges_.size(), 0);
            for(std::size_t i=0; i<edges_.size(); ++i) {
                for(auto edge_it=edges_[i].begin(); edge_it!=edges_[i].end(); ++edge_it) {
                    if(edge_it->head() > i) {
                        const auto head = edge_it->head();
                        auto* sister = &edges_(head, adjacency_list_count[head]++);
                        edge_it->sister_ = sister;
                        sister->sister_ = &(*edge_it);
                    }
                }
            }
        }

    // TODO: implement move operator
    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>& graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::operator=(const graph& o)
        {
            edges_ = o.edges_;
            set_sister_pointers();
            return *this;
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename LAMBDA>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::for_each_edge(LAMBDA f) const
        {
            for(std::size_t i=0; i<edges_.size(); ++i) {
                for(auto edge_it=edges_[i].begin(); edge_it!=edges_[i].end(); ++edge_it) {
                    const std::size_t j = edge_it->head();
                    if(i<j) {
                        f(i, edge_it->head(), edge_it->edge());
                    }
                }
            }

        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        bool graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::edge_present(const std::size_t i, const std::size_t j) const
        {
            assert(i != j && std::max(i,j) < no_nodes());
            auto edge_it = std::lower_bound(begin(i), end(i), j, edge_comparator{});
            return (edge_it != end(i) && edge_it->head() == j);
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        EDGE_INFORMATION& graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::edge(const std::size_t i, const std::size_t j)
        {
            assert(edge_present(i,j));
            auto edge_it = std::lower_bound(begin(i), end(i), j, edge_comparator{});
            assert(edge_it->head() == j);
            return edge_it->edge();
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        const EDGE_INFORMATION& graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::edge(const std::size_t i, const std::size_t j) const
        {
            assert(edge_present(i,j));
            auto edge_it = std::lower_bound(begin(i), end(i), j, edge_comparator{});
            assert(edge_it->head() == j);
            return edge_it->edge();
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        std::size_t graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::edge_no(const edge_type* e) const
        {
            const std::size_t i = e->tail();
            assert(std::distance(const_cast<const edge_type*>(edges_[i].begin()), e) < no_edges(i));
            return std::distance(const_cast<const edge_type*>(edges_[i].begin()), e);
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        std::size_t graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::edge_no(const edge_type& e) const
        {
            const std::size_t i = e.tail();
            assert(std::distance(const_cast<const edge_type*>(edges_[i].begin()), &e) < no_edges(i));
            return std::distance(const_cast<const edge_type*>(edges_[i].begin()), &e);
        }

    // enumerate all triangles and call f on each. f expects the three node indices (i<j<k) of triangles (sorted) and references to edges in lexicographical order (ij,ik,jk)
    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename LAMBDA>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::for_each_triangle(LAMBDA f) const
        {
            struct triangle_intersection_type {
                triangle_intersection_type(std::size_t _k, const edge_type* _ik, const edge_type* _jk) : k(_k), ik(_ik), jk(_jk) {}
                std::size_t k;
                const edge_type* ik;
                const edge_type* jk;
            };
            auto edge_intersection_merge = [](const edge_type& e1, const edge_type& e2) -> triangle_intersection_type {
                assert(e1.head() == e2.head());
                const auto k = e1.head();
                return triangle_intersection_type(k, &e1, &e2);
            };
            auto edge_intersection_sort = [](const edge_type& e1, const edge_type& e2) { return e1.head() < e2.head(); };

            std::vector<triangle_intersection_type> common_nodes;
            for(std::size_t i=0; i<no_nodes(); ++i) {
                for(auto edge_it=begin(i); edge_it!=end(i); ++edge_it) {
                    const auto j = edge_it->head();
                    if(i<j) {
                        // Now find all neighbors of both i and j to see where the triangles are
                        common_nodes.clear();
                        auto intersects_iter_end = set_intersection_merge(
                                edges_[i].begin(), edges_[i].end(), edges_[j].begin(), edges_[j].end(),
                                std::back_inserter(common_nodes), edge_intersection_sort, edge_intersection_merge);

                        for(const triangle_intersection_type t : common_nodes) {
                            const auto k = t.k;

                            // Since a triplet shows up three times as an edge plus a node, we only consider it for the case when i<j<k
                            if(!(j<k)) { continue; }

                            assert(t.ik->tail() == i && t.jk->tail() == j);
                            //f(i,j,k, *edge_it, *(t.ik), *(t.jk));
                            f(i,j,k, edge_it->edge(), t.ik->edge(), t.jk->edge()); // edge costs: 01, 02, 12
                        }
                    }
                }
            }
        }

    // we enumerate quadrangles with the method C4 from "Arboricity and subgraph listin algorithms" by Norishige Chiba and Takao Nishizeki
    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename LAMBDA>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::for_each_quadrangle(LAMBDA f)
        {
            struct node_degree { std::size_t i, degree; };
            std::vector<node_degree> node_degrees(no_nodes());
            for(std::size_t i=0; i<no_nodes(); ++i) {
                node_degrees[i].i = i;
                node_degrees[i].degree = no_edges(i);
            }
            std::sort(node_degrees.begin(), node_degrees.end(), [](const node_degree& n1, const node_degree& n2) { return n1.degree > n2.degree; });

            // structure for recording masked edges
            two_dim_variable_array<unsigned char> edge_mask(edges_.size_begin(), edges_.size_end());
            edge_mask.set(true);

            std::vector<std::vector<std::size_t>> U(no_nodes()); // TODO: more efficient datastructure?
            std::vector<std::size_t> relevant_U_entries;
            std::vector<std::size_t> nonempty_U_entries;

            for(const auto& n : node_degrees) {
                const std::size_t v = n.i;
                for(std::size_t c_v=0; c_v<no_edges(v); ++c_v) {
                    if(edge_mask(v,c_v) == false) { continue; }
                    const auto u = edges_(v,c_v).head();
                    assert(u != v);
                    for(std::size_t c_u=0; c_u<no_edges(u); ++c_u) {
                        if(edge_mask(u,c_u) == false) { continue; }
                        const std::size_t w = edges_(u,c_u).head();
                        if(w != v) {
                            U[w].push_back(u);
                            if(U[w].size() == 1) { nonempty_U_entries.push_back(w); }
                            if(U[w].size() == 2) { relevant_U_entries.push_back(w); }
                        }
                    }
                }
                for(const std::size_t w : relevant_U_entries) {
                    assert(U[w].size() >= 2);
                    for(std::size_t i=0; i<U[w].size(); ++i) {
                        for(std::size_t j=i+1; j<U[w].size(); ++j) {
                            f({v, U[w][i], w, U[w][j]});
                        }
                    }
                }
                relevant_U_entries.clear();
                for(const auto w : nonempty_U_entries)
                    U[w].clear();
                nonempty_U_entries.clear();
                for(const auto& l : U) { assert(l.empty()); }

                // mask all edges with v as its endpoint
                for(std::size_t c=0; c<no_edges(v); ++c) {
                    const auto u = edges_(v,c).head();
                    edge_mask(v,c) = false;
                    edge_mask(u, edge_no( edges_(v,c).sister() )) = false;
                }
            }
        }


    // return contracted graph and mapping from original nodes to contracted nodes
    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename EDGE_ITERATOR, typename MERGE_FUNC>
        std::tuple<graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>, std::vector<std::size_t>> graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::contract(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end, MERGE_FUNC merge_func) const
        {
            union_find uf(no_nodes());
            for(auto edge_it=edge_begin; edge_it!=edge_end; ++edge_it) {
                const auto i = (*edge_it)[0];
                const auto j = (*edge_it)[1];
                uf.merge(i,j);
            }
            auto contracted_ids = uf.get_contiguous_ids();
            struct contracted_edge_type : public std::array<std::size_t,2> {
                contracted_edge_type(const std::size_t i, const std::size_t j, EDGE_INFORMATION e) : std::array<std::size_t,2>({i,j}), edge(e) {}
                EDGE_INFORMATION edge;
            };
            std::vector<contracted_edge_type> contracted_edges;

            for_each_edge( [&](const std::size_t i, const std::size_t j, const auto& edge) {
                    const auto contracted_i = contracted_ids[ uf.find(i) ];
                    const auto contracted_j = contracted_ids[ uf.find(j) ];
                    if(contracted_i != contracted_j) {
                        contracted_edge_type contracted_edge(contracted_i, contracted_j, edge);
                        contracted_edges.push_back(contracted_edge);
                    }
            });

            // TODO: faster sorting by first sorting for tail ids via bucket sort, and then separately sorting for head ids.
            // merge contracted edges
            auto edge_sort = [](const auto e1, const auto e2) { return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end()); };
            std::sort(contracted_edges.begin(), contracted_edges.end(), edge_sort);
            std::vector<contracted_edge_type> unique_contracted_edges;
            for(std::size_t i=0; i<contracted_edges.size()-1; ++i) {
                if(contracted_edges[i][0] == contracted_edges[i+1][0] && contracted_edges[i][0] == contracted_edges[i+1][0]) {
                    contracted_edges[i+1].edge = merge_func(contracted_edges[i].edge, contracted_edges[i+1].edge);
                } else {
                    unique_contracted_edges.push_back(contracted_edges[i]);
                }
            }
            unique_contracted_edges.push_back(contracted_edges.back());

            type contracted_graph(unique_contracted_edges.begin(), unique_contracted_edges.end());

            std::vector<std::size_t> contraction_mapping(no_nodes());
            for(std::size_t i=0; i<no_nodes(); ++i) {
                contraction_mapping[i] = contracted_ids[ uf.find(i) ];
            }

            return std::make_tuple(contracted_graph, contraction_mapping);
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        template<typename EDGE_ITERATOR>
        std::tuple<graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>, std::vector<std::size_t>> graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::contract(EDGE_ITERATOR edge_begin, EDGE_ITERATOR edge_end) const
        {
            auto standard_merge_func = [](EDGE_INFORMATION& e1, EDGE_INFORMATION& e2) { return e1; };
            return contract(edge_begin, edge_end, standard_merge_func);
        }

    template<typename EDGE_INFORMATION, bool SUPPORT_SISTER, bool SUPPORT_MASKING>
        void graph<EDGE_INFORMATION, SUPPORT_SISTER, SUPPORT_MASKING>::check_graph() const
        {
            for(std::size_t i=0; i<edges_.size(); ++i) {
                // check that graph is simple
                assert(std::is_sorted(edges_[i].begin(), edges_[i].end()));

                // check sister pointers have been set correctly
                for(auto edge_it=edges_[i].begin(); edge_it!=edges_[i].end(); ++edge_it) {
                    assert(&edge_it->sister().sister() == &(*edge_it));
                    assert(edge_it->tail() == i);
                    if(edge_it+1 != edges_[i].end()) {
                        assert(edge_it->head() < (edge_it+1)->head());
                    }
                }

                // check access with indices
                for(std::size_t e=0; e<edges_[i].size(); ++e) {
                    assert(edges_(i,e).tail() == i);
                    assert(edges_(i,e).sister_->sister_ == &edges_(i,e));
                    assert(edge_present(i,edges_(i,e).head()));
                    assert(edge_no(&(edges_(i,e))) == e);
                }
            }
        }

    template<typename EDGE_TYPE>
        std::vector<std::size_t> nodes_of_path(const std::vector<EDGE_TYPE*>& edges)
        {
            assert(edges.size() > 0);
            std::vector<std::size_t> path;
            path.reserve(edges.size()+1);
            path.push_back(edges[0]->tail());
            for(std::size_t i=0; i<edges.size()-1; ++i) { assert(edges[i]->head() == edges[i+1]->tail()); }
            for(auto* e : edges) {
                path.push_back(e->head());
            }
            return path;
        }

    // bfs search
    template<typename GRAPH>
        void bfs_data<GRAPH>::update_graph()
        {
            d.resize(g.no_nodes());
            for(auto& i : d) { i.flag = 0; }
            flag1 = 0;
            flag2 = 1;

        }
    template<typename GRAPH>
        bfs_data<GRAPH>::bfs_data(const GRAPH& _g) : g(_g)
    {
        g.check_graph();
        update_graph();
    }
    template<typename GRAPH>
        void bfs_data<GRAPH>::reset()
        {
            visit.clear();
            flag1 += 2;
            flag2 += 2;
        }

    template<typename GRAPH>
        template<typename EDGE_OP>
        std::vector<std::size_t> bfs_data<GRAPH>::trace_path(const std::size_t i1, const std::size_t i2, EDGE_OP edge_op, const edge_information& edge_info) const
        {
            assert(i1 != i2);

            std::vector<std::size_t> path;
            path.push_back(i1);
            std::size_t j=i1;
            edge_op(i1,i2,edge_info);
            while(parent(j) != j) {
                edge_op(j, parent(j), d[j].e);
                j = parent(j);
                path.push_back(j);
            }
            std::reverse(path.begin(),path.end());
            path.push_back(i2);
            j=i2;
            while(parent(j) != j) {
                edge_op(j, parent(j), d[j].e);
                j = parent(j);
                path.push_back(j);
            }

            return path;
        }

    template<typename GRAPH>
        auto bfs_data<GRAPH>::find_path(const std::size_t start_node, const std::size_t end_node)
        {
            return find_path(start_node, end_node, no_mask_op, no_edge_op);
        }

    template<typename GRAPH>
        template<typename MASK_OP, typename EDGE_OP>
        std::vector<std::size_t> bfs_data<GRAPH>::find_path(const std::size_t start_node, const std::size_t end_node, MASK_OP mask_op, EDGE_OP edge_op)
        {
            assert(start_node != end_node);
            assert(start_node < g.no_nodes() && end_node < g.no_nodes());
            reset();
            visit.push_back({start_node, 0});
            label1(start_node);
            parent(start_node) = start_node;
            visit.push_back({end_node, 0});
            label2(end_node);
            parent(end_node) = end_node;

            while(!visit.empty()) {
                const std::size_t i = visit.front()[0];
                const std::size_t distance = visit.front()[1];
                visit.pop_front();

                assert(g.begin(i) <= g.end(i));

                if(labelled1(i)) {
                    for(auto a_it=g.begin(i); a_it!=g.end(i); ++a_it) {
                        const std::size_t j = a_it->head();

                        if(mask_op(i,j,a_it->edge(), distance)) {

                            if(!labelled(j)) {
                                visit.push_back({j, distance+1});
                                d[j].e = a_it->edge();
                                parent(j) = i;
                                label1(j);
                            } else if(labelled2(j)) { // shortest path found
                                // trace back path from j to end_node and from i to start_node
                                return trace_path(i,j, edge_op, a_it->edge());
                            }

                        }
                    }
                } else {
                    assert(labelled2(i));
                    for(auto a_it=g.begin(i); a_it!=g.end(i); ++a_it) {
                        const std::size_t j = a_it->head();

                        if(mask_op(i,j,a_it->edge(), distance)) {

                            if(!labelled(j)) {
                                visit.push_back({j, distance+1});
                                d[j].e = a_it->edge();
                                parent(j) = i;
                                label2(j);
                            } else if(labelled1(j)) { // shortest path found
                                // trace back path from j to end_node and from i to start_node
                                return trace_path(i,j, edge_op, a_it->edge());
                            }

                        }
                    }
                }


            }
            return std::vector<std::size_t>({});
        }

    // traverse all edges that have at least one endpoint in current component
    template<typename GRAPH>
        template<typename CUT_EDGE_OP>
        void bfs_data<GRAPH>::traverse_component(const std::size_t start_node, CUT_EDGE_OP cut_edge_op)
        {
            reset();
            visit.push_back({start_node, 0});
            label1(start_node);
            parent(start_node) = start_node;

            while(!visit.empty()) {
                const std::size_t i = visit.front()[0];
                const std::size_t distance = visit.front()[1];
                visit.pop_front();

                assert(g.begin(i) < g.end(i));

                if(labelled1(i)) {
                    for(auto a_it=g.begin(i); a_it!=g.end(i); ++a_it) {
                        const std::size_t j = a_it->head();

                        if(!cut_edge_op(i,j,*a_it)) {
                            if(!labelled(j)) {
                                visit.push_back({j, distance+1});
                                parent(j) = i;
                                label1(j);
                            }
                        }
                    }
                }
            }
        }

} // namespace LPMP
