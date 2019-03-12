#pragma once

namespace LPMP {

    template<class BASE_CONSTRUCTOR, 
        typename QUADRUPLET_FACTOR,
        typename TRIPLET_QUADRUPLET_MESSAGE_012,
        typename TRIPLET_QUADRUPLET_MESSAGE_013,
        typename TRIPLET_QUADRUPLET_MESSAGE_023,
        typename TRIPLET_QUADRUPLET_MESSAGE_123>
            class cut_base_quadruplet_constructor : public BASE_CONSTRUCTOR {
                public:
                    using FMC = typename BASE_CONSTRUCTOR::FMC;
                    using base_constructor = BASE_CONSTRUCTOR;

                    using quadruplet_factor_container = QUADRUPLET_FACTOR;
                    using triplet_quadruplet_message_012_container = TRIPLET_QUADRUPLET_MESSAGE_012;
                    using triplet_quadruplet_message_013_container = TRIPLET_QUADRUPLET_MESSAGE_013;
                    using triplet_quadruplet_message_023_container = TRIPLET_QUADRUPLET_MESSAGE_023;
                    using triplet_quadruplet_message_123_container = TRIPLET_QUADRUPLET_MESSAGE_123;

                    using base_constructor::base_constructor;

                    // add triplet indices additionally to tripletIndices_
                    //typename base_constructor::edge_factor* add_unary_factor(const std::size_t i1, const std::size_t i2, const double cost);
                    // add triplet indices additionally to tripletIndices_
                    struct triplet_connection {
                        std::array<std::size_t,2> nodes;
                        typename base_constructor::triplet_factor* f;

                        bool operator<(const triplet_connection& b) const
                        {
                            if(this->nodes[0] != b.nodes[0]) {
                                return this->nodes[0] < b.nodes[0];
                            } else {
                                return this->nodes[1] < b.nodes[1];
                            }
                        }
                    };

                    using triplet_connections = two_dim_variable_array<triplet_connection>;
                    triplet_connections compute_connected_triplets() const;
                    quadruplet_factor_container* add_quadruplet_factor(const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3);
                    bool has_quadruplet_factor(const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3) const;
                    quadruplet_factor_container* get_quadruplet_factor(const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3) const;

                    template<typename INSTANCE>
                        INSTANCE export_quadruplets() const;

                    template<typename ITERATOR>
                        void triangulate_odd_wheel(const std::size_t center_node, ITERATOR path_begin, ITERATOR path_end, const double cost);

                    std::size_t enumerate_quadruplets();

                    std::tuple<triplet_quadruplet_message_012_container*,triplet_quadruplet_message_013_container*,triplet_quadruplet_message_023_container*,triplet_quadruplet_message_123_container*>
                        get_triplet_quadruplet_messages(quadruplet_factor_container* f) const;
                    void send_messages_to_triplets(quadruplet_factor_container* f);
                    void send_messages_to_triplets();

                    auto& quadruplet_factors() { return quadruplet_factor_vector_; }
                    const auto& quadruplet_factors() const { return quadruplet_factor_vector_; }

                protected:
                    std::unordered_map<std::array<std::size_t,4>, quadruplet_factor_container*> quadruplet_factors_;
                    std::vector<std::pair<std::array<std::size_t,4>, quadruplet_factor_container*>> quadruplet_factor_vector_;
            };


    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        QUADRUPLET_FACTOR* cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::add_quadruplet_factor(const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3)
        {
            assert(!has_quadruplet_factor(i0,i1,i2,i3));
            assert(quadruplet_factors_.size() == quadruplet_factor_vector_.size());

            auto* f = this->lp_->template add_factor<quadruplet_factor_container>();
            std::array<std::size_t,4> idx{i0,i1,i2,i3};
            quadruplet_factors_.insert(std::make_pair(idx, f));
            quadruplet_factor_vector_.push_back(std::make_pair(idx, f));

            if(!this->has_triplet_factor(i0,i1,i2)) {
                this->add_triplet_factor(i0,i1,i2);
            }
            auto* t012 = this->get_triplet_factor(i0,i1,i2);
            auto* m012 = this->lp_->template add_message<triplet_quadruplet_message_012_container>(t012, f);
            this->lp_->add_factor_relation(t012, f);

            if(!this->has_triplet_factor(i0,i1,i3)) {
                this->add_triplet_factor(i0,i1,i3);
            }
            auto* t013 = this->get_triplet_factor(i0,i1,i3);
            auto* m013 = this->lp_->template add_message<triplet_quadruplet_message_013_container>(t013, f);

            if(!this->has_triplet_factor(i0,i2,i3)) {
                this->add_triplet_factor(i0,i2,i3);
            }
            auto* t023 = this->get_triplet_factor(i0,i2,i3);
            auto* m023 = this->lp_->template add_message<triplet_quadruplet_message_023_container>(t023, f);

            if(!this->has_triplet_factor(i1,i2,i3)) {
                this->add_triplet_factor(i1,i2,i3);
            }
            auto* t123 = this->get_triplet_factor(i1,i2,i3);
            auto* m123 = this->lp_->template add_message<triplet_quadruplet_message_123_container>(t123, f);
            this->lp_->add_factor_relation(f, t123);

            return f;
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        bool cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::has_quadruplet_factor(const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3) const
        {
            assert(i0 < i1 && i1 < i2 && i2 < i3 && i3 < this->no_nodes());
            return quadruplet_factors_.find(std::array<std::size_t,4>({i0,i1,i2,i3})) != quadruplet_factors_.end();
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        QUADRUPLET_FACTOR* cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::get_quadruplet_factor(const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3) const
        {
            assert(has_quadruplet_factor(i0,i1,i2,i3));
            return quadruplet_factors_.find(std::array<std::size_t,4>({i0,i1,i2,i3}))->second; 
        }

//    template<typename ITERATOR>
//        template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
//        void cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::triangulate_odd_wheel(const std::size_t i, const double cost, ITERATOR path_begin, ITERATOR path_end, std::vector<quadruplet_candidate>& candidates)
//        {
//            assert(std::distance(path_begin, path_end) >= 3);
//            this->cycle_normal_form(path_begin, path_end);
//            const std::size_t first_node = *path_begin;
//            for(auto it=path_begin+1; it+1!=path_end; ++it) {
//                std::array<std::size_t,4> nodes({i,first_node, *it, *(it+1)});
//                std::sort(nodes.begin(), nodes.end());
//                assert(HasUniqueValues(nodes));
//                candidates.push_back({nodes, cost});
//            } 
//        }
//
    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        template<typename INSTANCE>
        INSTANCE cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::export_quadruplets() const
        {
            INSTANCE output = this->template export_triplets<INSTANCE>(); // TODO: we can implement export_triplets more efficiently by iterating over triplet vector
            for(const auto& q : quadruplet_factors())
                output.add_quadruplet({q.first[0], q.first[1], q.first[2], q.first[3]}, *q.second->get_factor());
            return output;
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        template<typename ITERATOR>
        void cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::triangulate_odd_wheel(const std::size_t center_node, ITERATOR path_begin, ITERATOR path_end, const double cost)
        {
            using edge_factor_type = typename base_constructor::triplet_factor::FactorType;
            auto get_edge_func = [&,this](const std::array<std::size_t,2> edge_nodes) -> edge_factor_type& {
                std::array<std::size_t,3> nodes{edge_nodes[0], edge_nodes[1], center_node};
                std::sort(nodes.begin(), nodes.end());
                if(!this->has_triplet_factor(nodes[0], nodes[1], nodes[2]))
                    this->add_triplet_factor(nodes[0], nodes[1], nodes[2]);
                return *(this->get_triplet_factor(nodes[0], nodes[1], nodes[2])->get_factor());
            };

            using triplet_factor_type = typename quadruplet_factor_container::FactorType;
            auto get_triplet_func = [&,this](const std::array<std::size_t,3> triplet_nodes) -> triplet_factor_type& {
                std::array<std::size_t,4> nodes{triplet_nodes[0], triplet_nodes[1], triplet_nodes[2], center_node};
                std::sort(nodes.begin(), nodes.end());
                if(!this->has_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3]))
                    this->add_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3]);
                return *(this->get_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3])->get_factor());
            };

            using msg_type = std::variant<typename TRIPLET_QUADRUPLET_MESSAGE_012::MessageType, typename TRIPLET_QUADRUPLET_MESSAGE_013::MessageType, typename TRIPLET_QUADRUPLET_MESSAGE_023::MessageType, typename TRIPLET_QUADRUPLET_MESSAGE_123::MessageType>;
            auto get_msg_func = [&](const std::array<std::size_t,2> e, const std::array<std::size_t,3> t) -> msg_type {
                std::array<std::size_t,3> e_nodes {e[0], e[1], center_node};
                std::sort(e_nodes.begin(), e_nodes.end());
                std::array<std::size_t,4> t_nodes {t[0], t[1], t[2], center_node};
                std::sort(t_nodes.begin(), t_nodes.end());

                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[2])
                    return typename TRIPLET_QUADRUPLET_MESSAGE_012::MessageType();
                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[3])
                    return typename TRIPLET_QUADRUPLET_MESSAGE_013::MessageType();
                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[3])
                    return typename TRIPLET_QUADRUPLET_MESSAGE_023::MessageType();
                assert(e_nodes[0] == t_nodes[1] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[3]);
                return typename TRIPLET_QUADRUPLET_MESSAGE_123::MessageType();
            };

            LPMP::triangulate_cycle(path_begin, path_end, get_edge_func, get_triplet_func, get_msg_func, cost);
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        std::size_t cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::enumerate_quadruplets()
        {
            return 0;
            /*
            std::vector<std::size_t> adjacency_list_count(this->noNodes_,0);
            // first determine size for adjacency_list
            for(auto& it : this->edge_factors_vector_) {
                const std::size_t i = std::get<0>(it.first);
                const std::size_t j = std::get<1>(it.first);
                adjacency_list_count[i]++;
                adjacency_list_count[j]++; 
            }
            two_dim_variable_array<std::size_t> adjacency_list(adjacency_list_count);
            std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
            for(auto& it : this->edge_factors_vector_) {
                const std::size_t i = std::get<0>(it.first);
                const std::size_t j = std::get<1>(it.first);
                assert(i<j);
                adjacency_list[i][adjacency_list_count[i]] = j;
                adjacency_list_count[i]++;
                adjacency_list[j][adjacency_list_count[j]] = i;
                adjacency_list_count[j]++;
            }

            // Sort the adjacency list, for fast intersections later
#pragma omp parallel for  schedule(guided)
            for(int i=0; i < adjacency_list.size(); i++) {
                std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
            } 

            std::vector<quadruplet_candidate> quadruplet_candidates;

#pragma omp parallel
            {
                std::vector<quadruplet_candidate> quadruplet_candidates_local;
                typename quadruplet_factor_container::FactorType test_quadruplet_factor;

                std::vector<std::size_t> commonNodes(this->noNodes_);
#pragma omp for schedule(guided) nowait
                for(std::size_t i=0; i<triplet_vector_.size(); ++i) {
                    const std::size_t i1 = std::get<0>(triplet_vector_[i]);
                    const std::size_t i2 = std::get<1>(triplet_vector_[i]);
                    const std::size_t i3 = std::get<2>(triplet_vector_[i]);
                    auto* f = std::get<3>(triplet_vector_[i])->get_factor();
                    // search for node node k such that there exists an edge from i1, i2 and i3 to it.
                    // For this, intersect adjacency lists coming from i1,i2 and i3
                    // to do: do not intersect twice but do it at once
                    auto intersects_iter_end = std::set_intersection(adjacency_list[i1].begin(), adjacency_list[i1].end(), adjacency_list[i2].begin(), adjacency_list[i2].end(), commonNodes.begin());
                    intersects_iter_end = std::set_intersection(commonNodes.begin(), intersects_iter_end, adjacency_list[i3].begin(), adjacency_list[i3].end(), commonNodes.begin());
                    for(auto it=commonNodes.begin(); it!=intersects_iter_end; ++it) {
                        std::fill(test_quadruplet_factor.begin(), test_quadruplet_factor.end(), 0.0);
                        const std::size_t k = *it; 
                        std::array<std::size_t,4> nodes({i1,i2,i3,k});
                        std::sort(nodes.begin(), nodes.end());
                        if(!has_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3])
                                && this->has_triplet_factor(nodes[0], nodes[1], nodes[2])
                                && this->has_triplet_factor(nodes[0], nodes[1], nodes[3])
                                && this->has_triplet_factor(nodes[0], nodes[2], nodes[3])
                                && this->has_triplet_factor(nodes[1], nodes[2], nodes[3]) 
                                && k == nodes[3] // otherwise we will iterate over the same set of nodes 4 times
                          ) { 

                            double lb = 0.0;
                            // go over all 3-subset of nodes containing k
                            {
                                auto* sub_f = this->get_triplet_factor(nodes[0], nodes[1], nodes[2])->get_factor();
                                lb += sub_f->LowerBound();
                                typename triplet_quadruplet_message_012_container::MessageType m;
                                m.RepamRight(test_quadruplet_factor, *sub_f); 
                            }
                            {
                                auto* sub_f = this->get_triplet_factor(nodes[0], nodes[1], nodes[3])->get_factor();
                                lb += sub_f->LowerBound();
                                typename triplet_quadruplet_message_013_container::MessageType m;
                                m.RepamRight(test_quadruplet_factor, *sub_f); 
                            }
                            { 
                                auto* sub_f = this->get_triplet_factor(nodes[0], nodes[2], nodes[3])->get_factor();
                                lb += sub_f->LowerBound();
                                typename triplet_quadruplet_message_023_container::MessageType m;
                                m.RepamRight(test_quadruplet_factor, *sub_f); 
                            }
                            { 
                                auto* sub_f = this->get_triplet_factor(nodes[1], nodes[2], nodes[3])->get_factor();
                                lb += sub_f->LowerBound();
                                typename triplet_quadruplet_message_123_container::MessageType m;
                                m.RepamRight(test_quadruplet_factor, *sub_f); 
                            } 
                            const double best_labeling_cost = test_quadruplet_factor.LowerBound();
                            assert(lb <= best_labeling_cost + eps);
                            if(lb < best_labeling_cost - eps) {
                                quadruplet_candidates_local.push_back({nodes, best_labeling_cost - lb});
                            }
                        }
                    }
                }
#pragma omp critical
                quadruplet_candidates.insert(quadruplet_candidates.end(), quadruplet_candidates_local.begin(), quadruplet_candidates_local.end());
            }

            std::sort(quadruplet_candidates.begin(), quadruplet_candidates.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });
            if(quadruplet_candidates.size() > 2) {
                for(std::size_t i=0; i<quadruplet_candidates.size()-1; ++i) {
                    assert(quadruplet_candidates[i].cost > quadruplet_candidates[i+1].cost);
                }
            }

            std::size_t factors_added = 0;
            for(std::size_t i=0; i<quadruplet_candidates.size(); ++i) {
                auto& nodes = quadruplet_candidates[i].nodes;
                if(!has_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3])) {
                    add_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3]);
                    ++factors_added;
                }
                if(factors_added > max_factors_to_add) {
                    break;
                }
            }

            if(quadruplet_candidates.size() > 0 && diagnostics()) {
                std::cout << "added " << factors_added << " by local odd 3 wheel search with guaranteed dual improvement " << quadruplet_candidates[0].cost << "\n";
            }

            return factors_added;
            */
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
       std::tuple<TRIPLET_QUADRUPLET_MESSAGE_012*,TRIPLET_QUADRUPLET_MESSAGE_013*,TRIPLET_QUADRUPLET_MESSAGE_023*,TRIPLET_QUADRUPLET_MESSAGE_123*>
        cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::get_triplet_quadruplet_messages(quadruplet_factor_container* f) const
        {
            auto msg_012 = f->template get_messages<triplet_quadruplet_message_012_container>();
            auto msg_013 = f->template get_messages<triplet_quadruplet_message_013_container>();
            auto msg_023 = f->template get_messages<triplet_quadruplet_message_023_container>();
            auto msg_123 = f->template get_messages<triplet_quadruplet_message_123_container>();
            assert(msg_012.size() == 1 && msg_013.size() == 1 && msg_023.size() == 1 && msg_123.size() == 1);
            return std::make_tuple(msg_012[0], msg_013[0], msg_023[0], msg_123[0]);
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        void cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::send_messages_to_triplets(quadruplet_factor_container* f)
        {
            auto [msg_012, msg_013, msg_023, msg_123] = get_triplet_quadruplet_messages(f);

            msg_012->send_message_to_right();
            msg_013->send_message_to_right();
            msg_023->send_message_to_right();
            msg_123->send_message_to_right();

            msg_012->send_message_to_left(0.25);
            msg_013->send_message_to_left(1.0/3.0);
            msg_023->send_message_to_left(0.5);
            msg_123->send_message_to_left(1.0);

            msg_012->send_message_to_left(1.0/3.0);
            msg_013->send_message_to_left(0.5);
            msg_023->send_message_to_left(1.0);

            msg_012->send_message_to_left(0.5);
            msg_013->send_message_to_left(1.0);

            msg_012->send_message_to_left(1.0); 
        }

    template<class BASE_CONSTRUCTOR, typename QUADRUPLET_FACTOR, typename TRIPLET_QUADRUPLET_MESSAGE_012, typename TRIPLET_QUADRUPLET_MESSAGE_013, typename TRIPLET_QUADRUPLET_MESSAGE_023, typename TRIPLET_QUADRUPLET_MESSAGE_123>
        void cut_base_quadruplet_constructor<BASE_CONSTRUCTOR, QUADRUPLET_FACTOR, TRIPLET_QUADRUPLET_MESSAGE_012, TRIPLET_QUADRUPLET_MESSAGE_013, TRIPLET_QUADRUPLET_MESSAGE_023, TRIPLET_QUADRUPLET_MESSAGE_123>::send_messages_to_triplets()
        {
            for(auto f : quadruplet_factors())
                send_messages_to_triplets(f.second);
        }

} // namespace LPMP
