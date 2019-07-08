#pragma once

#include <array>
#include <vector>
#include <map> // TODO: change to unordered_map
#include "cut_base_apply_packing.hxx"
#include "LP.h"

namespace LPMP {

    // or additionally add CRTP?
    template<class FACTOR_MESSAGE_CONNECTION,
        typename EDGE_FACTOR, typename TRIPLET_FACTOR,
        typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2>
            class cut_base_triplet_constructor {
                public:
                    using FMC = FACTOR_MESSAGE_CONNECTION;

                    using edge_factor = EDGE_FACTOR;
                    using triplet_factor = TRIPLET_FACTOR;
                    using edge_triplet_message_0 = EDGE_TRIPLET_MESSAGE_0;
                    using edge_triplet_message_1 = EDGE_TRIPLET_MESSAGE_1;
                    using edge_triplet_message_2 = EDGE_TRIPLET_MESSAGE_2;

                    template<typename SOLVER>
                        cut_base_triplet_constructor(SOLVER& pd)
                        : lp_(&pd.GetLP())
                        {}

                    bool get_edge_label(const std::size_t i0, const std::size_t i1) const;
                    edge_factor* add_edge_factor(const std::size_t i1, const std::size_t i2, const double cost);
                    edge_factor* get_edge_factor(const std::size_t i1, const std::size_t i2) const;
                    std::size_t number_of_edges() const { return edge_factors_.size(); }
                    std::size_t number_of_triplets() const { return triplet_factors_.size(); }

                    template<typename MESSAGE_CONTAINER>
                        MESSAGE_CONTAINER* link_unary_triplet_factor(edge_factor* u, triplet_factor* t);

                    triplet_factor* add_triplet_factor(const std::size_t i1, const std::size_t i2, const std::size_t i3);
                    bool has_edge_factor(const std::array<std::size_t,2> e) const;
                    bool has_edge_factor(const std::size_t i1, const std::size_t i2) const;
                    bool has_triplet_factor(const std::size_t i1, const std::size_t i2, const std::size_t i3) const;
                    triplet_factor* get_triplet_factor(const std::size_t i1, const std::size_t i2, const std::size_t i3) const;
                    std::array<std::size_t,2> get_edge(const std::size_t i1, const std::size_t i2) const;
                    double get_edge_cost(const std::size_t i1, const std::size_t i2) const;

                    std::size_t no_nodes() const { return no_nodes_; }
                    std::size_t no_edges() const { return unary_factors_vector_.size(); }

                    std::array<std::size_t,2> get_edge(const std::size_t e) const;
                    double get_edge_cost(const std::size_t e) const; 

                    std::tuple<edge_triplet_message_0*, edge_triplet_message_1*, edge_triplet_message_2*>
                        get_edge_triplet_messages(triplet_factor* t);
                    void send_messages_to_triplets();
                    void send_messages_to_edges(triplet_factor* t);
                    void send_messages_to_edges();

                    template<typename ITERATOR>
                    void triangulate_cycle(ITERATOR cycle_begin, ITERATOR cycle_end, const double cycle_weight, const bool pack_cycle = true);

                    template<typename INSTANCE>
                        INSTANCE export_edges() const;
                    template<typename INSTANCE>
                        INSTANCE export_triplets() const;
                    template<typename LABELING>
                        LABELING export_edge_labeling() const;


                    template<typename LABELING>
                    void write_labeling_into_factors(const LABELING& labeling);

                    template<typename STREAM>
                        void WritePrimal(STREAM& s);

                        void Begin() {}
                        void End() {}
                    /*
                       struct triplet_candidate {
                       std::array<std::size_t,3> nodes;
                       double cost;
                       bool operator<(const triplet_candidate& o) const {
                       return this->cost > o.cost;
                       }
                       };

                       template<typename ITERATOR>
                       void cycle_normal_form(ITERATOR cycle_begin, ITERATOR cycle_end) const
                       {
                       assert(std::distance(cycle_begin, cycle_end) >= 3);
                    //assert(std::distance(cycle_begin, cycle_end)%2 == 1);
                    // first search for smallest entry and make it first
                    std::rotate(cycle_begin, std::min_element(cycle_begin, cycle_end), cycle_end);
                    // now two choices left: we can traverse cycle in forward or backward direction. Choose direction such that second entry is smaller than in reverse directoin.
                    if(*(cycle_begin+1) > *(cycle_end - 1)) {
                    std::reverse(cycle_begin+1, cycle_end);
                    }
                    } 

                    template<typename ITERATOR>
                    void triangulate_cycle(const double cost, ITERATOR path_begin, ITERATOR path_end, std::vector<triplet_candidate>& candidates)
                    {
                    assert(std::distance(path_begin, path_end) >= 3);
                    cycle_normal_form(path_begin, path_end);
                    const std::size_t first_node = *path_begin;
                    for(auto it=path_begin+1; it+1!=path_end; ++it) {
                    if(first_node != *it && first_node != *(it+1)) {
                    std::array<std::size_t,3> nodes({first_node, *it, *(it+1)});
                    std::sort(nodes.begin(), nodes.end());
                    assert(HasUniqueValues(nodes));
                    candidates.push_back({nodes, cost});
                    }
                    } 
                    }
                    */

                        auto& triplet_factors() { return triplet_factor_vector_; }
                        const auto& triplet_factors() const { return triplet_factor_vector_; }

                protected:
                    std::map<std::array<std::size_t,2>, edge_factor*> edge_factors_; // actually unary factors in multicut are defined on edges. assume first index < second one
                    std::size_t no_original_edges_ = std::numeric_limits<std::size_t>::max();
                    //std::unordered_map<std::array<std::size_t,2>, edge_factor*> edge_factors_; // actually unary factors in multicut are defined on edges. assume first index < second one
                    std::vector<std::pair<std::array<std::size_t,2>, edge_factor*>> unary_factors_vector_; // we store a second copy of unary factors for faster iterating
                    // sort triplet factors as follows: Let indices be i=(i1,i2,i3) and j=(j1,j2,j3). Then i<j iff i1+i2+i3 < j1+j2+j3 or for ties sort lexicographically
                    std::unordered_map<std::array<std::size_t,3>, triplet_factor*> triplet_factors_;
                    std::vector<std::pair<std::array<std::size_t,3>, triplet_factor*>> triplet_factor_vector_;
                    std::size_t no_nodes_ = 0;

                    LP<FMC>* lp_;
            };



    // implementation

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2>
        bool cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge_label(const std::size_t i0, const std::size_t i1) const
        {
            assert(i0 < i1 && i1 < no_nodes_);
            assert(has_edge_factor(i0,i1));
            auto* f = get_edge_factor(i0,i1);
            return f->get_factor()->primal()[0];
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2>
        EDGE_FACTOR* cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::add_edge_factor(const std::size_t i1, const std::size_t i2, const double cost)
        {
            assert(i1 < i2);
            no_nodes_ = std::max(no_nodes_,i2+1);
            assert(!has_edge_factor(i1,i2));

            auto* u = lp_->template add_factor<edge_factor>();
            (*u->get_factor())[0] = cost;
            auto it = edge_factors_.insert(std::make_pair(std::array<std::size_t,2>{i1,i2}, u)).first;
            unary_factors_vector_.push_back(std::make_pair(std::array<std::size_t,2>{i1,i2}, u));

            if(it != edge_factors_.begin()) {
                auto prevIt = it;
                --prevIt;
                assert(prevIt->second != u);
                lp_->add_factor_relation(prevIt->second, u);
            }
            auto nextIt = it;
            ++nextIt;
            if(nextIt != edge_factors_.end()) {
                assert(nextIt->second != u);
                lp_->add_factor_relation(u, nextIt->second);
            }

            return u;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2>
        EDGE_FACTOR* cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge_factor(const std::size_t i1, const std::size_t i2) const
        {
            assert(has_edge_factor(i1,i2));
            return edge_factors_.find(std::array<std::size_t,2>{i1,i2})->second;
        }

        template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2>
    template<typename MESSAGE_CONTAINER>
        MESSAGE_CONTAINER* cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::link_unary_triplet_factor(EDGE_FACTOR* u, TRIPLET_FACTOR* t) 
        {
            return lp_->template add_message<MESSAGE_CONTAINER>(u, t);
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        TRIPLET_FACTOR* cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::add_triplet_factor(const std::size_t i1, const std::size_t i2, const std::size_t i3) // declared virtual so that derived constructor notices when triplet factor is added
        {
            assert(i1 < i2 && i2 < i3);
            assert(!has_triplet_factor(i1,i2,i3));
            assert(triplet_factors_.size() == triplet_factor_vector_.size());
            if(!has_edge_factor(i1,i2))
                add_edge_factor(i1,i2,0.0);
            if(!has_edge_factor(i1,i3))
                add_edge_factor(i1,i3,0.0);
            if(!has_edge_factor(i2,i3))
                add_edge_factor(i2,i3,0.0);
            assert(has_edge_factor(i1,i2) && has_edge_factor(i1,i3) && has_edge_factor(i2,i3));
            auto* t = lp_->template add_factor<triplet_factor>();
            std::array<std::size_t,3> idx{i1,i2,i3};
            triplet_factors_.insert(std::make_pair( idx, t ));
            triplet_factor_vector_.push_back(std::make_pair(idx,t));
            // use following ordering of unary and triplet factors: triplet comes after edge factor (i1,i2) and before (i2,i3)
            auto* before = get_edge_factor(i1,i2);
            lp_->add_factor_relation(before,t);
            auto* middle = get_edge_factor(i1,i3);
            lp_->add_factor_relation(middle,t);
            auto* after = get_edge_factor(i2,i3);
            lp_->add_factor_relation(t,after);
            // link with all three unary factors
            link_unary_triplet_factor<edge_triplet_message_0>(before, t);
            link_unary_triplet_factor<edge_triplet_message_1>(middle, t);
            link_unary_triplet_factor<edge_triplet_message_2>(after, t);
            return t;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        bool cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::has_edge_factor(const std::array<std::size_t,2> e) const 
        {
            return has_edge_factor(std::get<0>(e), std::get<1>(e));
            assert(std::get<0>(e) < std::get<1>(e));
            return edge_factors_.find(std::array<std::size_t,2>{std::get<0>(e),std::get<1>(e)}) != edge_factors_.end();
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        bool cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::has_edge_factor(const std::size_t i1, const std::size_t i2) const 
        {
            assert(i1 < i2 && i2 < no_nodes_);
            return (edge_factors_.find(std::array<std::size_t,2>{i1,i2}) != edge_factors_.end());
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        bool cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::has_triplet_factor(const std::size_t i1, const std::size_t i2, const std::size_t i3) const 
        {
            assert(i1 < i2 && i2 < i3 && i3 < no_nodes_);
            return (triplet_factors_.find(std::array<std::size_t,3>{i1,i2,i3}) != triplet_factors_.end());
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        TRIPLET_FACTOR* cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_triplet_factor(const std::size_t i1, const std::size_t i2, const std::size_t i3) const 
        {
            assert(has_triplet_factor(i1,i2,i3));
            assert(triplet_factors_.size() == triplet_factor_vector_.size());
            return triplet_factors_.find(std::array<std::size_t,3>{i1,i2,i3})->second;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        std::array<std::size_t,2> cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge(const std::size_t i1, const std::size_t i2) const
        {
            return {std::min(i1,i2), std::max(i1,i2)};
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        double cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge_cost(const std::size_t i1, const std::size_t i2) const
        {
            assert(has_edge_factor(i1,i2));
            return *(edge_factors_.find(std::array<std::size_t,2>{i1,i2})->second->get_factor());
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        std::array<std::size_t,2> cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge(const std::size_t e) const 
        {
            return unary_factors_vector_[e].first;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        double cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge_cost(const std::size_t e) const {
            return *(unary_factors_vector_[e].second->get_factor());
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        std::tuple<EDGE_TRIPLET_MESSAGE_0*,EDGE_TRIPLET_MESSAGE_1*,EDGE_TRIPLET_MESSAGE_2*>
        cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::get_edge_triplet_messages(triplet_factor* t)
        {
            auto msg_0 = t->template get_messages<edge_triplet_message_0>();
            auto msg_1 = t->template get_messages<edge_triplet_message_1>();
            auto msg_2 = t->template get_messages<edge_triplet_message_2>();
            assert(msg_0.size() == 1 && msg_1.size() == 1 && msg_2.size() == 1);
            return {msg_0[0], msg_1[0], msg_2[0]};
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        void cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::send_messages_to_triplets()
        {
            for(const auto& e : unary_factors_vector_)
                e->template send_messages_uniform<EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>();
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        void cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::send_messages_to_edges(triplet_factor* t)
        {
            // TODO: use method send_messages_to_right_residual
            auto [msg_0, msg_1, msg_2] = get_edge_triplet_messages(t); 
            assert(msg_0->GetRightFactor() == t);
            assert(msg_1->GetRightFactor() == t);
            assert(msg_2->GetRightFactor() == t);
            assert(msg_0->GetLeftFactor() != msg_1->GetLeftFactor());
            assert(msg_0->GetLeftFactor() != msg_2->GetLeftFactor());
            assert(msg_1->GetLeftFactor() != msg_2->GetLeftFactor());

            msg_0->send_message_to_right();
            msg_1->send_message_to_right();
            msg_2->send_message_to_right();

            msg_0->send_message_to_left(1.0/3.0);
            msg_1->send_message_to_left(1.0/2.0);
            msg_2->send_message_to_left(1.0/1.0);

            msg_0->send_message_to_left(1.0/2.0);
            msg_1->send_message_to_left(1.0/1.0);

            msg_0->send_message_to_left(1.0/1.0);
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        void cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::send_messages_to_edges()
        {
            for(std::size_t i=0; i<triplet_factors().size(); ++i)
                send_messages_to_edges(triplet_factors()[i].second);
            //for(auto& t : triplet_factors())
            //    send_messages_to_edges(t.second); 
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        template<typename ITERATOR>
        void cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::triangulate_cycle(ITERATOR cycle_begin, ITERATOR cycle_end, const double cycle_weight, const bool pack_cycle)
        {
            using edge_factor_type = typename EDGE_FACTOR::FactorType;
            auto get_edge_func = [&,this](const std::array<std::size_t,2> nodes) -> edge_factor_type& {
                if(!this->has_edge_factor(nodes[0], nodes[1]))
                    this->add_edge_factor(nodes[0], nodes[1], 0.0);
                return *(this->get_edge_factor(nodes[0], nodes[1])->get_factor());
            };

            using triplet_factor_type = typename TRIPLET_FACTOR::FactorType;
            auto get_triplet_func = [&,this](const std::array<std::size_t,3> nodes) -> triplet_factor_type& {
                if(!this->has_triplet_factor(nodes[0], nodes[1], nodes[2]))
                    this->add_triplet_factor(nodes[0], nodes[1], nodes[2]);
                return *(this->get_triplet_factor(nodes[0], nodes[1], nodes[2])->get_factor());
            };

            using msg_type = std::variant<typename EDGE_TRIPLET_MESSAGE_0::MessageType, typename EDGE_TRIPLET_MESSAGE_1::MessageType, typename EDGE_TRIPLET_MESSAGE_2::MessageType>;
            auto get_msg_func = [&](const std::array<std::size_t,2> e, const std::array<std::size_t,3> t) -> msg_type {
                if(t[0] == e[0] && t[1] == e[1]) {
                    return typename EDGE_TRIPLET_MESSAGE_0::MessageType();
                } else if(t[0] == e[0] && t[2] == e[1]) {
                    return typename EDGE_TRIPLET_MESSAGE_1::MessageType();
                } else {
                    assert(t[1] == e[0] && t[2] == e[1]);
                    return typename EDGE_TRIPLET_MESSAGE_2::MessageType();
                } 
            };

            LPMP::triangulate_cycle(cycle_begin, cycle_end, get_edge_func, get_triplet_func, get_msg_func, cycle_weight, pack_cycle);
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        template<typename INSTANCE>
        INSTANCE cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::export_edges() const
        {
            INSTANCE output;
            for(const auto& e : unary_factors_vector_)
                output.add_edge(e.first[0], e.first[1], (*e.second->get_factor())[0]);
            return output;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        template<typename INSTANCE>
        INSTANCE cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::export_triplets() const
        {
            INSTANCE output = export_edges<INSTANCE>();
            for(const auto& t : triplet_factors())
                output.add_triplet({t.first[0], t.first[1], t.first[2]}, *t.second->get_factor());
            return output;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        template<typename LABELING>
        LABELING cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::export_edge_labeling() const
        {
            LABELING output;
            output.reserve(this->no_edges());
            for(const auto& e : unary_factors_vector_) {
                output.push_back(e.second->get_factor()->primal()[0]);
            }
            return output;
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        template<typename LABELING>
        void cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::write_labeling_into_factors(const LABELING& labeling) 
        {
            std::cout << "write labeling into factors\n";
            assert(labeling.size() <= unary_factors_vector_.size());
            for(std::size_t c=0; c<labeling.size(); ++c) {
                auto* f = unary_factors_vector_[c].second;
                f->get_factor()->primal()[0] = labeling[c];
                f->propagate_primal_through_messages();
            }

            // possibly, additional edges have been added because of tightening. infer labeling of those from union find datastructure
            if(labeling.size() < unary_factors_vector_.size()) {
                union_find uf(no_nodes_);
                for(std::size_t c=0; c<labeling.size(); ++c) {
                    edge_factor* f = unary_factors_vector_[c].second; 
                    if(f->get_factor()->primal()[0] == false) {
                        // connect components 
                        const std::size_t i = std::get<0>(unary_factors_vector_[c].first);
                        const std::size_t j = std::get<1>(unary_factors_vector_[c].first);
                        uf.merge(i,j);
                    }
                }
                if(debug()) {
                    std::cout << "built union find structure, propagate information now\n";
                }
                for(std::size_t c=labeling.size(); c<unary_factors_vector_.size(); ++c) {
                    edge_factor* f = unary_factors_vector_[c].second; 
                    const std::size_t i = std::get<0>(unary_factors_vector_[c].first);
                    const std::size_t j = std::get<1>(unary_factors_vector_[c].first);
                    if(uf.connected(i,j)) {
                        f->get_factor()->primal()[0] = false;
                    } else {
                        f->get_factor()->primal()[0] = true;
                    }
                    f->propagate_primal_through_messages();
                }
            }
        }

    template<class FACTOR_MESSAGE_CONNECTION, typename EDGE_FACTOR, typename TRIPLET_FACTOR, typename EDGE_TRIPLET_MESSAGE_0, typename EDGE_TRIPLET_MESSAGE_1, typename EDGE_TRIPLET_MESSAGE_2> 
        template<typename STREAM>
        void cut_base_triplet_constructor<FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, EDGE_TRIPLET_MESSAGE_0, EDGE_TRIPLET_MESSAGE_1, EDGE_TRIPLET_MESSAGE_2>::WritePrimal(STREAM& s)
        {
            assert(no_original_edges_ <= unary_factors_vector_.size());
            for(std::size_t e=0; e<std::min(no_original_edges_, std::size_t(unary_factors_vector_.size())); ++e) {
                const std::size_t i = unary_factors_vector_[e].first[0];
                const std::size_t j = unary_factors_vector_[e].first[1];
                auto* t = unary_factors_vector_[e].second->get_factor();
                const bool cut = t->primal()[0];
                s << i << " " << j << " " << cut << "\n";
            } 
        }

} // namespace LPMP
