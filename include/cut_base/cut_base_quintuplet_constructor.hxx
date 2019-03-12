#pragma once 
#include <array>
#include <unordered_map>
#include <tuple>
#include <vector>

#include "send_messages_residual.hxx"

namespace LPMP {

    template<typename QUINTUPLET_FACTOR>
        class cut_base_quintuplet_storage {
            public:
                template<typename LP>
                QUINTUPLET_FACTOR* add_quintuplet_factor(const std::array<std::size_t,5> idx, LP* lp);
                bool has_quintuplet_factor(const std::array<std::size_t,5>& idx) const;
                QUINTUPLET_FACTOR* get_quintuplet_factor(const std::array<std::size_t,5>& idx) const;

                template<typename INSTANCE>
                    void export_quintuplets(INSTANCE& instance) const;

                auto& quintuplet_factors() { return quintuplet_factor_vector_; }
                const auto& quintuplet_factors() const { return quintuplet_factor_vector_; }
            private:
                std::unordered_map<std::array<std::size_t,5>, QUINTUPLET_FACTOR*> quintuplet_factors_;
                std::vector<std::pair<std::array<std::size_t,5>, QUINTUPLET_FACTOR*>> quintuplet_factor_vector_;
        };

    template<typename QUADRUPLET_CONSTRUCTOR, 
        typename QUINTUPLET_FACTOR,
        typename QUADRUPLET_QUINTUPLET_MESSAGE_0123,
        typename QUADRUPLET_QUINTUPLET_MESSAGE_0124,
        typename QUADRUPLET_QUINTUPLET_MESSAGE_0134,
        typename QUADRUPLET_QUINTUPLET_MESSAGE_0234,
        typename QUADRUPLET_QUINTUPLET_MESSAGE_1234
            >
            class cut_base_quadruplet_quintuplet_constructor : public QUADRUPLET_CONSTRUCTOR, public cut_base_quintuplet_storage<QUINTUPLET_FACTOR> {
                public:
                    using base_constructor = QUADRUPLET_CONSTRUCTOR;
                    using FMC = typename base_constructor::FMC;
                    using quintuplet_storage = cut_base_quintuplet_storage<QUINTUPLET_FACTOR>;
                    using quadruplet_factor_container = typename QUADRUPLET_CONSTRUCTOR::quadruplet_factor_container;
                    using quintuplet_factor_container = QUINTUPLET_FACTOR;
                    using quadruplet_quintuplet_message_0123_container = QUADRUPLET_QUINTUPLET_MESSAGE_0123;
                    using quadruplet_quintuplet_message_0124_container = QUADRUPLET_QUINTUPLET_MESSAGE_0124;
                    using quadruplet_quintuplet_message_0134_container = QUADRUPLET_QUINTUPLET_MESSAGE_0134;
                    using quadruplet_quintuplet_message_0234_container = QUADRUPLET_QUINTUPLET_MESSAGE_0234;
                    using quadruplet_quintuplet_message_1234_container = QUADRUPLET_QUINTUPLET_MESSAGE_1234;

                    using base_constructor::base_constructor;

                    template<typename MSG_TYPE>
                        quadruplet_factor_container* connect_quintuplet_factor(quintuplet_factor_container* f, const std::array<std::size_t,4> idx);

                    quintuplet_factor_container* add_quintuplet_factor(const std::array<std::size_t,5> idx);

                    template<typename INSTANCE>
                        INSTANCE export_quintuplets() const;

                    template<typename ITERATOR>
                        void triangulate_odd_bicycle_wheel(std::array<std::size_t,2> axle, ITERATOR path_begin, ITERATOR path_end, const double cost);

                    std::tuple<quadruplet_quintuplet_message_0123_container*, quadruplet_quintuplet_message_0124_container*, quadruplet_quintuplet_message_0134_container*, quadruplet_quintuplet_message_0234_container*, quadruplet_quintuplet_message_1234_container*>
                        get_quadruplet_quintuplet_messages(quintuplet_factor_container* f) const;
                    void send_messages_to_quadruplets(quintuplet_factor_container* f);
                    void send_messages_to_quadruplets();
            };

    template<typename TRIPLET_CONSTRUCTOR,
        typename QUINTUPLET_FACTOR,
        typename TRIPLET_QUINTUPLET_MESSAGE_012,
        typename TRIPLET_QUINTUPLET_MESSAGE_013,
        typename TRIPLET_QUINTUPLET_MESSAGE_014,
        typename TRIPLET_QUINTUPLET_MESSAGE_023,
        typename TRIPLET_QUINTUPLET_MESSAGE_024,
        typename TRIPLET_QUINTUPLET_MESSAGE_034,
        typename TRIPLET_QUINTUPLET_MESSAGE_123,
        typename TRIPLET_QUINTUPLET_MESSAGE_124,
        typename TRIPLET_QUINTUPLET_MESSAGE_134,
        typename TRIPLET_QUINTUPLET_MESSAGE_234
            >
            class cut_base_triplet_quintuplet_constructor : public TRIPLET_CONSTRUCTOR, public cut_base_quintuplet_storage<QUINTUPLET_FACTOR> {
                public:
                    using base_constructor = TRIPLET_CONSTRUCTOR;
                    using quintuplet_storage = cut_base_quintuplet_storage<QUINTUPLET_FACTOR>;
                    using triplet_factor_container = typename TRIPLET_CONSTRUCTOR::triplet_factor; //TODO: name triplet_factor_container everywhere
                    using quintuplet_factor_container = QUINTUPLET_FACTOR;

                    using base_constructor::base_constructor;

                    template<typename MSG_TYPE>
                        triplet_factor_container* connect_quintuplet_factor(quintuplet_factor_container* f, const std::array<std::size_t,3> idx);

                    quintuplet_factor_container* add_quintuplet_factor(const std::array<std::size_t,5> idx);

                    template<typename INSTANCE>
                        INSTANCE export_quintuplets() const;

                    template<typename ITERATOR>
                        void triangulate_odd_bicycle_wheel(std::array<std::size_t,2> axle, ITERATOR path_begin, ITERATOR path_end, const double cost);

                    std::tuple<TRIPLET_QUINTUPLET_MESSAGE_012*, TRIPLET_QUINTUPLET_MESSAGE_013*, TRIPLET_QUINTUPLET_MESSAGE_014*, TRIPLET_QUINTUPLET_MESSAGE_023*, TRIPLET_QUINTUPLET_MESSAGE_024*, TRIPLET_QUINTUPLET_MESSAGE_034*, TRIPLET_QUINTUPLET_MESSAGE_123*, TRIPLET_QUINTUPLET_MESSAGE_124*, TRIPLET_QUINTUPLET_MESSAGE_134*, TRIPLET_QUINTUPLET_MESSAGE_234*>
                        get_triplet_quintuplet_messages(quintuplet_factor_container* f) const;

                    void send_messages_to_triplets();
            };

    // implementation cut_base_quintuplet_storage

    template<typename QUINTUPLET_FACTOR>
        template<typename LP>
        QUINTUPLET_FACTOR* cut_base_quintuplet_storage<QUINTUPLET_FACTOR>::add_quintuplet_factor(const std::array<std::size_t,5> idx, LP* lp)
        {
            assert(!has_quintuplet_factor(idx));
            assert(quintuplet_factors_.size() == quintuplet_factor_vector_.size());
            auto* f = lp->template add_factor<QUINTUPLET_FACTOR>();
            quintuplet_factors_.insert(std::make_pair(idx, f));
            quintuplet_factor_vector_.push_back(std::make_pair(idx, f));

            return f;
        }

    template<typename QUINTUPLET_FACTOR>
        bool cut_base_quintuplet_storage<QUINTUPLET_FACTOR>::has_quintuplet_factor(const std::array<std::size_t,5>& idx) const
        {
            assert(idx[0] < idx[1] && idx[1] < idx[2] && idx[2] < idx[3] && idx[3] < idx[4]);
            return quintuplet_factors_.find(idx) != quintuplet_factors_.end(); 
        }

    template<typename QUINTUPLET_FACTOR>
        QUINTUPLET_FACTOR* cut_base_quintuplet_storage<QUINTUPLET_FACTOR>::get_quintuplet_factor(const std::array<std::size_t,5>& idx) const
        {
            assert(has_quintuplet_factor(idx));
            return quintuplet_factors_.find(idx)->second;
        }

    template<typename QUINTUPLET_FACTOR>
        template<typename INSTANCE>
        void cut_base_quintuplet_storage<QUINTUPLET_FACTOR>::export_quintuplets(INSTANCE& instance) const
        {
            for(const auto& q : quintuplet_factors())
                instance.add_quintuplet({q.first[0], q.first[1], q.first[2], q.first[3], q.first[4]}, *q.second);
        }

    // implementation cut_base_quadruplet_quintuplet_storage

    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        template<typename MSG_TYPE>
        typename QUADRUPLET_CONSTRUCTOR::quadruplet_factor_container* cut_base_quadruplet_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::connect_quintuplet_factor(QUINTUPLET_FACTOR* f, const std::array<std::size_t,4> idx)
        {
            assert(idx[0] < idx[1] && idx[1] < idx[2] && idx[2] < idx[3]);
            auto* o = this->get_quadruplet_factor(idx[0], idx[1], idx[2], idx[3]);
            auto* m = this->lp_->template add_message<MSG_TYPE>(o, f);
            return o;
        }


    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        QUINTUPLET_FACTOR* cut_base_quadruplet_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::add_quintuplet_factor(const std::array<std::size_t,5> idx)
        {
            assert(!this->has_quintuplet_factor(idx));
            auto* f = this->add_quintuplet_factor(idx, this->lp_);

            if(!this->has_quadruplet_factor(idx[0], idx[1], idx[2], idx[3]))
                this->add_quadruplet_factor(idx[0], idx[1], idx[2], idx[3]);
            if(!this->has_quadruplet_factor(idx[0], idx[1], idx[2], idx[4]))
                this->add_quadruplet_factor(idx[0], idx[1], idx[2], idx[4]);
            if(!this->has_quadruplet_factor(idx[0], idx[1], idx[3], idx[4]))
                this->add_quadruplet_factor(idx[0], idx[1], idx[3], idx[4]);
            if(!this->has_quadruplet_factor(idx[0], idx[2], idx[3], idx[4]))
                this->add_quadruplet_factor(idx[0], idx[2], idx[3], idx[4]);
            if(!this->has_quadruplet_factor(idx[1], idx[2], idx[3], idx[4]))
                this->add_quadruplet_factor(idx[1], idx[2], idx[3], idx[4]);

            auto* o_0123 = connect_quintuplet_factor<quadruplet_quintuplet_message_0123_container>(f, {idx[0], idx[1], idx[2], idx[3]});
            this->lp_->add_factor_relation(o_0123, f);
            connect_quintuplet_factor<quadruplet_quintuplet_message_0124_container>(f, {idx[0], idx[1], idx[2], idx[4]});
            connect_quintuplet_factor<quadruplet_quintuplet_message_0134_container>(f, {idx[0], idx[1], idx[3], idx[4]});
            connect_quintuplet_factor<quadruplet_quintuplet_message_0234_container>(f, {idx[0], idx[2], idx[3], idx[4]});
            auto* o_1234 = connect_quintuplet_factor<quadruplet_quintuplet_message_1234_container>(f, {idx[1], idx[2], idx[3], idx[4]});
            this->lp_->add_factor_relation(o_1234, f);

            return f;
        }

    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        template<typename INSTANCE>
        INSTANCE cut_base_quadruplet_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::export_quintuplets() const
        {
            INSTANCE output = this->export_quadruplets();
            this->export_quintuplets(output);
            for(const auto& q : this->quintuplet_factors())
                output.add_quintuplet(q.first[0], q.first[1], q.first[2], q.first[3], q.first[4], *q.second);
            return output;
        }


    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        template<typename ITERATOR>
        void cut_base_quadruplet_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::triangulate_odd_bicycle_wheel(std::array<std::size_t,2> axle, ITERATOR path_begin, ITERATOR path_end, const double cost)
        {
            using edge_factor_type = typename base_constructor::quadruplet_factor_container::FactorType;
            auto get_edge_func = [&,this](const std::array<std::size_t,2> edge_nodes) -> edge_factor_type& {
                std::array<std::size_t,4> nodes{edge_nodes[0], edge_nodes[1], axle[0], axle[1]};
                std::sort(nodes.begin(), nodes.end());
                if(!this->has_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3]))
                    this->add_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3]);
                return *(this->get_quadruplet_factor(nodes[0], nodes[1], nodes[2], nodes[3])->get_factor());
            };

            using triplet_factor_type = typename quintuplet_factor_container::FactorType;
            auto get_triplet_func = [&,this](const std::array<std::size_t,3> triplet_nodes) -> triplet_factor_type& {
                std::array<std::size_t,5> nodes{triplet_nodes[0], triplet_nodes[1], triplet_nodes[2], axle[0], axle[1]};
                std::sort(nodes.begin(), nodes.end());
                if(!this->has_quintuplet_factor(nodes))
                    static_cast<cut_base_quintuplet_storage<QUINTUPLET_FACTOR>*>(this)->add_quintuplet_factor(nodes, this->lp_);
                return *(this->get_quintuplet_factor(nodes)->get_factor());
            };

            using msg_type = std::variant<typename QUADRUPLET_QUINTUPLET_MESSAGE_0123::MessageType, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124::MessageType, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134::MessageType, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234::MessageType, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234::MessageType>;
                
            auto get_msg_func = [&](const std::array<std::size_t,2> e, const std::array<std::size_t,3> t) -> msg_type {
                std::array<std::size_t,4> e_nodes {e[0], e[1], axle[0], axle[1]};
                std::sort(e_nodes.begin(), e_nodes.end());
                std::array<std::size_t,5> t_nodes {t[0], t[1], t[2], axle[0], axle[1]};
                std::sort(t_nodes.begin(), t_nodes.end());

                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[2] && e_nodes[3] == t_nodes[3])
                    return typename QUADRUPLET_QUINTUPLET_MESSAGE_0123::MessageType();
                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[2] && e_nodes[3] == t_nodes[4])
                    return typename QUADRUPLET_QUINTUPLET_MESSAGE_0124::MessageType();
                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[3] && e_nodes[3] == t_nodes[4])
                    return typename QUADRUPLET_QUINTUPLET_MESSAGE_0134::MessageType();
                if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[3] && e_nodes[3] == t_nodes[4])
                    return typename QUADRUPLET_QUINTUPLET_MESSAGE_0234::MessageType();
                assert(e_nodes[0] == t_nodes[1] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[3] && e_nodes[3] == t_nodes[4]);
                return typename QUADRUPLET_QUINTUPLET_MESSAGE_1234::MessageType();
            };

            LPMP::triangulate_cycle(path_begin, path_end, get_edge_func, get_triplet_func, get_msg_func, cost);

        }

    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        std::tuple<QUADRUPLET_QUINTUPLET_MESSAGE_0123*, QUADRUPLET_QUINTUPLET_MESSAGE_0124*, QUADRUPLET_QUINTUPLET_MESSAGE_0134*, QUADRUPLET_QUINTUPLET_MESSAGE_0234*, QUADRUPLET_QUINTUPLET_MESSAGE_1234*>
        cut_base_quadruplet_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::get_quadruplet_quintuplet_messages(QUINTUPLET_FACTOR* f) const
        {
            auto msg_0123 = f->template get_messages<quadruplet_quintuplet_message_0123_container >();
            auto msg_0124 = f->template get_messages<quadruplet_quintuplet_message_0124_container>();
            auto msg_0134 = f->template get_messages<quadruplet_quintuplet_message_0134_container>();
            auto msg_0234 = f->template get_messages<quadruplet_quintuplet_message_0234_container>();
            auto msg_1234 = f->template get_messages<quadruplet_quintuplet_message_1234_container>();
            assert(msg_0123.size() == 1 && msg_0124.size() == 1 && msg_0134.size() == 1 && msg_0234.size() == 1 && msg_1234.size() == 1);
            return std::make_tuple(msg_0123[0], msg_0124[0], msg_0134[0], msg_0234[0], msg_1234[0]);
        }

    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        void cut_base_quadruplet_quintuplet_constructor<QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::send_messages_to_quadruplets(QUINTUPLET_FACTOR* f)
        {
            auto [msg_0123, msg_0124, msg_0134, msg_0234, msg_1234] = get_quadruplet_quintuplet_messages(f);

            msg_0123->send_message_to_right();
            msg_0124->send_message_to_right();
            msg_0134->send_message_to_right();
            msg_0234->send_message_to_right();
            msg_1234->send_message_to_right();

            msg_0123->send_message_to_left(1.0/5.0);
            msg_0124->send_message_to_left(1.0/4.0);
            msg_0134->send_message_to_left(1.0/3.0);
            msg_0234->send_message_to_left(1.0/2.0);
            msg_1234->send_message_to_left(1.0/1.0);

            msg_0123->send_message_to_left(1.0/4.0);
            msg_0124->send_message_to_left(1.0/3.0);
            msg_0134->send_message_to_left(1.0/2.0);
            msg_0234->send_message_to_left(1.0/1.0);

            msg_0123->send_message_to_left(1.0/3.0);
            msg_0124->send_message_to_left(1.0/2.0);
            msg_0134->send_message_to_left(1.0/1.0);

            msg_0123->send_message_to_left(1.0/2.0);
            msg_0124->send_message_to_left(1.0/1.0);

            msg_0123->send_message_to_left(1.0/1.0);
        }

    template<  typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234 >
        void cut_base_quadruplet_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::send_messages_to_quadruplets()
        {
            // TODO: use send_messages_to_left_residual
            for(auto f : this->quintuplet_factors())
                send_messages_to_quadruplets(f.second);
        }

    // implementation cut_base_triplet_quintuplet_storage

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename TRIPLET_QUINTUPLET_MESSAGE_012, typename TRIPLET_QUINTUPLET_MESSAGE_013, typename TRIPLET_QUINTUPLET_MESSAGE_014, typename TRIPLET_QUINTUPLET_MESSAGE_023, typename TRIPLET_QUINTUPLET_MESSAGE_024, typename TRIPLET_QUINTUPLET_MESSAGE_034, typename TRIPLET_QUINTUPLET_MESSAGE_123, typename TRIPLET_QUINTUPLET_MESSAGE_124, typename TRIPLET_QUINTUPLET_MESSAGE_134, typename TRIPLET_QUINTUPLET_MESSAGE_234>
        template<typename MSG_TYPE>
        typename TRIPLET_CONSTRUCTOR::triplet_factor* cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, TRIPLET_QUINTUPLET_MESSAGE_012, TRIPLET_QUINTUPLET_MESSAGE_013, TRIPLET_QUINTUPLET_MESSAGE_014, TRIPLET_QUINTUPLET_MESSAGE_023, TRIPLET_QUINTUPLET_MESSAGE_024, TRIPLET_QUINTUPLET_MESSAGE_034, TRIPLET_QUINTUPLET_MESSAGE_123, TRIPLET_QUINTUPLET_MESSAGE_124, TRIPLET_QUINTUPLET_MESSAGE_134, TRIPLET_QUINTUPLET_MESSAGE_234>::connect_quintuplet_factor(QUINTUPLET_FACTOR* f, const std::array<std::size_t,3> idx)
        {
            auto* t = [&]() { 
                if(!this->has_triplet_factor(idx[0], idx[1], idx[2]))
                    return this->add_triplet_factor(idx[0], idx[1], idx[2]);
                else
                    return this->get_triplet_factor(idx[0], idx[1], idx[2]);
            }();

            auto* m = this->lp_->template add_message<MSG_TYPE>(t, f);
            return t; 
        }


    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename TRIPLET_QUINTUPLET_MESSAGE_012, typename TRIPLET_QUINTUPLET_MESSAGE_013, typename TRIPLET_QUINTUPLET_MESSAGE_014, typename TRIPLET_QUINTUPLET_MESSAGE_023, typename TRIPLET_QUINTUPLET_MESSAGE_024, typename TRIPLET_QUINTUPLET_MESSAGE_034, typename TRIPLET_QUINTUPLET_MESSAGE_123, typename TRIPLET_QUINTUPLET_MESSAGE_124, typename TRIPLET_QUINTUPLET_MESSAGE_134, typename TRIPLET_QUINTUPLET_MESSAGE_234>
            QUINTUPLET_FACTOR* cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, TRIPLET_QUINTUPLET_MESSAGE_012, TRIPLET_QUINTUPLET_MESSAGE_013, TRIPLET_QUINTUPLET_MESSAGE_014, TRIPLET_QUINTUPLET_MESSAGE_023, TRIPLET_QUINTUPLET_MESSAGE_024, TRIPLET_QUINTUPLET_MESSAGE_034, TRIPLET_QUINTUPLET_MESSAGE_123, TRIPLET_QUINTUPLET_MESSAGE_124, TRIPLET_QUINTUPLET_MESSAGE_134, TRIPLET_QUINTUPLET_MESSAGE_234>::add_quintuplet_factor(const std::array<std::size_t,5> idx)
            {
                auto* f = quintuplet_storage::add_quintuplet_factor(idx, this->lp_); 
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_012>(f, {idx[0], idx[1], idx[2]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_013>(f, {idx[0], idx[1], idx[3]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_014>(f, {idx[0], idx[1], idx[4]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_023>(f, {idx[0], idx[2], idx[3]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_024>(f, {idx[0], idx[2], idx[4]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_034>(f, {idx[0], idx[3], idx[4]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_123>(f, {idx[1], idx[2], idx[3]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_124>(f, {idx[1], idx[2], idx[4]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_134>(f, {idx[1], idx[3], idx[4]});
                connect_quintuplet_factor<TRIPLET_QUINTUPLET_MESSAGE_234>(f, {idx[2], idx[3], idx[4]});
                return f;
            }

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename TRIPLET_QUINTUPLET_MESSAGE_012, typename TRIPLET_QUINTUPLET_MESSAGE_013, typename TRIPLET_QUINTUPLET_MESSAGE_014, typename TRIPLET_QUINTUPLET_MESSAGE_023, typename TRIPLET_QUINTUPLET_MESSAGE_024, typename TRIPLET_QUINTUPLET_MESSAGE_034, typename TRIPLET_QUINTUPLET_MESSAGE_123, typename TRIPLET_QUINTUPLET_MESSAGE_124, typename TRIPLET_QUINTUPLET_MESSAGE_134, typename TRIPLET_QUINTUPLET_MESSAGE_234>
                    template<typename INSTANCE>
            INSTANCE cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, TRIPLET_QUINTUPLET_MESSAGE_012, TRIPLET_QUINTUPLET_MESSAGE_013, TRIPLET_QUINTUPLET_MESSAGE_014, TRIPLET_QUINTUPLET_MESSAGE_023, TRIPLET_QUINTUPLET_MESSAGE_024, TRIPLET_QUINTUPLET_MESSAGE_034, TRIPLET_QUINTUPLET_MESSAGE_123, TRIPLET_QUINTUPLET_MESSAGE_124, TRIPLET_QUINTUPLET_MESSAGE_134, TRIPLET_QUINTUPLET_MESSAGE_234>::export_quintuplets() const
            {
                auto output = this->template export_triplets<INSTANCE>();
                this->export_quintuplets(output);
                return output;
            }

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename TRIPLET_QUINTUPLET_MESSAGE_012, typename TRIPLET_QUINTUPLET_MESSAGE_013, typename TRIPLET_QUINTUPLET_MESSAGE_014, typename TRIPLET_QUINTUPLET_MESSAGE_023, typename TRIPLET_QUINTUPLET_MESSAGE_024, typename TRIPLET_QUINTUPLET_MESSAGE_034, typename TRIPLET_QUINTUPLET_MESSAGE_123, typename TRIPLET_QUINTUPLET_MESSAGE_124, typename TRIPLET_QUINTUPLET_MESSAGE_134, typename TRIPLET_QUINTUPLET_MESSAGE_234>
                    template<typename ITERATOR>
            void cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, TRIPLET_QUINTUPLET_MESSAGE_012, TRIPLET_QUINTUPLET_MESSAGE_013, TRIPLET_QUINTUPLET_MESSAGE_014, TRIPLET_QUINTUPLET_MESSAGE_023, TRIPLET_QUINTUPLET_MESSAGE_024, TRIPLET_QUINTUPLET_MESSAGE_034, TRIPLET_QUINTUPLET_MESSAGE_123, TRIPLET_QUINTUPLET_MESSAGE_124, TRIPLET_QUINTUPLET_MESSAGE_134, TRIPLET_QUINTUPLET_MESSAGE_234>::triangulate_odd_bicycle_wheel(std::array<std::size_t,2> axle, ITERATOR path_begin, ITERATOR path_end, const double cost)
            {
                using edge_factor_type = typename std::array<typename triplet_factor_container::FactorType*, 2>;
                auto get_edge_func = [&,this](const std::array<std::size_t,2> edge_nodes) -> edge_factor_type {
                    std::array<std::size_t,3> nodes_1{edge_nodes[0], edge_nodes[1], axle[0]};
                    std::sort(nodes_1.begin(), nodes_1.end());
                    if(!this->has_triplet_factor(nodes_1[0], nodes_1[1], nodes_1[2]))
                        this->add_triplet_factor(nodes_1[0], nodes_1[1], nodes_1[2]);
                    auto* f1 = this->get_triplet_factor(nodes_1[0], nodes_1[1], nodes_1[2])->get_factor();

                    std::array<std::size_t,3> nodes_2{edge_nodes[0], edge_nodes[1], axle[1]};
                    std::sort(nodes_2.begin(), nodes_2.end());
                    if(!this->has_triplet_factor(nodes_2[0], nodes_2[1], nodes_2[2]))
                        this->add_triplet_factor(nodes_2[0], nodes_2[1], nodes_2[2]);
                    auto* f2 = this->get_triplet_factor(nodes_2[0], nodes_2[1], nodes_2[2])->get_factor();

                    return {f1,f2};
                }; 

                using triplet_factor_type = typename quintuplet_factor_container::FactorType;
                auto get_triplet_func = [&,this](const std::array<std::size_t,3> triplet_nodes) -> triplet_factor_type& {
                    std::array<std::size_t,5> nodes{triplet_nodes[0], triplet_nodes[1], triplet_nodes[2], axle[0], axle[1]};
                    std::sort(nodes.begin(), nodes.end());
                    if(!this->has_quintuplet_factor(nodes))
                        this->add_quintuplet_factor(nodes);
                    return *(this->get_quintuplet_factor(nodes)->get_factor());
                };

                using single_msg_type = std::variant<typename TRIPLET_QUINTUPLET_MESSAGE_012::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_013::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_014::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_023::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_024::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_034::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_123::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_124::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_134::MessageType,
                      typename TRIPLET_QUINTUPLET_MESSAGE_234::MessageType
                    >;

                auto get_single_msg_func = [&](const std::array<std::size_t,2> e, const std::array<std::size_t,3> t, const std::size_t axle_node) -> single_msg_type {
                    std::array<std::size_t,3> e_nodes {e[0], e[1], axle_node};
                    std::sort(e_nodes.begin(), e_nodes.end());
                    std::array<std::size_t,5> t_nodes {t[0], t[1], t[2], axle[0], axle[1]};
                    std::sort(t_nodes.begin(), t_nodes.end());

                    if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[2])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_012::MessageType();

                    if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[3])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_013::MessageType();
                    
                    if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[1] && e_nodes[2] == t_nodes[4])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_014::MessageType();

                    if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[3])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_023::MessageType();
                    
                    if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[4])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_024::MessageType();
                    
                    if(e_nodes[0] == t_nodes[0] && e_nodes[1] == t_nodes[3] && e_nodes[2] == t_nodes[4])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_034::MessageType();

                    if(e_nodes[0] == t_nodes[1] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[3])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_123::MessageType();

                    if(e_nodes[0] == t_nodes[1] && e_nodes[1] == t_nodes[2] && e_nodes[2] == t_nodes[4])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_124::MessageType();

                    if(e_nodes[0] == t_nodes[1] && e_nodes[1] == t_nodes[3] && e_nodes[2] == t_nodes[4])
                        return typename TRIPLET_QUINTUPLET_MESSAGE_134::MessageType();

                    assert(e_nodes[0] == t_nodes[2] && e_nodes[1] == t_nodes[3] && e_nodes[2] == t_nodes[4]);
                    return typename TRIPLET_QUINTUPLET_MESSAGE_234::MessageType();
                };

                struct msg_type_impl : public  std::array<single_msg_type,2> {
                    using msg_val_type = std::array<array<double,triplet_factor_container::FactorType::size()>,2>;

                    void RepamRight(triplet_factor_type& t, const msg_val_type& msg_val)
                    { 
                        std::visit([&](auto& msg) { msg.RepamRight(t, msg_val[0]); }, (*this)[0]);
                        std::visit([&](auto& msg) { msg.RepamRight(t, msg_val[1]); }, (*this)[1]);
                    }

                    void RepamLeft(edge_factor_type e, const msg_val_type& msg_val)
                    {
                        std::visit([&](auto& msg) { msg.RepamLeft(*e[0], msg_val[0]); }, (*this)[0]);
                        std::visit([&](auto& msg) { msg.RepamLeft(*e[1], msg_val[1]); }, (*this)[1]);
                    }

                    void send_message_to_left(triplet_factor_type& t, msg_val_type& msg_val, const double scaling)
                    {
                        std::visit([&](auto& msg) { msg.send_message_to_left(t, msg_val[0], 0.5*scaling); }, (*this)[0]);
                        std::visit([&](auto& msg) { msg.send_message_to_left(t, msg_val[1], 0.5*scaling); }, (*this)[1]);
                    }

                    void send_message_to_right(edge_factor_type e, msg_val_type& msg_val, const double scaling)
                    {
                        std::visit([&](auto& msg) { msg.send_message_to_right(*e[0], msg_val[0], scaling); }, (*this)[0]);
                        std::visit([&](auto& msg) { msg.send_message_to_right(*e[1], msg_val[1], scaling); }, (*this)[1]);
                    }
                };

                using msg_type = std::variant<msg_type_impl>; // This is not nice, but needed fpr the triangulation code

                auto get_msg_func = [&](const std::array<std::size_t,2> e, const std::array<std::size_t,3> t) -> msg_type {
                    auto m1 = get_single_msg_func(e, t, axle[0]);
                    auto m2 = get_single_msg_func(e, t, axle[1]);
                    return msg_type_impl{m1,m2}; 
                };

                LPMP::triangulate_cycle(path_begin, path_end, get_edge_func, get_triplet_func, get_msg_func, cost);
            }

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename TRIPLET_QUINTUPLET_MESSAGE_012, typename TRIPLET_QUINTUPLET_MESSAGE_013, typename TRIPLET_QUINTUPLET_MESSAGE_014, typename TRIPLET_QUINTUPLET_MESSAGE_023, typename TRIPLET_QUINTUPLET_MESSAGE_024, typename TRIPLET_QUINTUPLET_MESSAGE_034, typename TRIPLET_QUINTUPLET_MESSAGE_123, typename TRIPLET_QUINTUPLET_MESSAGE_124, typename TRIPLET_QUINTUPLET_MESSAGE_134, typename TRIPLET_QUINTUPLET_MESSAGE_234>
            std::tuple<TRIPLET_QUINTUPLET_MESSAGE_012*, TRIPLET_QUINTUPLET_MESSAGE_013*, TRIPLET_QUINTUPLET_MESSAGE_014*, TRIPLET_QUINTUPLET_MESSAGE_023*, TRIPLET_QUINTUPLET_MESSAGE_024*, TRIPLET_QUINTUPLET_MESSAGE_034*, TRIPLET_QUINTUPLET_MESSAGE_123*, TRIPLET_QUINTUPLET_MESSAGE_124*, TRIPLET_QUINTUPLET_MESSAGE_134*, TRIPLET_QUINTUPLET_MESSAGE_234*>
            cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, TRIPLET_QUINTUPLET_MESSAGE_012, TRIPLET_QUINTUPLET_MESSAGE_013, TRIPLET_QUINTUPLET_MESSAGE_014, TRIPLET_QUINTUPLET_MESSAGE_023, TRIPLET_QUINTUPLET_MESSAGE_024, TRIPLET_QUINTUPLET_MESSAGE_034, TRIPLET_QUINTUPLET_MESSAGE_123, TRIPLET_QUINTUPLET_MESSAGE_124, TRIPLET_QUINTUPLET_MESSAGE_134, TRIPLET_QUINTUPLET_MESSAGE_234>::get_triplet_quintuplet_messages(quintuplet_factor_container* f) const
            {
                auto msg_012 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_012>();
                auto msg_013 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_013>();
                auto msg_014 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_014>();
                auto msg_023 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_023>();
                auto msg_024 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_024>();
                auto msg_034 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_034>();
                auto msg_123 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_123>();
                auto msg_124 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_124>();
                auto msg_134 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_134>();
                auto msg_234 = f->template get_messages<TRIPLET_QUINTUPLET_MESSAGE_234>();
                assert(msg_012.size() == 1);
                assert(msg_013.size() == 1);
                assert(msg_014.size() == 1);
                assert(msg_023.size() == 1);
                assert(msg_024.size() == 1);
                assert(msg_034.size() == 1);
                assert(msg_123.size() == 1);
                assert(msg_124.size() == 1);
                assert(msg_134.size() == 1);
                assert(msg_234.size() == 1);
                return std::make_tuple(msg_012[0], msg_013[0], msg_014[0], msg_023[0], msg_024[0], msg_034[0], msg_123[0], msg_124[0], msg_134[0], msg_234[0]); 
            }

    template<typename TRIPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename TRIPLET_QUINTUPLET_MESSAGE_012, typename TRIPLET_QUINTUPLET_MESSAGE_013, typename TRIPLET_QUINTUPLET_MESSAGE_014, typename TRIPLET_QUINTUPLET_MESSAGE_023, typename TRIPLET_QUINTUPLET_MESSAGE_024, typename TRIPLET_QUINTUPLET_MESSAGE_034, typename TRIPLET_QUINTUPLET_MESSAGE_123, typename TRIPLET_QUINTUPLET_MESSAGE_124, typename TRIPLET_QUINTUPLET_MESSAGE_134, typename TRIPLET_QUINTUPLET_MESSAGE_234>
            void cut_base_triplet_quintuplet_constructor<TRIPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, TRIPLET_QUINTUPLET_MESSAGE_012, TRIPLET_QUINTUPLET_MESSAGE_013, TRIPLET_QUINTUPLET_MESSAGE_014, TRIPLET_QUINTUPLET_MESSAGE_023, TRIPLET_QUINTUPLET_MESSAGE_024, TRIPLET_QUINTUPLET_MESSAGE_034, TRIPLET_QUINTUPLET_MESSAGE_123, TRIPLET_QUINTUPLET_MESSAGE_124, TRIPLET_QUINTUPLET_MESSAGE_134, TRIPLET_QUINTUPLET_MESSAGE_234>::send_messages_to_triplets()
            {
                for(const auto& f : this->unary_factors_vector_) {
                    f.second->template send_messages_uniform<typename TRIPLET_CONSTRUCTOR::edge_triplet_message_0, typename TRIPLET_CONSTRUCTOR::edge_triplet_message_1, typename TRIPLET_CONSTRUCTOR::edge_triplet_message_2>();
                }

                for(const auto& f : this->quintuplet_factors()) {
                    auto msgs = get_triplet_quintuplet_messages(f.second); 
                    send_messages_to_left_residual(msgs); 
                }
            } 

} // namespace LPMP
