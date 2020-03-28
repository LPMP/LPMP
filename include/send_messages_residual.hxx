#pragma once

#include <functional>

namespace LPMP {

    template<int N, int M, typename MSG_TUPLE>
        void send_message_to_left_residual_impl(const MSG_TUPLE& msgs)
        {
            static_assert(M < N);
            static_assert(N > 0);
            auto* msg = std::get<M>(msgs);
            msg->send_message_to_left(double(M+1)/double(N));

            if constexpr(M+1 < N)
                send_message_to_left_residual_impl<N,M+1>(msgs);
        }
    template<int N, typename MSG_TUPLE>
        void send_messages_to_left_residual_impl(const MSG_TUPLE& msgs)
        {
            send_message_to_left_residual_impl<N,0>(msgs);
            if constexpr(N > 1)
                send_messages_to_left_residual_impl<N-1>(msgs);

        }
    template<typename MSG_TUPLE>
        void send_messages_to_left_residual(const MSG_TUPLE& msgs)
        {
            static_assert(std::tuple_size_v<MSG_TUPLE> > 0);
            std::apply([](auto ...m){(..., m->send_message_to_right(1.0)); }, msgs);
            //std::get<0>(msgs)->GetRightFactor()->receive_all_messages();
            constexpr int N = std::tuple_size_v<MSG_TUPLE>;
            send_messages_to_left_residual_impl<N>(msgs);
        }

} // namespace LPMP
