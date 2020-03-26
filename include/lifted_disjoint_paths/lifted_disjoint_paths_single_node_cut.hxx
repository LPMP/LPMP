#pragma once

namespace LPMP {

    class ldp_single_node_cut_factor
    {
        public:
            constexpr static std::size_t no_edge_active = std::numeric_limits<std::size_t>::infinity();

            ldp_single_node_cut_factor(const std::size_t nr_outgoing_base_edges, const std::size_t nr_outgoing_lifted_edges);
            double LowerBound() const;
            double EvaluatePrimal() const;

            template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
            template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

            auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }

            void init_primal() { primal_ = no_edge_active; }

        private:
            std::size_t primal_; // the incoming resp. outgoing edge that is active.
    };

    class ldp_mcf_single_node_cut_message
    {
        template<typename SINGLE_NODE_CUT_FACTOR>
            void RepamLeft(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
            {
            } 

        template<typename MCF_FACTOR>
            void RepamRight(MCF_FACTOR& r, const double msg, const std::size_t msg_dim) const
            {
            } 

         template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
             void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
             {
             }

         template<typename MCF_FACTOR, typename MSG_ARRAY> 
             static void SendMessagesToRight(const MCF_FACTOR& leftRepam, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
             {
             } 
    };

}
