#pragma once

#include <lifted_disjoint_paths/ldp_single_node_cut.hxx>
#include "lifted_disjoint_paths/lifted_disjoint_paths_constructor.hxx"
#include "factors_messages.hxx"
#include "ldp_instance.hxx"

namespace LPMP {

    struct lifted_disjoint_paths_FMC {
        constexpr static const char* name = "lifted disjoint paths";

        using single_node_cut_factor_container = FactorContainer<ldp_single_node_cut_factor<lifted_disjoint_paths::LdpInstance>, lifted_disjoint_paths_FMC, 0>;

        using single_node_cut_lifted_edge_message_container = MessageContainer<ldp_snc_lifted_message, 0, 0, message_passing_schedule::left, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 0>;

        using FactorList = meta::list<single_node_cut_factor_container>;
        using MessageList = meta::list<single_node_cut_lifted_edge_message_container>;

        using problem_constructor = lifted_disjoint_paths_constructor<lifted_disjoint_paths_FMC, single_node_cut_factor_container, single_node_cut_lifted_edge_message_container>;
    };

}
