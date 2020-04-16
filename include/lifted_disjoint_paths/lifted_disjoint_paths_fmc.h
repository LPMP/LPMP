#pragma once

#include <lifted_disjoint_paths/ldp_single_node_cut.hxx>
#include "lifted_disjoint_paths_mcf.hxx"
#include "lifted_disjoint_paths/lifted_disjoint_paths_constructor.hxx"
#include "factors_messages.hxx"
#include "ldp_instance.hxx"

namespace LPMP {

    struct lifted_disjoint_paths_FMC {
        constexpr static const char* name = "lifted disjoint paths";

        using mcf_factor_container = FactorContainer<lifted_disjoint_paths_mcf_factor, lifted_disjoint_paths_FMC, 0>;
        using single_node_cut_factor_container = FactorContainer<ldp_single_node_cut_factor<lifted_disjoint_paths::LdpInstance>, lifted_disjoint_paths_FMC, 1>;

        using mcf_single_node_cut_message_container = MessageContainer<ldp_mcf_single_node_cut_message, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, lifted_disjoint_paths_FMC, 0>;

        using FactorList = meta::list<mcf_factor_container, single_node_cut_factor_container>;
        using MessageList = meta::list<mcf_single_node_cut_message_container>;

        using problem_constructor = lifted_disjoint_paths_constructor<lifted_disjoint_paths_FMC, mcf_factor_container, single_node_cut_factor_container, mcf_single_node_cut_message_container>;
    };

}
