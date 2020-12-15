#pragma once

#include <lifted_disjoint_paths/ldp_single_node_cut.hxx>
//#include <lifted_disjoint_paths/ldp_triangle_factor.hxx>
#include <lifted_disjoint_paths/ldp_cut_factor.hxx>
#include "lifted_disjoint_paths/lifted_disjoint_paths_constructor.hxx"
#include "factors_messages.hxx"
#include "ldp_instance.hxx"
#include "ldp_path_factor.hxx"

namespace LPMP {

    struct lifted_disjoint_paths_FMC {
        constexpr static const char* name = "lifted disjoint paths";

        using single_node_cut_factor_container = FactorContainer<ldp_single_node_cut_factor<lifted_disjoint_paths::LdpInstance>, lifted_disjoint_paths_FMC, 0>;

      //  using triangle_factor_container = FactorContainer<ldp_triangle_factor, lifted_disjoint_paths_FMC, 1>;
          using cut_factor_container = FactorContainer<ldp_cut_factor, lifted_disjoint_paths_FMC, 1>;

          using path_factor_container = FactorContainer<ldp_path_factor, lifted_disjoint_paths_FMC, 2>;


        using single_node_cut_lifted_edge_message_container = MessageContainer<ldp_snc_lifted_message, 0, 0, message_passing_schedule::only_send, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 0>;

       // using single_node_cut_node_message_container = MessageContainer<ldp_snc_node_message, 0, 0, message_passing_schedule::only_send, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 1>;

        // TODO: replace with atMostFourMessages
       // using snc_triangle_message_container=MessageContainer<ldp_snc_triangle_message, 1, 0, message_passing_schedule::only_send, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 1>;
         using snc_cut_message_container=MessageContainer<ldp_snc_cut_message, 1, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 1>;


         using snc_path_message_container=MessageContainer<ldp_snc_path_message, 2, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 2>;
     //    using snc_triangel_message_container=MessageContainer<ldp_snc_triangle_message, 1, 0, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 1>;

        //using FactorList = meta::list<single_node_cut_factor_container,cut_factor_container>;
        using FactorList = meta::list<single_node_cut_factor_container,cut_factor_container,path_factor_container>;
        //using MessageList = meta::list<single_node_cut_lifted_edge_message_container,snc_triangle_message_container>;
        using MessageList = meta::list<single_node_cut_lifted_edge_message_container,snc_cut_message_container,snc_path_message_container>;

        //using problem_constructor = lifted_disjoint_paths_constructor<lifted_disjoint_paths_FMC, single_node_cut_factor_container,cut_factor_container, single_node_cut_lifted_edge_message_container,snc_triangle_message_container>;
        using problem_constructor = lifted_disjoint_paths_constructor<lifted_disjoint_paths_FMC, single_node_cut_factor_container,cut_factor_container, single_node_cut_lifted_edge_message_container,snc_cut_message_container,path_factor_container,snc_path_message_container>;
    };

}


//#pragma once

//#include <lifted_disjoint_paths/ldp_single_node_cut.hxx>
////#include <lifted_disjoint_paths/ldp_triangle_factor.hxx>
//#include "lifted_disjoint_paths/lifted_disjoint_paths_constructor.hxx"
//#include "factors_messages.hxx"
//#include "ldp_instance.hxx"

//namespace LPMP {

//    struct lifted_disjoint_paths_FMC {
//        constexpr static const char* name = "lifted disjoint paths";

//        using single_node_cut_factor_container = FactorContainer<ldp_single_node_cut_factor<lifted_disjoint_paths::LdpInstance>, lifted_disjoint_paths_FMC, 0>;
//      //  using triangle_factor_container = FactorContainer<ldp_triangle_factor, lifted_disjoint_paths_FMC, 1>;

//        using single_node_cut_lifted_edge_message_container = MessageContainer<ldp_snc_lifted_message, 0, 0, message_passing_schedule::only_send, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 0>;

//        using single_node_cut_node_message_container = MessageContainer<ldp_snc_node_message, 0, 0, message_passing_schedule::only_send, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 1>;

//        //TODO implement class ldp_snc_lifted_triangle_message
//    //    using snc_lifted_edge_to_triangel_message_container=MessageContainer<ldp_snc_lifted_triangle_message, 0, 1, message_passing_schedule::full, variableMessageNumber, variableMessageNumber, lifted_disjoint_paths_FMC, 2>;

////        using FactorList = meta::list<single_node_cut_factor_container, triangle_factor_container>;
//        //using MessageList = meta::list<single_node_cut_lifted_edge_message_container,single_node_cut_node_message_container>;
//        using FactorList = meta::list<single_node_cut_factor_container>;
//        using MessageList = meta::list<single_node_cut_lifted_edge_message_container,single_node_cut_node_message_container>;

//        using problem_constructor = lifted_disjoint_paths_constructor<lifted_disjoint_paths_FMC, single_node_cut_factor_container, single_node_cut_lifted_edge_message_container,single_node_cut_node_message_container>;
//    };

//}
