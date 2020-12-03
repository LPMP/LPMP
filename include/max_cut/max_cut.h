#pragma once

#include "factors_messages.hxx"
#include "max_cut_factors_messages.h"
#include "max_cut_triplet_constructor.hxx"
#include "max_cut_quintuplet_constructor.hxx"
#include "solver.hxx"
#include "LP.h"

namespace LPMP {

struct FMC_MAX_CUT {
   constexpr static const char* name = "max-cut with cycle constraints";

   using edge_factor_container = FactorContainer<max_cut_edge_factor, FMC_MAX_CUT, 0, true>;
   using triplet_factor_container = FactorContainer<max_cut_triplet_factor, FMC_MAX_CUT, 1>;

   using edge_triplet_message_0_container = MessageContainer<max_cut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MAX_CUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<max_cut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MAX_CUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<max_cut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_MAX_CUT, 2 >;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container>;
   using MessageList = meta::list<edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;

   using max_cut_c = max_cut_triplet_constructor<FMC_MAX_CUT,edge_factor_container,triplet_factor_container,edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;
   using problem_constructor = max_cut_c;
};

struct FMC_ODD_BICYCLE_WHEEL_MAX_CUT {
   constexpr static const char* name = "max-cut with cycle and odd bicycle wheel constraints";

   using edge_factor_container = FactorContainer<max_cut_edge_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 0>;
   using triplet_factor_container = FactorContainer<max_cut_triplet_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 1>;
   using quintuplet_factor_container = FactorContainer<max_cut_quintuplet_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 2>;
      
   using edge_triplet_message_0_container = MessageContainer<max_cut_edge_triplet_message_0, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<max_cut_edge_triplet_message_1, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<max_cut_edge_triplet_message_2, 0, 1, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 2 >;

   using triplet_quintuplet__message_012_container = MessageContainer<max_cut_triplet_quintuplet_message_012, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 3 >;
   using triplet_quintuplet__message_013_container = MessageContainer<max_cut_triplet_quintuplet_message_013, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 4 >;
   using triplet_quintuplet__message_014_container = MessageContainer<max_cut_triplet_quintuplet_message_014, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 5 >;
   using triplet_quintuplet__message_023_container = MessageContainer<max_cut_triplet_quintuplet_message_023, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 6 >;
   using triplet_quintuplet__message_024_container = MessageContainer<max_cut_triplet_quintuplet_message_024, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 7 >;
   using triplet_quintuplet__message_034_container = MessageContainer<max_cut_triplet_quintuplet_message_034, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 8 >;
   using triplet_quintuplet__message_123_container = MessageContainer<max_cut_triplet_quintuplet_message_123, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 9 >;
   using triplet_quintuplet__message_124_container = MessageContainer<max_cut_triplet_quintuplet_message_124, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 10 >;
   using triplet_quintuplet__message_134_container = MessageContainer<max_cut_triplet_quintuplet_message_134, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 11 >;
   using triplet_quintuplet__message_234_container = MessageContainer<max_cut_triplet_quintuplet_message_234, 1, 2, message_passing_schedule::left, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 12 >;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container, quintuplet_factor_container>;
   using MessageList = meta::list<
      edge_triplet_message_0_container,
      edge_triplet_message_1_container,
      edge_triplet_message_2_container,  

      triplet_quintuplet__message_012_container,
      triplet_quintuplet__message_013_container,
      triplet_quintuplet__message_014_container,
      triplet_quintuplet__message_023_container,
      triplet_quintuplet__message_024_container,
      triplet_quintuplet__message_034_container,
      triplet_quintuplet__message_123_container,
      triplet_quintuplet__message_124_container,
      triplet_quintuplet__message_134_container,
      triplet_quintuplet__message_234_container
      >;

   using max_cut_c = max_cut_triplet_constructor<FMC_ODD_BICYCLE_WHEEL_MAX_CUT,edge_factor_container,triplet_factor_container,edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;
   using max_cut_cobw = max_cut_quintuplet_constructor<max_cut_c, quintuplet_factor_container, triplet_quintuplet__message_012_container, triplet_quintuplet__message_013_container, triplet_quintuplet__message_014_container, triplet_quintuplet__message_023_container, triplet_quintuplet__message_024_container, triplet_quintuplet__message_034_container, triplet_quintuplet__message_123_container, triplet_quintuplet__message_124_container, triplet_quintuplet__message_134_container, triplet_quintuplet__message_234_container>;
   using problem_constructor = max_cut_cobw;
};

} // namespace LPMP
