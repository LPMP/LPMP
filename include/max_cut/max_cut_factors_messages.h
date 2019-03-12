#pragma once

#include "factors/labeling_list_factor.hxx"

namespace LPMP {

using max_cut_edge_labelings = labelings< labeling<1> >;
using max_cut_edge_factor = labeling_factor< max_cut_edge_labelings, true >;

using max_cut_triplet_labelings = labelings<
   labeling<1,1,0>, // 100
   labeling<1,0,1>, // 010
   labeling<0,1,1>  // 001
      >;
using max_cut_triplet_factor = labeling_factor< max_cut_triplet_labelings, true >;

using max_cut_edge_triplet_message_0 = labeling_message< max_cut_edge_labelings, max_cut_triplet_labelings, 0 >;
using max_cut_edge_triplet_message_1 = labeling_message< max_cut_edge_labelings, max_cut_triplet_labelings, 1 >;
using max_cut_edge_triplet_message_2 = labeling_message< max_cut_edge_labelings, max_cut_triplet_labelings, 2 >;

/*

   ______i1_____
  /      |      \
 /       |       \
/   _____i3____   \
|  /           \  |
| /             \ |
i0_______________i2

   */

// edges 01,02,12,03,13,23
using max_cut_quadruplet_labelings = labelings<
   // two components
   // one node separated from the rest
   labeling<1,1,0,1,0,0>, // 1000
   labeling<1,0,1,0,1,0>, // 0100
   labeling<0,1,1,0,0,1>, // 0010
   labeling<0,0,0,1,1,1>, // 0001

   // two components of two nodes each
   labeling<0,1,1,1,1,0>, // 1100 // 4
   labeling<1,0,1,1,0,1>, // 0101 // 5
   labeling<1,1,0,0,1,1> // 1001 // 6
   >;
using max_cut_quadruplet_factor = labeling_factor< max_cut_quadruplet_labelings, true >;

using max_cut_triplet_quadruplet_message_012 = labeling_message< max_cut_triplet_labelings, max_cut_quadruplet_labelings, 0,1,2>;
using max_cut_triplet_quadruplet_message_013 = labeling_message< max_cut_triplet_labelings, max_cut_quadruplet_labelings, 0,3,4>;
using max_cut_triplet_quadruplet_message_023 = labeling_message< max_cut_triplet_labelings, max_cut_quadruplet_labelings, 1,3,5>;
using max_cut_triplet_quadruplet_message_123 = labeling_message< max_cut_triplet_labelings, max_cut_quadruplet_labelings, 2,4,5>;

using max_cut_edge_quadruplet_message_0 = labeling_message< max_cut_edge_labelings, max_cut_quadruplet_labelings, 0>;
using max_cut_edge_quadruplet_message_1 = labeling_message< max_cut_edge_labelings, max_cut_quadruplet_labelings, 1>;
using max_cut_edge_quadruplet_message_2 = labeling_message< max_cut_edge_labelings, max_cut_quadruplet_labelings, 2>;
using max_cut_edge_quadruplet_message_3 = labeling_message< max_cut_edge_labelings, max_cut_quadruplet_labelings, 3>;
using max_cut_edge_quadruplet_message_4 = labeling_message< max_cut_edge_labelings, max_cut_quadruplet_labelings, 4>;
using max_cut_edge_quadruplet_message_5 = labeling_message< max_cut_edge_labelings, max_cut_quadruplet_labelings, 5>;


/*
                       i1



               
                    i3_
                       \__
                          i4




   i0                                       i2
   */
// is also the full graph on 5 nodes
// edges: 01,02,12, 03,04,13,14,23,24, 34
using max_cut_quintuplet_labelings = labelings<
// two components
// one node separated from the rest
labeling<1,1,0,1,1,0,0,0,0,0>, // 10000
labeling<1,0,1,0,0,1,1,0,0,0>, // 01000
labeling<0,1,1,0,0,0,0,1,1,0>, // 00100
labeling<0,0,0,1,0,1,0,1,0,1>, // 00010
labeling<0,0,0,0,1,0,1,0,1,1>, // 00001
// two nodes, three nodes
labeling<0,1,1,1,1,1,1,0,0,0>, // 11000
labeling<1,0,1,1,1,0,0,1,1,0>, // 10100
labeling<1,1,0,0,1,1,0,1,0,1>, // 10010
labeling<1,1,0,1,0,0,1,0,1,1>, // 10001
labeling<1,1,0,0,0,1,1,1,1,0>, // 01100
labeling<1,0,1,1,0,0,1,1,0,1>, // 01010
labeling<1,0,1,0,1,1,0,0,1,1>, // 01001
labeling<0,1,1,1,0,1,0,0,1,1>, // 00110
labeling<0,1,1,0,1,0,1,1,0,1>, // 00101
labeling<0,0,0,1,1,1,1,1,1,0>  // 00011
   >;
using max_cut_quintuplet_factor = labeling_factor< max_cut_quintuplet_labelings, true >;

// edges 01,02,12, 03,04,13,14,23,24, 34
//        0  1  2   3  4  5  6  7  8   9
using max_cut_quadruplet_quintuplet_message_0123 = labeling_message< max_cut_quadruplet_labelings, max_cut_quintuplet_labelings, 0,1,2,3,5,7 >; // 01->01, 02->02, 12->12, 03->03, 13->13, 23->23
using max_cut_quadruplet_quintuplet_message_0124 = labeling_message< max_cut_quadruplet_labelings, max_cut_quintuplet_labelings, 0,1,2,4,6,8 >; // 01->01, 02->02, 12->12, 03->04, 13->14, 23->24
using max_cut_quadruplet_quintuplet_message_0134 = labeling_message< max_cut_quadruplet_labelings, max_cut_quintuplet_labelings, 0,3,5,4,6,9 >; // 01->01, 02->03, 13->13, 03->04, 13->14, 23->34
using max_cut_quadruplet_quintuplet_message_0234 = labeling_message< max_cut_quadruplet_labelings, max_cut_quintuplet_labelings, 1,3,7,4,8,9 >; // 01->02, 02->03, 12->23, 03->04, 13->24, 23->34
using max_cut_quadruplet_quintuplet_message_1234 = labeling_message< max_cut_quadruplet_labelings, max_cut_quintuplet_labelings, 2,5,7,6,8,9 >; // 01->12, 02->13, 12->23, 03->14, 13->24, 23->34


// edges 01,02,12, 03,04,13,14,23,24, 34
//        0  1  2   3  4  5  6  7  8   9
using max_cut_triplet_quintuplet_message_012 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 0,1,2>; // 01->01, 02->02, 12->12
using max_cut_triplet_quintuplet_message_013 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 0,3,5>; // 01->01, 02->03, 12->03
using max_cut_triplet_quintuplet_message_014 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 0,4,6>; // 01->01, 02->04, 12->14
using max_cut_triplet_quintuplet_message_023 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 1,3,7>; // 01->02, 02->03, 12->23
using max_cut_triplet_quintuplet_message_024 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 1,4,8>; // 01->02, 02->04, 12->24
using max_cut_triplet_quintuplet_message_034 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 3,4,9>; // 01->01, 02->02, 12->12
using max_cut_triplet_quintuplet_message_123 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 2,5,7>; // 01->12, 02->13, 12->23
using max_cut_triplet_quintuplet_message_124 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 2,6,8>; // 01->12, 02->14, 12->24
using max_cut_triplet_quintuplet_message_134 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 5,6,9>; // 01->13, 02->14, 12->34
using max_cut_triplet_quintuplet_message_234 = labeling_message< max_cut_triplet_labelings, max_cut_quintuplet_labelings, 7,8,9>; // 01->23, 02->24, 12->34 

} // end namespace LPMP 
