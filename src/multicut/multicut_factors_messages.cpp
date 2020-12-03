#include "multicut/multicut_factors_messages.h"

namespace LPMP {

template class labeling_factor< multicut_edge_labelings, true >;
template class labeling_factor< multicut_triplet_labelings, true >;

template class labeling_message< multicut_edge_labelings, multicut_triplet_labelings, 0 >;
template class labeling_message< multicut_edge_labelings, multicut_triplet_labelings, 1 >;
template class labeling_message< multicut_edge_labelings, multicut_triplet_labelings, 2 >;

template class labeling_factor< multicut_odd_3_wheel_labelings, true >;

template class labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 0,1,2>;
template class labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 0,3,4>;
template class labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 1,3,5>;
template class labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 2,4,5>;

template class labeling_factor< multicut_odd_bicycle_3_wheel_labelings, true >;

template class labeling_message< multicut_odd_3_wheel_labelings, multicut_odd_bicycle_3_wheel_labelings, 0,1,2,3,5,7 >; // 01->01, 02->02, 12->12, 03->03, 13->13, 23->23
template class labeling_message< multicut_odd_3_wheel_labelings, multicut_odd_bicycle_3_wheel_labelings, 0,1,2,4,6,8 >; // 01->01, 02->02, 12->12, 03->04, 13->14, 23->24
template class labeling_message< multicut_odd_3_wheel_labelings, multicut_odd_bicycle_3_wheel_labelings, 0,3,5,4,6,9 >; // 01->01, 02->03, 13->13, 03->04, 13->14, 23->34
template class labeling_message< multicut_odd_3_wheel_labelings, multicut_odd_bicycle_3_wheel_labelings, 1,3,7,4,8,9 >; // 01->02, 02->03, 12->23, 03->04, 13->24, 23->34
template class labeling_message< multicut_odd_3_wheel_labelings, multicut_odd_bicycle_3_wheel_labelings, 2,5,7,6,8,9 >; // 01->12, 02->13, 12->23, 03->14, 13->24, 23->34

} // end namespace LPMP 
