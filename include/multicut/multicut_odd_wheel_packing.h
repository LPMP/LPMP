#pragma once

#include "multicut/multicut_odd_wheel_packing.h"
#include "multicut/multicut_cycle_packing.h"
#include "multicut_instance.h"
#include "cut_base/cut_base_packing.h"
#include <vector>
#include <numeric>

namespace LPMP {

void multicut_odd_wheel_packing(const triplet_multicut_instance& input);
odd_wheel_packing compute_multicut_odd_wheel_packing(const triplet_multicut_instance& input);

quadruplet_multicut_instance pack_multicut_instance(const triplet_multicut_instance& input, const odd_wheel_packing& owp); 


} // namepsace LPMP
