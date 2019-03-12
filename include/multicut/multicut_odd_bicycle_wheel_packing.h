#pragma once

#include "cut_base/cut_base_packing.h"
#include "multicut/multicut_instance.h"

namespace LPMP {

    void multicut_odd_bicycle_wheel_packing(const quadruplet_multicut_instance& input);
    odd_bicycle_wheel_packing compute_multicut_odd_bicycle_wheel_packing(const quadruplet_multicut_instance& input); 

}
