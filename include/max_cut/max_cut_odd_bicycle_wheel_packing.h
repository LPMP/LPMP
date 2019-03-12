#pragma once

#include "cut_base/cut_base_packing.h"
#include "max_cut/max_cut_instance.hxx"

namespace LPMP {

    void max_cut_odd_bicycle_wheel_packing(const triplet_max_cut_instance& input);
    odd_bicycle_wheel_packing compute_max_cut_odd_bicycle_wheel_packing(const triplet_max_cut_instance& input); 

}
