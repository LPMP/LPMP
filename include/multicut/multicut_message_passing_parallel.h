#pragma once

#include "multicut_instance.h"

namespace LPMP {

multicut_edge_labeling multicut_message_passing_parallel(const multicut_instance& input, const bool record_cycles, const int nr_threads);


} // end namespace LPMP
