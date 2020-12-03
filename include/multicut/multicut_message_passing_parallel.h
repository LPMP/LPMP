#pragma once

#include "multicut_instance.h"

namespace LPMP {

std::pair<multicut_instance, double> multicut_message_passing_parallel(const multicut_instance& input, const bool record_cycles, const int nr_threads);

} // end namespace LPMP
