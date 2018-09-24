#include "mrf/qpbo_factor.hxx"
#include "mrf/dimacs_max_flow_input.h"
#include "mrf/transform_max_flow_instance.hxx"
#include <iostream>

using namespace LPMP;

int main(int argc, char** argv)
{
   if(argc == 1) throw std::runtime_error("input file argument needed");
   auto max_flow_input = dimacs_max_flow_input::parse_file(argv[1]);
   auto mrf_input = transform_QPBO_max_flow_to_binary_MRF(max_flow_input);
   qpbo_factor q(mrf_input);
   const std::size_t no_persistent_labels = q.no_persistent_labels();
   std::cout << no_persistent_labels << " are persistent out of " << q.no_variables() << " nodes\n"; 
} 
