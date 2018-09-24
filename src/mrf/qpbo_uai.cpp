#include "mrf/qpbo_factor.hxx"
#include "mrf/mrf_uai_input.h"
#include <iostream>

using namespace LPMP;

int main(int argc, char** argv)
{
   if(argc == 1) throw std::runtime_error("input file argument needed");
   auto input = mrf_uai_input::parse_file(argv[1]);
   for(std::size_t i=0; i<input.no_variables(); ++i)
      if(input.cardinality(i) != 2) 
         throw std::runtime_error("binary mrf expected");
   qpbo_factor q(input);
   const std::size_t no_persistent_labels = q.no_persistent_labels();
   std::cout << no_persistent_labels << " are persistent out of " << q.no_variables() << " nodes\n"; 
}
