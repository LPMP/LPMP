#include "multicut/multicut_cycle_packing.h"
#include "test.h"

using namespace LPMP;

int main()
{
   multicut_instance test_instance;
   test_instance.add_edge(0,1,-1);
   test_instance.add_edge(0,2,+1);
   test_instance.add_edge(1,2,+1);

   test(test_instance.no_edges() == 3);
   test(std::abs(test_instance.lower_bound() - -1.0) <= 1e-8);
   auto cp = compute_multicut_cycle_packing(test_instance);
   auto triplet_instance = pack_multicut_instance(test_instance, cp);
   test(triplet_instance.edges().size() == 3);
   test(triplet_instance.triplets().size() == 1);
   test(std::abs(triplet_instance.lower_bound()) <= 1e-8);
} 
