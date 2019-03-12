#include "multicut/multicut_odd_bicycle_wheel_packing.h"

using namespace LPMP;

int main()
{
   // construct minimal multicut instances with 3 triplets such that packing them into a quadruplet increases the lower bound
   quadruplet_multicut_instance test_instance;
   for(std::size_t i=0; i<5; ++i)
      for(std::size_t j=i+1; j<5; ++j)
         test_instance.add_edge(i,j,0.0);
   
   // quadruplet edge order: 01,02,03,12,13,23
   
   test_instance.add_quadruplet({0,1,2,3}, multicut_quadruplet_factor());
   test_instance.add_quadruplet({0,1,2,4}, multicut_quadruplet_factor());
   test_instance.add_quadruplet({0,1,3,4}, multicut_quadruplet_factor());
   test_instance.add_quadruplet({0,2,3,4}, multicut_quadruplet_factor());
   test_instance.add_quadruplet({1,2,3,4}, multicut_quadruplet_factor());


   auto obwp = compute_multicut_odd_bicycle_wheel_packing(test_instance);
} 

