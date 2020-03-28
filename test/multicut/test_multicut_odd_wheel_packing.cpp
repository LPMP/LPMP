#include "multicut/multicut_odd_wheel_packing.h"
#include "test.h"

using namespace LPMP;

int main()
{
   // construct minimal multicut instances with 3 triplets such that packing them into a quadruplet increases the lower bound
   triplet_multicut_instance test_instance;
   for(std::size_t i=0; i<4; ++i)
      for(std::size_t j=i+1; j<4; ++j)
         test_instance.add_edge(i,j,0.0);
   // minimal example that is non tight from Nowozin's dissertation
   const double e01 = 0.9;
   const double e02 = -0.9;
   const double e03 = -0.7;
   const double e12 = 0.8;
   const double e13 = 0.2;
   const double e23 = -0.9;

   multicut_triplet_factor t012;
   multicut_triplet_factor t013;
   multicut_triplet_factor t023;
   multicut_triplet_factor t123;

   multicut_edge_triplet_message_0 msg01;
   multicut_edge_triplet_message_1 msg02;
   multicut_edge_triplet_message_2 msg12;

   msg01.RepamRight(t012, std::array<double,1>{e01/2.0});
   msg01.RepamRight(t013, std::array<double,1>{e01/2.0});

   msg02.RepamRight(t012, std::array<double,1>{e02/2.0});
   msg01.RepamRight(t023, std::array<double,1>{e02/2.0});

   msg02.RepamRight(t013, std::array<double,1>{e03/2.0});
   msg02.RepamRight(t023, std::array<double,1>{e03/2.0});

   msg12.RepamRight(t012, std::array<double,1>{e12/2.0});
   msg01.RepamRight(t123, std::array<double,1>{e12/2.0});

   msg12.RepamRight(t013, std::array<double,1>{e13/2.0});
   msg02.RepamRight(t123, std::array<double,1>{e13/2.0});

   msg12.RepamRight(t023, std::array<double,1>{e23/2.0});
   msg12.RepamRight(t123, std::array<double,1>{e23/2.0});

   test_instance.add_triplet({0,1,2}, t012);
   test_instance.add_triplet({0,1,3}, t013);
   test_instance.add_triplet({0,2,3}, t023);
   test_instance.add_triplet({1,2,3}, t123);

   auto owp = compute_multicut_odd_wheel_packing(test_instance);
   test(owp.no_odd_wheels() >= 1);
   // check whether the bound obtained by the odd wheel packing is greater than the initial lower bound
   test(1e-8 < owp.get_odd_wheel_weight(0));
   const quadruplet_multicut_instance quadruplet_instance = pack_multicut_instance(test_instance, owp);
   test(test_instance.lower_bound() + owp.get_odd_wheel_weight(0) <= quadruplet_instance.lower_bound() + 1e-8);
   std::cout << owp.no_odd_wheels() << " odd wheels found. First weight = " << owp.get_odd_wheel_weight(0) << "\n";
   std::cout << "quadruplet instance lower bound = " << quadruplet_instance.lower_bound() << "\n";
} 
