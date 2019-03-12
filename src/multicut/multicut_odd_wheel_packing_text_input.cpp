#include "multicut/multicut_cycle_packing.h"
#include "multicut/multicut_odd_wheel_packing.h"
#include "multicut/multicut_text_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
   if(argc != 2) 
      throw std::runtime_error("input file expected as argument");
   auto input = LPMP::multicut_text_input::parse_file(argv[1]);
   auto cp = compute_multicut_cycle_packing(input);
   const triplet_multicut_instance tmi = pack_multicut_instance(input, cp);
   auto owp = compute_multicut_odd_wheel_packing(tmi);
} 
