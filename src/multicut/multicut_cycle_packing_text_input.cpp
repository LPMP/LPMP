#include "multicut/multicut_cycle_packing.h"
#include "multicut/multicut_text_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
   if(argc != 2) 
      throw std::runtime_error("input file expected as argument");
   auto input = LPMP::multicut_text_input::parse_file(argv[1]);
   multicut_cycle_packing(input);
} 
