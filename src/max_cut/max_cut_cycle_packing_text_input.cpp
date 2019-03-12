#include "max_cut/max_cut_instance.hxx"
#include "max_cut/max_cut_cycle_packing.h"
#include "max_cut/max_cut_text_input.h"

using namespace LPMP;
int main(int argc, char** argv) {
   if(argc != 2) 
      throw std::runtime_error("input file expected as argument");
   auto input = LPMP::max_cut_text_input::parse_file(argv[1]);
   max_cut_cycle_packing(input);
} 
