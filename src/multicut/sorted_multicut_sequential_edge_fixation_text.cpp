#include "multicut/multicut_sequential_edge_fixation.hxx"
#include "multicut/multicut_instance.hxx"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
   if(argc <= 1) throw std::runtime_error("at least input file argument must be present");
   auto input = multicut_text_input::parse_file(argv[1]);
   multicut_sequential_edge_fixation sef(input);
   if(!sef.multicut_determined()) throw std::runtime_error("no unique multicut determined. Is graph connected?");

   if(argc >= 3) {
      std::fstream f;
      f.open(argv[2]);
      sef.write_node_labeling(f);
   }
}
