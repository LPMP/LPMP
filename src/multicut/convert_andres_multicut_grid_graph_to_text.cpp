#include "multicut/andres_graph_utility.hxx"

int main(int argc, char** argv)
{
   assert(argc == 3); // first arg input in opengm format, second is output in text format
   const std::string input_file = argv[1];
   const std::string output_file = argv[2];

   read_graph<andres::graph::GridGraph<2>>(input_file, output_file); 
}
