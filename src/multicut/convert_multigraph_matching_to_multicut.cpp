#include <iostream>
#include <fstream>
#include "graph_matching/multigraph_matching_input.h"
#include "multicut/transform_multigraph_matching.h"

using namespace LPMP;

int main(int argc, char** argv)
{
   if(argc != 3)
      throw std::runtime_error("two input arguments expected: input file, output file.");

   const std::string input_file = argv[1];
   auto mgm_input = Torresani_et_al_multigraph_matching_input::parse_file(input_file);
   auto multicut_input = transform_multigraph_matching_to_correlation_clustering(mgm_input);
   multicut_input.transform_to_multicut();

   const std::string output_file = argv[2];
   std::ofstream output_file_stream(output_file, std::ofstream::out);
   multicut_input.write_problem(output_file_stream);
}
