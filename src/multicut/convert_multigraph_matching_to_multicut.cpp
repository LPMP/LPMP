#include <iostream>
#include <fstream>
#include "multicut/transform_multigraph_matching.h"
#include "multigraph_matching/multigraph_matching_input.h"
#include "multicut/transform_multigraph_matching.h"
#include "multicut/correlation_clustering_instance.h"

using namespace LPMP;

int main(int argc, char** argv)
{
   if(argc != 3)
      throw std::runtime_error("two input arguments expected: input file, output file.");

   const std::string input_file = argv[1];
   auto mgm_input = std::make_shared<multigraph_matching_input>(Torresani_et_al_multigraph_matching_input::parse_file(input_file));
   multigraph_matching_correlation_clustering_transform t(mgm_input);
   correlation_clustering_instance cc_input = t.get_correlatino_clustering_instance();

   const std::string output_file = argv[2];
   std::ofstream output_file_stream(output_file, std::ofstream::out);
   cc_input.write_problem(output_file_stream);
}
