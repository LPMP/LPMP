#include <iostream>
#include "multigraph_matching/multigraph_matching_input.h"

using namespace LPMP;

int main(int argc, char** argv) {
   if(argc != 3)
      throw std::runtime_error("two arguments expected: problem instance file, solution file");
   auto problem = Torresani_et_al_multigraph_matching_input::parse_file(argv[1]);
   auto solution = parse_multigraph_matching_result_file(argv[2]);

   std::cout << "objective = " << problem.evaluate(solution) << "\n";
}
