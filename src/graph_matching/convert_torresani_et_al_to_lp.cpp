#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching_input.h"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc <= 3)
        throw std::runtime_error("Two arguments expected: input and output");
    bool no_matching_allowed = true;
    if(argc == 4)
    {
        std::cout << "Additional argument \'" << argv[3] << "\'\n";
       if(std::string(argv[3]).compare("--bijective_assignment") || std::string(argv[3]).compare("-bijective_assignment"))
       {
           no_matching_allowed = false;
           std::cout << "Produce bijective assignment\n";
       }
    }
        

    const graph_matching_input instance = TorresaniEtAlInput::parse_file(argv[1]);

    std::ofstream out(argv[2]);
    instance.write_lp(out, no_matching_allowed); 
}
