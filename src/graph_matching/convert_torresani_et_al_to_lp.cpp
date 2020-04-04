#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching_input.h"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3)
        throw std::runtime_error("Two arguments expected: input and output");

    const graph_matching_input instance = TorresaniEtAlInput::parse_file(argv[1]);

    std::ofstream out(argv[2]);
    instance.write_lp(out); 
}