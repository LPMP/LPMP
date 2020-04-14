#include "cell-tracking/cell_tracking_input.h"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3)
        throw std::runtime_error("Two arguments expected: input and output");

    const cell_tracking_instance instance = cell_tracking_parser_2d::parse_file(argv[1]);

    std::ofstream out(argv[2]); 
    instance.write_to_lp(out);
}