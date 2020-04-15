#include "discrete_tomography/discrete_tomography_input.h"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3)
        throw std::runtime_error("Two arguments expected: input and output");

    const discrete_tomography_instance instance = discrete_tomography_UAI_input::parse_file(argv[1]);

    std::ofstream out(argv[2]); 
    instance.write_to_lp(out);
}
