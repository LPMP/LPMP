#include "mrf/mrf_input.h"
#include "mrf/mrf_opengm_input.h"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3)
        throw std::runtime_error("2 arguments expected: opengm input file, output file.");

    const auto mrf = mrf_opengm_input::parse_file(argv[1]);

    std::ofstream out;
    out.open(argv[2], std::ios::out);
    mrf.write(out);
    out.close();
}
