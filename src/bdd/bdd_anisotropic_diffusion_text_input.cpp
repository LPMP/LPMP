#include "bdd/bdd_anisotropic_diffusion.h"
#include "bdd/ILP_parser.h"
#include "cuddObj.hh"

#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("input filename must be present as argument");

    const ILP_input input = ILP_parser::parse_file(std::string(argv[1]));
    std::cout << "parsed input\n";
    std::ofstream test_out(std::string(argv[1]) + ".lp");
    input.write(test_out);
    test_out.close();

    Cudd bdd_mgr;
    bdd_anisotropic_diffusion bdds;
    std::cout << "Initializing bdds\n";
    bdds.init(input);
    std::cout << "Initialized bdds\n";

    std::cout << std::setprecision(10);
    const double initial_lb = bdds.lower_bound();
    std::cout << "initial lower bound = " << initial_lb << "\n";

    for(std::size_t iter=0; iter<1000; ++iter) {
        bdds.iteration();
        const double lb = bdds.lower_bound();
        std::cout << "iteration " << iter << ": lower bound = " << lb << "\n";
    } 
}

