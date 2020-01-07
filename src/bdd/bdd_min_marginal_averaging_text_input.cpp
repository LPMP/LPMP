#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/ILP_parser.h"
#include "cuddObj.hh"

#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("input filename must be present as argument");

    const ILP_input input = ILP_parser::parse_file(std::string(argv[1]));

    bdd_min_marginal_averaging bdds;
    bdds.init(input);

    std::cout << std::setprecision(10);
    const double initial_lb = bdds.lower_bound_backward_run();
    std::cout << "initial lower bound = " << initial_lb << "\n";

    for(std::size_t iter=0; iter<1000; ++iter) {
        const double lb = bdds.min_marginal_averaging_iteration();
        std::cout << "iteration " << iter << ": lower bound = " << lb << "\n";
    } 
}
