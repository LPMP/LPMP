#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/ILP_parser.h"
#include "cuddObj.hh"

#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc < 2)
        throw std::runtime_error("input filename must be present as argument");

    const double min_progress = 1e-01;
    const int max_iter = 1000;

    ILP_input input(ILP_parser::parse_file(std::string(argv[1])));
    //input.reorder_bfs();
    //input.reorder_Cuthill_McKee();

    bdd_min_marginal_averaging_options options(argc-1, argv+1);
    bdd_min_marginal_averaging bdds;
    bdds.set_options(options);
    bdds.init(input);

    std::cout << std::setprecision(10);
    const double initial_lb = bdds.compute_lower_bound();
    std::cout << "initial lower bound = " << initial_lb << "\n";

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<max_iter; ++iter) {
        std::cout << "iteration " << iter << ": " << std::flush;
        bdds.iteration();
        const double new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << "\n";
        if (new_lb - old_lb < min_progress)
            break;
        old_lb = new_lb;
    } 
}
