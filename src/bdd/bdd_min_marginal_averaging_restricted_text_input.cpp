#include "bdd/bdd_min_marginal_averaging_restricted.h"
#include "bdd/ILP_parser.h"
#include "bdd.h"

#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc < 2)
        throw std::runtime_error("input filename must be present as argument");

    const double min_progress = 1e-01;
    const int max_iter = 100;

    ILP_input input(ILP_parser::parse_file(std::string(argv[1])));
    input.reorder_bfs();
    //input.reorder_Cuthill_McKee();
    //input.reorder_minimum_degree_averaging();

    bdd_min_marginal_averaging_options options(argc-1, argv+1);

    bdd_min_marginal_averaging_restricted solver;
    solver.set_options(options);
    solver.init(input);

    std::cout << std::setprecision(10);
    const double initial_lb = solver.lower_bound();
    std::cout << "initial lower bound = " << initial_lb << "\n";

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<max_iter; ++iter) {
        solver.iteration();
        const double new_lb = solver.lower_bound();
        std::cout << "iteration " << iter << ": " << std::flush;
        std::cout << "lower bound = " << new_lb << "\n";
        if (new_lb - old_lb < min_progress)
            break;
        old_lb = new_lb;
    }
}
