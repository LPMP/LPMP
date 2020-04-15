#include "bdd/bdd_min_marginal_averaging_smoothed.h"
#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/ILP_parser.h"
#include "cuddObj.hh"

#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::cout << std::setprecision(10);

    if(argc < 2)
        throw std::runtime_error("input filename must be present as argument");

    const double min_progress = 1e-01;
    const int max_iter = 200;

    ILP_input input(ILP_parser::parse_file(std::string(argv[1])));
    //input.reorder_bfs();
    //input.reorder_Cuthill_McKee();
    input.reorder_minimum_degree_averaging();

    bdd_min_marginal_averaging_smoothed solver;
    solver.init(input);
    solver.set_cost_scaling(0.01);

    solver.smooth_forward_run();
    const double initial_lb_forward = solver.compute_smooth_lower_bound_forward();
    solver.smooth_backward_run();
    const double initial_lb_backward = solver.compute_smooth_lower_bound();
    std::cout << "initial lower bounds: " << initial_lb_backward << ", " << initial_lb_forward << "\n";
    //const double initial_lb = solver.lower_bound();
    //std::cout << "initial lower bound = " << initial_lb << "\n";

    //double old_lb = initial_lb;
    //double new_lb = old_lb;

    for(std::size_t iter=0; iter<max_iter; ++iter) {
        solver.smooth_iteration();
        const double new_lb = solver.compute_smooth_lower_bound(); // TODO: change, if available
        std::cout << "iteration " << iter << ": " << std::flush;
        std::cout << "lower bound = " << new_lb << "\n";
        //if (new_lb - old_lb < min_progress)
        //    break;
        //old_lb = new_lb;
    }

    solver.set_cost_scaling(1.0);
    solver.backward_run();
    double orig_lb = solver.compute_lower_bound();
    std::cout << "original final lower bound = " << orig_lb << "\n";
    for(std::size_t iter=0; iter<20; ++iter) {
        solver.iteration(); 
    }
    orig_lb = solver.compute_lower_bound();
    std::cout << "original final lower bound = " << orig_lb << "\n";
}

