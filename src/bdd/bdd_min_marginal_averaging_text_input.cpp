#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/bdd_primal_fixing.h"
#include "bdd/ILP_parser.h"
#include "cuddObj.hh"

#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc < 2)
        throw std::runtime_error("input filename must be present as argument");

    const double min_progress = 1e-06; // relative to objective function
    const int max_iter = 10000;

    const auto start_time = std::chrono::steady_clock::now();

    ILP_input input(ILP_parser::parse_file(std::string(argv[1])));
    bdd_min_marginal_averaging_options options(argc-1, argv+1);

    if (options.variable_order == bdd_min_marginal_averaging_options::variable_order::bfs)
        input.reorder_bfs();
    else if (options.variable_order == bdd_min_marginal_averaging_options::variable_order::cuthill)
        input.reorder_Cuthill_McKee();
    else if (options.variable_order == bdd_min_marginal_averaging_options::variable_order::mindegree)
        input.reorder_minimum_degree_averaging();

    bdd_mma_fixing solver;
    solver.set_options(options);
    solver.init(input);

    std::cout << "\#variables: " << solver.nr_variables() << std::endl;
    std::cout << "\#constraints: " << solver.nr_bdds() << std::endl;

    std::cout << std::setprecision(10);
    const double initial_lb = solver.compute_lower_bound();
    std::cout << "initial lower bound = " << initial_lb << std::flush;
    auto time = std::chrono::steady_clock::now();
    std::cout << ", time = " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(time - start_time).count() / 1000 << " s" << std::endl;

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<max_iter; ++iter) {
        std::cout << "iteration " << iter << ": " << std::flush;
        solver.iteration();
        const double new_lb = solver.lower_bound();
        std::cout << "lower bound = " << new_lb << std::flush;
        time = std::chrono::steady_clock::now();
        std::cout << ", time = " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(time - start_time).count() / 1000 << " s" << std::endl;
        if (std::abs((new_lb - old_lb) / old_lb) < min_progress)
        {
            std::cout << "Improvement less than " << min_progress*100 << "\%." << std::endl;
            break;
        }
        old_lb = new_lb;
        if (iter+1==max_iter)
            std::cout << "Maximum number of iterations reached." << std::endl;
    }
    std::cout << "Final lower bound: " << solver.lower_bound() << std::endl;

    // Primal heuristic
    double ub = std::numeric_limits<double>::infinity();
    std::cout << "Primal heuristic.. " << std::endl;
    if (solver.fix_variables()) {
        ub = solver.compute_upper_bound();
        std::cout << "\nPrimal solution value: " << ub << std::flush;
    } else
        std::cout << "\nNo primal solution found." << std::flush;
    time = std::chrono::steady_clock::now();
    std::cout << ", time = " << (double) std::chrono::duration_cast<std::chrono::milliseconds>(time - start_time).count() / 1000 << " s" << std::endl;
}
