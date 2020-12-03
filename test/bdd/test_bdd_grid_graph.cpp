#include "bdd/ILP_parser.h"
#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/bdd_primal_fixing.h"
#include "bdd/bdd_anisotropic_diffusion.h"
#include "bdd/convert_pb_to_bdd.h"
#include "ILP_sample_problems.h"
#include <vector>
#include <string>
#include <iostream>
#include "test.h"

using namespace LPMP;

void test_problem_anisotropic(const std::string input_string, const double expected_lb)
{
    const ILP_input input = ILP_parser::parse_string(input_string);
    bdd_anisotropic_diffusion bdds;
    bdds.init(input);

    bdds.backward_run(); 
    bdds.compute_lower_bound(); 
    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    std::cout << std::setprecision(10);
    std::cout << "initial lower bound = " << initial_lb << std::endl;

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<100; ++iter) {
        bdds.iteration();
        new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << std::endl;
        if (new_lb - old_lb < 1e-09)
            break;
        old_lb = new_lb;
    }

    bdds.compute_lower_bound(); 
    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}


void test_problem_fixing(const std::string input_string, const double expected_lb)
{
    ILP_input input = ILP_parser::parse_string(input_string);
    input.reorder_bfs();
    // input.reorder_Cuthill_McKee();

    bdd_mma_fixing bdds;
    bdds.init(input); 

    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    std::cout << std::setprecision(10);
    std::cout << "initial lower bound = " << initial_lb << std::endl;

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<100; ++iter) {
        std::cout << "iteration " << iter << ": " << std::flush;
        bdds.iteration();
        new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << std::endl;
        // const double ub = bdds.compute_upper_bound();
        // std::cout << "upper bound = " << ub << std::endl;
        if (new_lb - old_lb < 1e-09)
            break;
        old_lb = new_lb;
    }

    double ub = std::numeric_limits<double>::infinity();
    if (bdds.fix_variables())
        ub = bdds.compute_upper_bound();

    std::cout << "Primal solution value: " << ub << std::endl;

    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}


template<typename BDD_SOLVER>
void test_problem(const std::string input_string, const double expected_lb)
{
    ILP_input input = ILP_parser::parse_string(input_string);
    input.reorder_bfs();
    // input.reorder_Cuthill_McKee();

    BDD_SOLVER bdds;
    bdds.init(input); 

    bdds.backward_run(); 
    bdds.compute_lower_bound(); 
    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    std::cout << std::setprecision(10);
    std::cout << "initial lower bound = " << initial_lb << std::endl;

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<100; ++iter) {
        std::cout << "iteration " << iter << ": " << std::flush;
        bdds.min_marginal_averaging_iteration_SRMP();
        new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << std::endl;
        // const double ub = bdds.compute_upper_bound();
        // std::cout << "upper bound = " << ub << std::endl;
        if (new_lb - old_lb < 1e-09)
            break;
        old_lb = new_lb;
    }

    bdds.backward_run(); 
    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}


int main(int argc, char** arv)
{
    // test_problem_anisotropic(grid_graph_3x3, -9.0);

    test_problem<bdd_min_marginal_averaging>(grid_graph_3x3, -9.0);

    test_problem_fixing(grid_graph_3x3, -9.0);
}
