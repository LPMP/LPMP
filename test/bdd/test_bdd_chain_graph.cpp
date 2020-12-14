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

    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    std::cout << std::setprecision(10);
    std::cout << "initial lower bound = " << initial_lb << std::endl;

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<1; ++iter) {
        bdds.iteration();
        new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << std::endl;
        if (new_lb - old_lb < 1e-08)
            break;
        old_lb = new_lb;
    }

    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}


void test_problem_fixing(const std::string input_string, const double expected_lb)
{
    ILP_input input = ILP_parser::parse_string(input_string);
    // input.reorder_bfs();
    // input.reorder_Cuthill_McKee();

    bdd_mma_fixing bdds;
    bdd_min_marginal_averaging_options bdd_options;
    bdd_options.averaging_type = bdd_min_marginal_averaging_options::averaging_type::SRMP;
    bdds.set_options(bdd_options);
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
    // input.reorder_bfs();
    // input.reorder_Cuthill_McKee();

    BDD_SOLVER bdds;
    bdds.init(input);

    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    std::cout << std::setprecision(10);
    std::cout << "initial lower bound = " << initial_lb << std::endl;

    double old_lb = initial_lb;
    double new_lb = old_lb;

    for(std::size_t iter=0; iter<100; ++iter) {
        std::cout << "iteration " << iter << ": " << std::flush;
        bdds.min_marginal_averaging_iteration_SRMP();
        // bdds.min_marginal_averaging_iteration();
        new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << std::endl;
        if (new_lb - old_lb < 1e-09)
            break;
        old_lb = new_lb;
    }

    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}

void test_fixing(const std::string input_string)
{
    const ILP_input input = ILP_parser::parse_string(input_string);

    bdd_mma_fixing bdds_test;
    bdds_test.init(input);

    std::vector<size_t> indices{0,1,2,3,4,5,6,7};

    std::vector<char> feasible{1,0,1,0,1,0,0,0};
    test(bdds_test.fix_variables(indices, feasible));

}

int main(int argc, char** arv)
{
    // test_problem_anisotropic(small_chain, 1.0);
    // test_problem_anisotropic(small_chain_reversed, 1.0);

    test_problem<bdd_min_marginal_averaging>(small_chain, 1.0);

    test_fixing(small_chain);

    test_problem_fixing(small_chain, 1.0);
}
