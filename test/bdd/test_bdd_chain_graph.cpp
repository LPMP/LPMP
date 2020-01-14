#include "bdd/ILP_parser.h"
#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/bdd_anisotropic_diffusion.h"
#include "bdd/convert_pb_to_bdd.h"
#include <vector>
#include <string>
#include <iostream>
#include "test.h"

using namespace LPMP;

void test_problem_anisotropic(const std::string input_string, const double expected_lb)
{
    const ILP_input input = ILP_parser::parse_string(input_string);
    Cudd bdd_mgr;
    bdd_anisotropic_diffusion bdds;
    bdds.init(input);

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
        if (new_lb - old_lb < 1e-08)
            break;
        old_lb = new_lb;
    }

    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}


const char * small_chain = 
R"(Minimize
2 mu_1_0 + 1 mu_1_1 - 1 mu_2_0 + 0 mu_2_1
+ 1 mu_00 + 2 mu_10 + 1 mu_01 + 0 mu_11
Subject To
mu_1_0 + mu_1_1 = 1
mu_2_0 + mu_2_1 = 1
mu_00 + mu_10 + mu_01 + mu_11 = 1
mu_1_0 - mu_00 - mu_01 = 0
mu_1_1 - mu_10 - mu_11 = 0
mu_2_0 - mu_00 - mu_10 = 0
mu_2_1 - mu_01 - mu_11 = 0
End)";



int main(int argc, char** arv)
{
    test_problem_anisotropic(small_chain, 1.0);
}
