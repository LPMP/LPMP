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
    bdd_anisotropic_diffusion bdds;
    bdds.init(input); 

    //bdds.export_dot(std::cout);

    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    for(std::size_t iter=0; iter<100; ++iter) {
        bdds.iteration();
    }

    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}

template<typename BDD_SOLVER>
void test_problem(const std::string input_string, const double expected_lb)
{
    const ILP_input input = ILP_parser::parse_string(input_string);
    BDD_SOLVER bdds;
    bdds.init(input); 

    const double initial_lb = bdds.lower_bound();
    test(initial_lb <= expected_lb + 1e-8);

    for(std::size_t iter=0; iter<100; ++iter) {
        bdds.iteration();
    }

    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}

const char * matching_3x3_diag = 
R"(Minimize
-2 x_11 - 1 x_12 - 1 x_13
-1 x_21 - 2 x_22 - 1 x_23
-1 x_31 - 1 x_32 - 2 x_33
Subject To
x_11 + x_12 + x_13 = 1
x_21 + x_22 + x_23 = 1
x_31 + x_32 + x_33 = 1
x_11 + x_21 + x_31 = 1
x_12 + x_22 + x_32 = 1
x_13 + x_23 + x_33 = 1
End)";

const char * matching_3x3_first_row = 
R"(Minimize
-2 x_11 - 1 x_12 - 1 x_13
-2 x_21 - 1 x_22 - 1 x_23
-2 x_31 - 1 x_32 - 1 x_33
Subject To
x_11 + x_12 + x_13 = 1
x_21 + x_22 + x_23 = 1
x_31 + x_32 + x_33 = 1
x_11 + x_21 + x_31 = 1
x_12 + x_22 + x_32 = 1
x_13 + x_23 + x_33 = 1
End)";


int main(int argc, char** arv)
{
    test_problem<bdd_min_marginal_averaging>(matching_3x3_diag, -6.0);
    test_problem<bdd_min_marginal_averaging>(matching_3x3_first_row, -4.0);

    test_problem_anisotropic(matching_3x3_diag, -6.0);
    test_problem_anisotropic(matching_3x3_first_row, -4.0);
}
