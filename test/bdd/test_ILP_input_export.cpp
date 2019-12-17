#include "bdd/ILP_parser.h"
#include "bdd/bdd.h"
#include <string>
#include "test.h"

using namespace LPMP;

void test_problem(const std::string input_string, const double expected_lb)
{
    const ILP_input input = ILP_parser::parse_string(input_string);
    Cudd bdd_mgr;
    bdd_base bdds(bdd_mgr);
    bdds.init(input); 

    const double initial_lb = bdds.lower_bound_backward_run();
    test(initial_lb <= expected_lb + 1e-8);

    for(std::size_t iter=0; iter<100; ++iter) {
        bdds.min_marginal_averaging_iteration();
    }

    const double lb = bdds.min_marginal_averaging_iteration();

    test(std::abs(lb - expected_lb) <= 1e-8);
}

const char * matching_3x3 = 
R"(Minimize
-2 x_11 - 1 x_12 - 1 x_13
-1 x_21 - 2 x_22 - 1 x_23
-1 x_31 - 1 x_32 - 2 x_33
Subject To
+ 1 x_11 + 1 x_12 + 1 x_13 = 1
x_21 + x_22 + x_23 = 1
x_31 + x_32 + x_33 = 1
x_11 + x_21 + x_31 = 1
x_12 + x_22 + x_32 = 1
- x_13 - x_23 - x_33 = -1
End)";

int main(int argc, char** arv)
{
    const ILP_input input_orig = ILP_parser::parse_string(matching_3x3);
    std::stringstream exported;
    input_orig.write(exported);
    const ILP_input input_exported = ILP_parser::parse_string(exported.str());


    Cudd bdd_mgr_1;
    bdd_base bdds_1(bdd_mgr_1);
    bdds_1.init(input_orig); 
    const double initial_lb_1 = bdds_1.lower_bound_backward_run();

    Cudd bdd_mgr_2;
    bdd_base bdds_2(bdd_mgr_2);
    bdds_2.init(input_exported); 
    const double initial_lb_2 = bdds_2.lower_bound_backward_run();

    test(std::abs(initial_lb_1 - initial_lb_2) <= 1e-8);
}

