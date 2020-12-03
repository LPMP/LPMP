#include "bdd/ILP_parser.h"
#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/bdd_anisotropic_diffusion.h"
#include <string>
#include "test.h"

using namespace LPMP;

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

template<typename BDD_SOLVER>
void test_input_export()
{
    const ILP_input input_orig = ILP_parser::parse_string(matching_3x3);
    std::stringstream exported;
    input_orig.write(exported);
    const ILP_input input_exported = ILP_parser::parse_string(exported.str());


    BDD_SOLVER bdds_1;
    bdds_1.init(input_orig); 
    bdds_1.backward_run();
    bdds_1.compute_lower_bound();
    const double initial_lb_1 = bdds_1.lower_bound();

    bdd_min_marginal_averaging bdds_2;
    bdds_2.init(input_exported); 
    bdds_2.backward_run();
    bdds_2.compute_lower_bound();
    const double initial_lb_2 = bdds_2.lower_bound();

    test(std::abs(initial_lb_1 - initial_lb_2) <= 1e-8);
}

int main(int argc, char** arv)
{
    test_input_export<bdd_min_marginal_averaging>();
    test_input_export<bdd_anisotropic_diffusion>();
} 
