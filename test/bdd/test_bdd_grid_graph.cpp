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


template<typename BDD_SOLVER>
void test_problem(const std::string input_string, const double expected_lb)
{
    const ILP_input input = ILP_parser::parse_string(input_string);

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
        bdds.iteration();
        new_lb = bdds.lower_bound();
        std::cout << "lower bound = " << new_lb << std::endl;
        if (new_lb - old_lb < 1e-09)
            break;
        old_lb = new_lb;
    }

    bdds.backward_run(); 
    const double lb = bdds.lower_bound();

    test(std::abs(lb - expected_lb) <= 1e-8);
}


const char * grid_graph_3x3 = 
R"(Minimize
2 mu_0_0 - 1 mu_0_1 + 3 mu_1_0 - 1 mu_1_1
+ 3 mu_2_0 + 2 mu_2_1 - 1 mu_3_0 - 2 mu_3_1
- 2 mu_4_0 - 1 mu_4_1 + 3 mu_5_0 - 1 mu_5_1
+ 1 mu_6_0 + 1 mu_6_1 - 3 mu_7_0 + 2 mu_7_1
+ 0 mu_8_0 + 2 mu_8_1
+ 1 mu_01_00 - 2 mu_01_01 + 2 mu_01_10 - 1 mu_01_11
+ 0 mu_12_00 + 1 mu_12_01 + 1 mu_12_10 + 0 mu_12_11
- 1 mu_03_00 + 2 mu_03_01 + 0 mu_03_10 - 2 mu_03_11
+ 2 mu_14_00 + 0 mu_14_01 + 2 mu_14_10 + 2 mu_14_11
+ 1 mu_25_00 - 2 mu_25_01 - 3 mu_25_10 - 1 mu_25_11
+ 0 mu_34_00 + 1 mu_34_01 + 1 mu_34_10 + 1 mu_34_11
- 1 mu_45_00 - 2 mu_45_01 + 4 mu_45_10 - 2 mu_45_11
- 2 mu_36_00 + 0 mu_36_01 + 1 mu_36_10 + 3 mu_36_11
+ 3 mu_47_00 - 2 mu_47_01 - 2 mu_47_10 - 1 mu_47_11
+ 0 mu_58_00 + 1 mu_58_01 + 1 mu_58_10 + 1 mu_58_11
- 1 mu_67_00 + 2 mu_67_01 - 1 mu_67_10 - 1 mu_67_11
+ 2 mu_78_00 + 0 mu_78_01 + 2 mu_78_10 + 2 mu_78_11
Subject To
mu_0_0 + mu_0_1 = 1
mu_1_0 + mu_1_1 = 1
mu_2_0 + mu_2_1 = 1
mu_3_0 + mu_3_1 = 1
mu_4_0 + mu_4_1 = 1
mu_5_0 + mu_5_1 = 1
mu_6_0 + mu_6_1 = 1
mu_7_0 + mu_7_1 = 1
mu_8_0 + mu_8_1 = 1
mu_01_00 + mu_01_10 + mu_01_01 + mu_01_11 = 1
mu_12_00 + mu_12_10 + mu_12_01 + mu_12_11 = 1
mu_03_00 + mu_03_10 + mu_03_01 + mu_03_11 = 1
mu_14_00 + mu_14_10 + mu_14_01 + mu_14_11 = 1
mu_25_00 + mu_25_10 + mu_25_01 + mu_25_11 = 1
mu_34_00 + mu_34_10 + mu_34_01 + mu_34_11 = 1
mu_45_00 + mu_45_10 + mu_45_01 + mu_45_11 = 1
mu_36_00 + mu_36_10 + mu_36_01 + mu_36_11 = 1
mu_47_00 + mu_47_10 + mu_47_01 + mu_47_11 = 1
mu_58_00 + mu_58_10 + mu_58_01 + mu_58_11 = 1
mu_67_00 + mu_67_10 + mu_67_01 + mu_67_11 = 1
mu_78_00 + mu_78_10 + mu_78_01 + mu_78_11 = 1
mu_0_0 - mu_01_00 - mu_01_01 = 0
mu_0_1 - mu_01_10 - mu_01_11 = 0
mu_0_0 - mu_03_00 - mu_03_01 = 0
mu_0_1 - mu_03_10 - mu_03_11 = 0
mu_1_0 - mu_01_00 - mu_01_10 = 0
mu_1_1 - mu_01_01 - mu_01_11 = 0
mu_1_0 - mu_12_00 - mu_12_01 = 0
mu_1_1 - mu_12_10 - mu_12_11 = 0
mu_1_0 - mu_14_00 - mu_14_01 = 0
mu_1_1 - mu_14_10 - mu_14_11 = 0
mu_2_0 - mu_12_00 - mu_12_10 = 0
mu_2_1 - mu_12_01 - mu_12_11 = 0
mu_2_0 - mu_24_00 - mu_24_01 = 0
mu_2_1 - mu_24_10 - mu_24_11 = 0
mu_3_0 - mu_03_00 - mu_03_10 = 0
mu_3_1 - mu_03_01 - mu_03_11 = 0
mu_3_0 - mu_34_00 - mu_34_01 = 0
mu_3_1 - mu_34_10 - mu_34_11 = 0
mu_3_0 - mu_36_00 - mu_36_01 = 0
mu_3_1 - mu_36_10 - mu_36_11 = 0
mu_4_0 - mu_14_00 - mu_14_10 = 0
mu_4_1 - mu_14_01 - mu_14_11 = 0
mu_4_0 - mu_34_00 - mu_34_10 = 0
mu_4_1 - mu_34_01 - mu_34_11 = 0
mu_4_0 - mu_45_00 - mu_45_01 = 0
mu_4_1 - mu_45_10 - mu_45_11 = 0
mu_4_0 - mu_47_00 - mu_47_01 = 0
mu_4_1 - mu_47_10 - mu_47_11 = 0
mu_5_0 - mu_25_00 - mu_25_10 = 0
mu_5_1 - mu_25_01 - mu_25_11 = 0
mu_5_0 - mu_45_00 - mu_45_10 = 0
mu_5_1 - mu_45_01 - mu_45_11 = 0
mu_5_0 - mu_58_00 - mu_58_01 = 0
mu_5_1 - mu_58_10 - mu_58_11 = 0
mu_6_0 - mu_36_00 - mu_36_10 = 0
mu_6_1 - mu_36_01 - mu_36_11 = 0
mu_6_0 - mu_67_00 - mu_67_01 = 0
mu_6_1 - mu_67_10 - mu_67_11 = 0
mu_7_0 - mu_47_00 - mu_47_10 = 0
mu_7_1 - mu_47_01 - mu_47_11 = 0
mu_7_0 - mu_67_00 - mu_67_10 = 0
mu_7_1 - mu_67_01 - mu_67_11 = 0
mu_7_0 - mu_78_00 - mu_78_01 = 0
mu_7_1 - mu_78_10 - mu_78_11 = 0
mu_8_0 - mu_58_00 - mu_58_10 = 0
mu_8_1 - mu_58_01 - mu_58_11 = 0
mu_8_0 - mu_78_00 - mu_78_10 = 0
mu_8_1 - mu_78_01 - mu_78_11 = 0
End)";


int main(int argc, char** arv)
{
    test_problem_anisotropic(grid_graph_3x3, -9.0);

    test_problem<bdd_min_marginal_averaging>(grid_graph_3x3, -9.0);
}
