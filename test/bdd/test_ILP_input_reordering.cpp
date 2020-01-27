#include "bdd/ILP_input.h"
#include "bdd/ILP_parser.h"
#include "test.h"

using namespace LPMP;

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

int main(int argc, char** argv)
{
    ILP_input input = ILP_parser::parse_string(small_chain);
    std::vector<char> sol = {1,0, 0,1, 0,0,1,0};
    const double cost = input.evaluate(sol.begin(), sol.end());
    test(std::abs(cost - (2 + 0 + 1)) <= 1e-8);

    const auto new_order = input.reorder_bfs();
    const auto new_sol = new_order.permute(sol.begin(), sol.end());
    const double new_cost = input.evaluate(new_sol.begin(), new_sol.end());
    test(std::abs(new_cost - (2 + 0 + 1)) <= 1e-8); 
}
