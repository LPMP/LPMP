#include "bdd/ILP_parser.h"
#include <string>
#include "test.h"
#include <iostream>

using namespace LPMP;

const std::string ILP_example =
R"(Minimize
x1 + 2*x2 + 1.5 * x3 - 0.5*x4 - x5
Subject To
x1 + 2*x2 + 3 * x3 - 5*x4 - x5 >= 1
End)";

int main(int argc, char** argv)
{
    const ILP_input input = ILP_parser::parse_string(ILP_example);
    test(input.nr_variables() == 5);
    test(input.var_exists("x1"));
    test(input.var_exists("x2"));
    test(input.var_exists("x3"));
    test(input.var_exists("x4"));
    test(input.var_exists("x5"));

    test(input.objective("x1") == 1.0);;
    test(input.objective("x2") == 2.0);;
    test(input.objective("x3") == 1.5);;
    test(input.objective("x4") == -0.5);;
    test(input.objective("x5") == -1.0);;

    input.write(std::cout);
    test(input.nr_constraints() == 1);
}
