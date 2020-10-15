#include <string>
#include "graph_matching/graph_matching_koopmans_beckmann_input.h"
#include "test.h"

using namespace LPMP;

const std::string KBQ = 
R"(L
0, 1
2, 3

A
4, 5
6, 7

B
8, 9
10, 11
)";

int main(int argc, char** argv)
{
    const auto instance = graph_matching_koopmans_beckmann_input::parse_string(KBQ);
    test(instance.L.rows() == 2 && instance.L.cols() == 2);
    test(instance.L(0,0) == 0.0);
    test(instance.L(0,1) == 1.0);
    test(instance.L(1,0) == 2.0);
    test(instance.L(1,1) == 3.0);

    test(instance.A.rows() == 2 && instance.A.cols() == 2);
    test(instance.A(0,0) == 4.0);
    test(instance.A(0,1) == 5.0);
    test(instance.A(1,0) == 6.0);
    test(instance.A(1,1) == 7.0);

    test(instance.B.rows() == 2 && instance.B.cols() == 2);
    test(instance.B(0,0) == 8.0);
    test(instance.B(0,1) == 9.0);
    test(instance.B(1,0) == 10.0);
    test(instance.B(1,1) == 11.0);
}
