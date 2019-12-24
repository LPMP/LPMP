#include "cuthill-mckee.h"
#include "test.h"
#include <iostream>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::vector<std::vector<std::size_t>> a = {
        {1, 2, 3, 4, 5, 11},
        {0, 2, 3, 4, 5},
        {0, 1, 3, 4, 5},
        {0, 1, 2, 4, 5},
        {0, 1, 2, 3, 5},
        {0, 1, 2, 3, 4, 6},
        {5, 7, 8, 9, 10, 11},
        {6, 8, 9, 10, 11},
        {6, 7, 9, 10, 11},
        {6, 7, 8, 10, 11},
        {6, 7, 8, 9, 11},
        {6, 7, 8, 9, 10, 0}};

    two_dim_variable_array<std::size_t> adjacency(a);

    const auto permutation = Cuthill_McKee(adjacency);
    for(auto x : permutation)
        std::cout << x << " ";
    std::cout << "\n";

    test(permutation[0] != 0); 
}
