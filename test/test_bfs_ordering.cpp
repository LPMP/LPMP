#include "test.h"
#include "bfs_ordering.hxx"

using namespace LPMP;

int main(int arc, char** argv)
{
    std::vector<std::vector<std::size_t>> graph_adj = {
        {1,2},
        {0,3,4,5},
        {0,3,4,5},
        {1,2,4,5},
        {1,2,3,6},
        {1,2,3,6},
        {4,5} 
    };

    const auto order = bfs_ordering(graph_adj);
    test(order[0] == 0 || order[6] == 0);
    test(order[1] == 1 || order[1] == 2 || order[1] == 4 || order[1] == 5);
    test(order[2] == 1 || order[2] == 2 || order[2] == 4 || order[2] == 5);
    test(order[3] == 3);
    test(order[4] == 1 || order[4] == 2 || order[4] == 4 || order[4] == 5);
    test(order[5] == 1 || order[5] == 2 || order[5] == 4 || order[5] == 5);
    test(order[6] == 6 || order[6] == 6);
}
