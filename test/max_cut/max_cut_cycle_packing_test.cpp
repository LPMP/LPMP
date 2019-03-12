#include "max_cut/max_cut_instance.hxx"
#include "max_cut/max_cut_cycle_packing.h"
#include "test.h"

using namespace LPMP;

int main()
{
    max_cut_instance input;
    input.add_edge(0, 1, -1.0);
    input.add_edge(0, 2, -1.0);
    input.add_edge(0, 3, -1.0);
    input.add_edge(1, 2, -1.0);
    input.add_edge(1, 3, -1.0);
    input.add_edge(2, 3, -1.0);

    const auto cp = compute_max_cut_cycle_packing(input);

    test(cp.no_cycles() > 0);
}
