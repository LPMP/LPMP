#include "max_cut/max_cut_instance.hxx"
#include "max_cut/max_cut_cycle_packing.h"
#include "max_cut/max_cut_odd_bicycle_wheel_packing.h"
#include "test.h"

using namespace LPMP;

int main()
{
    triplet_max_cut_instance input;
    input.add_edge(0, 1, -1.0);
    input.add_edge(0, 2, -10);
    input.add_edge(0, 3, -1.0);
    input.add_edge(0, 4, -1.0);
    input.add_edge(1, 2, -1.0);
    input.add_edge(1, 3, -1.0);
    input.add_edge(1, 4, -1.0);
    input.add_edge(2, 3, -1.0);
    input.add_edge(2, 4, -1.0);
    input.add_edge(3, 4, -1.0);

    max_cut_triplet_factor f;
    test(f.size() == 3);
    f[0] = -1.0;
    f[1] = -1.0;
    f[2] = -1.0;
    input.add_triplet({0,1,2}, f);
    input.add_triplet({0,1,3}, f);
    input.add_triplet({0,1,4}, f);
    input.add_triplet({0,2,3}, f);
    input.add_triplet({0,2,4}, f);
    input.add_triplet({0,3,4}, f);
    input.add_triplet({1,2,3}, f);
    input.add_triplet({1,2,4}, f);
    input.add_triplet({1,3,4}, f);
    input.add_triplet({2,3,4}, f);

    const auto cp = compute_max_cut_odd_bicycle_wheel_packing(input);

    test(cp.no_odd_bicycle_wheels() > 0);
}
