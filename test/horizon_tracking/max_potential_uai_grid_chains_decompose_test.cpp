#include "LP.h"
#include "data_max_potential_grids_test.hxx"
#include "horizon_tracking_test_helper.hxx"

using namespace LPMP;

int main()
{
    // Test on 2x2 Grid:
    {
        TestUAIChains(solver_options_small, grid_uai_input_small, grid_uai_input_small_lb);
    }

    // Test on 3x3 Grid:
    {
        TestUAIChains(solver_options_small, grid_uai_input_3x3, grid_uai_input_3x3_lb);
    }

    // Test on 5x5 Grid:
    {
       TestUAIChains(solver_options_medium, grid_uai_input_medium, grid_uai_input_medium_lb);
    }
}
