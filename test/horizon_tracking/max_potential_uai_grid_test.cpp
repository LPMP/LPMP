#include "test_helper.hxx"
#include "data_max_potential_grids_test.hxx"

using namespace LPMP;

int main()
{
    // Test on 2x2 Grid:
    {
        TestUAI(solver_options_small, grid_uai_input_small, 0.242203);
    }

    // Test on 3x3 Grid:
    {
        TestUAI(solver_options_small, grid_uai_input_3x3, 0.590735);
    }

    // Test on 5x5 Grid:
    {
        TestUAI(solver_options_medium, grid_uai_input_medium, 4.307381);
    }

    // Test on 2x6 Grid having primal issues:
    {
        TestUAI(solver_options_medium, grid_uai_input_primal_issue, 4.307381);
    }
}
