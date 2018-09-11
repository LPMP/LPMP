#include "test.h"
#include "horizon_tracking/horizon_tracking_factors.hxx"

using namespace LPMP;

int main()
{
    // problem with zero linear potentials
    const int numNodes = 3;
    std::vector<INDEX> numLabels = {3, 2, 3};
    std::vector<std::array<INDEX,2>> potential_size{{3,2}, {2,3}};

    three_dimensional_variable_array<REAL> LinearPairwisePotentials(potential_size.begin(), potential_size.end());
    three_dimensional_variable_array<REAL> MaxPairwisePotentials(potential_size.begin(), potential_size.end());
    std::vector<std::array<INDEX, 2>> messagePassingSchedule;

    messagePassingSchedule.push_back({0, 1});
    messagePassingSchedule.push_back({2, 1});

    MaxPairwisePotentials(0,0,0) = 1;
    MaxPairwisePotentials(0,1,0) = 2;
    MaxPairwisePotentials(0,2,0) = 1;
    MaxPairwisePotentials(0,0,1) = 3;
    MaxPairwisePotentials(0,1,1) = 4;
    MaxPairwisePotentials(0,2,1) = 2;

    MaxPairwisePotentials(1,0,0) = 1.5;
    MaxPairwisePotentials(1,0,1) = 5;
    MaxPairwisePotentials(1,0,2) = 1;
    MaxPairwisePotentials(1,1,0) = 3;
    MaxPairwisePotentials(1,1,1) = 1;
    MaxPairwisePotentials(1,1,2) = 6;

    {
        max_potential_on_tree_iterative tree = max_potential_on_tree_iterative(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, messagePassingSchedule);
        tree.MaximizePotentialAndComputePrimal();
        auto best_marginal = tree.max_potential_marginal(tree.max_potential_index());
        test(best_marginal.MaxCost + best_marginal.LinearCost,  1, 0);
    }

    // now add linear potentials
    LinearPairwisePotentials(0,0,0) = 100;
    LinearPairwisePotentials(0,1,0) = 100;
    LinearPairwisePotentials(0,2,0) = 100;
    LinearPairwisePotentials(0,0,1) = 100;
    LinearPairwisePotentials(0,1,1) = 100;
    LinearPairwisePotentials(0,2,1) = 0;

    LinearPairwisePotentials(1,0,0) = 100;
    LinearPairwisePotentials(1,0,1) = 100;
    LinearPairwisePotentials(1,0,2) = 100;
    LinearPairwisePotentials(1,1,0) = 100;
    LinearPairwisePotentials(1,1,1) = 100;
    LinearPairwisePotentials(1,1,2) = 0;

    {
        max_potential_on_tree_iterative tree = max_potential_on_tree_iterative(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, messagePassingSchedule);
        tree.MaximizePotentialAndComputePrimal();
        auto best_marginal = tree.max_potential_marginal(tree.max_potential_index());
        test(best_marginal.MaxCost + best_marginal.LinearCost,  6, 0);
    }

    LinearPairwisePotentials(0,2,0) = 1.5;
    LinearPairwisePotentials(1,0,2) = 1.5;

    {
        max_potential_on_tree_iterative tree = max_potential_on_tree_iterative(MaxPairwisePotentials, LinearPairwisePotentials, numLabels, messagePassingSchedule);
        tree.MaximizePotentialAndComputePrimal();
        auto best_marginal = tree.max_potential_marginal(tree.max_potential_index());
        test(best_marginal.MaxCost + best_marginal.LinearCost,  4, 0);
    }
}
