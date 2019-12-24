#include "config.hxx"
#include "bdd/bdd.h"
#include "bdd/convert_pb_to_bdd.h"
#include <vector>
#include "test.h"

using namespace LPMP;

int main(int argc, char** arv)
{
    Cudd bdd_mgr;
    bdd_converter converter(bdd_mgr);

    // construct simplex constraint
    std::vector<int> simplex_weights = {1,1,1,1,1};
    auto simplex_bdd = converter.convert_to_bdd(simplex_weights.begin(), simplex_weights.end(), inequality_type::equal, 1);

    bdd_base bdds(bdd_mgr);

    std::vector<std::size_t> simplex_nodes1 = {0,1,2,3,4};
    bdds.add_bdd(simplex_bdd, simplex_nodes1.begin(), simplex_nodes1.end());

    std::vector<std::size_t> simplex_nodes2 = {4,5,6,7,8};
    bdds.add_bdd(simplex_bdd, simplex_nodes2.begin(), simplex_nodes2.end());

    bdds.init(); 

    std::vector<double> simplex_costs = {-2.0, 0.0, 0.0, 0.0, -3.5, 1.0, 1.0, 1.0, 1.0};
    bdds.set_costs(simplex_costs.begin(), simplex_costs.end());

    {
        bdds.forward_run();
        const double backward_lb = bdds.lower_bound_backward_run();
        test(backward_lb == -2.0 - 1.75);
    }

    {
        bdds.min_marginal_averaging_forward();
        const double backward_lb = bdds.lower_bound_backward_run();
        test(backward_lb == -3.5);
    }

    // now test backward averaging
    bdd_base bdds2(bdd_mgr);
    bdds2.add_bdd(simplex_bdd, simplex_nodes1.begin(), simplex_nodes1.end());
    bdds2.add_bdd(simplex_bdd, simplex_nodes2.begin(), simplex_nodes2.end());
    bdds2.init();
    bdds2.set_costs(simplex_costs.begin(), simplex_costs.end());
    test(bdds2.lower_bound_backward_run() == -3.75);
    bdds2.forward_run();
    std::cout << "\n\n\n";
    {
        const double lb = bdds2.min_marginal_averaging_backward();
        test(lb == -3.5);
    }
} 
