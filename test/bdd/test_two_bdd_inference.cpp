#include "config.hxx"
#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/bdd_anisotropic_diffusion.h"
#include "bdd/convert_pb_to_bdd.h"
#include <vector>
#include "test.h"

using namespace LPMP;

template<typename BDD_SOLVER>
void test_two_bdd_inference()
{
    Cudd bdd_mgr;
    bdd_converter converter(bdd_mgr);

    // construct simplex constraint
    std::vector<int> simplex_weights = {1,1,1,1,1};
    auto simplex_bdd = converter.convert_to_bdd(simplex_weights.begin(), simplex_weights.end(), inequality_type::equal, 1);

    BDD_SOLVER bdds;

    std::vector<std::size_t> simplex_nodes1 = {0,1,2,3,4};
    bdds.add_bdd(simplex_bdd, simplex_nodes1.begin(), simplex_nodes1.end(), bdd_mgr);

    std::vector<std::size_t> simplex_nodes2 = {4,5,6,7,8};
    bdds.add_bdd(simplex_bdd, simplex_nodes2.begin(), simplex_nodes2.end(), bdd_mgr);

    bdds.init(); 

    std::vector<double> simplex_costs = {-2.0, 0.0, 0.0, 0.0, -3.5, 1.0, 1.0, 1.0, 1.0};
    bdds.set_costs(simplex_costs.begin(), simplex_costs.end());

    {
        bdds.backward_run();
        bdds.compute_lower_bound();
        const double backward_lb = bdds.lower_bound();
        test(std::abs(backward_lb - -3.75) <= 1e-8);
    }

    {
        bdds.iteration();
        const double backward_lb = bdds.lower_bound();
        test(std::abs(backward_lb - -3.5) <= 1e-8);
    }
}

int main(int argc, char** arv)
{
    test_two_bdd_inference<bdd_min_marginal_averaging>();
    test_two_bdd_inference<bdd_anisotropic_diffusion>(); 
}
