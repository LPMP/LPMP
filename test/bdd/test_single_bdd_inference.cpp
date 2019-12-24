#include "config.hxx"
#include "bdd/bdd.h"
#include "bdd/convert_pb_to_bdd.h"
#include <vector>
#include <random>
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
    std::vector<std::size_t> simplex_nodes = {0,1,2,3,4};
    bdds.add_bdd(simplex_bdd, simplex_nodes.begin(), simplex_nodes.end());
    bdds.init(); 

    test(bdds.nr_variables() == 5);

    //std::vector<double> simplex_costs = {1.0, -1.0, 0.0, 2.0, -1.5};
    std::vector<double> simplex_costs = {1.0, -1.0, 0.0, 2.0, -0.5};
    bdds.set_costs(simplex_costs.begin(), simplex_costs.end());
    //bdds.backward_run(); 
    bdds.forward_run(); 

    const double backward_lb = bdds.lower_bound_backward_run();
    test(backward_lb == -1.0);

    {
        std::array<char,5> primal = {0,1,0,0,0};
        const double primal_cost = bdds.evaluate(primal.begin(), primal.end());
        test(primal_cost == -1.0);
    }

    {
        std::array<char,5> primal = {1,1,0,0,0};
        test(!bdds.check_feasibility(primal.begin(), primal.end()));
    }

    // test with random costs
    std::uniform_int_distribution<> d(-10,10);
    std::mt19937 gen;
    for(std::size_t i=0; i<100; ++i) {
        for(auto& x : simplex_costs)
            x = d(gen);
        bdds.set_costs(simplex_costs.begin(), simplex_costs.end());
        const double backward_lb = bdds.lower_bound_backward_run();
        test(backward_lb == *std::min_element(simplex_costs.begin(), simplex_costs.end()));
    }
}
