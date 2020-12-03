#include <random>
#include "graph_matching/matching_problem_input.h"
#include "multigraph_matching/multigraph_matching_input.h"
#include "multicut/transform_multigraph_matching.h"
#include "multicut/multicut_kernighan_lin.h"
#include "../graph_matching/generate_random_graph_matching_problem.h"
#include "test.h"

using namespace LPMP;

const std::string minimal_synchronization_example = 
R"(gm 0 1
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 0
a 2 1 0 0
a 3 1 1 -1

gm 0 2
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 0
a 2 1 0 0
a 3 1 1 -1

gm 1 2
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 0
a 2 1 0 0
a 3 1 1 -1
)";

void test_transformation(std::shared_ptr<multigraph_matching_input> mgm_instance)
{
    multigraph_matching_correlation_clustering_transform t(mgm_instance);
    auto cc_instance = t.get_correlatino_clustering_instance();
    auto mc_instance = cc_instance.transform_to_multicut();
    auto mc_sol = compute_multicut_kernighan_lin(mc_instance);
    test(mc_sol.check_primal_consistency(mc_instance));
    auto cc_sol = mc_sol.transform_to_correlation_clustering();
    auto mgm_sol = t.transform(cc_sol);
    test(mgm_sol.check_primal_consistency());
    test(std::abs(mgm_instance->evaluate(mgm_sol) - cc_instance.evaluate(cc_sol)) <= 1e-6); 
}

int main(int argc, char** argv)
{
    { // hand-crafted problem
        auto mgm_instance = std::make_shared<multigraph_matching_input>(Torresani_et_al_multigraph_matching_input::parse_string(minimal_synchronization_example));
        multigraph_matching_correlation_clustering_transform t(mgm_instance);
        auto cc_instance = t.get_correlatino_clustering_instance();;
        auto mc_instance = cc_instance.transform_to_multicut();
        auto mc_sol = compute_multicut_kernighan_lin(mc_instance);
        test(mc_sol.check_primal_consistency(mc_instance));
        auto cc_sol = mc_sol.transform_to_correlation_clustering();
        test(cc_instance.evaluate(cc_sol) == -6);
        auto mgm_sol = t.transform(cc_sol);
        test(mgm_sol.check_primal_consistency());
        test(mgm_instance->evaluate(mgm_sol) == -6);
    }

    // random problems
    std::random_device rd{};
    for(std::size_t no_graphs=3; no_graphs<=10; ++no_graphs) {
        for(std::size_t no_nodes=2; no_nodes<20; no_nodes+=3) {
            auto mgm_instance = std::make_shared<multigraph_matching_input>(generate_random_multigraph_matching_problem(no_graphs, no_nodes, rd));
         test_transformation(mgm_instance); 
      }
   }
}
