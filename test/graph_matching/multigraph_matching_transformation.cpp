#include <random>
#include "graph_matching/matching_problem_input.h"
#include "graph_matching/multigraph_matching_input.h"
#include "multicut/transform_multigraph_matching.h"
#include "multicut/multicut_kernighan_lin.h"
#include "generate_random_graph_matching_problem.h"
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

void test_transformation(const multigraph_matching_input& mgm_input)
{
    auto cc_input = transform_multigraph_matching_to_correlation_clustering(mgm_input);
    auto cc_sol = compute_multicut_kernighan_lin(cc_input);
    test(cc_sol.check_primal_consistency(cc_input));
    auto mgm_sol = transform_correlation_clustering_to_multigraph_matching(mgm_input, cc_input, cc_sol);
    test(mgm_sol.check_primal_consistency());
    test(std::abs(mgm_input.evaluate(mgm_sol) - cc_input.evaluate(cc_sol)) <= 1e-6); 
}

int main(int argc, char** argv)
{
   { // hand-crafted problem
    auto mgm_input = Torresani_et_al_multigraph_matching_input::parse_string(minimal_synchronization_example);
    auto cc_input = transform_multigraph_matching_to_correlation_clustering(mgm_input);
    auto cc_sol = compute_multicut_kernighan_lin(cc_input);
    test(cc_sol.check_primal_consistency(cc_input));
    test(cc_input.evaluate(cc_sol) == -6);
    auto mgm_sol = transform_correlation_clustering_to_multigraph_matching(mgm_input, cc_input, cc_sol);
    test(mgm_sol.check_primal_consistency());
    test(mgm_input.evaluate(mgm_sol) == -6);
   }

   // random problems
   std::random_device rd{};
   for(std::size_t no_graphs=3; no_graphs<=10; ++no_graphs) {
      for(std::size_t no_nodes=2; no_nodes<20; no_nodes+=3) {
         auto mgm_input = generate_random_multigraph_matching_problem(no_graphs, no_nodes, rd);
         test_transformation(mgm_input); 
      }
   }
}
