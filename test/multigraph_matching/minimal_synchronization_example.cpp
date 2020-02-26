#include "graph_matching/matching_problem_input.h"
#include "multigraph_matching/multigraph_matching.hxx"
#include "multicut/transform_multigraph_matching.h"
#include "visitors/standard_visitor.hxx"
#include "solver.hxx"
#include <string>
#include <vector>
#include "test.h"

using namespace LPMP;

const std::string minimal_synchronization_example = 
R"(gm 0 1
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 -10
a 2 1 0 -10
a 3 1 1 -1

gm 0 2
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 -10
a 2 1 0 -10
a 3 1 1 -1

gm 1 2
p 2 2 0 0
a 0 0 0 -1
a 1 0 1 -10
a 2 1 0 -10
a 3 1 1 -1
)";

const std::vector<std::string> options = 
{
"",
"--standardReparametrization", "anisotropic",
"--roundingReparametrization", "uniform:0.5",
"--tightenIteration", "10",
"--tightenInterval", "5",
"--tightenReparametrization", "uniform:0.5",
"--tightenConstraintsMax", "5",
"--maxIter", "100",
"--tighten",
"-v", "2"
};

int main(int argc, char** argv)
{
    ProblemConstructorRoundingSolver<Solver<LP<FMC_MGM<true>>,StandardTighteningVisitor>> solver(options);
    auto mgm_instance = std::make_shared<multigraph_matching_input>(Torresani_et_al_multigraph_matching_input::parse_string(minimal_synchronization_example));
    auto& mgm_constructor = solver.GetProblemConstructor();
    mgm_constructor.construct(*mgm_instance);
    solver.Solve();
    std::cout << solver.GetLP().number_of_messages() << "\n";
    test( std::abs(-42 - solver.lower_bound()) <= 1e-6 ); 
    solver.get_primal().write_primal_matching(std::cout);

    mgm_constructor.WritePrimal(std::cout);
    mgm_constructor.ComputePrimal();
    const multigraph_matching_input::labeling mgm_sol = solver.get_primal();
    mgm_sol.write_primal_matching(std::cout);
    std::cout << solver.lower_bound() << " = " << mgm_instance->evaluate(mgm_sol) << "\n";
    test( std::abs(solver.lower_bound() - mgm_instance->evaluate(mgm_sol)) <= eps );

    // test transformation to correlation clustering
    multigraph_matching_correlation_clustering_transform mgm_cc_trafo(mgm_instance);
    const auto& cc = mgm_cc_trafo.get_correlatino_clustering_instance();
    auto cc_sol = mgm_cc_trafo.transform(mgm_sol);
    test(std::abs(cc.evaluate(cc_sol) - mgm_instance->evaluate(mgm_sol)) <= eps);
}
