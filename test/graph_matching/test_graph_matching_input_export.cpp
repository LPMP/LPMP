#include "test.h"
#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching.h"
#include "graph_matching/graph_matching_constructor.hxx"
#include "generate_random_graph_matching_problem.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;

void test_instance_export(const std::size_t no_nodes, const double density, std::random_device& rd, int argc, char** argv)
{
    const graph_matching_input instance = generate_random_graph_matching_problem(no_nodes, density, rd);
    ProblemConstructorRoundingSolver<Solver<LP<FMC_MP>,StandardVisitor>>solver(argc,argv);
    solver.GetProblemConstructor().construct(instance);
    const graph_matching_input exported_instance = solver.GetProblemConstructor().export_graph_matching_input();

    solver.GetProblemConstructor().ComputePrimal();
    const graph_matching_input::labeling l = solver.GetProblemConstructor().write_out_labeling();
    test(std::abs(instance.evaluate(l) - exported_instance.evaluate(l)) <= 1e-8); 

    solver.Solve();
    const graph_matching_input::labeling l_sol = solver.GetProblemConstructor().write_out_labeling();
    const graph_matching_input solved_instance = solver.GetProblemConstructor().export_graph_matching_input();
    test(std::abs(instance.evaluate(l_sol) - exported_instance.evaluate(l_sol)) <= 1e-8); 
    solver.GetProblemConstructor().read_in_labeling(l_sol);
    test(std::abs(instance.evaluate(l_sol) - solver.EvaluatePrimal()) <= 1e-8);
    test(std::abs(instance.evaluate(l_sol) - solved_instance.evaluate(l_sol)) <= 1e-8); 
    test(std::abs(instance.evaluate(l) - exported_instance.evaluate(l)) <= 1e-8); 
    test(std::abs(instance.evaluate(l) - solved_instance.evaluate(l)) <= 1e-8); 
}

int main(int argc, char** argv)
{
    std::random_device rd{};

    for(std::size_t i=10; i<20; ++i)
        test_instance_export(i,-1.0, rd, argc, argv);

    for(std::size_t i=10; i<20; ++i)
        test_instance_export(i,0.5, rd, argc, argv);

    for(std::size_t i=10; i<20; ++i)
        test_instance_export(i,2.0, rd, argc, argv);
}
