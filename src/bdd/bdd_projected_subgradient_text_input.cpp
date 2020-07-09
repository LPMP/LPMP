#include "bdd/bdd_gradient_based_optimization.hxx"

using namespace LPMP;

int main(int argc, char** argv)
{
    std::cout << std::setprecision(10);

    if(argc < 2)
        throw std::runtime_error("input filename must be present as argument");

    ILP_input input(ILP_parser::parse_file(std::string(argv[1])));
    bdd_gradient_opt bdd_solver;
    bdd_solver.init(input);

    const double initial_lb = bdd_solver.compute_lower_bound();
    std::cout << "initial lower bound = " << initial_lb << "\n";

    for (std::size_t iter = 0; iter < 100; ++iter)
    {
        const auto g = bdd_solver.gradient();
        bdd_solver.apply_and_project(g);
        const double lb = bdd_solver.compute_lower_bound();
        std::cout << "iteration " << iter << ", lower bound = " << lb << "\n";
    }
}
