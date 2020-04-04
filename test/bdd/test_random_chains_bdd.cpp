#include "test.h"
#include "bdd/bdd_anisotropic_diffusion.h"
#include "bdd/bdd_min_marginal_averaging.h"
#include "bdd/ILP_parser.h"
#include "mrf/mrf_input.h"
#include <random>
#include <sstream>
#include <fstream>

using namespace LPMP;

mrf_input generate_random_mrf_chain(const std::size_t n, const std::size_t l, const unsigned long int seed)
{
    mrf_input instance;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    std::cout << unif(gen) << "\n";
    
    std::vector<std::vector<double>> unaries;
    for(std::size_t i=0; i<n; ++i) {
        unaries.push_back({});
        for(std::size_t j=0; j<l; ++j) {
            unaries.back().push_back(unif(gen));
        } 
    }
    instance.unaries = two_dim_variable_array<double>(unaries);

    std::vector<std::array<std::size_t,2>> pairwise_size;
    for(std::size_t i=0; i+1<n; ++i) {
        pairwise_size.push_back({l,l});
    }
    three_dimensional_variable_array<double> binaries(pairwise_size.begin(), pairwise_size.end());

    for(std::size_t i=0; i+1<n; ++i) {
        for(std::size_t l1=0; l1<l; ++l1) {
            for(std::size_t l2=0; l2<l; ++l2) {
                binaries(i,l1,l2) = unif(gen);
            }
        }
    }

    instance.pairwise_values = binaries;

    for(std::size_t i=0; i+1<n; ++i) {
        instance.pairwise_indices.push_back({i, i+1});
    }

    return instance;
}

void test_randon_chains_problem(const std::size_t n, const std::size_t l, const unsigned long int seed)
{
    const mrf_input instance = generate_random_mrf_chain(n,l,seed);
    std::stringstream ss;
    instance.write(ss);
    //instance.write(std::cout);
    ILP_input input = ILP_parser::parse_string(ss.str());

    bdd_anisotropic_diffusion bdd_solver_bfs;
    bdd_anisotropic_diffusion_options bfs_options;
    bfs_options.order = bdd_anisotropic_diffusion_options::order::BFS;
    bdd_solver_bfs.set_options(bfs_options);
    bdd_solver_bfs.init(input);
    bdd_solver_bfs.compute_lower_bound();
    const double initial_lb_bfs = bdd_solver_bfs.lower_bound();
    bdd_solver_bfs.iteration();
    const double end_lb_bfs = bdd_solver_bfs.lower_bound();
    bdd_solver_bfs.iteration();
    const double end_lb_2_bfs = bdd_solver_bfs.lower_bound();

    bdd_anisotropic_diffusion bdd_solver_min_deg;
    bdd_anisotropic_diffusion_options min_deg_options;
    min_deg_options.order = bdd_anisotropic_diffusion_options::order::minimum_degree;
    bdd_solver_min_deg.set_options(min_deg_options);
    bdd_solver_min_deg.init(input);
    bdd_solver_min_deg.compute_lower_bound();
    const double initial_lb_min_deg = bdd_solver_min_deg.lower_bound();
    bdd_solver_min_deg.iteration();
    const double end_lb_min_deg = bdd_solver_min_deg.lower_bound();
    bdd_solver_min_deg.iteration();
    const double end_lb_2_min_deg = bdd_solver_min_deg.lower_bound();

    test(initial_lb_bfs <= end_lb_bfs);
    test(std::abs(end_lb_bfs - end_lb_2_bfs) <= 1e-8);

    test(initial_lb_min_deg <= end_lb_min_deg);
    test(std::abs(end_lb_min_deg - end_lb_2_min_deg) <= 1e-8);

    test(std::abs(end_lb_min_deg - end_lb_bfs) <= 1e-8);
}

void test_randon_chains_problem_mma(const std::size_t n, const std::size_t l, const unsigned long int seed)
{
    const mrf_input instance = generate_random_mrf_chain(n,l,seed);
    std::stringstream ss;
    instance.write(ss);
    // instance.write(std::cout);
    ILP_input input = ILP_parser::parse_string(ss.str());
    input.reorder_bfs();

    bdd_mma_fixing bdd_solver;
    bdd_min_marginal_averaging_options bdd_options;
    bdd_options.averaging_type = bdd_min_marginal_averaging_options::averaging_type::SRMP;
    bdd_solver.set_options(bdd_options);
    bdd_solver.init(input);
    const double initial_lb = bdd_solver.compute_lower_bound();
    double lb = initial_lb;
    double prev_lb = -std::numeric_limits<double>::infinity();
    size_t iterations = 0;
    while ((lb - prev_lb) > 1e-08)
    {
        iterations++;
        bdd_solver.iteration();
        prev_lb = lb;
        lb = bdd_solver.lower_bound();
    }
    double ub;

    std::cout << "\niterations = " << iterations << ", lower bound = " << lb << std::endl;
    if (bdd_solver.fix_variables())
    {
        ub = bdd_solver.compute_upper_bound();
        std::cout << "\nUpper bound = " << ub << std::endl;
    }
    else
        std::cout << "\nNo primal solution found." << std::endl;

    test(initial_lb <= lb);
    test(std::abs(ub - lb) <= 1e-8);
}

int main(int argc, char** argv)
{
    std::random_device rd;
    unsigned long int seed = 123;
    for(std::size_t n=2; n<20; ++n) {
        for(std::size_t l=2; l<20; ++l) {
            std::cout << "n = " << n << ", l = " << l << "\n";
            test_randon_chains_problem_mma(n,l,seed++);
        }
    }
}
