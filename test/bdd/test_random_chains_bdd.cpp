#include "test.h"
#include "bdd/bdd_anisotropic_diffusion.h"
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
    bdd_anisotropic_diffusion bdd_solver;
    bdd_solver.init(input);
    const double initial_lb = bdd_solver.lower_bound();
    bdd_solver.iteration();
    const double end_lb = bdd_solver.lower_bound();
    bdd_solver.iteration();
    const double end_lb_2 = bdd_solver.lower_bound();
    test(initial_lb <= end_lb);
    test(std::abs(end_lb - end_lb_2) <= 1e-8);
}

int main(int argc, char** argv)
{
    std::random_device rd;
    unsigned long int seed = 123;
    for(std::size_t n=2; n<20; ++n) {
        for(std::size_t l=2; l<20; ++l) {
            test_randon_chains_problem(n,l,seed++);
        }
    }
}
