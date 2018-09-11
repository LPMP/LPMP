#include "test.h"
#include "mrf/mrf_uai_input.h"
#include "mrf/binary_MRF_instance.hxx"
#include "mrf/qpbo_factor.hxx"
#include "../generate_random_graph.hxx"
#include <random>
#include <sstream>
#include <iostream>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::random_device rd{};

    // general pairwise potentials
    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=0; m<(n*(n-1))/2; m+=100) {
            const auto input = generate_random_binary_MRF_instance(n,m,rd);
            qpbo_factor mrf_factor(input);
            const auto mrf_lower_bound = mrf_factor.LowerBound();

            // write uai model into string, read from it afterwards and check if gives the same optimal value as original problem.
            std::ostringstream ss;
            input.write_uai(ss);
            const std::string uai_input = ss.str();
            //std::cout << uai_input << "\n\n";

            const auto read_input = mrf_uai_input::parse_string(uai_input);
            qpbo_factor mrf_factor_2(read_input);
            const auto mrf_lower_bound_2 = mrf_factor_2.LowerBound();

            test(std::abs(mrf_lower_bound - mrf_lower_bound) < 1e-8);
        }
    }

    // Potts pairwise potentials
    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=0; m<(n*(n-1))/2; m+=100) {
            const auto input = generate_random_binary_Potts_instance(n,m,rd);
            qpbo_factor mrf_factor(input);
            const auto mrf_lower_bound = mrf_factor.LowerBound();

            // write uai model into string, read from it afterwards and check if gives the same optimal value as original problem.
            std::ostringstream ss;
            input.write_uai(ss);
            const std::string uai_input = ss.str();
            //std::cout << uai_input << "\n\n";

            const auto read_input = mrf_uai_input::parse_string(uai_input);
            qpbo_factor mrf_factor_2(read_input);
            const auto mrf_lower_bound_2 = mrf_factor_2.LowerBound();

            test(std::abs(mrf_lower_bound - mrf_lower_bound) < 1e-8);
        }
    }
} 
