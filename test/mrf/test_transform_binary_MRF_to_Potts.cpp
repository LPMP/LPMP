#include "test.h"
#include "mrf/binary_MRF_instance.hxx"
#include "../generate_random_graph.hxx"
#include "mrf/transform_max_flow_instance.hxx"
#include <random>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::bernoulli_distribution bd(0.5);
    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=2*n; m<(n*(n-1))/2; m+=100) {
            const max_flow_instance mf = generate_random_qpbo_instance(n,m,rd);
            const binary_MRF_instance mrf = transform_QPBO_max_flow_to_binary_MRF(mf);

            const binary_Potts_instance mrf_Potts = transform_binary_MRF_to_Potts(mrf);

            binary_MRF_instance::labeling l;
            l.resize(n);
            for(auto& x : l)
                x = bd(gen);
            test(std::abs(mrf.evaluate(l) - mrf_Potts.evaluate(l)) <= 1e-8);
        }
    }
}
