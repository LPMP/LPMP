#include "test.h"
#include "mrf/max_flow_instance.hxx"
#include "mrf/binary_MRF_instance.hxx"
#include "mrf/transform_max_flow_instance.hxx"
#include "mrf/maxflow_factor.hxx"
#include "mrf/qpbo_factor.hxx"
#include "../generate_random_graph.hxx"
#include <random>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::bernoulli_distribution bd(0.5);

    // graph cut instances
    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=2*n; m<(n*(n-1))/2; m+=100) {
            const max_flow_instance mf = generate_random_graph_cut_instance(n,m,rd);
            const binary_Potts_instance mrf = transform_graph_cut_max_flow_to_binary_Potts(mf);

            maxflow_factor mf_factor(mf);
            const auto mf_lower_bound = mf_factor.LowerBound();

            qpbo_factor mrf_factor(mrf);
            const auto mrf_lower_bound = mrf_factor.LowerBound();

            test(std::abs(mf_lower_bound - mrf_lower_bound) <= 1e-8);
        }
    } 

    // qpbo instances
    // submodular
    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=0; m<(n*(n-1))/2; m+=100) {
            // TODO: balance unary and pairwise weights, so that non-trivial solutions are generated.
            const max_flow_instance mf = generate_random_qpbo_instance(n,m,rd,true);
            const binary_MRF_instance mrf = transform_QPBO_max_flow_to_binary_MRF(mf);
            const binary_Potts_instance mrf_Potts = transform_binary_MRF_to_Potts(mrf);

            for(const auto& p : mrf_Potts.pairwise_potentials) {
                test(p.cost >= 0.0, "potential must be submodular.");
            }

            maxflow_factor mf_factor(mf);
            const auto mf_lower_bound = mf_factor.LowerBound();

            maxflow_factor mf_factor_Potts(mrf_Potts);
            const auto mf_Potts_lower_bound = mf_factor_Potts.LowerBound();

            qpbo_factor mrf_factor(mrf);
            const auto mrf_lower_bound = mrf_factor.LowerBound();

            test(std::abs(mf_lower_bound - mf_Potts_lower_bound) <= 1e-8);
            test(std::abs(mf_lower_bound - mrf_lower_bound) <= 1e-8);
        }
    } 

    // general problems
    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=2*n; m<(n*(n-1))/2; m+=100) {
            const max_flow_instance mf = generate_random_qpbo_instance(n,m,rd);
            const binary_MRF_instance mrf = transform_QPBO_max_flow_to_binary_MRF(mf);

            maxflow_factor mf_factor(mf);
            const auto mf_lower_bound = mf_factor.LowerBound();

            qpbo_factor mrf_factor(mrf);
            const auto mrf_lower_bound = mrf_factor.LowerBound();

            test(std::abs(mf_lower_bound - mrf_lower_bound) <= 1e-8);
        }
    } 
}

