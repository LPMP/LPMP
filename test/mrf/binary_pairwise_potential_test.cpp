#include "test.h"
#include "mrf/binary_MRF_instance.hxx"
#include <random>
#include <cmath>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution nd(0.0, 1.0);

    for(std::size_t i=0; i<100; ++i) {
        binary_MRF_instance::binary_pairwise_potential p(0,1, {{ {nd(gen), nd(gen)}, {nd(gen), nd(gen)}}});
        auto [c, msg_1, msg_2] = p.make_potts();

        p.reparametrize(msg_1, msg_2);

        test(std::abs(p.cost[0][0]) < 1e-8);
        test(std::abs(p.cost[1][1]) < 1e-8);

        test(std::abs(p.cost[0][1] - c) < 1e-8);
        test(std::abs(p.cost[1][0] - c) < 1e-8);
    } 
}
