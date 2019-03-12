#include "test.h"
#include "max_cut/max_cut_instance.hxx"
#include "../generate_random_graph.hxx"
#include <random>

using namespace LPMP;

int main(int argc, char** argv)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::bernoulli_distribution bd(0.5);

    for(std::size_t n=10; n<100; n+=10) {
        for(std::size_t m=2*n; m<(n*(n-1))/2; m+=100) {
            const max_cut_instance instance = generate_random_max_cut_instance(n, m, rd);

            // generate random labeling
            for(std::size_t c=0; c<10; ++c) {
                max_cut_node_labeling l;
                for(std::size_t i=0; i<n; ++i)
                    l.push_back(bd(gen));

                max_cut_edge_labeling el(instance, l);
                test(el.check_primal_consistency(instance));
                test(std::abs(instance.evaluate(el) - instance.evaluate(l)) <= 1e-8);
                max_cut_node_labeling ll = el.transform_to_node_labeling(instance);
                test(std::abs(instance.evaluate(el) - instance.evaluate(ll)) <= 1e-8);
            } 
        }
    } 
}
