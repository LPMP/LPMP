#include "test.h"
#include "mrf/binary_MRF_instance.hxx"
#include "max_cut/max_cut_instance.hxx"
#include "max_cut/transform_binary_MRF.hxx"
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
            const binary_Potts_instance Potts = generate_random_binary_Potts_instance(n,m,rd);
            const max_cut_instance max_cut = transform_binary_Potts_to_max_cut(Potts);

            // generate random labeling
            for(std::size_t c=0; c<10; ++c) {
                binary_Potts_instance::labeling pl;
                for(std::size_t k=0; k<Potts.unaries.size(); ++k) {
                    pl.push_back(bd(gen));
                } 

                max_cut_edge_labeling mcl = transform_binary_Potts_labeling_to_max_cut(Potts,pl);

                test(std::abs(Potts.evaluate(pl) - max_cut.evaluate(mcl)) <= 1e-8);
            } 
        }
    } 
}
