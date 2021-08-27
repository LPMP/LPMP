#include "multicut/multicut_instance.h"
#include "multicut/multicut_balanced_edge_contraction.h"
#include "multicut/multicut_text_input.h"
#include <iostream>
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("input file expected as argument");

    const multicut_instance input = multicut_text_input::parse_file(argv[1]);

    const auto begin_time = std::chrono::steady_clock::now();
    multicut_edge_labeling sol = multicut_balanced_edge_contraction(input);
    const auto end_time = std::chrono::steady_clock::now();
    std::cout << "bec energy = " << input.evaluate(sol) << "\n";
    std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
}
