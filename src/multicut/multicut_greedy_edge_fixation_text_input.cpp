#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_edge_fixation.h"
#include "multicut/multicut_text_input.h"
#include "multicut/multicut_kernighan_lin.h"
#include <iostream>
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("input file expected as argument");

    const multicut_instance input = multicut_text_input::parse_file(argv[1]);

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const multicut_edge_labeling andres_labeling = compute_multicut_greedy_edge_fixation(input);
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Andres greedy edge fixation energy = " << input.evaluate(andres_labeling) << "\n";
        std::cout << "Andres optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    } 

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const multicut_edge_labeling sol = multicut_greedy_edge_fixation(input);
        std::cout << "greedy edge fixation energy = " << input.evaluate(sol) << "\n";
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }
}
