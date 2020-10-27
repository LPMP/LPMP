#include "max_cut/max_cut_instance.hxx"
#include "max_cut/max_cut_greedy_additive_edge_contraction.h"
#include "max_cut/max_cut_local_search.h"
#include "max_cut/max_cut_text_input.h"
#include <iostream>
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("input file expected as argument");

    const max_cut_instance input = max_cut_text_input::parse_file(argv[1]);

    const auto begin_time = std::chrono::steady_clock::now();
    const max_cut_edge_labeling sol = greedy_additive_edge_contraction(input);
    std::cout << "gaec energy = " << input.evaluate(sol) << "\n";
    const auto end_time = std::chrono::steady_clock::now();
    std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

    const auto ls_begin_time = std::chrono::steady_clock::now();
    const auto sol_node_labeling = sol.transform_to_node_labeling(input);
    max_cut_local_search ls(input, sol_node_labeling);
    const double improvement = ls.perform_swaps();
    std::cout << "swap improvement = " << improvement  << "\n";
    const auto improved_sol = ls.get_labeling();
    std::cout << "energy after local search = " << input.evaluate(improved_sol) << "\n";
    const auto ls_end_time = std::chrono::steady_clock::now();
    std::cout << "local search took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(ls_end_time - ls_begin_time).count() << " milliseconds\n";
}

