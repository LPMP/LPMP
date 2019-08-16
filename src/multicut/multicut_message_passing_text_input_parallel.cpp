#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "multicut/multicut_text_input.h"
#include "multicut/multicut_message_passing_parallel.h"
#include <iostream>
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3)
        throw std::runtime_error("[prog_name] [input_file] [nr_threads]");
    const int nr_thread = std::stoi(argv[2]);
    const multicut_instance input = multicut_text_input::parse_file(argv[1]);

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(input, nr_thread, "non-blocking");
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Parallel gaec energy = " << input.evaluate(sol) << "\n";
        std::cout << "Parallel optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const multicut_edge_labeling sol = multicut_message_passing_parallel(input, false, nr_thread);
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Parallel CP + GAEC energy = " << input.evaluate(sol) << "\n";
        std::cout << "CP + GAEC took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }
    return 0;

}