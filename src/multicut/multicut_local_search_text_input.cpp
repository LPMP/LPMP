#include "multicut/multicut_instance.h"
#include "multicut/multicut_text_input.h"
#include "multicut/multicut_greedy_additive_edge_contraction.h"
#include "multicut/multicut_local_search.h"
#include "multicut/multicut_kernighan_lin.h"
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("Expected filename as argument");

    multicut_instance instance = multicut_text_input::parse_file(argv[1]);

    const auto begin_time = std::chrono::steady_clock::now();
    multicut_edge_labeling base_sol= greedy_additive_edge_contraction(instance);
    multicut_node_labeling base_sol_n = base_sol.transform_to_node_labeling(instance);
    const auto end_time = std::chrono::steady_clock::now();
    std::cout << "base objective = " << instance.evaluate(base_sol_n) << "\n";
    std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    
    {
        const auto begin_time = std::chrono::steady_clock::now();
        multicut_local_search ls(instance, base_sol_n);
        ls.perform_swaps();
        auto improved_sol = ls.get_labeling();
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "improved objective own = " << instance.evaluate(improved_sol) << "\n";
        std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }

    {
        const auto begin_time = std::chrono::steady_clock::now();
        auto improved_sol = compute_multicut_kernighan_lin(instance, base_sol);
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "improved objective Andres = " << instance.evaluate(improved_sol) << "\n";
        std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }
}
