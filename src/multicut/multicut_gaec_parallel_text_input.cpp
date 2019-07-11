#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "multicut/multicut_greedy_additive_edge_contraction.h"
#include "multicut/multicut_text_input.h"
#include "multicut/multicut_kernighan_lin.h"
#include <iostream>
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3)
        throw std::runtime_error("usage: [prog_name] [input_file] [nr_of_threads]");

    const multicut_instance input = multicut_text_input::parse_file(argv[1]);
    const int nr_of_threads = std::stoi(argv[2]);

    {
	const auto begin_time = std::chrono::steady_clock::now();
	const multicut_edge_labeling sol = greedy_additive_edge_contraction(input);
	const auto end_time = std::chrono::steady_clock::now();
	std::cout << "gaec energy = " << input.evaluate(sol) << "\n";
	std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }


    {
        const auto begin_time = std::chrono::steady_clock::now();
        const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(input, nr_of_threads);
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Parallel gaec energy = " << input.evaluate(sol) << "\n";
        std::cout << "Parallel optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    }

    {
        //const auto begin_time = std::chrono::steady_clock::now();
        //const multicut_edge_labeling andres_gaec = compute_gaec(input);
        //const auto end_time = std::chrono::steady_clock::now();
        //std::cout << "Andres gaec energy = " << input.evaluate(andres_gaec) << "\n";
        //std::cout << "Andres optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    } 

    //{
    //    const multicut_edge_labeling sol = greedy_additive_edge_contraction(input);
    //    const auto begin_time = std::chrono::steady_clock::now();
    //    const multicut_edge_labeling kl = compute_multicut_kernighan_lin(input, sol);
    //    std::cout << "K&L energy = " << input.evaluate(kl) << "\n";
    //    const auto end_time = std::chrono::steady_clock::now();
    //    std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
    //}

}
