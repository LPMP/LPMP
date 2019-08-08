#include "multicut/multicut_instance.h"
#include "multicut/multicut_greedy_additive_edge_contraction_parallel.h"
#include "multicut/multicut_greedy_additive_edge_contraction.h"
#include "multicut/multicut_text_input.h"
#include "multicut/multicut_kernighan_lin.h"
#include <iostream>
#include <chrono>
#include <tclap/CmdLine.h>

using namespace LPMP;

int main(int argc, char** argv)
{
    try {
        TCLAP::CmdLine cmd("GAEC Parallel");
        //TCLAP: command line library. -i ${input file} --numThreads ${x} --edgeDistribution ${...}
        std::vector<std::string> allowed;
		allowed.push_back("round_robin_not_sorted");
		allowed.push_back("endnodes");
		allowed.push_back("chunk_sorted");
		allowed.push_back("round_robin_sorted");
		allowed.push_back("chunk_not_sorted");
        allowed.push_back("non-blocking");
		TCLAP::ValuesConstraint<std::string> allowedVals(allowed);

        TCLAP::ValueArg<std::string> nameArg("i","inputFile","Path to the input file.",true,"","string");
        TCLAP::ValueArg<std::string> threadArg("t","numThreads","Number of threads.",true,"1","int");
        TCLAP::ValueArg<std::string> optArg("x","edgeDistribution","Methods to distribute edges.",true,"round_robin_sorted", &allowedVals);

        cmd.add(nameArg);
        cmd.add(threadArg);
        cmd.add(optArg);
        cmd.parse(argc, argv);

        const std::string filename = nameArg.getValue();
        const multicut_instance input = multicut_text_input::parse_file(filename);
        const int nr_of_threads = std::stoi(threadArg.getValue());
        const std::string option = optArg.getValue();

        // {
        // 	const auto begin_time = std::chrono::steady_clock::now();
        // 	const multicut_edge_labeling sol = greedy_additive_edge_contraction(input);
        // 	const auto end_time = std::chrono::steady_clock::now();
        // 	std::cout << "gaec energy = " << input.evaluate(sol) << "\n";
        // 	std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
        // }


        {
            const auto begin_time = std::chrono::steady_clock::now();
            const multicut_edge_labeling sol = greedy_additive_edge_contraction_parallel(input, nr_of_threads, option);
            const auto end_time = std::chrono::steady_clock::now();
            std::cout << "Parallel gaec energy = " << input.evaluate(sol) << "\n";
            std::cout << "Parallel optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
        }

    } catch (TCLAP::ArgException &e) { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }


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
