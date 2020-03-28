#include "max_cut/max_cut_instance.hxx"
#include "max_cut/max_cut_sahni_gonzalez.h"
#include "max_cut/max_cut_text_input.h"
#include "max_cut/max_cut_local_search.h"
#include <iostream>
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("input file expected as argument");

    const max_cut_instance input = max_cut_text_input::parse_file(argv[1]);

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const max_cut_node_labeling sol = max_cut_sahni_gonzalez_1(input);
        std::cout << "sahni gonzalez 1 energy = " << input.evaluate(sol) << "\n";
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        const auto ls_begin_time = std::chrono::steady_clock::now();
        max_cut_local_search ls(input, sol);
        std::cout << "swap improvement = " << ls.perform_swaps() << "\n";
        const auto improved_sol = ls.get_labeling();
        std::cout << "energy after local search = " << input.evaluate(improved_sol) << "\n";
        const auto ls_end_time = std::chrono::steady_clock::now();
        std::cout << "local search took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(ls_end_time - ls_begin_time).count() << " milliseconds\n";
    }

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const max_cut_node_labeling sol = max_cut_sahni_gonzalez_2(input);
        std::cout << "sahni gonzalez 2 energy = " << input.evaluate(sol) << "\n";
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        const auto ls_begin_time = std::chrono::steady_clock::now();
        max_cut_local_search ls(input, sol);
        std::cout << "swap improvement = " << ls.perform_swaps() << "\n";
        const auto improved_sol = ls.get_labeling();
        std::cout << "energy after local search = " << input.evaluate(improved_sol) << "\n";
        const auto ls_end_time = std::chrono::steady_clock::now();
        std::cout << "local search took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(ls_end_time - ls_begin_time).count() << " milliseconds\n";
    }

    {
        const auto begin_time = std::chrono::steady_clock::now();
        const max_cut_node_labeling sol = max_cut_sahni_gonzalez_3(input);
        std::cout << "sahni gonzalez 3 energy = " << input.evaluate(sol) << "\n";
        const auto end_time = std::chrono::steady_clock::now();
        std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";

        const auto ls_begin_time = std::chrono::steady_clock::now();
        max_cut_local_search ls(input, sol);
        std::cout << "swap improvement = " << ls.perform_swaps() << "\n";
        const auto improved_sol = ls.get_labeling();
        std::cout << "energy after local search = " << input.evaluate(improved_sol) << "\n";
        const auto ls_end_time = std::chrono::steady_clock::now();
        std::cout << "local search took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(ls_end_time - ls_begin_time).count() << " milliseconds\n";
    }

}
