#include "graph_matching/graph_matching_input.h"
#include "graph_matching_frank_wolfe.cpp"
#include <chrono>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("filename expected as argument");

    const graph_matching_input instance = TorresaniEtAlInput::parse_file(argv[1]);

    const auto begin_time = std::chrono::steady_clock::now();
    graph_matching_frank_wolfe_sparse gm_fw(instance);
    gm_fw.solve();
    const auto end_time = std::chrono::steady_clock::now();
    std::cout << "Optimization took " <<  std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " milliseconds\n";
}

