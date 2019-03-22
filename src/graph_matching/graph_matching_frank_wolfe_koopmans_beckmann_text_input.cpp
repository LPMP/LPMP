#include "graph_matching/graph_matching_instance_koopmans_beckmann.h"
#include "graph_matching/graph_matching_koopmans_beckmann_input.h"
#include "graph_matching_frank_wolfe_koopmans_beckmann.cpp"

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("Expected filename as argument");

    const auto instance = graph_matching_koopmans_beckmann_input::parse_file(std::string(argv[1]));

    std::cout << instance.L << "\n";

    graph_matching_frank_wolfe_koopmans_beckmann s(instance);
    s.solve(); 
}
