#include "multigraph_matching_universe_frank_wolfe.cpp"
#include "graph_matching/multigraph_matching.hxx"
#include <iostream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("Expected filename as argument");

    auto input = Torresani_et_al_multigraph_matching_input::parse_file(argv[1]);

    multigraph_matching_frank_wolfe_universe s(input);
    s.solve();
}
