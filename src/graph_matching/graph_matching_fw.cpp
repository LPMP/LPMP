#include "graph_matching/graph_matching_input.h"
#include "graph_matching_frank_wolfe.cpp"

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("filename expected as argument");

    const graph_matching_input instance = TorresaniEtAlInput::parse_file(argv[1]);
    graph_matching_frank_wolfe gm_fw(instance);
    gm_fw.solve();
}
