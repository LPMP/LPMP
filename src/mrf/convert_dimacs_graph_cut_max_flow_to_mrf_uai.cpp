#include "mrf/max_flow_instance.hxx"
#include "mrf/binary_MRF_instance.hxx"
#include "mrf/transform_max_flow_instance.hxx"
#include "mrf/dimacs_max_flow_input.h"
#include <string>
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3) { throw std::runtime_error("expected three arguments"); }
    const auto max_flow_input = dimacs_max_flow_input::parse_file(argv[1]);
    const auto binary_MRF_input = transform_graph_cut_max_flow_to_binary_Potts(max_flow_input);

    std::fstream s(argv[2], s.out);
    if(!s.is_open()) throw std::runtime_error(std::string("output file ") + argv[2] + " could not be opened");
    binary_MRF_input.write_uai(s); 
}
