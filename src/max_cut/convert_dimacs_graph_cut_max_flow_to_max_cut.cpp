#include "mrf/max_flow_instance.hxx"
#include "mrf/binary_MRF_instance.hxx"
#include "mrf/transform_max_flow_instance.hxx"
#include "mrf/dimacs_max_flow_input.h"
#include "max_cut/max_cut_instance.hxx"
#include "max_cut/transform_binary_MRF.hxx"
#include <fstream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 3) { throw std::runtime_error("expected three arguments"); }
    const auto max_flow_input = dimacs_max_flow_input::parse_file(argv[1]);
    const auto binary_Potts_input = transform_graph_cut_max_flow_to_binary_Potts(max_flow_input);
    const auto max_cut_input = transform_binary_Potts_to_max_cut(binary_Potts_input);

    std::fstream s(argv[2], s.out);
    if(!s.is_open()) throw std::runtime_error(std::string("output file ") + argv[2] + " could not be opened");
    max_cut_input.write(s); 
}

