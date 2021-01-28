#include "asymmetric_multiway_cut/asymmetric_multiway_cut_gaec.h"
#include "asymmetric_multiway_cut/asymmetric_multiway_cut_parser.h"
#include <iostream>

using namespace LPMP;

int main(int argc, char** argv)
{
    if(argc != 2)
        throw std::runtime_error("file name expected as argument");
    auto instance = asymmetric_multiway_cut_parser::parse_file(argv[1]);
    std::cout << "put multicut in normal form\n";
    instance.edge_costs.normalize();

    std::cout << "Compute labeling\n";
    const auto labeling = asymmetric_multiway_cut_gaec(instance);
    std::cout << "labeling energy = " << instance.evaluate(labeling) << "\n";
}
