#include "asymmetric_multiway_cut/asymmetric_multiway_cut_parser.h"
#include "test.h"

using namespace LPMP; 

const std::string instance_3x3x3 = 
R"(ASYMMETRIC MULTIWAY CUT
MULTICUT
0 1 1.0
0 2 -1.0
1 2 -1.0
NODE COSTS
0 1 2.0
0 1 2.0
0 1 2.0)";

int main(int argc, char** argv)
{
    const asymmetric_multiway_cut_instance instance = asymmetric_multiway_cut_parser::parse_string(instance_3x3x3);

    test(instance.nr_nodes() == 3);
    test(instance.nr_edges() == 3);
    test(instance.nr_labels() == 3);

}

