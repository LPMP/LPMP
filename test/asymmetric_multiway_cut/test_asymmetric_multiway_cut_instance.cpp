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

    asymmetric_multiway_cut_labeling labeling;
    labeling.edge_labels.push_back(1);
    labeling.edge_labels.push_back(1);
    labeling.edge_labels.push_back(1);

    labeling.node_labels.push_back(0);
    labeling.node_labels.push_back(0);
    labeling.node_labels.push_back(0);

    test(instance.feasible(labeling));
    test(instance.evaluate(labeling) == -1.0);
}


