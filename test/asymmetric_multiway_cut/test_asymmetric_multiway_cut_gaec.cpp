#include "asymmetric_multiway_cut/asymmetric_multiway_cut_parser.h"
#include "asymmetric_multiway_cut/asymmetric_multiway_cut_gaec.h"
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

    const asymmetric_multiway_cut_labeling labeling = asymmetric_multiway_cut_gaec(instance);

    test(instance.feasible(labeling));
    test(instance.evaluate(labeling) == -2.0);

}
