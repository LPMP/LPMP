#include "max_cut/max_cut.h"
#include "solver.hxx"
#include "test.h"
#include "visitors/standard_visitor.hxx"
#include <random>

using namespace LPMP;

int main(int argc, char** argv)
{
    Solver<LP<FMC_ODD_BICYCLE_WHEEL_MAX_CUT>,StandardVisitor> s;
    auto& mc = s.GetProblemConstructor();
    mc.add_edge_factor(0, 1, 0.0);
    mc.add_edge_factor(0, 2, 0.0);
    mc.add_edge_factor(0, 3, 0.0);
    mc.add_edge_factor(0, 4, 0.0);
    mc.add_edge_factor(1, 2, 0.0);
    mc.add_edge_factor(1, 3, 0.0);
    mc.add_edge_factor(1, 4, 0.0);
    mc.add_edge_factor(2, 3, 0.0);
    mc.add_edge_factor(2, 4, 0.0);
    mc.add_edge_factor(3, 4, 0.0);
    
    mc.add_triplet_factor(0,1,2);
    mc.add_triplet_factor(0,1,3);
    mc.add_triplet_factor(0,1,4);
    mc.add_triplet_factor(0,2,3);
    mc.add_triplet_factor(0,2,4);
    mc.add_triplet_factor(0,3,4);
    mc.add_triplet_factor(1,2,3);
    mc.add_triplet_factor(1,2,4);
    mc.add_triplet_factor(1,3,4);
    mc.add_triplet_factor(2,3,4);

    auto set_labeling_and_test = [&](const std::array<char,10> labeling) {
        auto* f01 = mc.get_edge_factor(0,1);
        auto* f02 = mc.get_edge_factor(0,2);
        auto* f03 = mc.get_edge_factor(0,3);
        auto* f04 = mc.get_edge_factor(0,4);
        auto* f12 = mc.get_edge_factor(1,2);
        auto* f13 = mc.get_edge_factor(1,3);
        auto* f14 = mc.get_edge_factor(1,4);
        auto* f23 = mc.get_edge_factor(2,3);
        auto* f24 = mc.get_edge_factor(2,4);
        auto* f34 = mc.get_edge_factor(3,4);

        // edges of quintuplet: 01,02,12, 03,04,13,14,23,24, 34
        f01->get_factor()->primal() = labeling[0];
        f02->get_factor()->primal() = labeling[1];
        f12->get_factor()->primal() = labeling[2];
        f03->get_factor()->primal() = labeling[3];
        f04->get_factor()->primal() = labeling[4];
        f13->get_factor()->primal() = labeling[5];
        f14->get_factor()->primal() = labeling[6];
        f23->get_factor()->primal() = labeling[7];
        f24->get_factor()->primal() = labeling[8];
        f34->get_factor()->primal() = labeling[9];

        std::cout << f01->get_factor()->primal() << ",";
        std::cout << f02->get_factor()->primal() << ",";
        std::cout << f12->get_factor()->primal() << ",";
        std::cout << f03->get_factor()->primal() << ",";
        std::cout << f04->get_factor()->primal() << ",";
        std::cout << f13->get_factor()->primal() << ",";
        std::cout << f14->get_factor()->primal() << ",";
        std::cout << f23->get_factor()->primal() << ",";
        std::cout << f24->get_factor()->primal() << ",";
        std::cout << f34->get_factor()->primal() << "\n";

        f01->propagate_primal_through_messages();
        f02->propagate_primal_through_messages();
        f03->propagate_primal_through_messages();
        f04->propagate_primal_through_messages();
        f12->propagate_primal_through_messages();
        f13->propagate_primal_through_messages();
        f14->propagate_primal_through_messages();
        f23->propagate_primal_through_messages();
        f24->propagate_primal_through_messages();
        f34->propagate_primal_through_messages();

        std::cout << f01->get_factor()->primal() << ",";
        std::cout << f02->get_factor()->primal() << ",";
        std::cout << f12->get_factor()->primal() << ",";
        std::cout << f03->get_factor()->primal() << ",";
        std::cout << f04->get_factor()->primal() << ",";
        std::cout << f13->get_factor()->primal() << ",";
        std::cout << f14->get_factor()->primal() << ",";
        std::cout << f23->get_factor()->primal() << ",";
        std::cout << f24->get_factor()->primal() << ",";
        std::cout << f34->get_factor()->primal() << "\n"; 

        auto* f012 = mc.get_triplet_factor(0,1,2);
        auto* f013 = mc.get_triplet_factor(0,1,3);
        auto* f014 = mc.get_triplet_factor(0,1,4);
        auto* f023 = mc.get_triplet_factor(0,2,4);
        auto* f024 = mc.get_triplet_factor(0,3,4);
        auto* f034 = mc.get_triplet_factor(1,2,3);
        auto* f123 = mc.get_triplet_factor(1,2,4);
        auto* f124 = mc.get_triplet_factor(1,2,4);
        auto* f134 = mc.get_triplet_factor(1,3,4);
        auto* f234 = mc.get_triplet_factor(2,3,4);
        std::cout << "012: " << f012->get_factor()->primal() << "\n";
        std::cout << "013: " << f013->get_factor()->primal() << "\n";
        std::cout << "014: " << f014->get_factor()->primal() << "\n";
        std::cout << "023: " << f023->get_factor()->primal() << "\n";
        std::cout << "024: " << f024->get_factor()->primal() << "\n";
        std::cout << "034: " << f034->get_factor()->primal() << "\n";
        std::cout << "123: " << f123->get_factor()->primal() << "\n";
        std::cout << "124: " << f124->get_factor()->primal() << "\n";
        std::cout << "134: " << f134->get_factor()->primal() << "\n";
        std::cout << "234: " << f234->get_factor()->primal() << "\n";
        test(std::abs(s.GetLP().LowerBound()) < 1e-8);
        test(s.GetLP().CheckPrimalConsistency());
        test(mc.CheckPrimalConsistency());
        test(std::abs(s.GetLP().EvaluatePrimal()) < 1e-8);
    };

    // labelings as in max_cut_quintuplet_factor
    set_labeling_and_test({0,0,0,0,0,0,0,0,0,0});
    // one node, four nodes
    set_labeling_and_test({1,1,0,1,1,0,0,0,0,0}); // 10000
    set_labeling_and_test({1,0,1,0,0,1,1,0,0,0}); // 01000
    set_labeling_and_test({0,1,1,0,0,0,0,1,1,0}); // 00100
    set_labeling_and_test({0,0,0,1,0,1,0,1,0,1}); // 00010
    set_labeling_and_test({0,0,0,0,1,0,1,0,1,1}); // 00001
    // two nodes, three nodes
    set_labeling_and_test({0,1,1,1,1,1,1,0,0,0}); // 11000
    set_labeling_and_test({1,0,1,1,1,0,0,1,1,0}); // 10100
    set_labeling_and_test({1,1,0,0,1,1,0,1,0,1}); // 10010
    set_labeling_and_test({1,1,0,1,0,0,1,0,1,1}); // 10001
    set_labeling_and_test({1,1,0,0,0,1,1,1,1,0}); // 01100
    set_labeling_and_test({1,0,1,1,0,0,1,1,0,1}); // 01010
    set_labeling_and_test({1,0,1,0,1,1,0,0,1,1}); // 01001
    set_labeling_and_test({0,1,1,1,0,1,0,0,1,1}); // 00110
    set_labeling_and_test({0,1,1,0,1,0,1,1,0,1}); // 00101
    set_labeling_and_test({0,0,0,1,1,1,1,1,1,0}); // 00011

    mc.add_quintuplet_factor({0,1,2,3,4});

    // labelings as in max_cut_quintuplet_factor
    set_labeling_and_test({0,0,0,0,0,0,0,0,0,0});
    // one node, four nodes
    set_labeling_and_test({1,1,0,1,1,0,0,0,0,0}); // 10000
    set_labeling_and_test({1,0,1,0,0,1,1,0,0,0}); // 01000
    set_labeling_and_test({0,1,1,0,0,0,0,1,1,0}); // 00100
    set_labeling_and_test({0,0,0,1,0,1,0,1,0,1}); // 00010
    set_labeling_and_test({0,0,0,0,1,0,1,0,1,1}); // 00001
    // two nodes, three nodes
    set_labeling_and_test({0,1,1,1,1,1,1,0,0,0}); // 11000
    set_labeling_and_test({1,0,1,1,1,0,0,1,1,0}); // 10100
    set_labeling_and_test({1,1,0,0,1,1,0,1,0,1}); // 10010
    set_labeling_and_test({1,1,0,1,0,0,1,0,1,1}); // 10001
    set_labeling_and_test({1,1,0,0,0,1,1,1,1,0}); // 01100
    set_labeling_and_test({1,0,1,1,0,0,1,1,0,1}); // 01010
    set_labeling_and_test({1,0,1,0,1,1,0,0,1,1}); // 01001
    set_labeling_and_test({0,1,1,1,0,1,0,0,1,1}); // 00110
    set_labeling_and_test({0,1,1,0,1,0,1,1,0,1}); // 00101
    set_labeling_and_test({0,0,0,1,1,1,1,1,1,0}); // 00011
}
