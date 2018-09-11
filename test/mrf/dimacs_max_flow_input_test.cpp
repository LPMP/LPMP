#include "test.h"
#include "mrf/dimacs_max_flow_input.h"
#include "mrf/qpbo_factor.hxx"
#include "mrf/maxflow_factor.hxx"
#include "mrf/transform_max_flow_instance.hxx"
#include <cmath>
#include <iostream>

using namespace LPMP;


/*
1---------
|         |
3----4    5--
      \__   |
6----7    8 |
      \_____|
2----


*/
const std::string test_dimacs_qpbo_max_flow =
R"(p max 8 11
n 1 s
n 2 t
c unary potentials
a 1 3 1
a 6 2 1
a 4 2 1.5
a 1 7 1.5
a 1 5 2
a 8 2 2
c pairwise potentials
a 3 4 3
a 4 3 3
a 6 7 3
a 7 6 3
a 4 8 3
a 8 4 3
a 7 5 3
a 5 7 3
)";

const std::string test_dimacs_qpbo_max_flow_non_symmetric = 
R"(p max 8 11
n 1 s
n 2 t
c unary potentials
a 1 3 1
a 6 2 1
a 4 2 1.5
a 1 7 1.5
a 1 5 2
a 8 2 2
c pairwise potentials
a 3 4 1.5
a 4 3 3
a 6 7 3
a 7 6 1.5
a 4 8 3.5
a 8 4 3
a 7 5 3.5
a 5 7 3
)";
;

const std::string test_dimacs_graph_cut_max_flow =
R"(p max 5 7
n 1 s
n 2 t
a 1 3 2
a 1 5 3
a 4 2 4
a 3 4 3
a 4 3 3
a 4 5 3
a 5 4 3
)";



int main(int argc, char** argv)
{
    // submodular
    {
        auto input = dimacs_max_flow_input::parse_string(test_dimacs_graph_cut_max_flow);
        test(input.no_nodes = 5);
        test(input.arcs.size() == 7); 

        test(input.arcs[0][0] == 0 && input.arcs[0][1] == 2 && input.arcs[0].capacity == 2);
        test(input.arcs[1][0] == 0 && input.arcs[1][1] == 4 && input.arcs[1].capacity == 3);
        test(input.arcs[2][0] == 3 && input.arcs[2][1] == 1 && input.arcs[2].capacity == 4);

        test(input.arcs[3][0] == 2 && input.arcs[3][1] == 3 && input.arcs[3].capacity == 3);
        test(input.arcs[4][0] == 3 && input.arcs[4][1] == 2 && input.arcs[4].capacity == 3);
        test(input.arcs[5][0] == 3 && input.arcs[5][1] == 4 && input.arcs[5].capacity == 3);
        test(input.arcs[6][0] == 4 && input.arcs[6][1] == 3 && input.arcs[6].capacity == 3);

        maxflow_factor m(input);
        const auto maxflow_bound = m.LowerBound();

        const auto mrf = transform_graph_cut_max_flow_to_binary_Potts(input);
        qpbo_factor f(mrf);
        const auto mrf_bound = f.LowerBound();

        test(mrf.unaries.size() == 3);
        test(mrf.unaries[0][0] == 2 && mrf.unaries[0][1] == 0);
        test(mrf.unaries[1][0] == 0 && mrf.unaries[1][1] == 4);
        test(mrf.unaries[2][0] == 3 && mrf.unaries[2][1] == 0);

        test(mrf.pairwise_potentials.size() == 2);
        test(mrf.pairwise_potentials[0][0] == 0 && mrf.pairwise_potentials[0][1] == 1 && mrf.pairwise_potentials[0].cost == 3); 
        test(mrf.pairwise_potentials[1][0] == 1 && mrf.pairwise_potentials[1][1] == 2 && mrf.pairwise_potentials[0].cost == 3); 

        test(std::abs(maxflow_bound - mrf_bound) <= 1e-8);
    }

    // qpbo, nonsubmodular, but tight
    {
        auto input = dimacs_max_flow_input::parse_string(test_dimacs_qpbo_max_flow);
        test(input.no_nodes = 8);
        test(input.arcs.size() == 14); 

        // unaries
        test(input.arcs[0][0] == 0 && input.arcs[0][1] == 2 && input.arcs[0].capacity == 1);
        test(input.arcs[1][0] == 5 && input.arcs[1][1] == 1 && input.arcs[1].capacity == 1);
        test(input.arcs[2][0] == 3 && input.arcs[2][1] == 1 && input.arcs[2].capacity == 1.5);
        test(input.arcs[3][0] == 0 && input.arcs[3][1] == 6 && input.arcs[3].capacity == 1.5);
        test(input.arcs[4][0] == 0 && input.arcs[4][1] == 4 && input.arcs[4].capacity == 2);
        test(input.arcs[5][0] == 7 && input.arcs[5][1] == 1 && input.arcs[5].capacity == 2);

        // pairwise
        test(input.arcs[6][0] == 2 && input.arcs[6][1] == 3 && input.arcs[6].capacity == 3);
        test(input.arcs[7][0] == 3 && input.arcs[7][1] == 2 && input.arcs[7].capacity == 3);
        test(input.arcs[8][0] == 5 && input.arcs[8][1] == 6 && input.arcs[8].capacity == 3);
        test(input.arcs[9][0] == 6 && input.arcs[9][1] == 5 && input.arcs[9].capacity == 3);

        test(input.arcs[10][0] == 3 && input.arcs[10][1] == 7 && input.arcs[10].capacity == 3);
        test(input.arcs[11][0] == 7 && input.arcs[11][1] == 3 && input.arcs[11].capacity == 3);
        test(input.arcs[12][0] == 6 && input.arcs[12][1] == 4 && input.arcs[12].capacity == 3);
        test(input.arcs[13][0] == 4 && input.arcs[13][1] == 6 && input.arcs[13].capacity == 3);

        maxflow_factor m(input);
        const auto maxflow_bound = m.LowerBound();

        const auto mrf = transform_QPBO_max_flow_to_binary_MRF(input);
        test(mrf.unaries.size() == 3);
        test(mrf.unaries[0][0] == 0 && mrf.unaries[0][1] == 2);
        test(mrf.unaries[1][0] == 3 && mrf.unaries[1][1] == 0);
        test(mrf.unaries[2][0] == 0 && mrf.unaries[2][1] == 4);

        test(mrf.pairwise_potentials.size() == 2);
        test(mrf.pairwise_potentials[0].i == 0 && mrf.pairwise_potentials[0].j == 1); 
        test(mrf.pairwise_potentials[0].cost[0][0] == 0 && mrf.pairwise_potentials[0].cost[0][1] == 6 &&
             mrf.pairwise_potentials[0].cost[1][0] == 6 && mrf.pairwise_potentials[0].cost[1][1] == 0);

        test(mrf.pairwise_potentials[1].i == 1 && mrf.pairwise_potentials[1].j == 2); 
        test(mrf.pairwise_potentials[1].cost[0][0] == 6 && mrf.pairwise_potentials[1].cost[0][1] == 0 &&
             mrf.pairwise_potentials[1].cost[1][0] == 0 && mrf.pairwise_potentials[1].cost[1][1] == 6);
        
        qpbo_factor f(mrf);
        const auto mrf_bound = f.LowerBound();

        test(std::abs(maxflow_bound - mrf_bound) <= 1e-8);
    }

    // qpbo, nonsymmetric, nonsubmodular, but tight
    {
        auto input = dimacs_max_flow_input::parse_string(test_dimacs_qpbo_max_flow_non_symmetric);
        test(input.no_nodes = 8);
        test(input.arcs.size() == 14); 

        maxflow_factor m(input);
        const auto maxflow_bound = m.LowerBound();

        const auto mrf = transform_QPBO_max_flow_to_binary_MRF(input);
        test(mrf.unaries.size() == 3);
        test(mrf.unaries.size() == 3);

        test(mrf.pairwise_potentials.size() == 2);
        test(mrf.pairwise_potentials[0].i == 0 && mrf.pairwise_potentials[0].j == 1); 
        test(mrf.pairwise_potentials[0].cost[0][0] == 0 && mrf.pairwise_potentials[0].cost[0][1] == 3 &&
             mrf.pairwise_potentials[0].cost[1][0] == 6 && mrf.pairwise_potentials[0].cost[1][1] == 0);

        test(mrf.pairwise_potentials[1].i == 1 && mrf.pairwise_potentials[1].j == 2); 
        //test(mrf.pairwise_potentials[1].cost[0][0] == 3 && mrf.pairwise_potentials[1].cost[0][1] == 0 &&
        //     mrf.pairwise_potentials[1].cost[1][0] == 0 && mrf.pairwise_potentials[1].cost[1][1] == 6);
        
        qpbo_factor f(mrf);
        const auto mrf_bound = f.LowerBound();

        test(std::abs(maxflow_bound - mrf_bound) <= 1e-8);
    }
}
