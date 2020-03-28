#!/usr/bin/python

from collections import namedtuple
solver = namedtuple("solver_info", "preamble FMC LP filename")

preamble_mp = """
#include "graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"
"""

preamble_BCFW = """
#include "graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"
#include "LP_FWMAP.hxx"
"""

solvers = [
    solver(preamble_mp, 'FMC_MP', 'LP', "graph_matching_mp.cpp"),
    solver(preamble_mp, 'FMC_GM', 'LP', "graph_matching_gm.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP', 'LP', "graph_matching_hungarian_bp.cpp"),
    solver(preamble_mp, 'FMC_MP_T', 'LP', "graph_matching_mp_tightening.cpp"),
    solver(preamble_mp, 'FMC_GM_T', 'LP', "graph_matching_gm_tightening.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP_T', 'LP', "graph_matching_hungarian_bp_tightening.cpp"),
    solver(preamble_BCFW, 'FMC_GM', 'LP_tree_FWMAP', "graph_matching_gm_proximal_bundle.cpp"),
    solver(preamble_BCFW, 'FMC_MCF', 'LP_tree_FWMAP', "graph_matching_mcf_proximal_bundle.cpp"),
    solver(preamble_mp, 'FMC_MP_Q', 'LP', "graph_matching_mp_inter_quadratic_message.cpp"),
    solver(preamble_mp, 'FMC_MP_Q_T', 'LP', "graph_matching_mp_inter_quadratic_message_tightening.cpp")
    ]

for e in solvers:
   f = open(e.filename,'w')
   lp_type = e.LP + "<" + e.FMC + ">"
   solver_type = "Solver<" + lp_type + ",StandardTighteningVisitor>"
   f.write(e.preamble)
   f.write("\nusing namespace LPMP;\nint main(int argc, char** argv) {\n")
   f.write("ProblemConstructorRoundingSolver<" + solver_type + ">")
   f.write("solver(argc,argv);\n")
   f.write("auto input = LPMP::TorresaniEtAlInput::parse_file(solver.get_input_file());\n")
   f.write("solver.GetProblemConstructor().construct(input);\n")
   f.write("return solver.Solve();\n}")
   f.close()
