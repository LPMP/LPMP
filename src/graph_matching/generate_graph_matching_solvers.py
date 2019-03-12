#!/usr/bin/python

from collections import namedtuple
solver = namedtuple("solver_info", "preamble FMC LP parse_fun filename")

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
    solver(preamble_mp, 'FMC_MCF<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_mcf_left.cpp"),
    solver(preamble_mp, 'FMC_MCF<PairwiseConstruction::Right>', 'LP', 'parse_problem',"graph_matching_mcf_right.cpp"),
    solver(preamble_mp, 'FMC_MCF<PairwiseConstruction::BothSides>', 'LP', 'parse_problem',"graph_matching_mcf_both_sides.cpp"),
    solver(preamble_mp, 'FMC_MP<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_mp_left.cpp"),
    solver(preamble_mp, 'FMC_MP<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_mp_right.cpp"),
    solver(preamble_mp, 'FMC_MP<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_mp_both_sides.cpp"),
    solver(preamble_mp, 'FMC_GM<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_gm_left.cpp"),
    solver(preamble_mp, 'FMC_GM<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_gm_right.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_hungarian_bp_left.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_hungarian_bp_right.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP<PairwiseConstruction::BothSides>', 'LP', 'parse_problem', "graph_matching_hungarian_bp_both_sides.cpp"),
    solver(preamble_mp, 'FMC_MCF_T<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_mcf_left_tightening.cpp.cpp"),
    solver(preamble_mp, 'FMC_MCF_T<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_mcf_right_tightening.cpp.cpp"),
    solver(preamble_mp, 'FMC_MCF_T<PairwiseConstruction::BothSides>', 'LP', 'parse_problem', "graph_matching_mcf_both_sides_tightening.cpp"),
    solver(preamble_mp, 'FMC_MP_T<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_mp_left_tightening.cpp"),
    solver(preamble_mp, 'FMC_MP_T<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_mp_right_tightening.cpp"),
    solver(preamble_mp, 'FMC_MP_T<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_mp_both_sides_tightening.cpp"),
    solver(preamble_mp, 'FMC_GM_T<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_gm_left_tightening.cpp"),
    solver(preamble_mp, 'FMC_GM_T<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_gm_right_tightening.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP_T<PairwiseConstruction::Left>', 'LP', 'parse_problem', "graph_matching_hungarian_bp_left_tightening.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP_T<PairwiseConstruction::Right>', 'LP', 'parse_problem', "graph_matching_hungarian_bp_right_tightening.cpp"),
    solver(preamble_mp, 'FMC_HUNGARIAN_BP_T<PairwiseConstruction::BothSides>', 'LP', 'parse_problem', "graph_matching_hungarian_bp_both_sides_tightening.cpp"),
    solver(preamble_BCFW, 'FMC_MCF<PairwiseConstruction::Left>', 'LP_tree_FWMAP', 'ParseProblemMCF_trees', "graph_matching_mcf_proximal_bundle_left.cpp"),
    solver(preamble_BCFW, 'FMC_GM<PairwiseConstruction::Left>', 'LP_tree_FWMAP', 'ParseProblemGM_trees', "graph_matching_gm_proximal_bundle_left.cpp"),
    solver(preamble_BCFW, 'FMC_LOCAL_SUBPROBLEM<PairwiseConstruction::Left>', 'LP_tree_FWMAP', 'ParseProblemLocalSubproblems_trees', "graph_matching_local_subproblems_proximal_bundle_left.cpp"),

    solver(preamble_mp, 'FMC_MP_Q<PairwiseConstruction::BothSides>', 'LP', 'parse_problem', "graph_matching_mp_both_sides_inter_quadratic_message.cpp"),
    solver(preamble_mp, 'FMC_MP_Q_T<PairwiseConstruction::BothSides>', 'LP', 'parse_problem', "graph_matching_mp_both_sides_inter_quadratic_message_tightening.cpp")
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
