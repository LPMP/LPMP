#!/usr/bin/python

from collections import namedtuple
solver = namedtuple("solver_info", "preamble FMC LP parse_fun filename")

preamble = """
#include "horizon_tracking/horizon_tracking.h"
#include "horizon_tracking/horizon_tracking_primal_rounding.hxx"
#include "visitors/standard_visitor.hxx"
#include "LP.h"
"""

#include "LP_FWMAP.hxx"
#include "LP_conic_bundle.hxx"
#include "tree_decomposition.hxx"

solvers = [
    solver(preamble + "#include \"LP_FWMAP.hxx\"", 'FMC_HORIZON_TRACKING_CHAINS', 'LP_tree_FWMAP', 'horizon_tracking_uai_input::parse_file', "horizon_tracking_fwmap.cpp"),
    solver(preamble + "#include \"LP_conic_bundle.hxx\"", 'FMC_HORIZON_TRACKING_CHAINS', 'LP_conic_bundle', 'horizon_tracking_uai_input::parse_file', "horizon_tracking_conic_bundle.cpp"),
    solver(preamble + "#include \"tree_decomposition.hxx\"", 'FMC_HORIZON_TRACKING_CHAINS', 'LP_subgradient_ascent', 'horizon_tracking_uai_input::parse_file', "horizon_tracking_subgradient.cpp"),
    ]

for e in solvers:
   f = open(e.filename,'w')
   lp_type = e.LP + "<" + e.FMC + ">"
   solver_type = "Solver<" + lp_type + ",StandardVisitor>"
   f.write(e.preamble)
   f.write("\nusing namespace LPMP;\nint main(int argc, char** argv) {\n")
   f.write(solver_type + " solver(argc,argv);\n")
   f.write("auto input = " + e.parse_fun + "(solver.get_input_file());\n")
   f.write("construct_horizon_tracking_problem_on_grid_to_chains(input, solver, solver.template GetProblemConstructor<0>());\n")
   f.write("order_nodes_by_label_space_cadinality(solver.template GetProblemConstructor<0>());\n")
   f.write("solver.Solve();\n")
   f.write("solver.GetLP().write_back_reparametrization();\n")
   f.write("round_primal_solution(solver);\n")
   f.write("solver.WritePrimal();\n")
   f.write("std::cout<<\"\\n\\n Primal Cost: \"<<solver.primal_cost();\n")
   f.write("std::cout<<\"\\n Percentage duality gap: \"<<100.0 * (solver.primal_cost() - solver.lower_bound()) / solver.lower_bound() <<\"\%\\n\\n\";\n")
   f.write("}")
   f.close()

