#!/usr/bin/python

from collections import namedtuple
solver = namedtuple("solver_info", "preamble FMC LP filename")

preamble = """
#include "graph_matching/multigraph_matching.hxx"
#include "visitors/standard_visitor.hxx"
"""

solvers = [
    solver(preamble, 'FMC_MGM<true>', 'LP', "multigraph_matching_mp.cpp"),
    solver(preamble, 'FMC_MGM<false>', 'LP', "multigraph_matching_mcf.cpp"),
    solver(preamble, 'FMC_MGM_T<true>', 'LP', "multigraph_matching_tightening_mp.cpp"),
    solver(preamble, 'FMC_MGM_T<false>', 'LP', "multigraph_matching_tightening_mcf.cpp"),
    ]

for e in solvers:
   f = open(e.filename,'w')
   lp_type = e.LP + "<" + e.FMC + ">"
   solver_type = "Solver<" + lp_type + ",StandardTighteningVisitor>"
   f.write(e.preamble)
   f.write("\nusing namespace LPMP;\nint main(int argc, char** argv) {\n")
   f.write("ProblemConstructorRoundingSolver<" + solver_type + ">")
   f.write("solver(argc,argv);\n")
   f.write("auto input = Torresani_et_al_multigraph_matching_input::parse_file(solver.get_input_file());\n");
   f.write("solver.template GetProblemConstructor<0>().construct(input);\n")
   f.write("return solver.Solve();\n}")
   f.close()

