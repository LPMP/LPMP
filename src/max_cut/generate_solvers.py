#!/usr/bin/python
from collections import namedtuple

solver = namedtuple("solver_info", "preamble FMC LP parse_fun filename")

preamble = """
#include "max_cut/max_cut.h"
#include "visitors/standard_visitor.hxx"
"""

preamble_text_input = preamble + "\n#include \"max_cut/max_cut_text_input.h\"\n";
preamble_uai_input = preamble + "\n#include \"mrf/mrf_uai_input.h\"\n";

solvers = [
    solver(preamble_text_input, 'FMC_MAX_CUT', 'LP', 'max_cut_text_input::parse_file', "max_cut_cycle_text_input.cpp"),
    solver(preamble_text_input, 'FMC_ODD_BICYCLE_WHEEL_MAX_CUT', 'LP', 'max_cut_text_input::parse_file', "max_cut_odd_bicycle_wheel_text_input.cpp"),

    solver(preamble_uai_input, 'FMC_MAX_CUT', 'LP', 'binary_MRF_uai_input::parse_file', 'max_cut_cycle_qpbo_uai_input.cpp'),
    solver(preamble_uai_input, 'FMC_ODD_BICYCLE_WHEEL_MAX_CUT', 'LP', 'binary_MRF_uai_input::parse_file', 'max_cut_odd_bicycle_wheel_qpbo_uai_input.cpp'),
    ]


for e in solvers:
   f = open(e.filename,'w')
   lp_type = e.LP + "<" + e.FMC + ">"
   solver_type = "Solver<" + lp_type + ",StandardTighteningVisitor>"
   f.write(e.preamble)
   f.write("using namespace LPMP;\n")
   f.write("int main(int argc, char** argv) {\n")
   f.write("ProblemConstructorRoundingSolver<" + solver_type + "> solver(argc,argv);\n")
   f.write("auto input = LPMP::" + e.parse_fun + "(solver.get_input_file());\n")
   f.write("solver.GetProblemConstructor().construct(input);\n")
   f.write("return solver.Solve();\n")
   f.write("}")
   f.close()
