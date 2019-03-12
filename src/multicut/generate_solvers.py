#!/usr/bin/python
from collections import namedtuple

solver = namedtuple("solver_info", "preamble FMC LP parse_fun filename")

preamble = """
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
"""

preamble_text_input = preamble + "\n#include \"multicut/multicut_text_input.h\"\n";
preamble_opengm_input = preamble + "\n#include \"multicut/multicut_opengm_input.h\"\n";
preamble_andres_input = preamble + "\n#include \"multicut/multicut_andres_input.h\"\n";

solvers = [
    solver(preamble_text_input, 'FMC_MULTICUT<MessageSendingType::SRMP>', 'LP', 'multicut_text_input::parse_file', "multicut_cycle_text_input.cpp"),
    solver(preamble_text_input, 'FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>', 'LP', 'multicut_text_input::parse_file', "multicut_odd_wheel_text_input.cpp"),
    solver(preamble_text_input, 'FMC_ODD_BICYCLE_WHEEL_MULTICUT', 'LP', 'multicut_text_input::parse_file', "multicut_odd_bicycle_wheel_text_input.cpp"),

    solver(preamble_opengm_input, 'FMC_MULTICUT<MessageSendingType::SRMP>', 'LP', 'multicut_opengm_input::parse_file', "multicut_cycle_opengm_input.cpp"),
    solver(preamble_opengm_input, 'FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>', 'LP', 'multicut_opengm_input::parse_file', "multicut_odd_wheel_opengm_input.cpp"),
    solver(preamble_opengm_input, 'FMC_ODD_BICYCLE_WHEEL_MULTICUT', 'LP', 'multicut_opengm_input::parse_file', "multicut_odd_bicycle_wheel_opengm_input.cpp"),

    solver(preamble_andres_input, 'FMC_MULTICUT<MessageSendingType::SRMP>', 'LP', 'multicut_andres_input::parse_file', "multicut_cycle_andres_input.cpp"),
    solver(preamble_andres_input, 'FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>', 'LP', 'multicut_andres_input::parse_file', "multicut_odd_wheel_andres_input.cpp"),
    solver(preamble_andres_input, 'FMC_ODD_BICYCLE_WHEEL_MULTICUT', 'LP', 'multicut_andres_input::parse_file', "multicut_odd_bicycle_wheel_andres_input.cpp")
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
