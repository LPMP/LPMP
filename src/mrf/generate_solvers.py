#!/usr/bin/python

from collections import namedtuple
solver = namedtuple("solver_info", "preamble FMC LP parse_fun filename visitor construct_trees")

preamble = """
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_uai_input.h"
"""

opengm_preamble = """
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_opengm_input.h"
"""

combiLP_preamble = """
#include "mrf/graphical_model.h"
#include "mrf/combiLP.hxx"
#include "gurobi_interface.hxx"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_uai_input.h"
"""

opengm_combiLP_preamble = """
#include "mrf/graphical_model.h"
#include "combiLP.hxx"
#include "gurobi_interface.hxx"
#include "visitors/standard_visitor.hxx"
#include "mrf/mrf_opengm_input.h"
"""

dd_preamble = """
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "LP_FWMAP.hxx"
#include "mrf/mrf_uai_input.h"
"""

dd_preamble_opengm = """
#include "mrf/graphical_model.h" 
#include "visitors/standard_visitor.hxx" 
#include "LP_FWMAP.hxx"
#include "mrf/mrf_opengm_input.h"
"""

conic_bundle_preamble = """
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "LP_conic_bundle.hxx"
#include "mrf/mrf_uai_input.h"
"""

conic_bundle_preamble_opengm = """
#include "mrf/graphical_model.h"
#include "visitors/standard_visitor.hxx"
#include "LP_conic_bundle.hxx"
#include "mrf/mrf_opengm_input.h"
"""

solvers = [
    solver(opengm_preamble, 'FMC_SRMP', 'LP<FMC_SRMP>', 'mrf_opengm_input::parse_file', 'srmp_opengm.cpp', 'StandardVisitor', False),
    solver(opengm_preamble, 'FMC_SRMP_T', 'LP<FMC_SRMP_T>', 'mrf_opengm_input::parse_file', 'srmp_opengm_tightening.cpp', 'StandardTighteningVisitor', False),
    solver(opengm_preamble, 'FMC_MPLP', 'LP<FMC_MPLP>', 'mrf_opengm_input::parse_file', 'mplp_opengm.cpp', 'StandardVisitor', False),
    solver(preamble, 'FMC_SRMP', 'LP<FMC_SRMP>', 'mrf_uai_input::parse_file', 'srmp_uai.cpp', 'StandardVisitor', False),
    solver(preamble, 'FMC_SRMP_T', 'LP<FMC_SRMP_T>', 'mrf_uai_input::parse_file', 'srmp_uai_tightening.cpp', 'StandardTighteningVisitor', False),
    solver(preamble, 'FMC_MPLP', 'LP<FMC_MPLP>', 'mrf_uai_input::parse_file', 'mplp_uai.cpp', 'StandardVisitor', False),
    solver(combiLP_preamble, 'FMC_SRMP', 'combiLP<DD_ILP::gurobi_interface, LP<FMC_SRMP>>', 'mrf_uai_input::parse_file', 'srmp_uai_combiLP.cpp', 'StandardVisitor', False),
    solver(opengm_combiLP_preamble, 'FMC_SRMP', 'combiLP<DD_ILP::gurobi_interface, LP<FMC_SRMP>>', 'mrf_opengm_input::parse_file', 'srmp_opengm_combiLP.cpp', 'StandardVisitor', False),
    solver(dd_preamble, 'FMC_SRMP', 'LP_tree_FWMAP<FMC_SRMP>', 'mrf_uai_input::parse_file', 'FWMAP_uai.cpp', 'StandardVisitor', True),
    solver(dd_preamble_opengm, 'FMC_SRMP', 'LP_tree_FWMAP<FMC_SRMP>', 'mrf_opengm_input::parse_file', 'FWMAP_opengm.cpp', 'StandardVisitor', True),
    solver(conic_bundle_preamble, 'FMC_SRMP', 'LP_conic_bundle<FMC_SRMP>', 'mrf_uai_input::parse_file', 'conic_bundle_uai.cpp', 'StandardVisitor', True),
    solver(conic_bundle_preamble_opengm, 'FMC_SRMP', 'LP_conic_bundle<FMC_SRMP>', 'mrf_opengm_input::parse_file', 'conic_bundle_opengm.cpp', 'StandardVisitor', True)
    ]


for e in solvers:
   f = open(e.filename,'w')
   f.write(e.preamble)
   f.write("\nusing namespace LPMP;\n")
   f.write("int main(int argc, char** argv) {\n")
   f.write("MpRoundingSolver<Solver<" + e.LP + "," + e.visitor + ">> solver(argc,argv);\n")
   f.write("auto input = " + e.parse_fun + "(solver.get_input_file());\n") 
   f.write("solver.GetProblemConstructor().construct(input);\n")
   if e.construct_trees:
       f.write("auto trees = solver.GetProblemConstructor().compute_forest_cover();\n")
       f.write("for(auto& tree : trees) { solver.GetLP().add_tree(tree); }\n") 
   f.write("return solver.Solve();\n}")
