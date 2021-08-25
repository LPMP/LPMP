#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import lpmp_ldp

pathToFiles="./input_files/"

#Command line parameters of the solver
solverParameters=["solveFromFiles","-o","myOutputPython.txt","--maxIter","50"]

solver=lpmp_ldp.ldp.Solver(solverParameters)

#Create a parser for parameters
paramsParser=lpmp_ldp.ldp.ParametersParser()

#Parses parameters from file
paramsParser.init_from_file(pathToFiles+"params_sequence.ini")

#Initializes structure for holding solver parameters. It expects a string to string map (dictionary) as an input. ParametersParser.get_parsed_params() can be used for providing such map.
params=lpmp_ldp.ldp.LdpParams(paramsParser.get_parsed_params())
print("params read")

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=lpmp_ldp.ldp.TimeFramesToVertices()

#Initalizing the structure from a file
timeFrames.init_from_file(pathToFiles+"problemDesc_frames",params)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=lpmp_ldp.ldp.GraphStructure(timeFrames)

#Adding edges to graph structure from a file.
completeGraphStructure.add_edges_from_file(pathToFiles+"problemDesc",params,0)

instance=lpmp_ldp.ldp.LdpInstance(params,completeGraphStructure)



lpmp_ldp.ldp.construct(solver,instance)

solver.solve()
