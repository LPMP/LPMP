#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import ldpMessagePassingPy as ldpMP


#pathToFiles="/home/fuksova/codes/higher-order-disjoint-paths/data/newSolverInput/"
pathToFiles="/BS/Hornakova/nobackup/newSolverInput/"

#Command line parameters of the solver
solverParameters=["solveFromFiles","-o",pathToFiles+"myOutputPython.txt","--maxIter","15","--tighten","--tightenConstraintsMax","10"]



#Create a parser for parameters
paramsParser=ldpMP.ParametersParser()

#Parses parameters from file
paramsParser.init_from_file(pathToFiles+"params_sequence.ini")

#Initializes structure for holding solver parameters. It expects a string to string map (dictionary) as an input. ParametersParser.get_parsed_params() can be used for providing such map.
params=ldpMP.LdpParams(paramsParser.get_parsed_params())
print("params read")

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpMP.TimeFramesToVertices()

#Initalizing the structure from a file
timeFrames.init_from_file(pathToFiles+"problemDesc_frames",params)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpMP.GraphStructure(timeFrames)

#Adding edges to graph structure from a file.
completeGraphStructure.add_edges_from_file(pathToFiles+"problemDesc",params)

instance=ldpMP.LdpInstance(params,completeGraphStructure)

solver=ldpMP.Solver(solverParameters)

ldpMP.construct(solver,instance)

solver.solve()
