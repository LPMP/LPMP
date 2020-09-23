#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import ldpMessagePassing as ldpMP

pathToFiles="/home/fuksova/codes/higher-order-disjoint-paths/data/newSolverInput/"

#Initializes structure for holding solver parameters. It expects the path to the solver parameter file.
params=ldpMP.LdpParams("/home/fuksova/codes/higher-order-disjoint-paths/data/newSolverInput/params_sequence.ini")
print("params read")
#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpMP.TimeFramesToVertices()

#Initalizing the structure from a file
timeFrames.init_from_file(pathToFiles+"problemDesc_frames",params)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpMP.GraphStructure(timeFrames)

#Adding edges to graph structure from a file.
completeGraphStructure.add_edges_from_file(pathToFiles+"problemDesc",params)

instance=LdpInstance(params,completeGraphStructure)


ldpMP.construct(solver,instance)

solver=ldpMP.Solver()

