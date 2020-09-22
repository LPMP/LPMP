#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import ldpMessagePassing as ldpMP

#Initializes structure for holding solver parameters. It expects the path to the solver parameter file.
params=ldpMP.LdpParams("/BS/Hornakova/nobackup/newSolverInput/params_sequence.ini")

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpMP.TimeFramesToVertices()

#Initalizing the structure from a file
timeFrames.init_from_file("/BS/Hornakova/nobackup/newSolverInput/problemDesc_frames",params)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpMP.GraphStructure(timeFrames)

#Adding edges to graph structure from a file.
completeGraphStructure.add_edges_from_file("/BS/Hornakova/nobackup/newSolverInput/problemDesc",params)

instance=LdpInstance(params,&completeGraphStructure);

solver=ldpMP.Solver();

ldpMP.construct(solver,)


#The resulting paths are allways automatically saved into the file with suffix "-all-paths-FINAL" in the output directory. 
#You can optionally save the resulting paths into another file. 
ldpPy.write_output_to_file(paths,"../data/exampleSolverILP/my_python_output.txt")

for path in paths:
  for v in path:
    print(v, end =" ")
  print("")
