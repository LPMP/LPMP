#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import ldpMessagePassingPy as ldpMP
import numpy as np


paramsMap={}
paramsMap["INPUT_COST"]="0"
paramsMap["OUTPUT_COST"]="0"
paramsMap["SPARSIFY"]="1"
paramsMap["KNN_GAP"]="3"
paramsMap["KNN_K"]="3"
paramsMap["BASE_THRESHOLD"]="0.9"   #remove 90% of base edges
paramsMap["DENSE_TIMEGAP_LIFTED"]="12"
paramsMap["NEGATIVE_THRESHOLD_LIFTED"]="0.1"  #remove 10% of negative lifted edges
paramsMap["POSITIVE_THRESHOLD_LIFTED"]="0.1"  #remove 10% of positive lifted edges
paramsMap["LONGER_LIFTED_INTERVAL"]="4"
paramsMap["MAX_TIMEGAP_BASE"]="60"
paramsMap["MAX_TIMEGAP_LIFTED"]="60"
paramsMap["MAX_TIMEGAP_COMPLETE"]="60"
paramsMap["USE_ADAPTIVE_THRESHOLDS"]="1"
paramsMap["TIGHT_MIN_IMPROVEMENT"]="0.00001"
paramsMap["TIGHT_MAX_EDGE_USAGE"]="4"


pathToFiles="/home/fuksova/codes/higher-order-disjoint-paths/data/newSolverInput/"
#pathToFiles="/BS/Hornakova/nobackup/newSolverInput/"

#Command line parameters of LPMP
#solverParameters=["solveFromFiles","-o",pathToFiles+"myOutputPython.txt","--maxIter","5","-v","0"]
#Command line parameters of the solver with tightening enabled
solverParameters=["solveFromFiles","-o",pathToFiles+"myOutputPython.txt","--maxIter","15","--tighten","--tightenConstraintsPercentage","0.02","--tightenInterval","10","--tightenIteration","10","-v","1"]



#Construct the solver.
solver=ldpMP.Solver(solverParameters)


#Initializes structure for holding solver parameters. It expects a string to string map (dictionary) as an input. ParametersParser.get_parsed_params() can be alternatively used for providing such map.
params=ldpMP.LdpParams(paramsMap)


#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpMP.TimeFramesToVertices()

#This array contains information about how many detections are in each time frame. The array must contain non-negative integer values. 
#Numbering of graph vertices in time frames: Frame 1: 0,1,2; Frame 2: 3,4; Frame 3: 4,5,6
vect= np.array([3, 2, 3])
numberOfVertices=8;

#Initalizing the structure from a file
timeFrames.init_from_vector(vect)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpMP.GraphStructure(timeFrames)

#Edge vector. Each edge is represented by a pair of vertices. Vertices ar numbered from zero
edgeVector=np.array([[0,3],[0,4],[1,3],[1,4],[2,3],[2,4],[3,5],[3,6],[3,7],[4,5],[4,6],[4,7],[0,5],[0,6],[0,7],[1,5],[1,6],[1,7],[2,6],[2,7]])
costVector=np.array([-2.2,-1,-1,-2.2,-1,-1,-2.2,-1,-1,-1,-2.2,-1,  -2.2,-1,-1, -1,-2.2,-1, -1,-2.2])
verticesScore=np.array([0.0,-10,-20,0,-10,-2,0,0.5])

#Adding edges to graph structure from vectors. This method adds all edges regardles restrictions given in params.
completeGraphStructure.add_edges_from_vectors_all(edgeVector,costVector)
completeGraphStructure.set_score_of_vertices(verticesScore)

#Alternative method where either MAX_TIMEGAP_COMPLETE or the maximum from MAX_TIMEGAP_BASE and MAX_TIMEGAP_LIFTED applies for selection of edges. 
#completeGraphStructure.add_edges_from_vectors(edgeVector,costVector,params)

#Create the problem instance.
instance=ldpMP.LdpInstance(params,completeGraphStructure)



#Run problem constructor using the solver and problem instance.
ldpMP.construct(solver,instance)

#Run the message passing solver.
solver.solve()

#Extract the resulting paths from the best primal solution.
paths=solver.get_best_primal()

#Obtaining edge labels based on the obtained paths. Label 1 (True) is given iff endpoints of the respective edge belong to the same path.
edgeLabels=ldpMP.get_lifted_edge_labels(edgeVector,paths,numberOfVertices)

#Get lower bound 
lowerBound=solver.get_lower_bound()

#Get primal value
primalValue=solver.get_best_primal_value()

print(lowerBound)
print(primalValue)

#An alternative method gives labels to the edges actually used in the input graph.
#edgeLabels=completeGraphStructure.get_edge_labels(paths)

#for edge in edgeLabels:
#  print(edge)
  
#for edge in usedEdgesVector:
#   for v in edge:
#      print(v,end =" ")
#   print("")

#for path in paths:
#  for v in path:
#    print(v, end =" ")
#  print("")






