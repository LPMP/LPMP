#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import numpy as np
import lpmp_ldp

#Parameters of the problem instance and parameters related to LDP methods come here (LDP parameters)
paramsMap={}
paramsMap["INPUT_COST"]="0"
paramsMap["OUTPUT_COST"]="0"
paramsMap["SPARSIFY"]="1"
paramsMap["KNN_GAP"]="3"
paramsMap["KNN_K"]="3"
paramsMap["BASE_THRESHOLD"]="0.9"   #remove 90% of base edges
paramsMap["DENSE_TIMEGAP_LIFTED"]="12"
paramsMap["NEGATIVE_THRESHOLD_LIFTED"]="-0.1"  #remove 10% of negative lifted edges
paramsMap["POSITIVE_THRESHOLD_LIFTED"]="0.1"  #remove 10% of positive lifted edges
paramsMap["LONGER_LIFTED_INTERVAL"]="4"
paramsMap["MAX_TIMEGAP_BASE"]="60"
paramsMap["MAX_TIMEGAP_LIFTED"]="60"
paramsMap["MAX_TIMEGAP_COMPLETE"]="60"
#paramsMap["USE_ADAPTIVE_THRESHOLDS"]="0"
paramsMap["TIGHT_MIN_IMPROVEMENT"]="0.00001"
paramsMap["TIGHT_MAX_EDGE_USAGE"]="4"


#Command line parameters of LPMP
#solverParameters=["solveFromFiles","-o",pathToFiles+"myOutputPython.txt","--maxIter","5","-v","0"]
#Command line parameters of the solver with tightening enabled
solverParameters=["solveFromFiles","-o","myOutputPython.txt","--maxIter","15","--tighten","--tightenConstraintsPercentage","0.02","--tightenInterval","10","--tightenIteration","10","-v","1"]



#Construct the solver.
solver=lpmp_ldp.ldp.Solver(solverParameters)


#Initializes structure for holding LDP parameters. It expects a string to string map (dictionary) as an input. 
#Such map can be alternatively obtained from file params_sequence.ini
params=lpmp_ldp.ldp.LdpParams(paramsMap)


#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=lpmp_ldp.ldp.TimeFramesToVertices()

#This array contains information about how many detections are in each time frame. The array must contain non-negative integer values. 
#Numbering of graph vertices in time frames: Frame 1: 0,1,2; Frame 2: 3,4; Frame 3: 4,5,6
vect= np.array([3, 2, 3])
numberOfVertices=8;

#Initalizing the structure from an array
#This structure can be alternatively initialized from file problemDesc_frames
timeFrames.init_from_vector(vect)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=lpmp_ldp.ldp.GraphStructure(timeFrames)

#Edge vector. Each edge is represented by a pair of vertices. Vertices ar numbered from zero
edgeVector=np.array([[0,3],[0,4],[1,3],[1,4],[2,3],[2,4],[3,5],[3,6],[3,7],[4,5],[4,6],[4,7],[0,5],[0,6],[0,7],[1,5],[1,6],[1,7],[2,6],[2,7]])

#Cost vector. The cost values must correspond to the edge ordering in edgeVector.
costVector=np.array([-2.2,-1,-1,-2.2,-1,-1,-2.2,-1,-1,-1,-2.2,-1,  -2.2,-1,-1, -1,-2.2,-1, -1,-2.2])

#This vector can be used for setting cost of vertices. It is an optional parameter. If it is not used, all vertices
#have cost zero
verticesCost=np.array([0.0,-10,-20,0,-10,-2,0,0.5])

#Adding edges to graph structure from vectors. This method adds all edges regardles restrictions given in params.
completeGraphStructure.add_edges_from_vectors_all(edgeVector,costVector)
#Alternative method where either MAX_TIMEGAP_COMPLETE or the maximum from MAX_TIMEGAP_BASE and MAX_TIMEGAP_LIFTED applies for selection of edges. 
#completeGraphStructure.add_edges_from_vectors(edgeVector,costVector,params)

#Setting costs of vertices. Optional
completeGraphStructure.set_score_of_vertices(verticesCost)
print("score set")


#Create the problem instance from the parameters and the complete graph structure.
instance=lpmp_ldp.ldp.LdpInstance(params,completeGraphStructure)
print("instance created")


#Run problem constructor using the solver and problem instance.
lpmp_ldp.ldp.construct(solver,instance)
print("contructor called")


#Run the message passing solver.
solver.solve()

#Extract the resulting paths from the best primal solution.
paths=solver.get_best_primal()

#Obtaining edge labels based on the obtained paths. Label 1 (True) is given iff endpoints of the respective edge belong to the same path.
edgeLabels=lpmp_ldp.ldp.get_lifted_edge_labels(edgeVector,paths,numberOfVertices)
#An alternative method gives labels to the edges actually used in the input graph.
#edgeLabels=completeGraphStructure.get_edge_labels(paths)

#Get lower bound 
lowerBound=solver.get_lower_bound()

#Get primal value
primalValue=solver.get_best_primal_value()

print(lowerBound)
print(primalValue)


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






