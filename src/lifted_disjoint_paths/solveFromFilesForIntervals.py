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


#pathToFiles="/home/fuksova/codes/higher-order-disjoint-paths/data/newSolverInput/"
pathToFiles="/home/fuksova/codes/LPMPData/ICML_INPUTS/MAX_S-2IO0ID333/MOT17-09-DPM/"
outputPath="/home/fuksova/codes/myLPMP"


#Command line parameters of the solver
solverParameters=["solveFromFiles","-o",outputPath+"/myOutputPython.txt","--maxIter","2","--tighten","--tightenConstraintsPercentage","0.05","--tightenInterval","20","--tightenIteration","20","-v","1","--roundingReparametrization","uniform:0.5","--tightenReparametrization","uniform:0.5","--standardReparametrization","uniform:0.5"]


##-----------------One independent task-----------------------
solver1=ldpMP.Solver(solverParameters)

#DO NOT SET THIS! This is just for initialization from files. For Initializing timeFrames from vectors, this is not needed.
paramsMap["MAX_TIMEGAP"]="100"
paramsMap["MIN_TIME_FRAME"]="51"

params=ldpMP.LdpParams(paramsMap)
timeFrames=ldpMP.TimeFramesToVertices()
#You will use:
#init_from_vector(vector_with_time_frames). The vector must contain distribution of vertices in the given time interval
timeFrames.init_from_file(pathToFiles+"problemDesc_frames",params)

#We need ID of minimal vertex in the interval. I can extract this from timeFrames when I call initializatio from file. You need to figure out differently.
minimalVertexID=timeFrames.get_vertex_shift()

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpMP.GraphStructure(timeFrames)

 
#You will load edges from python variable as follows:
#completeGraphStructure.add_edges_from_vectors(vectorWithEdges,vectorWithCosts,minimalVertexID)
#It is the same as before, the only difference is adding the minimalVertexID. If not given,its value
#is zero. You can omit it for the initial interval.
completeGraphStructure.add_edges_from_file(pathToFiles+"problemDesc",params,minimalVertexID)

instance1=ldpMP.LdpInstance(params,completeGraphStructure)

ldpMP.construct(solver1,instance1)

solver1.solve()

#You need to store these paths. Their vertices are numbered from zero. Do not be confused by that.
paths1=solver1.get_best_primal()  


##-----------------Second independent task-----------------------


solver2=ldpMP.Solver(solverParameters)
#This is just for initialization from files. For Initializing timeFrames from vectors, this is not needed.
paramsMap["MAX_TIMEGAP"]="150"
paramsMap["MIN_TIME_FRAME"]="101"

params2=ldpMP.LdpParams(paramsMap)
timeFrames2=ldpMP.TimeFramesToVertices()
timeFrames2.init_from_file(pathToFiles+"problemDesc_frames",params2)
minimalVertexID2=timeFrames2.get_vertex_shift()

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure2=ldpMP.GraphStructure(timeFrames2)

#Adding edges to graph structure from a file.
completeGraphStructure2.add_edges_from_file(pathToFiles+"problemDesc",params2,minimalVertexID2)

instance2=ldpMP.LdpInstance(params2,completeGraphStructure2)

ldpMP.construct(solver2,instance2)

solver2.solve()

paths2=solver2.get_best_primal()  #vertices in these paths do not have global IDs, they are numbered from zero

##------------------Thired independent task: Computing connection of two intervals------------

solverConnection=ldpMP.Solver(solverParameters)

#You have to initialize timeFrames and timeFrames2 from vectors again
#You have to load resulting paths into variables paths1, and paths2
cutOff=10 #how many frames should be cut from the interval solution I suggest the edge time gap here!
isFirst=True  #Set up if the given interval is first or last in the dataset
isLast=False
pathsExtractor1=ldpMP.PathsExtractor(timeFrames,paths1,cutOff,isFirst,isLast,minimalVertexID)

isFirst2=False
isLast2=False
pathsExtractor2=ldpMP.PathsExtractor(timeFrames2,paths2,cutOff,isFirst2,isLast2,minimalVertexID2)


intervalConnection=ldpMP.IntervalConnection(pathsExtractor1,pathsExtractor2)
intervalConnection.init_from_file(pathToFiles+"problemDesc")

instanceConnection=ldpMP.LdpInstance(params,intervalConnection)
ldpMP.construct(solverConnection,instanceConnection)

solverConnection.solve()
connectedPaths=solverConnection.get_best_primal()

#intervalConnection.decode_paths(connectedPaths)

intervalConnection.create_result_structures(connectedPaths)

#The following variables need to be saved for the last part which is label LabelAssignment
middlePaths=intervalConnection.get_middle_paths()  #Paths computed for the connection between the first and the second interval
#These two variables contain connections from one set of paths to the other
firstToMiddle=intervalConnection.get_first_to_middle()
middleToSecond=intervalConnection.get_middle_to_second()

#You can either save theese extracted paths or initializae them from the old stored paths and create pathExtractors again in the last phase
extractedPaths1=pathsExtractor1.get_extracted_paths()
extractedPaths2=pathsExtractor2.get_extracted_paths()



##-----------------Last phase: assign labels to detections from output pats and connection vectors

labelAssignment=ldpMP.LabelAssignment()
labelAssignment.init(extractedPaths1)  #Init is called with the paths extracted from the first interval

 #You will call update with the new set of paths and with the mapping from the previous to the new set of paths
labelAssignment.update(middlePaths,firstToMiddle)

#Again, the same principle
labelAssignment.update(extractedPaths2,middleToSecond)

#At the end, you obtain the labeled vertices.
labels=labelAssignment.get_labels()

for gl in labels:
  for i in gl:
    print(i, end =" ")
  print("")





