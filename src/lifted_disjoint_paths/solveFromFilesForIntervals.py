#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

import lpmp_ldp
import numpy as np

#LDP parameters
paramsMap={}
paramsMap["INPUT_COST"]="0"
paramsMap["OUTPUT_COST"]="0"
paramsMap["SPARSIFY"]="1"
paramsMap["KNN_GAP"]="3"
paramsMap["KNN_K"]="3"
paramsMap["BASE_THRESHOLD"]="0.0"   
paramsMap["DENSE_TIMEGAP_LIFTED"]="12"
paramsMap["NEGATIVE_THRESHOLD_LIFTED"]="-0.1" 
paramsMap["POSITIVE_THRESHOLD_LIFTED"]="0.1"  
paramsMap["LONGER_LIFTED_INTERVAL"]="4"
paramsMap["MAX_TIMEGAP_BASE"]="60"
paramsMap["MAX_TIMEGAP_LIFTED"]="60"
paramsMap["MAX_TIMEGAP_COMPLETE"]="60"
paramsMap["USE_ADAPTIVE_THRESHOLDS"]="0"
paramsMap["TIGHT_MIN_IMPROVEMENT"]="0.00001"
paramsMap["TIGHT_MAX_EDGE_USAGE"]="4"
paramsMap["SAVE_INTERMEDIATE"]="0"



pathToFiles="./input_files/"



#Command line parameters of the solver
solverParameters=["solveFromFiles","-o","myOutputPython.txt","--maxIter","2","--tighten","--tightenConstraintsPercentage","0.05","--tightenInterval","20","--tightenIteration","20","-v","1","--roundingReparametrization","uniform:0.5","--tightenReparametrization","uniform:0.5","--standardReparametrization","uniform:0.5"]


##-----------------One independent task-----------------------
solver1=lpmp_ldp.ldp.Solver(solverParameters)

#This is just for initialization from files. For Initializing timeFrames from vectors, this is not needed.
paramsMap["MAX_TIMEGAP"]="180"
paramsMap["MIN_TIME_FRAME"]="1"

#Initialization of LDP parameters from paramsMap
params=lpmp_ldp.ldp.LdpParams(paramsMap)

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=lpmp_ldp.ldp.TimeFramesToVertices()

#Initialization of mapping between time frames and graph vertices from text file.
#In order to read the time frame data from python variables, you would use: 
#init_from_vector(vector_with_time_frames). The vector must contain distribution of vertices in the given time interval
timeFrames.init_from_file(pathToFiles+"problemDesc_frames",params)

#We need to provide the ID of the minimal vertex in the interval. I can extract this from timeFrames when 
#I call initializatio from file. When you use python variables for time frame input, you need to figure out differently.
minimalVertexID=timeFrames.get_vertex_shift()

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=lpmp_ldp.ldp.GraphStructure(timeFrames)

#Loading the edges from the text file. 
#In order to load edges from python variable, you use:
#completeGraphStructure.add_edges_from_vectors(vectorWithEdges,vectorWithCosts,minimalVertexID)
#It is the same as in solveFromFiles.py, the only difference is adding the minimalVertexID. If not given,its value
#is zero. Therefore, you can omit minimalVertexID for the initial interval.
completeGraphStructure.add_edges_from_file(pathToFiles+"problemDesc",params,minimalVertexID)

#LDP problem instance initialization
instance1=lpmp_ldp.ldp.LdpInstance(params,completeGraphStructure)

#Calling solver constructor using solver and instance
lpmp_ldp.ldp.construct(solver1,instance1)

#Running the solver
solver1.solve()

#You will typically need to store these paths if you solve many independent intervals. The vertices are numbered from zero 
#even if they do not belong to the initial interval. Do not be confused by that. The numbers are adjusted automatically during the 
#paths merging.
paths1=solver1.get_best_primal()  


##-----------------Second independent task-----------------------
#This is the second interval to solve. Since this task is completely independent on solving the first interval,
#you can run them in parallel as different processes.

#Another solver instace
solver2=lpmp_ldp.ldp.Solver(solverParameters)

#This is just for initialization from files. For Initializing timeFrames from vectors, this is not needed.
paramsMap["MAX_TIMEGAP"]="360"
paramsMap["MIN_TIME_FRAME"]="181"

#Initialization of parameters, time frame structure and minimalVertexID2.
params2=lpmp_ldp.ldp.LdpParams(paramsMap)
timeFrames2=lpmp_ldp.ldp.TimeFramesToVertices()
timeFrames2.init_from_file(pathToFiles+"problemDesc_frames",params2)
minimalVertexID2=timeFrames2.get_vertex_shift()

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure2=lpmp_ldp.ldp.GraphStructure(timeFrames2)

#Adding edges to graph structure from a file.
completeGraphStructure2.add_edges_from_file(pathToFiles+"problemDesc",params2,minimalVertexID2)

#Constructor of problem instance.
instance2=lpmp_ldp.ldp.LdpInstance(params2,completeGraphStructure2)

#Constructor of the solver
lpmp_ldp.ldp.construct(solver2,instance2)

#Running the solver
solver2.solve()

#Resulting paths to store.
#vertices in these paths do not have global IDs, they are numbered from zero. The vertex numbers are automatically adjusted
#during paths merging.
paths2=solver2.get_best_primal()  

##------------------Third independent task: Computing connection of two previous intervals------------
#Once the two previous intervals are solved, you can compute their connection here.
solverConnection=lpmp_ldp.ldp.Solver(solverParameters)

#how many frames should be cut from the interval solution. I suggest the edge time gap here!
cutOff=int(paramsMap["MAX_TIMEGAP_COMPLETE"]) 


#NOTE: You will typically run this from a separate script. Therefore, you have to initialize timeFrames and timeFrames2 again
#and you have to load resulting paths into variables paths1, and paths2

#Set up if the given interval is first or last in the sequence
isFirst=True  
isLast=False
#This is a structure that is used for extracting paths segments from the results on the separate intervals
pathsExtractor1=lpmp_ldp.ldp.PathsExtractor(timeFrames,paths1,cutOff,isFirst,isLast,minimalVertexID)

isFirst2=False
isLast2=False
pathsExtractor2=lpmp_ldp.ldp.PathsExtractor(timeFrames2,paths2,cutOff,isFirst2,isLast2,minimalVertexID2)

#This structure is necessary for creating the LDP instance for connecting results of two neighboring intervals
intervalConnection=lpmp_ldp.ldp.IntervalConnection(pathsExtractor1,pathsExtractor2)

#We need to load the edges between detections. Edges that are irrelevant for the present detections are ignored.
#Be sure that you provide also all edges between the paths segments and the unassigned detections.
#In order to init from edge vectors, you would use:
#intevalConnection.init_edges_from_vectors(edgeVector,costVector)
intervalConnection.init_from_file(pathToFiles+"problemDesc")

instanceConnection=lpmp_ldp.ldp.LdpInstance(params,intervalConnection)
lpmp_ldp.ldp.construct(solverConnection,instanceConnection)

solverConnection.solve()
connectedPaths=solverConnection.get_best_primal()

#intervalConnection.decode_paths(connectedPaths)

intervalConnection.create_result_structures(connectedPaths)

#The following variables need to be saved for the last part which is label LabelAssignment

#Paths computed for the connection between the first and the second interval
middlePaths=intervalConnection.get_middle_paths()  

#These two variables contain connections from one set of paths to the other
firstToMiddle=intervalConnection.get_first_to_middle()
middleToSecond=intervalConnection.get_middle_to_second()

#You can either save theese extracted paths or initializae them from the old stored paths and create pathExtractors again
#in the last phase
extractedPaths1=pathsExtractor1.get_extracted_paths()
extractedPaths2=pathsExtractor2.get_extracted_paths()



##-----------------Last phase: assign labels to detections from output pats and connection vectors
#Constructor of the structore for label assignment
labelAssignment=lpmp_ldp.ldp.LabelAssignment()

#Init is called with the paths extracted from the first interval
labelAssignment.init(extractedPaths1)  

 #You will call update with the new set of paths and with the mapping from the previous to the new set of paths
labelAssignment.update(middlePaths,firstToMiddle)

#Again, the same principle
labelAssignment.update(extractedPaths2,middleToSecond)

#At the end, you obtain the labeled vertices. It is a vector of size nx2 where n is the number of all input detections assigned
#to a trajectory. First element of each row is the detection ID, the second element is its label
labels=labelAssignment.get_labels()

##Uncomment for printing the result
#for gl in labels:
#  for i in gl:
#    print(i, end =" ")
#  print("")





