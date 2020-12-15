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
paramsMap["BASE_THRESHOLD"]="0"   
paramsMap["DENSE_TIMEGAP_LIFTED"]="60"
paramsMap["NEGATIVE_THRESHOLD_LIFTED"]="0"  
paramsMap["POSITIVE_THRESHOLD_LIFTED"]="0"  
paramsMap["LONGER_LIFTED_INTERVAL"]="4"
paramsMap["MAX_TIMEGAP_BASE"]="60"
paramsMap["MAX_TIMEGAP_LIFTED"]="60"
paramsMap["MAX_TIMEGAP_COMPLETE"]="60"
paramsMap["USE_ADAPTIVE_THRESHOLDS"]="0"
#!Do NOT set ALL_BASE_TO_ZERO= 0 in this settings!



pathToFiles="/home/fuksova/codes/higher-order-disjoint-paths/data/newSolverInput/"
#pathToFiles="/BS/Hornakova/nobackup/newSolverInput/"

#Command line parameters of the solver
solverParameters=["solveFromFiles","-o",pathToFiles+"myOutputPython.txt","--maxIter","2","-v","1"]

params=ldpMP.LdpParams(paramsMap)

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpMP.TimeFramesToVertices()

#Initalizing the structure from a file
timeFrames.init_from_file(pathToFiles+"problemDesc_frames",params)
maxTimeFrame=timeFrames.get_max_time()
#maxTimeFrame=580

maxTimeGapLifted=60 #value that corresponds to maximal length of a lifted edge
batchSize=300
lengthOfOldPathsToCut=maxTimeGapLifted #This can be reconsidered
lengthOfOldPathsToUse=maxTimeGapLifted #This makes sense, longer gap would not influence anything anyway, smaller would lack some information
minTime=1
maxTime=minTime+batchSize-1
maxUsedLabel=0
maxTimeForUsingLabels=0
oldLabels=[]
#globalLabels=[]
stop=False

while not stop:
  print(minTime," ",maxTime)
  batchProc=ldpMP.BatchProcess(timeFrames,oldLabels,maxUsedLabel,maxTimeForUsingLabels,minTime,maxTime)
  batchProc.init_from_file(pathToFiles+"problemDesc")


  firstSolver=ldpMP.Solver(solverParameters)
  #You will use instead:
  #batchProc.init_edges_from_vectors(arrayWithEdgeVertices,arrayWithCosts) #first argument n x 2 array with edge vertices, second n x 1 array of costs
  #batchProc.init_vertices_from_vectors(arrayWithVertices,arrayWithCostsOfVertices) #Not necessary, run AFTER initializing edges! State only those vertices that should have non-zero costs
  #Vertices and edges out of range can be present, they are ignored


  firstInstance=ldpMP.LdpInstance(params,batchProc)

  ldpMP.construct(firstSolver,firstInstance)
  firstSolver.solve()
  firstPaths=firstSolver.get_best_primal()
  batchProc.decode_solution(firstPaths)
  labels=batchProc.get_labels();
  
  if oldLabels:
    indexToDelete=batchProc.get_index_to_delete()
    if indexToDelete<len(oldLabels):
      del oldLabels [indexToDelete:len(oldLabels)]
    oldLabels=oldLabels+labels
  else:
    oldLabels=labels
  maxUsedLabel=batchProc.get_max_used_label()
  #globalLabels.append(labels)
  if maxTime==maxTimeFrame:
    stop=True
  else:
    minTime=maxTime-lengthOfOldPathsToCut-lengthOfOldPathsToUse+1;
    maxTimeForUsingLabels=minTime+lengthOfOldPathsToUse-1;
    maxTime=min(minTime+batchSize-1,maxTimeFrame)

for gl in oldLabels:
  for i in gl:
    print(i, end =" ")
  print("")
