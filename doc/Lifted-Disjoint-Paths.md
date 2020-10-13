## Lifted Disjoint Paths

### Problem

 define problem

### File format

### Parameters
   
  - `SPARSIFY = 1`  
    Expects value 0/1. If set to 0, no sparsification is done, so parameters related to sparsificatioin are not needed. If set to 1, sparsification is done. Default and recommended value is 1.
    
  - `MAX_TIMEGAP = 900`  
    Expects a positive integer. Maximum number of time frames to be used  for the input. Only frames 1,2,...,`MAX\_TIMEGAP` are used for the computation. If you want to use all time frames, either leave this parameter out or set it to zero or to a sufficiently large value.
    
  - `INPUT_COST = 0.0`  
    Expects a real value. Cost of starting a track. Default value is 0.

  - `OUTPUT_COST =0.0`  
    Expects a real value. Cost of terminating a track. Default value is 0.
    
  - `MAX_TIMEGAP_COMPLETE = 60`  
  Expects a positive integer. Influences the second step of the two-step procedure. If `COMPLETE_GAP_IN_TRACKLET` is set to one, this value is used for selection of base and lifted edges in creating the tracklet graph. Its default value is the maximum from `MAX_TIMEGAP_BASE` and `MAX_TIMEGAP_LIFTED`.
    
 
**Base graph sparsificatioin parameters**

  - `MAX_TIMEGAP_BASE = 50`  
    Expects a positive integer. Maximal time gap for the base edges (amount of time frames).
    
  - `KNN_K = 3`  
    Expects a positive integer. For every vertex and every time gap, up to `KNN_K` nearest neighbors can be connected by a base edge with that vertex.
    
  - `KNN_GAP = 6`  
    Expects a non-negative integer. Up to `KNN_GAP`, all found nearest neighbors are connected with a processed vertex. Beyond `KNN_GAP`, edges to nearest neighbors have to additionally satisfy other constraints in order to be used.
    
  - `BASE_THRESHOLD = 0.0`  
    Expects a real value. Upper threshold for cost of base edges with a bigger time gap than `KNN_GAP`. Edges longer than `KNN_GAP` have to have cost lower than this threshold in order to be used. Default value is 0.
    
    

**Lifted graph sparsification parameters**

  - `MAX_TIMEGAP_LIFTED = 60`  
    Expects a non-negative integer. Maximal time gap for the lifted edges.
    
  - `DENSE_TIMEGAP_LIFTED =4`  
    Expects a non-negative integer. Every lifted edge with time gap lower or equal to this value will be added if it satisfies the constraints on its cost value given by `NEGATIVE_THRESHOLD_LIFTED` and `POSITIVE_THRESHOLD_LIFTED`.
    
  - `NEGATIVE_THRESHOLD_LIFTED = -1.0`  
    Expects a non-positive real value. Lifted edges with negative cost will only be added if their cost value is lower or equal to this number.
    
  - `POSITIVE_THRESHOLD_LIFTED = 1.0`  
    Expects a non-negative real value. Lifted edges with positive cost will only be added if their cost value is higher or equal to this number.
    
  - `LONGER_LIFTED_INTERVAL = 4`  
    Expects a positive integer. Influences lifted edges with time gap between `DENSE_TIMEGAP_LIFTED` and `MAX_TIMEGAP_LIFTED`. If set to value \(n\), only every \(n\)-th time gap will be used for adding lifted edges.
    

```
code here
```

### Run from python

description

