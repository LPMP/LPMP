## Lifted Disjoint Paths

### Problem
Lifted Disjoint Paths problem was introduced in [1]. Its original implementation including example input files is available here [LifT_Solver](https://github.com/AndreaHor/LifT_Solver).



### Compilation and Running
After running `cmake` for the whole LPMP project, go to your build directory, than to `src/lifted-disjoint-paths` and run `make` here.  
Now, you can either run the solver from the command line and use input from specific input files 
```
./lifted_disjoint_paths_text_input -i /path/to/your/input/inputFile.txt 
``` 

Another possibility is to use python script for running the solver on an example instance in. Here, no input files are needed. The whole problem instance is specified directly in the python script.
```
python3 solveFromVectors.py
```
You can use other command line parameters that are applicable for other problems in LPMP. You can display their full list by running 
```
./lifted_disjoint_paths_text_input --help
```
In case of running the solver from python, you provide the parameters and their values in an array of strings. See `solveFromVectors.py` for an example.

### File format
In case of running from command line, you have to provide several input files. The main input file `inputFile.txt` passed to the solver as the command line argument has the following structure:

```
INPUT_GRAPH=/path/to/your/input/problemDesc
INPUT_FRAMES=/path/to/your/input/problemDesc_frames
INPUT_PARAMS=/path/to/your/input/params_sequence.ini
```
Example input files can be downloaded [here](https://github.com/AndreaHor/LifT_Solver/tree/master/data/exampleSolverILP). Most of the parameters listed in `params_sequence.ini` are not applicable for this solver. The list of relevant parameters is written below. Detailed description of the format of the other two files can be found [here](https://github.com/AndreaHor/LifT_Solver/tree/master/solverILP) in sections "Graph File" and "File with Time Frames".

### Parameters
Parameters of the problem instance are either passed to the solver in the file `params_sequence.ini` in case of running from the command line or are specified in a python dictionary in case of running from python (see `solveFromVectors.py` for an example). Do not forget the keyword `[SOLVER]` in your file `params_sequence.ini`.
   
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
    
### References
[1]: `A. Hornakova, R. Henschel, B. Rosenhahn, P. Swoboda. Lifted Disjoint Paths with Application in Multiple Object Tracking, ICML 2020`
