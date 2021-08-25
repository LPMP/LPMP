## Lifted Disjoint Paths

### Problem
Lifted Disjoint Paths problem was introduced in [1]. Its original implementation including example input files is available here [LifT_Solver](https://github.com/AndreaHor/LifT_Solver).



### Compilation and Running

#### For usage from the command line
In order to run the solver from the command line, follow these steps.
1. If you have not done it yet, clone the LPMP project and get its submodules:
   ``` 
   git clone https://github.com/LPMP/LPMP.git
   cd LPMP
   git submodule update --init --remote --recursive```
2. If you have not done it yet, create a build directory and run cmake
   ```
   mkdir LPMP-build
   cd LPMP-build
   cmake -D CMAKE_BUILD_TYPE=Release path/to/LPMP
   ```
3. Run make for the lifted disjoint paths project
   ```
   cd src/lifted-disjoint-paths
   make   
   ```

4. Run the solver from command line

   ```
   ./lifted_disjoint_paths_text_input -i /path/to/your/input/inputFile.txt -o /path/to/your/output/outputFile.txt
   ```

#### For usage from python
If you want to use the solver from python, you have two options:
1. Either: Download the source code with `git clone` (see Step 1 in the previous section) and run `pip install` on the `LPMP` directory. 
   ```    
    PACKAGES="ldp" python3 -m pip install path/to/LPMP
   ```
   
2. OR: Install the solver directly from the repository without downloading the source code.
   ```
   PACKAGES="ldp" python3 -m pip install git+https://github.com/LPMP/LPMP.git
   ```
   If you want to use the ldp solver together with other packages from this repository, follow the python installation instructions on the [`front page`](https://github.com/LPMP/LPMP).

You can test if the python installation was successfull by runnig an example python script
      ```python3 LPMP/src/lifted-disjoint-paths/solveFromVectors.py```

Another possibility is to use python script for running the solver on an example instance. Here, no input files are needed. The whole problem instance is specified directly in the python script. It is possible to use one graph structure as an input. The solver will extract the base graph and the lifted graph from it. You can run the respective example script by running
 ```
python3 solveFromVectors.py
 ```
Alternatively, you can directly specify the base graph and the lifted graph separately. You can run the respective example script by running
```
python3 solveFromVectorsTwoGraphs.py
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
  - `MIN_TIMEFRAME = 1`  
    Expects a positive integer. The index of the first time frame to be used for the input. Only frames `MIN_TIMEFRAME`,...,`MAX\_TIMEGAP` are used for the computation. Note that the frames are numbered from 1. If you want to use all time frames from the first one, either leave this parameter out or set it to 1. The default value is 1.
  - `INPUT_COST = 0.0`  
    Expects a real value. Cost of starting a track. Default value is 0.
  - `OUTPUT_COST =0.0`  
    Expects a real value. Cost of terminating a track. Default value is 0.
  - `MAX_TIMEGAP_COMPLETE = 60`  
    Expects a positive integer. Influences the second step of the two-step procedure. If `COMPLETE_GAP_IN_TRACKLET` is set to one, this value is used for selection of base and lifted edges in creating the tracklet graph. Its default value is the maximum from `MAX_TIMEGAP_BASE` and `MAX_TIMEGAP_LIFTED`.
  - `USE_ADAPTIVE_THRESHOLDS=1`
    Expects value 0/1. If set to 1, three specific parameters expecting thresholds on edge cost will expect real values within the interval [0,...,1) denoting the rate of edges to keep instead of exact threshold values. If set to 0, the three parameters expect exact threshold values. These parameters are `BASE_THRESHOLD, NEGATIVE_THRESHOLD_LIFTED, POSITIVE_THRESHOLD_LIFTED`. Default value is 0.   
  - `TIGHT_MAX_EDGE_USAGE=4`
    Expects a positive integer. The maximal number of newly added path or cut subproblems that can share one edge. Default value is 4.
  - `TIGHT_MIN_IMPROVEMENT=0.00001`
    Expects a positive real number. A new path or a new cut subproblem is added to the optimization only if the expected lower bound improvement is higher or equal to this threshold. The lower bound improvement is considered before rescaling for compensation of edge sharing among multiple new subproblems. The default value is 0.00001.
    
  - `PRIMAL_HEURISTIC_ITERATIONS=1`
    Expects value 0/1. If set to 1, local search heuristic is applied after the MCF heuristic for finding primal solutions. It is recommended because the local search improves the primal solution significantly. Default value is 1. 

  - `REPAM_COST_IN_HEURISTIC=0`
    Expects value 0/1. If set to 1, reparametrized cost is used in the local search heuristic instead of the original cost. Default value is 0.

  - `USE_PRE_ITERATE=0`
    Expects value 0/1. If set to 1, an extra MCF subproblem is created and an extra message passing step is used before every MP iteration where messages are exchanged between the MCF subproblem and the inflow and outflow subproblems. However, it does not seem to improve neither the lower bound, nor the upper bound. Therefore, its default value is 1. 

  - `MERGE_THRESHOLD=0.25`
    Expects a real value in the interval [0,...1). Applies in the local search heuristic for deciding if two paths can be merged. Two paths are merged only if the sum of positive edges cost between the two paths is lower or equal than the sum of absolute values of the negative edges cost between the paths multiplied by the ``MERGE_THRESHOLD``. Default value is 0.25.

  - `SAVE_INTERMEDIATE`
    Expects value 0/1. If set to 1, the solver stores intermediate primal solutions. That is, each primal solution having better cost than the previously found primal solutions is saved. The location for saving is derived from the command line parameter ``-o`` of the solver and the iteration number when the solution was obtained.  For instance, if the solver is run with the following command
    ``./lifted_disjoint_paths_text_input -i /path/to/your/input/inputFile.txt -o /path/to/your/output/outputFile.txt``
    and an intermediate primal solution is found in iteration 11, it is stored in file ``/path/to/your/output/outputFile-iter-11.txt``. Default value is 1. 

**Base graph sparsificatioin and cost adjustment parameters**

  - `MAX_TIMEGAP_BASE = 50`  
    Expects a positive integer. Maximal time gap for the base edges (amount of time frames).
    
  - `KNN_K = 3`  
    Expects a positive integer. For every vertex and every time gap, up to `KNN_K` nearest neighbors can be connected by a base edge with that vertex.
    
  - `KNN_GAP = 6`  
    Expects a non-negative integer. Up to `KNN_GAP`, all found nearest neighbors are connected with a processed vertex. Beyond `KNN_GAP`, edges to nearest neighbors have to additionally satisfy other constraints in order to be used.
    
  - `BASE_THRESHOLD = 0.0`  
    Expects a real value. Upper threshold for cost of base edges with a bigger time gap than `KNN_GAP`. Edges longer than `KNN_GAP` have to have cost lower than this threshold in order to be used. Default value is 0.
    
  - `ALL_BASE_TO_ZERO = 1`  
    Expects value 0/1. If set to 1 (default value),  the costs of base edges between non-consecutive frames are set to zero after the base graph sparsification. The purpose is to not duplicate the cost contribution if a base edge overlaps with a lifted edge. On the other hand, some lifted edges are removed during lifted graph sparsification and we want to ensure that the activation of each base edge contributes to the objective value. Parameter `COVER_BASE_WITH_LIFTED` offers two possibilities for dealing with such cases. If `ALL_BASE_TO_ZERO = 0`, the costs of base edges do not change after the base graph sparsification.
    
  - `COVER_BASE_WITH_LIFTED=1`

    Expects value 0/1. This parameter applies only if `ALL_BASE_TO_ZERO=1`. It ensures that an activation of each base edge contributes to the objective value. If set to 1, each zero-valued base edge is covered by a lifted edge. That is, if the cost of a base edge is set to zero, the lifted edge overlapping with this base edge is never removed during lifted graph sparsification. If set to 0, some of the zero-valued base edges get back their non-zero cost. That is, if the cost of a base edge is set to zero and the overlapping lifted edge is removed, the base edge gets its cost back. Default value is 1.

**Lifted graph sparsification and graph adjustment parameters**

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
    
  - `MISSING_AS_MUST_CUT=0` 

    Expects value 0/1. If set to `MISSING_AS_MUST_CUT=1`, it enables to save the input data size by omitting pairs of detections that are obviously non-matching. If the cost between a pair of detections within the maximal time distance is not provided in the input data, the cost of the lifted edge connecting the two detections is automatically set to a high positive value given by the parameter `MUST_CUT_PENALTY`. Like other lifted edges, these special lifted edges are added only between detections that are reachable in the base graph and undergo further sparsifcation.  The maximal time distance for which these lifted edges are created is determined by the minimum of the two values (i) `MAX_TIMEGAP_LIFTED` and (ii) the maximal time gap between two detections in the input data. If `MISSING_AS_MUST_CUT=0`, missing values in the input data for a pair of vertices means that no edge (neither base, nor lifted) is created between them.   

  - `MUST_CUT_PENALTY=10.0`
    Expects a real value. This parameter applies only if `MISSING_AS_MUST_CUT=1`. It gives the cost value of the extra added lifted edges. Default value is 10.

### References
* [1]: [`A. Hornakova, R. Henschel,  B. Rosenhahn,  P. Swoboda. Lifted Disjoint Paths with Application in Multiple Object Tracking, ICML 2020`](http://proceedings.mlr.press/v119/hornakova20a.html)
* [2]: [`A. Hornakova, T. Kaiser, P. Swoboda,  M. Rolinek,  B. Rosenhahn, R. Henschel. Making Higher Order MOT Scalable: An Efficient Approximate Solver for Lifted Disjoint Paths, ICCV 2021`](https://arxiv.org/abs/2108.10606)

