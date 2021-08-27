from .raw_solvers import mc_solver_gaec, mc_solver_gef, mc_solver_bec

try:
    import torch
except ImportError as e:
    import warnings
    warnings.warn("Missing Torch, skipping import of torch modules")
else:
    from .multicut import *
