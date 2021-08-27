LPMP
========

<!--- [![Build Status](https://travis-ci.org/LPMP/LPMP.svg?branch=master)](https://travis-ci.org/LPMP/LPMP) --->

LPMP is a C++ framework for developing scalable dual (Lagrangean) decomposition based solvers for a wide range of LP-relaxations to discrete optimization problems.
For a theoretical introduction to the techniques used and the class of problems that can be optimized see [1,2].

## Solvers
We provide a range of solvers for various discrete optimization problems, including
* **[Discrete graphical models](/include/mrf)**,
* **[Multicut](/include/multicut)**, [3]
* **[Max-cut](/include/max_cut)**, 
* **[Graph matching](include/graph_matching)**, [4]
* **[Multi-graph matching](/include/multigraph_matching)**, [5]
* **[Discrete graphical models with bottleneck terms](/include/horizon_tracking)**, [6]
* **[Lifted disjoint paths](include/lifted_disjoint_paths)**.

Benchmark problems for various solvers above can be found in [datasets](/datasets).

## Optimization techniques
Optimization techniques include
* **Messsage passing [1]**,
* **Subgradient ascent with a proximal bundle method based on the Frank-Wolfe algorithm [2]**, [Vladimir Kolmogorov's](http://http://pub.ist.ac.at/~vnk/) [original implementation](http://pub.ist.ac.at/~vnk/papers/FWMAP.html).
* An interface to **external solvers** is provided by [DD_ILP](https://github.com/pawelswoboda/DD_ILP).

## Differentiable wrappers

The solvers can be wrapped as differentiable PyTorch modules using the technique of [7]. Currently, wrappers are available for graph matching, multigraph matching and lifted disjoint paths solvers. For usage examples see an application to keypoint matching [8] ([code](https://github.com/martius-lab/blackbox-deep-graph-matching)) or the general [repository](https://github.com/martius-lab/blackbox-backprop) of [7].

All these can be `pip` installed with

```python3 -m pip install git+https://github.com/LPMP/LPMP.git```

In order to install only graph matching and multigraph matching solvers, type

`PACKAGES="gm" python3 -m pip install git+https://github.com/LPMP/LPMP.git`

In order to install only the lifted disjoint paths solver, type

`PACKAGES="ldp" python3 -m pip install git+https://github.com/LPMP/LPMP.git`

In order to install only the all variants of multicut solvers (multicut, asymmetric multiway cut [9], multiway cut [9]), type

`PACKAGES="mc" python3 -m pip install git+https://github.com/LPMP/LPMP.git`

In order to get the precise version of graph matching as used in [8], type

```python3 -m pip install git+https://github.com/lpmp/LPMP.git@keypiont_submission```

NOTE:
If you already installed with one of the above options and would like to change the option or if you have an older version installed, you must add `--upgrade` argument to `pip`. So, you type for instance

`PACKAGES="ldp" python3 -m pip install git+https://github.com/lpmp/LPMP.git --upgrade`

## Installation
Type `git clone https://github.com/LPMP/LPMP.git` for downloading, then `cd LPMP` and `git submodule update --init --remote --recursive` for downloading dependencies and finally `cmake .` for building.

Prerequisites:
* Clang 5.0 or GCC 8.0 upwards for C++17 compatibility (see [here](https://solarianprogrammer.com/2016/10/07/building-gcc-ubuntu-linux/) for installation instructions).
* HDF5 (install with `apt install libhdf5-serial-dev`)
* cmake (install with `apt install cmake`)


## Documentation

A tutorial on writing a new solver from scratch can be found [here](/doc/Getting-Started.md).
The documentation for using all the existing solvers is in the directory [doc](/doc).

## References
* [1]: [`P. Swoboda, J. Kuske and B. Savchynskyy. A Dual Ascent Framework for Lagrangean Decomposition of Combinatorial Problems. In CVPR 2017.`](http://openaccess.thecvf.com/content_cvpr_2017/html/Swoboda_A_Dual_Ascent_CVPR_2017_paper.html)
* [2]: [`P. Swoboda and V. Kolmogorov. MAP inference via Block-Coordinate Frank-Wolfe Algorithm. In CVPR 2019`](http://openaccess.thecvf.com/content_CVPR_2019/html/Swoboda_MAP_Inference_via_Block-Coordinate_Frank-Wolfe_Algorithm_CVPR_2019_paper.html)
* [3]: [`P. Swoboda and B. Andres. A Message Passing Algorithm for the Minimum Cost Multicut Problem. In CVPR 2017.`](http://openaccess.thecvf.com/content_cvpr_2017/html/Swoboda_A_Message_Passing_CVPR_2017_paper.html)
* [4]: [`P. Swoboda, C. Rother, H. A. Alhaija, D. Kainmuller, B. Savchynskyy. A Study of Lagrangean Decompositions and Dual Ascent Solvers for Graph Matching. In CVPR 2017.`](http://openaccess.thecvf.com/content_cvpr_2017/html/Swoboda_A_Study_of_CVPR_2017_paper.html)
* [5]: [`P. Swoboda, D. Kainmueller, A. Mokarian, C. Theobalt and F. Bernard. A convex approach to multi-graph matching. In CVPR 2019.`](http://openaccess.thecvf.com/content_CVPR_2019/html/Swoboda_A_Convex_Relaxation_for_Multi-Graph_Matching_CVPR_2019_paper.html)
* [6]: [`A. Abbas, P. Swoboda. Bottleneck Potentials in Markov Random Fields. In ICCV 2019.`](http://openaccess.thecvf.com/content_ICCV_2019/html/Abbas_Bottleneck_Potentials_in_Markov_Random_Fields_ICCV_2019_paper.html)
* [7]: [`M. Vlastelica, A. Paulus, V. Musil, G. Martius, M. Rolínek. Differentiation of Blackbox Combinatorial Solvers. In ICLR 2020.`](https://openreview.net/forum?id=BkevoJSYPB)
* [8]: [`M. Rolínek, P. Swoboda, D. Zietlow, A. Paulus, V. Musil, G. Martius. Deep Graph Matching via Blackbox Differentiation of Combinatorial Solvers.`](https://arxiv.org/abs/2003.11657)
* [9]: [`A. Abbas, P. Swoboda. Combinatorial Optimization for Panoptic Segmentation: An End-to-End Trainable Approach.`](https://arxiv.org/abs/2106.03188)