LPMP
========

[![Build Status](https://travis-ci.org/LPMP/LPMP.svg?branch=master)](https://travis-ci.org/LPMP/LPMP)

LPMP is a C++ framework for developing scalable dual (Lagrangean) decomposition based solvers for a wide range of LP-relaxations to discrete optimization problems.
For a theoretical introduction to the techniques used and the class of problems that can be optimized see [1,2].

## Solvers
We provide a range of solvers for various discrete optimization problems, including
* **[Discrete graphical models](/include/mrf)**,
* **[Multicut](/include/multicut)**, [3]
* **[Max-cut](/include/max_cut)**, 
* **[Graph matching](include/graph_matching)**, [4]
* **[Multi-graph matching](/include/multigraph_matching)**, [5]

Benchmark problems for various solvers above can be found in [Multi-graph matching](/datasets).

## Optimization techniques
Optimization techniques include
* **Messsage passing [1]**,
* **Subgradient ascent with a proximal bundle method based on the Frank-Wolfe algorithm [2]**, [Vladimir Kolmogorov's](http://http://pub.ist.ac.at/~vnk/) [original implementation](http://pub.ist.ac.at/~vnk/papers/FWMAP.html).
* An interface to **external solvers** is provided by [DD_ILP](https://github.com/pawelswoboda/DD_ILP).

## Installation
Type `git clone https://github.com/LPMP/LPMP.git` for downloading, then `cd LPMP` and `git submodule update --init` for downloading dependencies` and finally `cmake .` for building.

Prerequisites:
* Clang 5.0 or GCC 7.0 upwards for C++17 compatibility.

## References
* [1]: [`P. Swoboda, J. Kuske and B. Savchynskyy. A Dual Ascent Framework for Lagrangean Decomposition of Combinatorial Problems. In CVPR 2017.`](http://openaccess.thecvf.com/content_cvpr_2017/html/Swoboda_A_Dual_Ascent_CVPR_2017_paper.html)
* [2]: `P. Swoboda and V. Kolmogorov. MAP inference via Block-Coordinate Frank-Wolfe Algorithm. In CVPR 2019`
* [3]: [`P. Swoboda and B. Andres. Message passing for multicut. In CVPR 2017.`]()
* [4]: [`P. Swoboda and B. Savchynskyy. A study of. In CVPR 2017.`]()
* [5]: `P. Swoboda, D. Kainmueller, A. Mokarian, C. Theobalt and F. Bernard. A convex approach to multi-graph matching. In CVPR 2019.`
