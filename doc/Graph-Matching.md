## Graph matching

### Problem

The graph matching problem can be stated as

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{X&space;\in&space;\mathbb{P}_{n&space;\times&space;m}}&space;\text{vec}(X)^{\top}&space;W&space;\text{vec}(X)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{X&space;\in&space;\mathbb{R}^{n&space;\times&space;m}}&space;\text{vec}(X)^{\top}&space;W&space;\text{vec}(X)" title="\min_{X \in \mathbb{R}^{n \times m}} \text{vec}(X)^{\top} W \text{vec}(X)" /></a>

where the set of partial assignments is

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{P}_{n&space;\times&space;m}&space;=&space;\{&space;X&space;\in&space;\{0,1\}^{n&space;\times&space;m}&space;:&space;X&space;\mathbbmss{1}&space;\leq&space;\mathbbmss{1},&space;X^\top&space;\mathbbmss{1}&space;\leq&space;\mathbbmss{1}&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{P}_{n&space;\times&space;m}&space;=&space;\{&space;X&space;\in&space;\{0,1\}^{n&space;\times&space;m}&space;:&space;X&space;\mathbbmss{1}&space;\leq&space;\mathbbmss{1},&space;X^\top&space;\mathbbmss{1}&space;\leq&space;\mathbbmss{1}&space;\}" title="\mathbb{P}_{n \times m} = \{ X \in \{0,1\}^{n \times m} : X \mathbbmss{1} \leq \mathbbmss{1}, X^\top \mathbbmss{1} \leq \mathbbmss{1} \}" /></a>

We decompose the cost into the following terms

<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{ij&space;\in&space;A}&space;c_{ij}&space;X_{ij}&space;&plus;&space;\sum_{(ij,kl)&space;\in&space;E}&space;d_{ijkl}&space;X_{ij}&space;X_{kl}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{ij&space;\in&space;A}&space;c_{ij}&space;X_{ij}&space;&plus;&space;\sum_{(ij,kl)&space;\in&space;E}&space;d_{ijkl}&space;X_{ij}&space;X_{kl}" title="\sum_{ij \in A} c_{ij} X_{ij} + \sum_{(ij,kl) \in E} d_{ijkl} X_{ij} X_{kl}" /></a>

The elements in A correspond to non-zero entries on the diagonal of W and the elements in E to non-zero entries on the off-diagonal of W.

### File format

We use the file format used by the dual decomposition based graph matching solver by Torresani et al, [corresponding publication](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6197199).
Namely, it is structured as follows:

```
p <n> <m> <#A> <#E>   // n, m, # elements of A, # elements of E
a <a> <i> <j> {cost}  // specify assignment
...

e <a0> <a1> {cost}    // nodes ij and kl come from the assignments endpoints
...
```

### Datasets

* [C.elegans annotation dataset (worms)](https://datarep.app.ist.ac.at/57/1/wormMatchingProblems.zip), 
by Kainmueller et al, [corresponding publication](http://dx.doi.org/10.1007/978-3-319-10404-1_11).

* [GraphFlow – 6D Large Displacement Scene Flow via Graph Matching](https://datarep.app.ist.ac.at/id/eprint/82) by Alhaija et al, [corresponding publication](https://link.springer.com/chapter/10.1007/978-3-319-24947-6_23).

* [cars and motor dataset](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/car_motor_graph_matching.zip) by M. Leordeanu and M. Hebert, [corresponding publication](https://ieeexplore.ieee.org/document/5206533/).

