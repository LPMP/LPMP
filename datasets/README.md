# Datasets

## Markov Random Fields

### Problem

Given a Markov Random Field (MRF), i.e. a graph G=(V,E) with a label space 
<a href="https://www.codecogs.com/eqnedit.php?latex=X_i\&space;\forall&space;i&space;\in&space;V" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_i\&space;\forall&space;i&space;\in&space;V" title="X_i\ \forall i \in V" /></a>,
unary potentials
<a href="https://www.codecogs.com/eqnedit.php?latex=\theta_i&space;:&space;X_i&space;\rightarrow&space;\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_i&space;:&space;X_i&space;\rightarrow&space;\mathbb{R}" title="\theta_i : X_i \rightarrow \mathbb{R}" /></a>
and pairwise potentials 
<a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{ij}&space;:&space;X_i&space;\times&space;X_j&space;\rightarrow&space;\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{ij}&space;:&space;X_i&space;\times&space;X_j&space;\rightarrow&space;\mathbb{R}" title="\theta_{ij} : X_i \times X_j \rightarrow \mathbb{R}" /></a>,
the associated Maximum-A-Posteriori (MAP) inference problem is

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{x&space;\in&space;\prod_{i&space;\in&space;V}&space;X_i}&space;\sum_{i&space;\in&space;V}&space;\theta_i(x_i)&space;&plus;&space;\sum_{ij&space;\in&space;E}&space;\theta_{ij}(x_i,x_j)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{x&space;\in&space;\prod_{i&space;\in&space;V}&space;X_i}&space;\sum_{i&space;\in&space;V}&space;\theta_i(x_i)&space;&plus;&space;\sum_{ij&space;\in&space;E}&space;\theta_{ij}(x_i,x_j)" title="\min_{x \in \prod_{i \in V} X_i} \sum_{i \in V} \theta_i(x_i) + \sum_{ij \in E} \theta_{ij}(x_i,x_j)" /></a>

### File format

We use the [UAI file format](http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php).

### Datasets

* [QPBO models from various computer vision tasks](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/QPBO_CV_problems.zip), collected by [Dhruv Batra](https://ttic.uchicago.edu/~dbatra/research/mfcomp/)

---

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

---

## Multi-graph matching

The multigraph matching problem can be stated as a number of pairwise graph matching problems

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{(X^{[pq]}&space;\in&space;\mathbb{P}_{n_p&space;\times&space;n_q})_{p,q&space;\in&space;[d]]}&space;}\sum_{p,q&space;\in&space;[d]}&space;\text{vec}(X^{[pq]})^\top&space;W^{[pq]}&space;\text{vec}(X^{[pq]})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{(X^{[pq]}&space;\in&space;\mathbb{P}_{n_p&space;\times&space;n_q})_{p,q&space;\in&space;[d]]}&space;}\sum_{p,q&space;\in&space;[d]}&space;\text{vec}(X^{[pq]})^\top&space;W^{[pq]}&space;\text{vec}(X^{[pq]})" title="\min_{(X^{[pq]} \in \mathbb{P}_{n_p \times n_q})_{p,q \in [d]]} }\sum_{p,q \in [d]} \text{vec}(X^{[pq]})^\top W^{[pq]} \text{vec}(X^{[pq]})" /></a>

with the additional cycle consistency constraints

<a href="https://www.codecogs.com/eqnedit.php?latex=X^{[pq]}&space;X^{[qr]}&space;\leq&space;X^{[pr]}&space;\quad&space;\forall&space;p,q&space;\in&space;[d]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X^{[pq]}&space;X^{[qr]}&space;\leq&space;X^{[pr]}&space;\quad&space;\forall&space;p,q&space;\in&space;[d]" title="X^{[pq]} X^{[qr]} \leq X^{[pr]} \quad \forall p,q \in [d]" /></a>

### File format

We use a file format similar to the one used for graph matching.

```
gm <p> <q>                  // pairwise graph matching problem between graph p and q
p <n_p> <n_q> <#A> <#E>     // now comes the graph matching problem file format
a <a> <i> <j> {cost}
...
e <a0> <a1> {cost}
...
```

### Datasets

* [C.elegans annotation dataset (worms)](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/worms_mgm.zip)

* [Hotel/House/Synthetic](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/hotel_house_synthetic_mgm.zip) by Florian Bernard et al, [corresponding publication](https://arxiv.org/pdf/1711.10733.pdf).

---

## Multicut

The multicut problem is to find a decomposition of a given weighted graph G=(V,E), namely

<a href="https://www.codecogs.com/eqnedit.php?latex=\Pi&space;=&space;(\Pi_1,\ldots,\Pi_k)\quad&space;\text{s.t.&space;}&space;\Pi_i&space;\cap&space;\Pi_j&space;=&space;\varnothing&space;\text{&space;for&space;}&space;i&space;\neq&space;j,\&space;\Pi_1&space;\cup&space;\ldots&space;\cup&space;\Pi_k&space;=&space;V" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Pi&space;=&space;(\Pi_1,\ldots,\Pi_k)\quad&space;\text{s.t.&space;}&space;\Pi_i&space;\cap&space;\Pi_j&space;=&space;\varnothing&space;\text{&space;for&space;}&space;i&space;\neq&space;j,\&space;\Pi_1&space;\cup&space;\ldots&space;\cup&space;\Pi_k&space;=&space;V" title="\Pi = (\Pi_1,\ldots,\Pi_k)\quad \text{s.t. } \Pi_i \cap \Pi_j = \varnothing \text{ for } i \neq j,\ \Pi_1 \cup \ldots \cup \Pi_k = V" /></a>

The number of components is determined as part of the optimization process
The objective function is given by the weighted sum of edges straddling distinct components, i.e.

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{\Pi&space;=&space;(\Pi_1,\ldots,\Pi_k)}&space;\sum_{ij&space;\in&space;E,&space;ij&space;\notin&space;\Pi_l&space;\forall&space;l&space;\in&space;[k]}&space;c_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{\Pi&space;=&space;(\Pi_1,\ldots,\Pi_k)}&space;\sum_{ij&space;\in&space;E,&space;ij&space;\notin&space;\Pi_l&space;\forall&space;l&space;\in&space;[k]}&space;c_{ij}" title="\min_{\Pi = (\Pi_1,\ldots,\Pi_k)} \sum_{ij \in E, ij \notin \Pi_l \forall l \in [k]} c_{ij}" /></a>

### File format

The file format for the multicut problem is

```
MULTICUT
<i> <j> <cost>
...
```

### Datasets

* [multicut problems for circuit reconstruction from electron microscopy images from the CREMI challenge](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/CREMI_multicut_nature_methods.zip) by Thorsten Beier et al, [corresponding publication](https://www.nature.com/articles/nmeth.4151)

* [large scale fruit fly brain segmentation problems](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/fruit_fly_brain_segmentation_Pape.zip) by Constantin Pape et al, [corresponding publication](http://openaccess.thecvf.com/content_ICCV_2017_workshops/papers/w1/Pape_Solving_Large_Multicut_ICCV_2017_paper.pdf)

---

## Max-cut

The max-cut problem consists in finding a bipartition of a given graph G=(V,E) such that the weighted sum of edges having endpoints in distinct components is maximized.
The edge weights are arbitrary real numbers.
The problem can be stated formally as

<a href="https://www.codecogs.com/eqnedit.php?latex=\max_{U&space;\subset&space;V}&space;\sum_{ij&space;\in&space;E:&space;i&space;\in&space;V,&space;j&space;\notin&space;V}&space;c_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\max_{U&space;\subset&space;V}&space;\sum_{ij&space;\in&space;E:&space;i&space;\in&space;V,&space;j&space;\notin&space;V}&space;c_{ij}" title="\max_{U \subset V} \sum_{ij \in E: i \in V, j \notin V} c_{ij}" /></a>

### File format

```
MAX-CUT
<i> <j> <cost>
...
```

### Datasets

* [Ising models from physical applications [1]](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/Ising_models.zip)

* [Converted QPBO models from various computer vision tasks](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/QPBO_CV_problems.zip), collected by [Dhruv Batra](https://ttic.uchicago.edu/~dbatra/research/mfcomp/)

---

## Discrete tomography

### Problem

The discrete tomography problem is, given a graph G=(V,E) to find an assignment of natural numbers to nodes that (i) minimize an energy coming from an MRF on G

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{x&space;\in&space;\{0,1,\ldots,k\}^{\abs{V}}}&space;\sum_{i&space;\in&space;V}&space;\theta_i(x_i)&space;&plus;&space;\sum_{ij&space;\in&space;E}&space;\theta_{ij}(x_i,&space;x_j)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{x&space;\in&space;\{0,1,\ldots,k\}^{\abs{V}}}&space;\sum_{i&space;\in&space;V}&space;\theta_i(x_i)&space;&plus;&space;\sum_{ij&space;\in&space;E}&space;\theta_{ij}(x_i,&space;x_j)" title="\min_{x \in \{0,1,\ldots,k\}^{\abs{V}}} \sum_{i \in V} \theta_i(x_i) + \sum_{ij \in E} \theta_{ij}(x_i, x_j)" /></a>

and (ii) fulfills constraints that come from tomographic projections:

<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_{i&space;\in&space;V}&space;a_{ij}&space;x_i&space;=&space;b_i&space;\quad&space;\forall&space;j&space;\in&space;[m]&space;\text{&space;where&space;}&space;A&space;\in&space;\mathbb{R}_&plus;^{m&space;\times&space;\abs{V}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_{i&space;\in&space;V}&space;a_{ij}&space;x_i&space;=&space;b_i&space;\quad&space;\forall&space;j&space;\in&space;[m]&space;\text{&space;where&space;}&space;A&space;\in&space;\mathbb{R}_&plus;^{m&space;\times&space;\abs{V}}" title="\sum_{i \in V} a_{ij} x_i = b_i \quad \forall j \in [m] \text{ where } A \in \mathbb{R}_+^{m \times \abs{V}}" /></a>

or more generally a non-linear penalizer that depends on the sum of values as

<a href="https://www.codecogs.com/eqnedit.php?latex=f_j\left(\sum_{i&space;\in&space;V}&space;a_{ij}&space;x_i\right)&space;\quad&space;\forall&space;j&space;\in&space;[m]&space;\text{&space;where&space;}&space;A&space;\in&space;\mathbb{R}_&plus;^{m&space;\times&space;\abs{V}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_j\left(\sum_{i&space;\in&space;V}&space;a_{ij}&space;x_i\right)&space;\quad&space;\forall&space;j&space;\in&space;[m]&space;\text{&space;where&space;}&space;A&space;\in&space;\mathbb{R}_&plus;^{m&space;\times&space;\abs{V}}" title="f_j\left(\sum_{i \in V} a_{ij} x_i\right) \quad \forall j \in [m] \text{ where } A \in \mathbb{R}_+^{m \times \abs{V}}" /></a>

### File format

We extend the UAI file format ([see here](http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php)) by additional fields for the tomographic projection constraints.

```
MARKOV
${UAI File for MRF}

PROJECTIONS
<var_1> + ... + <var_k> = (<val sum=0>, <val sum=1>, ..., <val sum=max_sum>)
...
```

### Datasets

[synthetic multi-label problems](https://datarep.app.ist.ac.at/46/1/discrete_tomography_synthetic.zip)
by Kuske et al.


## References

[1]: `Liers, M. J ̈unger, G. Reinelt, and G. Rinaldi. Computing Exact Ground States of Hard Ising Spin Glass Problems by
Branch-and-Cut, chapter 4, pages 47–69. Wiley-Blackwell, 2005`
