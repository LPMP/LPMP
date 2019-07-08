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

