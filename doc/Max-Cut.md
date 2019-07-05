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

## References

[1]: `Liers, M. J ̈unger, G. Reinelt, and G. Rinaldi. Computing Exact Ground States of Hard Ising Spin Glass Problems by
Branch-and-Cut, chapter 4, pages 47–69. Wiley-Blackwell, 2005`
