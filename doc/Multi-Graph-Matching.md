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

