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

