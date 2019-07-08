# Markov Random Fields

* [Problem formulation](#problem-formulation)
* [File format](#file-format)
* [Datasets](#datasets)

## Problem

Given a Markov Random Field (MRF), i.e. a graph G=(V,E) with a label space 
<a href="https://www.codecogs.com/eqnedit.php?latex=X_i\&space;\forall&space;i&space;\in&space;V" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_i\&space;\forall&space;i&space;\in&space;V" title="X_i\ \forall i \in V" /></a>,
unary potentials
<a href="https://www.codecogs.com/eqnedit.php?latex=\theta_i&space;:&space;X_i&space;\rightarrow&space;\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_i&space;:&space;X_i&space;\rightarrow&space;\mathbb{R}" title="\theta_i : X_i \rightarrow \mathbb{R}" /></a>
and pairwise potentials 
<a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{ij}&space;:&space;X_i&space;\times&space;X_j&space;\rightarrow&space;\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{ij}&space;:&space;X_i&space;\times&space;X_j&space;\rightarrow&space;\mathbb{R}" title="\theta_{ij} : X_i \times X_j \rightarrow \mathbb{R}" /></a>,
the associated Maximum-A-Posteriori (MAP) inference problem is

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{x&space;\in&space;\prod_{i&space;\in&space;V}&space;X_i}&space;\sum_{i&space;\in&space;V}&space;\theta_i(x_i)&space;&plus;&space;\sum_{ij&space;\in&space;E}&space;\theta_{ij}(x_i,x_j)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{x&space;\in&space;\prod_{i&space;\in&space;V}&space;X_i}&space;\sum_{i&space;\in&space;V}&space;\theta_i(x_i)&space;&plus;&space;\sum_{ij&space;\in&space;E}&space;\theta_{ij}(x_i,x_j)" title="\min_{x \in \prod_{i \in V} X_i} \sum_{i \in V} \theta_i(x_i) + \sum_{ij \in E} \theta_{ij}(x_i,x_j)" /></a>

## File format

We use the [UAI file format](http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php).

## Datasets

* [QPBO (binary pairwise MRF) models from various computer vision tasks](https://datasets.d2.mpi-inf.mpg.de/discrete_cv_problems/QPBO_CV_problems.zip), collected by [Dhruv Batra](https://ttic.uchicago.edu/~dbatra/research/mfcomp/)

