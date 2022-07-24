The following repository contains an implementation of the [standard Lanczos Algorithm](https://en.wikipedia.org/wiki/Lanczos_algorithm).

# Installation

To play around with the algorithm yourself, `git clone` the repository onto your local machine. Then, in the cloned directory, initialize the conda environment by running
```
conda env create --file environment.yml
```
From there, you can activate the environment like so
```
conda activate LanczosAlgo
```
and run the notebook in `src/lanczos.ipynb`.

# Introduction

The Lanczos algorithm is a direct algorithm devised by Cornelius Lanczos that is an adaptation of power methods to find the $m$ "most useful" (tending towards extreme highest/lowest) eigenvalues and eigenvectors of an $n\times n$ Hermitian matrix, where $m$ is often but not necessarily much smaller than $n$.

# Algorithm

The algorithm proceeds as follows:

1. Given a Hermitian matrix $A$ of size $n\times n$, and an arbitrary vector $v_1$ with Euclidean norm 1, specify a default number of function calls $m=n$
2. Let $w_1'$=$Av_1$
3. Let $\alpha_1=w_1'^* v_1$
4. Let $w_1=w_1'^*-\alpha_1 v_1$

We refer to steps 2 - 4 as the first iteration steps. Subsequently, 

5. Let $\beta_j=||w_{j-1}||$
6. If $\beta_j\neq0$, let $v_j=\frac{w_{j-1}}{\beta_j}$. Else, let $v_j$ be an arbitrary vector with Euclidean norm 1 that is orthogonal to $v_1, ..., v_{j-1}$
7. Let $w_j'=Av_j$
8. Let $\alpha_j=w_j'v_j$
9. Let $w_j=w_j'-\alpha_j v_j-\beta_j v_{j-1}$

where $j$ denotes the iteration number, and must satisfy $2\leq j\leq m$. Finally, the output is a tridiagonalized matrix $T$ with $\alpha_1, ..., \alpha_m$ along the main diagonal, and $\beta_2, ..., \beta_m$ along the super- and subdiagonals.

# Future Work

Efforts are currently ongoing to extend the algorithm to include scattering-state systems. Star this repo if you'd like to be notified when that goes live!
