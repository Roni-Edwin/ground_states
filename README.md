# Abstract

 We're exploring the ground state configurations at high densities under the potentials $f_\alpha:x \mapsto \left(x^\alpha+1\right)^{-1}$, where $\alpha>2$ is an even integer.
For a configuration of points $X=\left(x_n\right)_{n=-\infty}^\infty$ with fixed density $\rho$, given by the following limit
```math
\rho=\lim_{r \to \infty}\frac{\left\{i:\left|x_{i}\right|\le r\right\}}{2r}.
```
We defined its (lower) $f_\alpha$-energy $E_{f_\alpha}(X)$ as the following limit inferior:  
```math
E_{f_\alpha}(X)=\liminf_{r \to \infty}\frac1{\left|X_r\right|}\sum_{\substack{i,j \in X_r \\ i\neq j}}f_\alpha\left(\left|x_i-x_j\right|\right) \text{ where } X_r=\left\{i:\ \left|x_{i}\right|\le r\right\}.
```
You can think of this as the average potential energy per particle of the points in $X$. The problem then is given a fixed density $\rho$, to find the configuration of points $X_\rho$ that minimises lower $f_\alpha$ energy, particularly for large values of $\rho$. One way to do analyse this, as $\rho \to \infty$, is to extend this idea of average energy to measure the average potential energy of a borel measure $\mu$ on $\mathbb{R}$. For such a measure $\mu$ with $\frac{\mu\left([-r,r]\right)}{2r} \to 1$ as $r \to \infty$, we define its continuous $f_\alpha$-energy, $\mathcal{E}_{f}\left(\mu\right)$, by 
```math
\mathcal{E}_{f_\alpha}(\mu)=\liminf_{r \to \infty}\frac1{\mu\left([-r,r]\right)}\int_{[-r,r]}\int_{[-r,r]}f_\alpha\left(|x-y|\right)d\mu(x)d\mu(y).
```
The goal then is to find a borel measure $\mu_{opt}$, with  $\frac{\mu_{opt}\left([-r,r]\right)}{2r} \to 1$ as $r \to \infty$, that minimises continuous $f_\alpha$-energy. The idea here is that the discrete problem 'discretizes' the continuous version, so then we get better approximations as the density $\rho$ tends to $\infty$. Experimentally, it seems that the optimal measure $\mu_{opt,\alpha}$ is the counting measure on the lattice $s_\alpha\mathbb{Z}$, normalise to have average density $1$. So
```math
\mu_{opt,\alpha}\left(A\right)=s_\alpha \sum_{n =-\infty}^\infty \delta\left(s_\alpha n\right),
```
where $\delta(x)$ is the dirac measure at $x$. This repository contains the necessary code used in proving this. Here is the original arXiv eprint for reference: https://doi.org/10.48550/arXiv.2405.11428


