# Simulations

A quick recap of what we're doing. We're exploring the ground state configurations at high densities under the potential $\phi:x \mapsto \frac1{x^4+1}$. For a configuration of points $X=\left(x_n\right)_{n=-\infty}^\infty$ with fixed density $\rho$, given by the following limit
```math
\rho=\lim_{r \to \infty}\frac{\left\{i:\left|x_{i}\right|\le r\right\}}{2r},
```
we defined its (lower) $\phi$-energy $E(X)$ as the following limit inferior:  
```math
E(X)=\liminf_{r \to \infty}\frac1{\left|X_r\right|}\sum_{\substack{i,j \in X_r \\ i\neq j}}\phi\left(\left|x_i-x_j\right|\right) \text{ where } X_r=\left\{i:\ \left|x_{i}\right|\le r\right\}.
```
You can think of this as the average potential energy per particle of the points in $X$. The problem then is given a fixed density $\rho$, to find the configuration of points $X_\rho$ that minimises $E$, particularly for large values of $\rho$. One way to do analyse this, as $\rho \to \infty$, is to extend this idea of average energy to measure the average potential energy of a borel measure $\mu$ on $\mathbb{R}$. For such a measure $\mu$ with $\frac{\mu\left([-r,r]\right)}{2r} \to 1$ as $r \to \infty$, we define its average energy (under $\phi$) $\mathcal{E}(\mu)$ by 
```math
\mathcal{E}(\mu)=\liminf_{r \to \infty}\frac1{\mu\left([-r,r]\right)}\int_{[-r,r]}\int_{[-r,r]}\phi\left(|x-y|\right)d\mu(x)d\mu(y).
```
The goal then is to find a borel measure $\mu_{opt}$, with  $\frac{\mu_{opt}\left([-r,r]\right)}{2r} \to 1$ as $r \to \infty$, that minimises $\mathcal{E}$. The idea here is that the discrete problem 'discretizes' the continuous version, so then we get better approximations as the density $\rho$ tends to $\infty$. Experimentally, it seemed that the optimal measure $\mu_{opt}$ is the counting measure on the lattice $\sqrt{2}\mathbb{Z}$, normalised by $\sqrt{2}$ to have average density $1$. Explicitly, for each set $A \subset \mathbb{R}$, we set 
```math
\mu_{opt}\left(A\right)=\sqrt{2}\cdot \left|\left\{n:\sqrt{2}n \in A\right\}\right|
```
One way to test this, experimentally, is to consider the compact version of the minimisation problem. That is, a given radius $r$ and number of points $N$, the goal is to find a configuration of $N$ points $x_1,x_2, . . . , x_N$ that minimise the average potential energy $E$ given by 
```math
E=\frac1{N}\sum_{\substack{i,j=1 \\ i\neq j}}^Np\left(\left|x_i-x_j\right|\right), \text{ with } \left|x_i\right|\le r.
```
The idea here is that for large values of $r$, and $\frac{N}{r}$, the distribution of points should approximate the 'mass distribution' of the optimal measure $\mu^*$. I implemented this using a basic gradient descent algorithm and periodically plotted the current distribution of points. The points tend to gather into equally spaced clusters, which would seem to support the idea that $\mu_{opt}$ is as described above. Another behavior observed is that this behavior seems to hold for potentials $p$ of the form $p(x)=\frac1{x^{4n}+1}$, for $n \in \mathbb{N}$. That is, the optimal configuration for large values of $r$ and $\frac{N}{r}$ is that the points tend to gather into equally spaced clusters.



# Interval Arithmetic
A recap and the code for Interval Arithmetic. We have the function $\mathcal{R}$ given by 
```math
\mathcal{R}(\xi):=\sum_{n=-\infty}^\infty H(\xi,n)
```
where $H(\xi,n)$ for $\xi\ge 0$, $n \in \mathbb{Z}$ is given by
```math
H(\xi,n)=\frac{2\xi^2\left(\widehat{\phi}(\xi)-\widehat{\phi}\left(\frac{n}{\sqrt{2}}\right)\right)}{\left(n-\sqrt{2}\xi\right)^2}.
```
Our goal in this part is to show that 
```math
\mathcal{R}(\xi)\ge -\sum_{n=-\infty}^\infty\widehat{\phi}\left(\frac{n}{\sqrt{2}}\right)=-\frac{\pi}{\sqrt{2}}\tanh\left(\frac{\pi}{\sqrt{2}}\right)
```
for $0\le \xi\le 8$. We defined the truncated series for $\mathcal{R}(\xi)$, $\mathcal{R}_p(\xi)$ given by 
```math
\mathcal{R}_p(\xi)=\sum_{|n|\le p}H(\xi,n).
```

We then showed that if $p>11$, then 
```math
\mathcal{R}(\xi)\ge \mathcal{R}_p(\xi)-B(\xi,p)
```
where 
```math
 B(\xi,p):=\frac{4p\xi^2\left|\widehat{\phi}(\xi)\right|}{p^2-2\xi^2}+2\xi^2 e^{-\pi p}\left(\frac{0.1}{p-\sqrt{2}\xi}+\frac{0.15}{\left(p+\sqrt{2}\xi\right)^2}\right)
```
So we will now show that 
```math
\mathcal{R}_p(\xi)-B(\xi,p)\ge -\frac{\pi}{\sqrt{2}}\tanh\left(\frac{\pi}{\sqrt{2}}\right)
```
for $0\le \xi\le 8$. To avoid large uncertainties in the computed value of $H(\xi,n)$ (the denominator might get really small), we introduced a lower bound for $H(\xi,n)$, $\widetilde{H}(\xi,n)$, given by 
```math
H_{l}(\xi,n)=-\frac{n^2\sqrt{2}\pi^3\left(-1\right)^{n}e^{-\pi n}}{2}-\left|\xi-\frac{n}{\sqrt{2}}\right|\left(\frac{\xi^2}{6}e^{\pi\left|\sqrt{2}\xi-n\right|}\cdot 8\pi^4e^{-\sqrt{2}\pi \xi}+2\pi^3e^{-\pi n}\left(\xi+\frac{n}{\sqrt{2}}\right)\right).
```
We showed that $H(\xi,n)\ge H_{l}(\xi,n)$ for $\xi\ge 0$ and $n\ge 0$ an integer. Choose some precision value $\varepsilon>0$, and define $\widetilde{H}(\xi,n)$ by 
```math
\widetilde{H}(\xi,n)=\begin{cases}
    H(\xi,n) & \text{if } \left|n-\sqrt{2}\xi\right|\ge \varepsilon \text{ or } n<0 \\
     H_{l}(\xi,n)& \text{if } \left|n-\sqrt{2}\xi\right|<\varepsilon,
\end{cases}
```
 so that $H(\xi,n)\ge \widetilde{H}(\xi,n)$, and consequently, 
 ```math
\mathcal{R}_p(\xi)\ge \sum_{|n|\le p}\widetilde{H}(\xi,n) 
```
So our new goal is to show that 
```math
\sum_{|n|\le p}\widetilde{H}(\xi,n)-B(\xi,p)\ge -\frac{\pi}{\sqrt{2}}\tanh\left(\frac{\pi}{\sqrt{2}}\right)
```
for $0\le \xi\le p$. In the code, we set $\varepsilon=2^{-6}$. 

To do this we use the IntervalArithemtic package. The way it works basically, is that once you've defined your function $f$, you pass in you interval $I$ and returns an interval $J$ such that $f\left(I\right)\subset J$. Here attention to detail is important, especially when working with numbers that aren't representable exactly with floating point arithmetic, such as $\pi$, $e$, and $\frac{4}{3}$. So when implenting the function $f$, it's important to make sure you represent these constants as degenrate intervals. So for example, you would use $[\pi,\pi]$ (instead of just $\pi$), which Julia would represent as a small interval containing the true value of $\pi$, rather than its floating point representation. With this implementation, the general outline of the code is as follows. We have our function $f$. Let $f^*(I)$ denote the interval returned by Julia, so we wish to show  and wish to show $f\left(I\right)> [l,l]$. For two intervals $[a,b]$, $[c,d]$, $[a,b]>[c,d]$ means $a>d$. The way the code works is like so. Suppose $I=[a,b]$. Intialize $d$ to $b-a$, and repeatedly check if $f\left([a,a+d]\right)>[l,l]$, and if not, reassign $l \to \frac{l}{2}$. Once we find some $l$ for which this works, we reassing $a \to l$, and repeat until we've covered the whole interval $[a,b]$. 
