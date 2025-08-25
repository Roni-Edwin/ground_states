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

The first order of business is verifying that $s_\alpha>1$. For $t>0$, define the energy $E_\alpha(t)$ by 
```math
    E_\alpha(t)=\sum_{n \in \mathbb{Z}}\frac{t}{1+t^\alpha n^\alpha},
```
and let $s_\alpha$ minimise $E_\alpha(t)$. Then $s_\alpha >1$. This is proven in https://doi.org/10.48550/arXiv.2405.11428 when $\alpha=4$ and $\alpha\ge 16$, so it remains to check this for $\alpha \in \{6,8,10,12\}$. For each such $\alpha$, the table below shows the range of the candidate value of $\alpha$: 



| $\alpha$ |Lower bound for $$\left(s_{\alpha}\right)^\alpha$$ | Upper bound for $$\left(s_\alpha\right)^\alpha$$ | Upper bound for $$E_\alpha(s_\alpha)$$ |
| ------ | ---- | ----- | ----- |
| $6$ | $7.83668$ | $7.83669$ | $1.73458$ |
|  $8$ | $11.813148$ | $11.813151$ | $1.57506$ |
|  $10$ | $15.830697$ | $15.830699$ | $1.47491$ |
| $12$ | $19.854857$ | $19.854859$ | $1.40585$ |
| $14$ | $23.876137$ | $23.876138$ | $1.355228$ |

We implement $E_\alpha(t)$ numerically by bounding the tail, and to do that, we start by showing $s_\alpha$ does not get arbitrarily small (so we can bound the tail independently of $t$). Note 
```math
    E_\alpha(t)=t+2\sum_{n=1}^\infty \frac{t}{1+t^\alpha n^\alpha}\ge t+2\int_t^\infty \frac1{1+x^\alpha }dx.
```
 The right-hand side above is decreasing in $x$ for $x\le 1$, so if $t<s$ for some parameter $s$ small, this implies 
 ```math
    E_\alpha(t)\ge s+2\int_s^\infty \frac1{1+x^\alpha}dx\ge s+2\sum_{n=2}^\infty \frac{s}{1+s^\alpha n^\alpha}.
```
 When $s$ is small, this is larger than the calculated upper bound for the energy in the table above. To bound the tail of $E_\alpha(t)$, note 
 ```math
    \sum_{|n|\ge p+1}\frac{t}{1+t^{\alpha}n^\alpha}\le \frac1{t^{\alpha-1}}\sum_{n=p+1}^\infty \frac2{n^\alpha}\le \frac2{t^{\alpha-1}(p+1)^\alpha}\left(\frac{\alpha+p}{\alpha-1}\right),
```

For $\alpha\ge 14$, we use the bounds in Proposition 4.4 in https://doi.org/10.48550/arXiv.2405.11428 which says 
```math
s_\alpha^\alpha=\alpha-2+\sqrt{\left(\alpha-2\right)^2-3}+\mathfrak{G}(\alpha),
```
where 
```math
    \left|\mathfrak{G}(\alpha)\right|\le \frac{2\alpha}{0.99^2}\cdot \frac{\left(\alpha+1\right)^2}{2^{\alpha-1}\left(\alpha-1\right)}.
```
 Using this we are able to verify the needed inequalities. We start with those for $\widehat{\psi_\alpha}$:

 <br>
 <br>

 For $\alpha=4,6,8,10$, define a constant $\mathcal{T}(\alpha)$ by 
 ```math
        \mathcal{T}(\alpha)\coloneqq\frac1{2}\left(1-\sum_{n=1}^\infty 2F_\alpha(n)\right)-\frac1{\pi}\sum_{n=2}^\infty \left|F_\alpha'(n)\right|.
```
 Then $\mathcal{T}(\alpha)\ge 0$. 

We just need to bound the tails of the series, which is not too hard. Pick some cutoff positive integer $p$. We know $n\ge p$, $\left|F_\alpha'(n)\right|\le \alpha F_\alpha(n)/p$, and $F_\alpha(n)\le n^{-\alpha}$,  and so 
```math
    \sum_{n=p+1}^\infty\left(F_\alpha(n)+\frac1{\pi}\left|F_\alpha'(n)\right|\right)\le \left(1+\frac{\alpha}{p\pi}\right)\sum_{n=p+1}^\infty F_\alpha(n)\le \sum_{n=p+1}^\infty \frac{1+\frac{\alpha}{\pi p}}{n^\alpha} \le \left(1+\frac{\alpha}{\pi p}\right)\left(\frac{\alpha+p}{\alpha-1}\right)\frac1{(p+1)^\alpha},
```

<br>

The next inequality to verify is as follows:

 Let $\mathcal{R}(x)$ be given by 
 ```math
        \mathcal{R}(x)=\frac{\sin x-x}{x^3}+\frac1{6}
 ```
 for $x\neq 0$, and $\mathcal{R}(0)=0$. For $\alpha =6,8,10$, define a constant $\mathfrak{L}(\alpha)$ by 
 ```math
        \mathfrak{L}(\alpha):=\sum_{n=-\infty}^\infty n^3F_\alpha'(n)\left(-\frac{2}{3}+4\mathcal{R}\left(\pi n \right)\right)-\sum_{n=-\infty}^\infty 2n^2F_\alpha(n).
```
 Then $\mathfrak{L}(\alpha)\ge 0$.

  We need to the bound the tail of the series. We write the series for $\mathfrak{L}(\alpha)$ as 
 ```math
\mathfrak{L}(\alpha)=-\sum_{n=1}^\infty 2n^3F_\alpha'(n)\left(\frac2{3}-4\mathcal{R}(\pi n)\right)-\sum_{n=1}^\infty 4n^2F_\alpha(n).
```
We showed showed in https://doi.org/10.48550/arXiv.2405.11428 that $\mathcal{R}(x)$ is increasing in $|x|$, and $\mathcal{R}(0)=0$, so
```math
    \left|-\sum_{n=p+1}^\infty 2n^3F_\alpha'(n)\left(\frac2{3}-4\mathcal{R}(\pi n)\right)-\sum_{n=p+1}^\infty 4n^2F_\alpha(n)\right|\le \sum_{n=p+1}^\infty \left(-\frac4{3}n^3F_\alpha'(n)+4n^2F_\alpha(n)\right).
```
We also know $-nF_\alpha'(n)\le \alpha F_\alpha(n)$, so this implies 
```math
     \left|\sum_{n=p+1}^\infty \left(-2n^3F_\alpha'(n)\left(\frac2{3}-4\mathcal{R}(\pi n)\right)- 4n^2F_\alpha(n)\right)\right|\le \sum_{n=p+1}^\infty 4n^2F_\alpha(n)\left(1+\frac{\alpha}{3}\right) \le \left(1+\frac{\alpha}{3}\right)\sum_{n=p+1}^\infty \frac4{n^{\alpha-2}}, 
```
and so
```math
     \left|\sum_{n=p+1}^\infty \left(-2n^3F_\alpha'(n)\left(\frac2{3}-4\mathcal{R}(\pi n)\right)- 4n^2F_\alpha(n)\right)\right|\le \frac{4\left(1+\alpha/3\right)}{(p+1)^{\alpha-2}}\left(\frac{\alpha-2+p}{\alpha-3}\right).
```
<br>

We now move on to finishing the proof of $\psi_4(x)\le F_4(x)$. We want to show $\psi_4(x)\le F_4(x)$ for $0\le x\le 9$, which we showed in https://doi.org/10.48550/arXiv.2405.11428 is equivalent to showing the following:
  Let $F_4(x)=\left(4x^4+1\right)^{-1}$. Then for $0\le x\le 9$, 
  ```math
        \sum_{n \in \mathbb{Z}}\left(\frac{F_4(x)-F_4(x)-F_4'(n)(x-n)}{(x-n)^2}\right)\ge 0.
```
We start by bounding the tail of this series. Let $\eta(x)$ be $x$ rounded up to the nearest integer, and let
```math
    H(x,n)=\frac{F_4(x)-F_4(n)-F'_4(n)(x-n)}{(x-n)^2},
```
 so we want to show $\sum_{n \in \mathbb{Z}}H(x,n)\ge 0$. Note since $F_4$ is even 
 ```math
    \frac{H(x,n)+H(x,-n)}{2}=\frac1{2}\left(F_4(x)-F_4(n)\right)\left(\frac1{(x-n)^2}+\frac1{(x+n)^2}\right)-\frac1{2}F_4'(n)\left(\frac1{x-n}-\frac1{x+n}\right)=\left(F_4(x)-F_4(n)\right)\left(\frac{x^2+n^2}{\left(x-n\right)^2(x+n)^2}\right)-\frac{nF_4'(n)}{(x-n)(x+n)}.
```
 We choose a cutoff  $p\ge \eta(x)+1$, and consider $n\ge p+1$, in which case we have 
 ```math
    \left|\frac{H(x,n)+H(x,-n)}{2}\right|\le \frac{\left|F_4(x)-F_4(n)\right|}{|x-n|^2}-\frac{nF_4'(n)}{\left|n^2-x^2\right|}\le \frac{\left|F_4(x)-F_4(n)\right|-nF_4'(n)}{\left(n-x\right)^2}.
```
We have the inequality, $-nF'(n)\le 4 F_4(n)$, so this implies 
```math
    \left|\frac{H(x,n)+H(x,-n)}{2}\right|\le \frac{\left|F_4(x)-F_4(n)\right|+4 F_4(n)}{(n-x)^2}.
```
 $F_4$ is decreasing, so this implies for $n\ge p\ge \eta(x)+2$, 
 ```math
    \sum_{|n|\ge p}\left|H(x,n)\right|\le \sum_{n=p+1}^\infty \frac{2\left(F_4(x)+5F_4(p)\right)}{(n-x)^2}\le  \int_p^\infty \frac{2\left(F_4(x)+5F_4(p)\right)}{(n-x)^2}dn= \frac{2\left(F_4(x)+5F_4(p)\right)}{p-x}.
```
Of course there is the problem of the $\eta(x)$ term being unstable when $n$ is close to an integer. When that happens, we use the mean-value theorem which implies 
```math
    \frac{F_4(x)-F_4(\eta(x))-F_4'(\eta(x))\left(x-\eta(x)\right)}{\left(x-\eta(x)\right)^2}\ge \frac1{2}\min_{\substack{s \\ s \in [x,\eta(x)]}}F_4''(s),
```
with the second derivative satisfying
```math
    F_4''(s)=4F_4(s)\left(\frac{1-F_4(s)}{s^2}\right)\left(4\left(1-2F_4(s)\right)+1\right).
```
When $x$ is close to an integer, we replace $H_4\left(x,\eta(x)\right)$ with the minimum of the right-hand side above over $s \in [x,\eta(x)]$.
This takes care of the $\alpha=4$ case. We now move on to $\alpha\ge 6$.

<br>

Finishing the proof of $\psi_\alpha(x)\le F_\alpha(x)$ for $\alpha\ge 6$.
First we need to verify the following inequality:
Let $\alpha=6,8,10$. Then 
```math
        4\left(F_\alpha\left(\frac1{2}\right)-1\right)+\sum_{n \in \mathbb{Z}\setminus\{0\}}\frac{F_\alpha\left(\frac1{2}\right)-F_\alpha(n)}{n^2}\ge \sum_{n \in \mathbb{Z}\setminus\{0\}}\frac{nF_\alpha'(n)}{\frac1{4}-n^2}.
```
We start by symmetrising the sum, so we want to show 
```math
    2\left(F_\alpha\left(\frac1{2}\right)-1\right)+\sum_{n=1}^\infty \frac{F_\alpha\left(\frac1{2}\right)-F_\alpha(n)}{n^2}-\sum_{n=1}^\infty \frac{nF_\alpha'(n)}{\frac1{4}-n^2}\ge 0,
```
 and evaluating the sum with $F_\alpha\left(\frac1{2}\right)$, we want to show 
 ```math
     2\left(F_\alpha\left(\frac1{2}\right)-1\right)+\frac{\pi^2F_\alpha\left(\frac1{2}\right)}{6}-\sum_{n=1}^\infty \left(\frac{F_\alpha(n)}{n^2}+\frac{-nF_\alpha'(n)}{n^2-\frac1{4}}\right) \ge 0.
```
 We know that $-nF_\alpha'(n)=\alpha F_\alpha(n)\left(1-F_\alpha(n)\right)$, so we want to show 
 ```math
     2\left(F_\alpha\left(\frac1{2}\right)-1\right)+\frac{\pi^2F_\alpha\left(\frac1{2}\right)}{6}-\sum_{n=1}^\infty F_\alpha\left(n\right)\left(\frac1{n^2}+\frac{\alpha\left(1-F_\alpha(n)\right)}{n^2-\frac1{4}}\right)\ge 0.
```
Now for bounding the tail. Pick a cutoff integer $p$, and using the fact that $s_\alpha>1$,  
```math
    \sum_{n=p+1}^\infty F_\alpha\left(n\right)\left(\frac1{n^2}+\frac{\alpha\left(1-F_\alpha(n)\right)}{n^2-\frac1{4}}\right)\le \frac{4\left(\alpha+1\right)}{3}\sum_{n=p+1}^\infty F_\alpha(n)\le \frac{4(\alpha+1)}{3}\sum_{n=p+1}^\infty \frac1{n^\alpha},
```
 so this implies 
 ```math
   \sum_{n=p+1}^\infty F_\alpha\left(n\right)\left(\frac1{n^2}+\frac{\alpha\left(1-F_\alpha(n)\right)}{n^2-\frac1{4}}\right) \le \frac{4(\alpha+1)}{3 (p+1)^\alpha}\left(\frac{\alpha+p}{\alpha-1}\right).
```

<br>

The next thing to do is to prove the following:
 Let $6\le \alpha\le 998$, $\alpha$ an even integer, and let 
 ```math
       \mathcal{B}(\alpha,t)=\sum_{\substack{n \in \mathbb{Z} \\ |n|\ge 2}}\left(\frac{F_\alpha(n)}{(1+t-n)^2}+\frac{F_\alpha'(n)}{1+t-n}\right).
```
   Then 
```math
        \frac{F_\alpha(1+t)-F_\alpha(1)-tF_\alpha'(1)}{t^2}+\sum_{n \in \mathbb{Z}\setminus\{0\}}\frac{F_\alpha(1+t)}{\left(n-t\right)^2}\ge \frac1{(1+t)^2}+\frac{F_\alpha(1)}{(2+t)^2}-\frac{F_\alpha'(1)}{(2+t)}+\mathcal{B}(\alpha,t)
```
 for all $t \in \left[-\frac1{2},\frac1{2}\right]$.

 We start by bounding the tail of $\mathcal{B}(\alpha,t)$. Choosing $p\ge 3$ means $n-(1+|t|)\ge 1$ for $n\ge p$. We know $\left|F_\alpha'(n)\right|\le \alpha F_\alpha(n)$, so 
 ```math
    \sum_{n=p+1}^\infty \left|\frac{F_\alpha(n)}{(1+t-n)^2}+\frac{F_\alpha'(n)}{1+t-n}\right|\le \sum_{n=p+1}^\infty \left(F_\alpha(n)+\left|F_\alpha'(n)\right|\right)\le (\alpha+1)\sum_{n=p+1}^\infty F_\alpha(n)\le \sum_{n=p+1}^\infty \frac{\alpha+1}{n^\alpha}.
``` 
This implies 
```math
    \sum_{n=p+1}^\infty \left|\frac{F_\alpha(n)}{(1+t-n)^2}+\frac{F_\alpha'(n)}{1+t-n}\right|\le \frac{\alpha+1}{(p+1)^\alpha}\cdot \left(\frac{\alpha+p}{\alpha-1}\right).
```
Similarly, for $t \in \left[-\frac1{2},\frac1{2}\right]$, 
```math
    \sum_{|n|\ge p+1}\frac1{(n-t)^2}\le 2\sum_{n=p+1}^\infty \frac1{(n-|t|)^2}\le 2\int_p^\infty \frac1{(n-|t|)^2}dn=\frac2{p-|t|}\le \frac4{2p-1}
```
 When $t$ is small, we use the mean value theorem, so 
 ```math
    \frac{F_\alpha(1+t)-F_\alpha(1)-tF_\alpha'(1)}{t^2}\ge \frac1{2}\min_{|s|\le t}F_\alpha''(1+s), 
```
and evaluating the second derivative, we get
```math
 F_\alpha''(x)=\alpha F_\alpha(x)\left(\frac{1-F_\alpha(x)}{x^2}\right)\left(\alpha\left(1-2F_\alpha(x)\right)+1\right).
```
This implies 
```math
     \frac{F_\alpha(1+t)-F_\alpha(1)-tF_\alpha'(1)}{t^2}\ge \frac{\alpha F_\alpha(1+|t|)\left(1-F_\alpha(1-|t|)\right)}{2(1+|t|)^2}\left(\alpha\left(1-2F_\alpha\left(1-|t|\right)\right)+1\right),
```
so from small $t$, we instead show 
```math
     \frac{\alpha F_\alpha(1+|t|)\left(1-F_\alpha(1-|t|)\right)}{2(1+|t|)^2}\left(\alpha\left(1-2F_\alpha\left(1-|t|\right)\right)+1\right)+\sum_{n \in \mathbb{Z}\setminus\{0\}}\frac{F_\alpha(1+t)}{(n-t)^2}\ge  
    \frac1{(1+t)^2}+\frac{F_\alpha(1)}{(2+t)^2}-\frac{F_\alpha'(1)}{2+t}+\mathcal{B}(\alpha,t).
```

<br>
<br>
<br>
<br>

The penultimate thing to do is to show the following. Let $\eta(x)$ be $x$ rounded up to the nearest integer. Then when $6\le \alpha\le 14$, and $1.5\le x\le 10$,

```math
 \sum_{n \in \mathbb{Z}\setminus\{\eta(x)\}}\left(\frac{F_\alpha(n)}{\left(x-n\right)^2}+\frac{F_\alpha'(n)}{x-n}\right) \le 0.
```
Once again, we need bound the tail. Choosing a cutoff integer $p\ge 12$, 
```math
\sum_{\substack{n \in \mathbb{Z} \\ |n|\ge p}}\left(\frac{F_\alpha(n)}{\left(x-n\right)^2}+\frac{F_\alpha'(n)}{x-n}\right)=
\sum_{\substack{n \in \mathbb{Z} \\ |n|\ge p}}\left(\frac{F_\alpha(n)}{\left(x-n\right)^2}+\frac{nF_\alpha'(n)}{x^2-n^2}\right),
```
and if $n\neq \eta(x)$ with $x\ge 1.5$, then $(n-x)^2,\left|x^2-n^2\right|\ge \frac1{4}$, so this implies 
```math
\left|\sum_{\substack{n \in \mathbb{Z} \\ |n|\ge p}}\left(\frac{F_\alpha(n)}{\left(x-n\right)^2}+\frac{F_\alpha'(n)}{x-n}\right)\right|\le 2\sum_{n=p}^\infty 4\left(F_\alpha(n)+\left|nF_\alpha(n)\right|\right).
```
Using the derivative formula 
```math
F_\alpha'(x)=-\frac{\alpha F_\alpha(x)\left(1-F_\alpha(x)\right)}{x},
```
we get 
```math 
\left|\sum_{\substack{n \in \mathbb{Z} \\ |n|\ge p}}\left(\frac{F_\alpha(n)}{\left(x-n\right)^2}+\frac{F_\alpha'(n)}{x-n}\right)\right|\le 8\sum_{n=p}^\infty(\alpha+1)F_\alpha(n)\le 8(\alpha+1)\sum_{n=p}^\infty \frac1{n^\alpha}\le 8(\alpha+1)\cdot \frac1{p^\alpha}\left(\frac{\alpha+p-1}{\alpha-1}\right).
```

<br>
<br>

Finally we have to check that 
```math
 \sum_{\substack{n \in \mathbb{Z}}}\left(3n^2F_\alpha(n)+n^3F_\alpha'(n)\right)+\frac{16}{100s_\alpha^\alpha}+\frac{10F_\alpha(1)+2F_\alpha'(1)}{99}+\sum_{n=2}^\infty\frac{\left|10n^4F_\alpha(n)+2n^5F_\alpha'(n)\right|}{5}+\frac{8\alpha+2}{10^{\alpha-2}}\le 0.
```
From the derivative formula
```math
F_\alpha'(x)=-\frac{\alpha F_\alpha(x)\left(1-F_\alpha(x)\right)}{x},
```
we see that $3n^2F_\alpha(n)+n^3F_\alpha'(n)\le 0$ if $|n|\ge 2$. What is left is to bound the tail of the other series. We have 
```math
\sum_{n=p}^\infty\frac{\left|10n^4F_\alpha(n)+2n^5F_\alpha'(n)\right|}{5}\le \sum_{n=p}^\infty \frac{(2\alpha+10)n^4F_\alpha(n)}{5}\le \frac{2\alpha+10}{5}\sum_{n=p}^\infty \frac1{n^{\alpha-4}}\le \frac{2\alpha+10}{5}\cdot \frac1{p^{\alpha-4}}\left(\frac{\alpha-5+p}{\alpha-5}\right).
```







