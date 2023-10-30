# Interval Arithmetic Code

A recap and the code for Interval Arithmetic.
We define a function $\widetilde{R}$ as follows: Start by defining $\widehat{\phi}:\mathbb{R} \to \mathbb{R}$ by $\widehat{\phi}(x)=\pi e^{-\sqrt{2}\pi |x|}\cos\left(\sqrt{2}\pi|x|-\frac{\pi}{4}\right)$. We then define $T_n,\mathcal{T}_n:\mathbb{R} \to \mathbb{R}$ by
```math
T_n(x)=\frac{n^2+2x^2}{\left(n+\sqrt{2}x\right)^2}\cdot \frac{2x^2\left(\widehat{\phi}(x)-\widehat{\phi}\left(\frac{n}{\sqrt{2}}\right)\right)}{\left(x-\frac{n}{\sqrt{2}}\right)^2}, \ \mathcal{T}_n(x)=\frac{2x^2\left(n^2+2x^2\right)}{\left(n+\sqrt{2}x\right)^2 }\sum_{k=1}^{9}\frac{\widehat{\phi}^{(k)}\left(\frac{n}{\sqrt{2}}\right)}{k!}\left(x-\frac{n}{\sqrt{2}}\right)^{k-2}
```
Here $\widehat{\phi}^{(k)}$ is the $k$-th derivative of $\widehat{\phi}$. Explicitly,
```math
\widehat{\phi}^{(k)}\left(\frac{n}{\sqrt{2}}\right)=\pi\left(-2\pi\right)^{k}\cos\left(\frac{\pi\left(k+1\right)}{4}\right)\left(-e^{-\pi}\right)^{n}\ .
```
 You can think of $\mathcal{T_n}$ as giving a taylor expansion for $T_n$ when $x$ is close to $\frac{n}{\sqrt{2}}$ for some $n$. Then define $\mathcal{I}_n$ by, for some chosen $\epsilon<10^{-2}$, 
```math
\mathcal{I}_n(x)=\begin{cases}
			\mathcal{T}_n(x) & \text{if $\left|x-\frac{n}{\sqrt{2}}\right|<\epsilon$}\\
            T_n(x) & \text{otherwise}
		 \end{cases}
```
We then set 
```math
\widetilde{R}(x)=\widehat{\phi}\left(x\right)-\widehat{\phi}(0)+\sum_{n=1}^{5\cdot 10^4} \mathcal{I}_n(x)
```
Witht his definition, we aim to show $\widetilde{R}\ge  -\frac{\pi}{\sqrt{2}}\tanh\left(\frac{\pi}{2}\right)+\frac{681}{2^{20}}+10^{-10}$ on $[0,14]$. To do this, we introduce an auxialiary function $\widetilde{R_2}$ given by
```math
\widetilde{R}_2(x):=\widehat{\phi}(x)-\widehat{\phi}(0)+T_1\left(\frac1{\sqrt{2}}\right)+\sum_{n=2}^{5\cdot 10^4}\mathcal{I}_n(x)=\widehat{\phi}(x)-\widehat{\phi}(0)+\frac{\pi^3 e^{-\pi}}{\sqrt{2}}+\sum_{n=2}^{5\cdot 10^4}\mathcal{I}_n(x)\ .
```
We've seen(inlcude link) that for $\left|x-\frac1{\sqrt{2}}\right|\le 0.015$, $\left|\mathcal{R}(x)-\widetilde{R}_2(x)\right|\le \frac{681}{2^{20}}+0.2+10^{-10}$. So we the way we show the desired inequality is to first check that 
```math
\widetilde{R}_2(x)\ge -\frac{\pi}{\sqrt{2}}\tanh\left(\frac{\pi}{2}\right)+\frac{681}{2^{20}}+0.2+10^{-10} \text{ for } x \in \left[\frac1{\sqrt{2}}-0.015,\frac1{\sqrt{2}}+0.015\right]
```
then show 
```math
\widetilde{R}(x)\ge  -\frac{\pi}{\sqrt{2}}\tanh\left(\frac{\pi}{2}\right)+\frac{681}{2^{20}}+10^{-10} \text{ on } \left[0,14\right]\setminus \left[\frac1{\sqrt{2}}-0.015,\frac1{\sqrt{2}}+0.015\right]\ .
```
To do this we use the IntervalArithemtic package. The way it works basically, is that once you've defined your function $f$, you pass in you interval $I$ and returns an interval $J$ such that $f\left(I\right)\subset J$. Here attention to detail is important, especially when working with numbers that aren't representable exactly with floating point arithmetic, such as $\pi$, $e$, and $\frac{4}{3}$. So when implenting the function $f$, it's important to make sure you represent these constants as degenrate intervals. So for example, you would use $[\pi,\pi]$ (instead of just $\pi$), which Julia would represent as a small interval containing the true value of $\pi$, rather than its floating point representation. With this implementation, the general outline of the code is as follows. We have our function $f$. Let $f^*(I)$ denote the interval returned by Julia, so we wish to show  and wish to show $f\left(I\right)> [l,l]$. For two intervals $[a,b]$, $[c,d]$, $[a,b]>[c,d]$ means $a>d$. The way the code works is like so. Suppose $I=[a,b]$. Intialize $d$ to $b-a$, and repeatedly check if $f\left([a,a+d]\right)>[l,l]$, and if not, reassign $l \to \frac{l}{2}$. Once we find some $l$ for which this works, we reassing $a \to l$, and repeat until we've covered the whole interval $[a,b]$. This is why I defined $\mathcal{R}_2$, because I was running into an infinite loop.
```

