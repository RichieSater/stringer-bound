# An all-sample-size one-cap no-uplift theorem

## Operational statement

Let `t_1 <= ... <= t_n` be the observed taints and let `SB_B` be the
binomial-factor Stringer value at tail probability `alpha`.  This note proves

```text
0 < alpha <= exp(1-e),  n >= 1,  and  SB_B >= t_n
    ==> SB_B >= G_alpha(t_1,...,t_n),                 (1)
```

where `G_alpha` is the valid one-sided Gaffke bounded-mean upper limit.
Consequently, whenever nominal confidence is at least
`1-exp(1-e) = 82.063...%`, the pre-specified Stringer--Gaffke safeguard has
zero uplift on every sample satisfying the displayed, directly checkable
condition.  This includes the conventional 90%, 95%, and 99% levels.  The
result is true for both factor conventions: throughout this range Poisson
Stringer is pointwise no smaller than binomial Stringer.

The result has no upper bound on sample size.  It is a sample-wise comparison,
not a conditional-coverage theorem for ordinary Stringer.  The coverage
guarantee continues to belong to the pre-specified maximum rule; (1) proves
that the maximum returns the familiar Stringer value on one complete region
of the sample space.

## 1. A dimension-free one-cap lemma

Let

\[
  x_0\le x_1\le\cdots\le x_{n-1}\le x_n=1,
  \qquad c_i>0,\qquad \sum_{i=0}^n c_i=1,
\]

and put

\[
  s=\sum_{i=0}^n c_i x_i,
  \qquad C_r=\sum_{i=0}^{r-1}c_i\quad(1\le r\le n).
\]

Let `(D_0,...,D_n) ~ Dirichlet(1,...,1)`.  If `s >= x_(n-1)`, then

\[
 \Pr\!\left\{\sum_{i=0}^n x_iD_i>s\right\}
 \le
 \max_{1\le r\le n}
 \left(\frac{C_r}{C_r+c_n}\right)^r .                 \tag{2}
\]

This lemma is analytic and has no finite-`n` restriction.

### Proof

Put `y_i=1-x_i` for `i<n` and `S=1-s`.  The one-upper-knot cap formula is

\[
 V:=\Pr\!\left\{\sum_{i=0}^n x_iD_i>s\right\}
   =\frac{S^n}{\prod_{i=0}^{n-1}y_i}.                 \tag{3}
\]

If `S=0`, then `V=0` and there is nothing to prove.  If one of the `y_i`
vanishes, the region condition forces `S=0`; equivalently, the remaining
case follows by continuity.  Hence assume `S>0` and every `y_i>0`.

Because `s>=x_(n-1)`, we have `S<=y_(n-1)`.  Set

\[
 z_i=\frac{y_i}{S},\qquad u_i=z_i-1,
 \qquad i=0,\ldots,n-1.
\]

Then

\[
 u_0\ge u_1\ge\cdots\ge u_{n-1}\ge0,
 \qquad \sum_{i=0}^{n-1}c_i u_i=c_n,                 \tag{4}
\]

and (3) becomes `V=1/prod_i(1+u_i)`.

Let `v_r` be the length-`n` vector whose first `r` entries are one and whose
remaining entries are zero.  With `u_n=0`, every vector in (4) has the cone
decomposition

\[
 u=\sum_{r=1}^n(u_{r-1}-u_r)v_r.
\]

The budget equation in (4) says

\[
 \lambda_r:=\frac{C_r(u_{r-1}-u_r)}{c_n}\ge0,
 \qquad \sum_{r=1}^n\lambda_r=1,
\]

and therefore

\[
 u=\sum_{r=1}^n\lambda_r\frac{c_n}{C_r}v_r.          \tag{5}
\]

Thus the feasible budget section is the convex hull of the `n` prefix-step
vectors `(c_n/C_r)v_r`.  The function

\[
 u\longmapsto\sum_{i=0}^{n-1}\log(1+u_i)
\]

is concave.  Its value at any convex combination is at least the smallest
value at a vertex.  By (5),

\[
 \prod_{i=0}^{n-1}(1+u_i)
 \ge \min_{1\le r\le n}
       \left(1+\frac{c_n}{C_r}\right)^r.
\]

Taking reciprocals proves (2).

## 2. Specialization to the Stringer weights

Write `p_j=p_n(j)` for the binomial Clopper--Pearson factors.  When the
sample taints are the ascending knots `x_i=t_(i+1)` for `i<n` and `x_n=1`,
the Stringer threshold is `s=sum_i c_i x_i`, where

\[
\begin{aligned}
c_0&=1-p_{n-1},\\
c_j&=p_{n-j}-p_{n-j-1},\qquad 1\le j\le n-1,\\
c_n&=p_0.
\end{aligned}                                          \tag{6}
\]

The factors are strictly increasing, so these weights are positive and sum
to one.  Their prefix sums telescope:

\[
 C_r=1-p_{n-r}.                                        \tag{7}
\]

It follows from (2) that the one-cap comparison is proved once

\[
 \left(\frac{1-p_{n-r}}{1-p_{n-r}+p_0}\right)^r
 \le\alpha,\qquad r=1,\ldots,n.                       \tag{8}
\]

The terminal term is an identity.  Indeed, `C_n=1-p_0`, while the defining
Clopper--Pearson equation for zero observed errors is

\[
 (1-p_0)^n=\alpha.                                     \tag{9}
\]

The next section proves every nonterminal inequality in (8) analytically
when `alpha<=exp(1-e)`.

## 3. The all-`n` Clopper--Pearson comparison

Fix `1<=r<n`, write

\[
 \alpha=e^{-x},\qquad
 a=e^{-x/n}=\alpha^{1/n},\qquad
 b=e^{-x/r}=\alpha^{1/r},
\]

and put

\[
 q_{n,r}:=\frac{b(1-a)}{1-b}
          =\frac{1-e^{-x/n}}{e^{x/r}-1}.              \tag{10}
\]

By binomial symmetry, `q=1-p_(n-r)` is the unique solution of

\[
 \Pr\{\operatorname{Bin}(n,q)\ge r\}=\alpha.         \tag{11}
\]

The inequality in (8) is equivalent to `q<=q_(n,r)`.  Since the upper
binomial tail is increasing in its success probability, it is enough to
prove

\[
 \Pr\{\operatorname{Bin}(n,q_{n,r})\ge r\}\ge\alpha. \tag{12}
\]

We prove (12) for every `n` and `r` when `0<alpha<=exp(1-e)`.

### 3.1 Reduction to `n=r+1`

For fixed `r`, let

\[
 \mu_n=nq_{n,r}
 =\frac{b}{1-b}\,n(1-e^{-x/n}).                       \tag{13}
\]

The function `z -> z(1-e^(-x/z))` is strictly increasing because its
derivative is `1-(1+x/z)e^(-x/z)>0`.  Hence `mu_n` increases with `n`.
Also

\[
 q_{n,r}=\frac{1-e^{-x/n}}{e^{x/r}-1}<\frac{x/n}{x/r}=\frac rn,          \tag{14}
\]

so `mu_n<r`.

Suppose first that `r>=2` and `mu_n>=r-1`.  If
`X~Bin(n,q_(n,r))`, then `r-1<=E[X]<r`, so
`{X>=r}={X>E[X]}`.  Moreover `nq_(n,r)>=1>ln(4/3)`.  Pinelis's sharp
binomial-mean inequality therefore gives

\[
 \Pr(X\ge r)=\Pr(X>\mathbb EX)\ge\frac14>\alpha.      \tag{15}
\]

Now suppose `mu_n<r-1`, and set `mu_0=mu_(r+1)`.  Then
`mu_0<=mu_n<r-1`.  Monotonicity of the binomial tail in its success
probability gives

\[
 \Pr\{\operatorname{Bin}(n,\mu_n/n)\ge r\}
 \ge
 \Pr\{\operatorname{Bin}(n,\mu_0/n)\ge r\}.          \tag{16}
\]

Anderson and Samuels's fixed-mean monotonicity theorem says that, when
`k>lambda`, the lower tail
`Pr{Bin(m,lambda/m)<=k}` decreases as `m` increases.  Apply it with
`k=r-1` and `lambda=mu_0`.  Taking complements in (16) yields

\[
 \Pr\{\operatorname{Bin}(n,q_{n,r})\ge r\}
 \ge
 \Pr\{\operatorname{Bin}(r+1,q_{r+1,r})\ge r\}.      \tag{17}
\]

Thus, apart from `r=1`, every case either already has probability at least
`1/4` or reduces to the two-term boundary tail at `n=r+1`.

For `r=1`, (12) follows directly.  Here

\[
 q_{n,1}=\frac{\alpha}{1-\alpha}(1-\alpha^{1/n}),
\]

and `mu_n=nq_(n,1)` is increasing, so

\[
 \mu_n\ge\mu_2=\frac{2\alpha}{1+\sqrt\alpha}
 \ge\frac{\alpha}{1-\alpha}
 \ge-\log(1-\alpha).                                  \tag{18}
\]

The middle inequality holds for `alpha<=1/4`.  Therefore

\[
 \Pr\{\operatorname{Bin}(n,q_{n,1})\ge1\}
 =1-(1-q_{n,1})^n
 \ge1-e^{-\mu_n}\ge\alpha.                           \tag{19}
\]

### 3.2 The boundary tail

It remains to bound the right side of (17).  Put

\[
 u=\frac{x}{r+1},\qquad v=\frac{x}{r},\qquad
 A=\frac{1-e^{-u}}{1-e^{-v}},\qquad q_{r+1,r}=bA.
\]

Because `alpha=b^r`, the exact two-term tail is

\[
\begin{aligned}
 \Pr\{\operatorname{Bin}(r+1,bA)\ge r\}
 &= (bA)^r(r+1-rbA)\\
 &= \alpha A^r(r+1-rbA).                              \tag{20}
\end{aligned}
\]

Let `f(t)=log(1-e^(-t))`.  Its derivative `1/(e^t-1)` is decreasing.
Since `0<u<v`,

\[
 r\log A
 =-r\int_u^v\frac{dt}{e^t-1}
 \ge-\frac{u}{e^u-1}.                                 \tag{21}
\]

Also `0<A<1`, so `r+1-rbA>=1+r(1-b)`.  It follows that the logarithm
of the multiplier of `alpha` in (20) is at least

\[
 H_r(x):=log\!\left(1+r(1-e^{-x/r})\right)
          -\frac{x/(r+1)}{e^{x/(r+1)}-1}.              \tag{22}
\]

The first term in (22) increases with `x`; the subtracted term decreases.
We prove that `H_r(x)>0` already when `x>=e-1`.

First suppose `r>=3`, put `x_0=e-1`, and set `t=1/r`.  For positive
`z,u`,

\[
 1-e^{-z}>\frac{z}{1+z/2+z^2/12},
 \qquad
 \frac{u}{e^u-1}<\frac1{1+u/2}.                       \tag{23}
\]

For the first inequality, the derivative of

\[
 \log\!\left(\frac{z^2+6z+12}{z^2-6z+12}\right)-z
\]

is
`-z^4/((z^2-6z+12)(z^2+6z+12))<0`; the second follows from
`e^u>1+u+u^2/2`.  Hence `H_r(x_0)` is strictly larger than

\[
 F(t)=
 \log\!\left(1+\frac{x_0}{1+x_0t/2+x_0^2t^2/12}\right)
 -\frac1{1+x_0t/(2(1+t))}.                            \tag{24}
\]

Define `F(0)=0` by continuity.  Direct differentiation shows that `F'(t)`
has the sign of

\[
\begin{aligned}
 P_x(t)={}&t^4x^4-12t^3x^4-36t^3x^3-48t^3x^2
 -72t^2x^3-180t^2x^2\\
 &-144t^2x-120tx^2-144tx+144
\end{aligned}                                       \tag{25}
\]

at `x=x_0`.  Since `0<x_0<2` and `0<t<=1/3`, the sole positive term in
`P_x'(t)` is dominated by `-36t^2x^4`; every other term is negative.
Thus `P_x` is strictly decreasing and `F'` changes sign at most once, from
positive to negative.  The minimum of `F` on `[0,1/3]` is at an endpoint.

At the nonzero endpoint, put

\[
 z_0=\frac{54(e-1)}{e^2+70e+37}.
\]

Then

\[
 F(1/3)=2\operatorname{arctanh}(z_0)-\frac8{e+7}
 >2\left(z_0+\frac{z_0^3}{3}\right)-\frac8{e+7}.     \tag{26}
\]

The last expression increases with `e` on `[65/24,3]`: `z_0` increases
there, while `8/(e+7)` decreases.  The exponential series gives
`65/24<e<3`, and exact rational simplification at `65/24` makes the last
difference greater than `1/200`.  Thus both endpoints of (24) are
nonnegative, and its nonzero endpoint is positive.  This proves
`H_r(e-1)>0` for `r>=3`.

For `r=2`, monotonicity permits the smaller rational value
`x=41/24<e-1`.  Put `z=x/2` and `u=x/3`.  The alternating and positive
exponential series give

\[
 1-e^{-z}>z-\frac{z^2}{2}+\frac{z^3}{6}-\frac{z^4}{24},
 \qquad
 \frac{u}{e^u-1}<\frac1{1+u/2+u^2/6}.                \tag{27}
\]

Apply
`log y=2 arctanh((y-1)/(y+1))>2(v+v^3/3)` to the rational lower bound for
the first term in `H_2(x)`.  Exact simplification leaves a margin greater
than `1/100`.  Hence `H_r(x)>0` for every `r>=2` and `x>=e-1`.

Finally, `e>8/3` gives both `e-1>3/2` and `e^(3/2)>4`, so
`exp(1-e)<1/4`; the separate `r=1` argument in (18)--(19) therefore
continues to apply.  Equations (15)--(27) prove (12) in full and establish
(8) for every sample size whenever `alpha<=exp(1-e)`.

## 4. Reproducibility and independent finite regression

The analytic proof above is the theorem; it does not depend on a finite
search.  The source
[`one_cap_all_n_check.py`](../computations/python/one_cap_all_n_check.py)
checks its algebraic identities, the exact rational constant comparisons,
and a broad numerical regression grid.  Regenerate that check with:

```bash
make one-cap-all-n-check
```

The earlier exact finite-range program
[`one_cap_certificate.py`](../computations/python/one_cap_certificate.py)
remains as an arithmetically independent regression.  It encloses every
Clopper--Pearson root on a dyadic grid, evaluates endpoint signs with integer
arithmetic, and checks all `59,700` nonterminal vertices through `n=200` at
`alpha` equal to `0.10`, `0.05`, and `0.01`.  Its committed output is
[`one-cap-certificate.json`](../computations/certificates/one-cap-certificate.json).
It is no longer the logical basis for the all-sample-size result.

## 5. Scope and next mathematical target

The theorem removes the entire one-upper-knot region from the unresolved
comparison at every sample size and every nominal confidence level of at
least `1-exp(1-e)=82.063...%`.  It does **not** prove any of the following:

- unrestricted pointwise Stringer--Gaffke domination for `n>=6`;
- ordinary Stringer coverage at arbitrary sample size;
- the one-cap factor inequalities when nominal confidence is below
  `1-exp(1-e)`; or
- validity under sampling designs outside the paper's independent
  `[0,1]`-taint model.

The remaining pointwise problem is to control the additional alternating
simplex-cap terms when the Stringer threshold lies below the largest observed
taint.  That is now the only geometric region left by the one-cap theorem.

## Sources used in the all-`n` reduction

- T. W. Anderson and S. M. Samuels, “Some inequalities among binomial and
  Poisson probabilities,” *Proceedings of the Fifth Berkeley Symposium on
  Mathematical Statistics and Probability*, Vol. 1, 1967, pp. 1--12,
  Theorem 2.1.
- I. Pinelis, “Best lower bound on the probability of a binomial exceeding
  its expectation,” *Statistics & Probability Letters* **179** (2021),
  109224, <https://doi.org/10.1016/j.spl.2021.109224>.
