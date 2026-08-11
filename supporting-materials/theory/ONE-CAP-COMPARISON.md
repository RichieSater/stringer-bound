# An all-sample-size one-cap no-uplift theorem

## Operational statement

Let `t_1 <= ... <= t_n` be the observed taints and let `SB_B` be the
binomial-factor Stringer value at tail probability `alpha`.  This note proves

```text
0 < alpha <= 0.10,  n >= 1,  and  SB_B >= t_n
    ==> SB_B >= G_alpha(t_1,...,t_n),                 (1)
```

where `G_alpha` is the valid one-sided Gaffke bounded-mean upper limit.
Consequently, whenever nominal confidence is at least 90%, the
pre-specified Stringer--Gaffke safeguard has zero uplift on every sample
satisfying the displayed, directly checkable condition.  This is true for
both factor conventions: at these levels Poisson Stringer is pointwise no
smaller than binomial Stringer.

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
when `alpha<=0.10`.

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

We prove (12) for every `n` and `r` when `0<alpha<=0.10`.

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
It therefore suffices to take `x=log 10`.

For `r>=4`, the function `r(1-e^(-x/r))` increases with `r`, and

\[
 4(1-10^{-1/4})>\frac74>e-1.
\]

Here `10^(1/4)>16/9` and `e<11/4` give the two strict rational
comparisons.  The logarithm in (22) is consequently greater than one,
whereas `u/(e^u-1)<1`.

For `r=2`, use `sqrt(10)>3` and `log 10>2` to get

\[
 \log(3-2/\sqrt{10})>\log(7/3)>3/4,
 \qquad \frac{u}{e^u-1}<\frac{1}{1+u/2}<3/4.
\]

For `r=3`, similarly,

\[
 \log(4-3/\sqrt[3]{10})>\log(5/2)>4/5,
 \qquad \frac{u}{e^u-1}<\frac{1}{1+u/2}<4/5.
\]

The logarithmic comparisons follow, for example, from
`(11/4)^2<10`, `(11/4)^3<(7/3)^4`, and
`(11/4)^4<(5/2)^5`.  Thus `H_r(x)>0` for every
`r>=2` and `x>=log 10`.  Equations (20)--(22) prove the boundary case,
and (15)--(19) prove (12) in full.  This establishes (8) for every sample
size whenever `alpha<=0.10`.

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
least 90%.  It does **not** prove any of the following:

- unrestricted pointwise Stringer--Gaffke domination for `n>=6`;
- ordinary Stringer coverage at arbitrary sample size;
- the one-cap factor inequalities when nominal confidence is below 90%; or
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
