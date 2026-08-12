# An all-sample-size one-cap no-uplift theorem

## Operational statement

Let `t_1 <= ... <= t_n` be the observed taints and let `SB_B` be the
binomial-factor Stringer value at tail probability `alpha`.  This note proves

```text
0 < alpha <= 1/4,  n >= 1,  and  SB_B >= t_n
    ==> SB_B >= G_alpha(t_1,...,t_n),                 (1)
```

where `G_alpha` is the valid one-sided Gaffke bounded-mean upper limit.
Consequently, whenever nominal confidence is at least
`75%`, the pre-specified Stringer--Gaffke safeguard has
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
when `alpha<=1/4`.

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

We prove (12) for every `n` and `r` when `0<alpha<=1/4`.

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
 \Pr(X\ge r)=\Pr(X>\mathbb EX)\ge\frac14\ge\alpha.  \tag{15}
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

Write

\[
 M_r(x)=A^r(r+1-rbA),                                  \tag{21}
\]

the multiplier of `alpha` in (20).  This multiplier is strictly increasing
in `x`.  Indeed, with `psi(z)=z/(e^z-1)`,

\[
 \frac{d}{dx}\log A=\frac{\psi(u)-\psi(v)}x>0.         \tag{22}
\]

Moreover, `psi` is decreasing and

\[
 \psi'(z)+1=\frac{e^z(e^z-1-z)}{(e^z-1)^2}>0.
\]

Consequently,

\[
 \frac{d}{dx}\log(bA)
 =-\frac1r+\frac{\psi(u)-\psi(v)}x
 <-\frac1r+\frac{v-u}{x}=-\frac1{r+1}<0.              \tag{23}
\]

Thus both factors defining `M_r(x)` increase with `x`.

It remains to prove `M_r(1)>1`.  Put `t=1/r` and `d=t^2/(1+t)`.  At
`x=1`,

\[
 1-A=\frac{e^d-1}{e^t-1}=: \delta.                    \tag{24}
\]

For `0<z<2`, the elementary bounds

\[
 z+\frac{z^2}{2}<e^z-1<\frac{z}{1-z/2},
 \qquad
 1-z<e^{-z}<1-z+\frac{z^2}{2}                         \tag{25}
\]

follow from the exponential series and
`2 arctanh(z/2)>z`.  They give

\[
 \delta<\frac{d/(1-d/2)}{t+t^2/2},
 \qquad
 \delta>\frac{d(1-t/2)}t.                             \tag{26}
\]

Since `(e^z-1)/z` is increasing, also
`delta<d/t=1/(r+1)`, and hence `A>r/(r+1)`.  The second bound in (26),
together with `b=e^(-t)>1-t` and `1-b>t-t^2/2`, gives

\[
\begin{aligned}
 r+1-rbA
 &=1+r(1-b)+rb\delta\\
 &>1+\left(1-\frac t2\right)
   +\frac{(1-t)(1-t/2)}{1+t}
 =\frac3{1+t}=\frac{3r}{r+1}.                         \tag{27}
\end{aligned}
\]

Therefore, for `r>=5`,

\[
 M_r(1)>3\left(\frac r{r+1}\right)^{r+1}
 \ge3\left(\frac56\right)^6>1,                       \tag{28}
\]

because the sequence `(r/(r+1))^(r+1)` is increasing.  For `r=2,3,4`, the first bound
in (26) gives these exact rational lower margins:

| `r` | lower bound for `A` | lower bound for `(3r/(r+1)) A^r - 1` |
|---:|---:|---:|
| 2 | `39/55` | `17/3025` |
| 3 | `125/161` | `885001/16693124` |
| 4 | `287/351` | `1841131309/25297477335` |

Thus `M_r(1)>1` for every `r>=2`.  Because `x>=log(4)>1` throughout the
stated range, the boundary probability in (17) is at least `alpha`.
Together with the separate `r=1` argument in (18)--(19), this proves (12)
and establishes (8) for every sample size whenever `alpha<=1/4`.

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
least `75%`.  It does **not** prove any of the following:

- unrestricted pointwise Stringer--Gaffke domination outside the separately
  certified sample-size/level pairs (`n=3,4,5` at 90%, 95%, and 99%, and
  `n=6` at 95%);
- ordinary Stringer coverage at arbitrary sample size;
- the one-cap factor inequalities when nominal confidence is below
  `75%`; or
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
