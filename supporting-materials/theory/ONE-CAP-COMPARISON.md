# A certified one-cap no-uplift region through `n = 200`

## Operational statement

Let `t_1 <= ... <= t_n` be the observed taints and let
`SB_B` be the binomial-factor Stringer value at tail probability `alpha`.
For

```text
alpha in {0.10, 0.05, 0.01},       1 <= n <= 200,
```

this note proves the following pointwise implication:

```text
SB_B >= t_n    ==>    SB_B >= G_alpha(t_1,...,t_n),       (1)
```

where `G_alpha` is the valid one-sided Gaffke bounded-mean upper limit.
Consequently, on every sample satisfying the condition in (1), the
pre-specified Stringer--Gaffke safeguard has zero uplift.  This holds for
both factor conventions: Poisson Stringer is pointwise no smaller than
binomial Stringer at the three listed levels.

The implication is directly checkable from an observed sample.  It is not a
conditional-coverage theorem for ordinary Stringer.  The coverage guarantee
continues to belong to the pre-specified maximum rule; (1) proves that the
maximum returns the familiar Stringer value on this part of the sample
space.

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
 \left(\frac{C_r}{C_r+c_n}\right)^r .                  \tag{2}
\]

This lemma is analytic and has no finite-`n` restriction.

### Proof

Put `y_i=1-x_i` for `i<n` and `S=1-s`.  The one-upper-knot cap formula is

\[
 V:=\Pr\!\left\{\sum_{i=0}^n x_iD_i>s\right\}
   =\frac{S^n}{\prod_{i=0}^{n-1}y_i}.                  \tag{3}
\]

If `S=0`, then `V=0` and there is nothing to prove.  If one of the `y_i`
vanishes, the region condition forces `S=0`; equivalently, the remaining
case follows by continuity.  Hence assume `S>0` and every `y_i>0`.

Because `s>=x_(n-1)`, we have `S<=y_(n-1)`.  Set

\[
 z_i=\frac{y_i}{S},\qquad u_i=z_i-1,qquad i=0,\ldots,n-1.
\]

Then

\[
 u_0\ge u_1\ge\cdots\ge u_{n-1}\ge0,
 \qquad \sum_{i=0}^{n-1}c_i u_i=c_n,                  \tag{4}
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
 u=\sum_{r=1}^n\lambda_r\frac{c_n}{C_r}v_r.           \tag{5}
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
\end{aligned}                                           \tag{6}
\]

The factors are strictly increasing, so these weights are positive and sum
to one.  Their prefix sums telescope:

\[
 C_r=1-p_{n-r}.                                         \tag{7}
\]

It follows from (2) that the one-cap comparison is proved once

\[
 \left(\frac{1-p_{n-r}}{1-p_{n-r}+p_0}\right)^r
 \le\alpha,qquad r=1,\ldots,n.                        \tag{8}
\]

The terminal term is an identity.  Indeed, `C_n=1-p_0`, while the defining
Clopper--Pearson equation for zero observed errors is

\[
 (1-p_0)^n=\alpha.                                      \tag{9}
\]

Only the `r<n` inequalities in (8) require certification.

## 3. Exact finite-range certificate

The source
[`one_cap_certificate.py`](../computations/python/one_cap_certificate.py)
checks every nonterminal inequality in (8) for `1<=n<=200` at the three
listed levels.  For each factor it:

1. locates the root between dyadic rationals on a `2^-64` grid;
2. evaluates the binomial-CDF sign at both endpoints with integer
   arithmetic;
3. uses the upper enclosure

   \[
   C_r\le1-p_{n-r}^{\rm lower},\qquad
   p_0\ge p_0^{\rm lower};
   \]

4. raises the resulting rational ratio to the integer power `r`; and
5. compares that exact rational strictly with the decimal design value of
   `alpha`, represented as a `Fraction`.

There are `19,900` strict nonterminal comparisons at each level and `200`
analytic terminal equalities.  The largest certified nonterminal upper
bounds are:

| `alpha` | nominal confidence | `n` | `r` | upper bound | margin below `alpha` |
|---:|---:|---:|---:|---:|---:|
| 0.01 | 99% | 200 | 199 | 0.009766418939443551 | 0.000233581060556449 |
| 0.05 | 95% | 200 | 199 | 0.049444579494823319 | 0.000555420505176681 |
| 0.10 | 90% | 200 | 199 | 0.099327536847359528 | 0.000672463152640472 |

The displayed decimals are summaries.  The committed
[`one-cap-certificate.json`](../computations/certificates/one-cap-certificate.json)
contains the exact numerator and denominator of every worst-case upper
bound and margin, plus a SHA-256 digest over all `59,700` nonterminal
rational bounds.  Regenerate it byte-for-byte with:

```bash
make one-cap-certificate-check
```

No floating-point value decides a root sign or an inequality sign.

## 4. Scope and next mathematical target

This result materially enlarges the region in which the all-sample-size
safeguard is known to leave Stringer unchanged, but it does not prove any of
the following:

- unrestricted pointwise Stringer--Gaffke domination for `6<=n<=200`;
- ordinary Stringer coverage for those sample sizes;
- the inequalities in (8) for `n>200`; or
- validity under sampling designs outside the paper's independent
  `[0,1]`-taint model.

The upper limit `200` is the endpoint of the exact certificate, not a
failure boundary.  The analytic lemma (2) holds in every dimension.  Two
natural extensions are therefore sharply separated:

1. prove (8) analytically at conventional confidence levels for every `n`;
2. control the additional alternating simplex-cap terms when the Stringer
   threshold lies below the largest observed taint.

The second problem is the remaining obstacle to a dimension-free pointwise
comparison.  The present lemma removes the entire one-upper-knot region from
that obstacle.
