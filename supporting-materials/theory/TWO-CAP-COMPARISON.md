# A two-cap collar for the Stringer--Gaffke comparison

## Operational statement

Let `t_1<=...<=t_n` be the observed taints, let `s=SB_B` be the
binomial-factor Stringer value, and let `G_alpha` be the valid Gaffke
bounded-mean upper limit.  The one-cap theorem already proves `s>=G_alpha`
whenever `s>=t_n` at confidence at least 75%.

This note moves strictly below that boundary.  If

\[
 t_{n-1}\le s\le t_n<1,
 \qquad q:=\frac{1-t_n}{1-s},                         \tag{1}
\]

then the following certified collars imply `s>=G_alpha`:

| nominal confidence | sample sizes | sufficient relative-gap condition |
|---:|---:|---:|
| 90% | `2<=n<=200` | `q>=5/12` |
| 95% | `2<=n<=200` | `q>=1/3` |
| 99% | `2<=n<=200` | `q>=1/4` |

Consequently the pre-specified Stringer--Gaffke safeguard has zero uplift on
this complete two-upper-knot collar.  The result is sample-wise, not a
conditional-coverage theorem for ordinary Stringer.  It does not assert that
the rest of the two-cap region is unsafe; that residual pointwise comparison
remains open beyond the separately proved complete dimensions.

## 1. Exact two-cap formula

Use the general ordered-simplex notation

\[
 x_0\le\cdots\le x_{n-2}\le s\le x_{n-1}<x_n=1,
 \qquad s=\sum_{i=0}^n c_i x_i,
\]

where the positive `c_i` sum to one.  Put

\[
 S=1-s,\qquad q=\frac{1-x_{n-1}}S,\qquad
 z_i=\frac{1-x_i}S\quad(0\le i\le n-2).              \tag{2}
\]

Then `0<q<=1`, the `z_i` are decreasing and at least one, and

\[
 sum_{i=0}^{n-2}c_i(z_i-1)
 =K(q):=c_n+c_{n-1}(1-q).                            \tag{3}
\]

The uniform-simplex spline formula has exactly two active knots.  For
distinct knots it gives

\[
 \Phi_c(x)=\frac1q\left{
 \frac1{\prod_{i=0}^{n-2}z_i}
 -\frac{(1-q)^n}{\prod_{i=0}^{n-2}(z_i-q)}
 \right}.                                           \tag{4}
\]

Coincident-knot and endpoint cases follow by continuity.  Equation (4) is
also an exact closed-form test for any particular sample in the two-cap
region.

## 2. A dimension-free uniform upper bound

Define prefix weights

\[
 C_r=\sum_{i=0}^{r-1}c_i,\qquad 1\le r\le n-1,
\]

and

\[
 g_n(q)=\frac{1-(1-q)^n}{q}
       =\sum_{k=0}^{n-1}(1-q)^k.                     \tag{5}
\]

**Two-cap lemma.**  Under (2)--(3),

\[
 \Phi_c(x)\le g_n(q)
 \max_{1\le r\le n-1}
 \left(\frac{C_r}{C_r+c_n+c_{n-1}(1-q)}\right)^r.   \tag{6}
\]

**Proof.**  Since every `z_i>=1`, the second product in (4) satisfies

\[
 \frac{(1-q)^n}{\prod_i(z_i-q)}
 \ge \frac{(1-q)^n}{\prod_i z_i}.
\]

Thus (4) is at most `g_n(q)/prod_i z_i`.  Put `u_i=z_i-1` and append
`u_(n-1)=0`.  The decreasing nonnegative vector `u` has the prefix-step
decomposition

\[
 u=\sum_{r=1}^{n-1}\lambda_r\frac{K(q)}{C_r}v_r,
 \qquad \lambda_r\ge0,\quad\sum_r\lambda_r=1,
\]

where `v_r` has `r` initial ones.  Concavity of
`sum_i log(1+u_i)` therefore gives

\[
 \prod_i z_i\ge
 \min_r\left(1+\frac{K(q)}{C_r}\right)^r.
\]

Taking reciprocals proves (6).  \(\square\)

The lemma is analytic, applies to arbitrary positive barycentric weights,
and has no sample-size cutoff.

## 3. Why one endpoint certificate covers an entire collar

Put `x=1-q`, `A=C_r+c_n`, and `b=c_(n-1)`.  The `r`-th expression on the
right of (6), apart from the constant `C_r^r`, is

\[
 B_r(x)=\frac{1+x+\cdots+x^{n-1}}{(A+bx)^r}.         \tag{7}
\]

Its derivative has the sign of

\[
 P_r(x)=(A+bx)(1+x+\cdots+x^{n-1})'
        -rb(1+x+\cdots+x^{n-1}).                    \tag{8}
\]

For `0<=j<=n-2`, the coefficient of `x^j` in (8) is

\[
 A(j+1)+b(j-r),                                      \tag{9}
\]

and consecutive coefficients increase by `A+b`.  The last coefficient is
`b(n-1-r)>=0` (and vanishes when `r=n-1`).  After trailing zeros are removed,
the coefficient signs therefore have at most one change, from negative to
positive.  Descartes' rule shows that `P_r` has at most one positive root;
if it has one, `B_r` decreases and then increases.  Hence `B_r` has no
interior maximum on any positive interval.  Its maximum on
`0<=x<=x_0` occurs at `x=0` or `x=x_0`.

At `x=0` (the one-cap boundary `q=1`), the analytic all-sample-size theorem
in `ONE-CAP-COMPARISON.md` bounds every term by `alpha` whenever
`alpha<=1/4`.  It remains only to check the other endpoint.  Exact dyadic
Clopper--Pearson brackets and rational arithmetic certify that endpoint for
every `2<=n<=200` and every `r<n` at

\[
 (\alpha,q_0)=(0.10,5/12),\ (0.05,1/3),\ (0.01,1/4). \tag{10}
\]

Equations (6)--(10) prove the operational statement.

## 4. Verification boundary

`computations/python/two_cap_certificate.py` regenerates
`computations/certificates/two-cap-certificate.json`.  It encloses every
Clopper--Pearson factor on a 64-bit dyadic grid, checks the defining binomial
CDF signs with integer arithmetic, and evaluates all 59,700 collar-endpoint
inequalities as exact rational numbers.  The run is serial; on the reference
machine it used about 80 MB peak resident memory.

The exact certificate proves only `n<=200` at the three stated levels.  The
dimension-free lemma and the endpoint-maximization argument are ordinary
mathematics, but an all-`n` proof of the remaining factor endpoint inequality
has not been supplied.  Nor does this result settle the part of the two-cap
region with smaller `q`, cap regions with three or more upper knots, or
general Stringer coverage.
