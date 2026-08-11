# Program for an all-`n` Poisson-Stringer theorem

> **Status.** This note separates proved reductions from two open analytic
> inequalities.  It is a research roadmap, not a coverage claim.  The paper
> continues to state that ordinary Stringer coverage is unresolved for
> arbitrary `n`.

A separate corrected simultaneous-band argument now proves ordinary Poisson
Stringer coverage for every `n<=8` at 90%, every `n<=11` at 95%, and every
`n<=20` at 99% confidence; see
[`POISSON-SIMULTANEOUS-BAND.md`](POISSON-SIMULTANEOUS-BAND.md). The sufficient
band event is already below nominal at the next sample size in each row, so
that finite-range theorem does not replace the dimension-free program here.

An all-sample-size theorem for the Poisson factors used in audit practice
would have substantially greater practice relevance than extending the exact
Bernstein certificates one dimension at a time.  The reductions below turn
that target into two distributional-quantile questions with no audit-specific
notation left in them.

## 1. An exact ordered-weight identity

Fix tail probability `alpha`.  Let `lambda_j` be the upper-tail gamma
quantile

\[
 \Pr\{\Gamma(j+1,1)>\lambda_j\}=\alpha,
 \qquad j=0,1,\ldots .                       \tag{1}
\]

These are exactly the Poisson count limits, so
`p_n^P(j)=lambda_j/n`.

Given taints `t_(1)>=...>=t_(n)`, put

\[
a=(a_1,\ldots,a_{n+1})=(1,t_{(1)},\ldots,t_{(n)}),
\qquad a_{n+2}=0.
\]

Summation by parts gives the exact identity

\[
n\,\operatorname{SB}_{\rm P}
=\sum_{k=1}^{n+1}(a_k-a_{k+1})\lambda_{k-1}.       \tag{2}
\]

Thus Poisson Stringer is the ordered-weight linear interpolation of the
gamma quantiles in (1).

## 2. Reduction to weighted exponentials and a Dirichlet average

Let `E_1,...,E_m` be iid unit exponentials, where `m=n+1`, and define

\[
q_\alpha(a)=Q_{1-\alpha}\!\left(\sum_{i=1}^m a_iE_i\right).       \tag{3}
\]

The ordered vector has the vertex decomposition

\[
a=\sum_{k=1}^m(a_k-a_{k+1})\,\mathbf 1_{{1,\ldots,k\}}.          \tag{4}
\]

Positive homogeneity and (1) give

\[
q_\alpha(\mathbf 1_{{1,\ldots,k\}})=\lambda_{k-1}.              \tag{5}
\]

Consequently, convexity of `q_alpha` on the nonnegative orthant would imply

\[
q_\alpha(a)
\le\sum_{k=1}^m(a_k-a_{k+1})\lambda_{k-1}
=n\,\operatorname{SB}_{\rm P}.                                  \tag{6}
\]

Now put

\[
S=\sum_{i=1}^mE_i\sim\Gamma(n+1,1),
\qquad D_i=E_i/S,
\qquad T_a=\sum_{i=1}^m a_iD_i.
\]

The Dirichlet vector `D` is independent of `S`, and

\[
Z_a:=\sum_i a_iE_i=S T_a.                                        \tag{7}
\]

The Gaffke endpoint for the observed sample is `Q_(1-alpha)(T_a)`.
Therefore the following two statements would combine with (2)--(7) to prove
that Poisson Stringer pointwise dominates the valid Gaffke endpoint for every
sample size:

**A. Exponential-quantile convexity**

\[
a\longmapsto q_\alpha(a)\quad\text{is convex on }\mathbb R_+^m.  \tag{A}
\]

**B. Dirichlet Poissonization**

\[
Q_{1-\alpha}(T_a)\le \frac1n Q_{1-\alpha}(Z_a).                  \tag{B}
\]

Both inequalities are strongly supported at conventional levels, but neither
is proved here.

## 3. A sharp obstruction for route A

The numerical evidence suggests that (A) holds at least when

\[
\alpha\le \alpha_\star:=4e^{-3}\approx0.1991482735.              \tag{8}
\]

The constant cannot be increased in any dimension-free convexity theorem.
This necessity is exact already in dimension two.  Let

\[
q(t)=Q_{1-\alpha}(E_1+tE_2)
\]

and let `q_0=q(1)`, so `(1+q_0)e^{-q_0}=alpha`.  Implicit
differentiation of the two-exponential tail at `t=1` gives

\[
q'(1)=\frac{q_0}{2},
\qquad
q''(1)=\frac{q_0(q_0-3)}{12}.                                  \tag{9}
\]

Hence local convexity at equal weights requires `q_0>=3`, which is exactly
`alpha<=4e^-3`.  Just above that tail probability, (A) fails locally.

### A differential target

There is also an exact Hessian identity that isolates what remains to prove.
For `X=a dot E`, a direction `z`, `V=z dot E`, and a smooth quantile
`q(t)=Q_(1-alpha)(X+tV)`, differentiation of
`Pr{X+tV<=q(t)}=1-alpha` gives

\[
q'(0)=\mathbb E[V\mid X=q],                                     \tag{10}
\]

and

\[
q''(0)
=-\frac{1}{f_X(q)}
  \left.\frac{d}{dx}\left{
  f_X(x)\operatorname{Var}(V\mid X=x)\right}\right|_{x=q}.     \tag{11}
\]

Thus (A) is equivalent to showing that the density-weighted conditional
variance in (11) is nonincreasing at the relevant high quantile for every
coefficient vector and every direction.  At equal weights it is proportional
to `x^(m+1)e^(-x/a)`, whose turning point yields (9).  A successful proof must
show that unequal weights cannot create a worse high-tail turning point.

Random Jensen tests conducted during exploration found no violation of (A)
at or below (8) in dimensions through 12 and found the predicted failures
just above (8) in dimension two.  These are heuristic checks only.

## 4. A kernel and localization target for route B

For `y>0`, set

\[
r_n(y)=\Pr\{\Gamma(n+1,1)>n/y\}
=e^{-n/y}\sum_{k=0}^n\frac{(n/y)^k}{k!}.                         \tag{12}
\]

After scaling a candidate threshold to one, the difference relevant to (B)
is

\[
\Pr\{ST_a>n\}-\Pr\{T_a>1\}
=\mathbb E\{r_n(T_a)-\mathbf1_{\{T_a>1\}}\}.                   \tag{13}
\]

The sign-changing kernel in (13) has the exact antiderivative identity

\[
\frac1{n!}\frac{d^n}{dy^n}
\left[y^ne^{-n/y}-(y-1)_+^n\right]
=r_n(y)-\mathbf1_{\{y>1\}}.                                   \tag{14}
\]

A nondegenerate projection of the uniform `n`-simplex has a
`1/(n-1)`-concave density; equivalently its probability law is
`1/n`-concave.  A natural route is therefore a one-dimensional
`s`-concave localization theorem.  It should reduce the sign of (13), under
the relevant tail constraint, to densities proportional to

\[
(A+By)^{n-1}\mathbf1_{[\ell,r]}(y).                              \tag{15}
\]

For (15), integrating (14) by parts `n` times leaves only endpoint terms.
Exploratory checks of these endpoint expressions and random Dirichlet
profiles support (B) for `alpha<=1/2`.  Turning that observation into a
complete localization statement, including every endpoint and degenerate
case, is the remaining task.

## 5. Shortcuts that do not work

The valid SymPol upper interval is easier to compute than the Gaffke
quantile, but Poisson Stringer does not dominate it for all `n`, even at
conventional levels.  For example, at `alpha=0.05`, `n=15`, and a sample
containing one taint equal to one and fourteen zeros, the numerical values
are

\[
U_{\rm SymPol}\approx0.3187079309,
\qquad
\operatorname{SB}_{\rm P}\approx0.3162576346.                   \tag{16}
\]

Thus SymPol cannot supply the desired pointwise comparison.  A
mean-only Bernoulli-KL confidence bound also exceeds Poisson Stringer on
some samples at 95% confidence.  Generic Chernoff and simultaneous-band
bounds are still wider.  These checks prevent replacing (A)--(B) with a
seemingly simpler but false domination claim.

## 6. Next proof tasks

1. Prove the high-tail monotonicity in (11), or find the first exact
   counterexample to (A) below `4e^-3`.
2. State and verify the precise one-dimensional localization theorem needed
   for (13), including whether the tail constraint is preserved under
   localization.
3. Reduce the `s`-affine endpoint expression from (14)--(15) to a single
   sign inequality and identify its sharp tail threshold.
4. Only after both analytic steps pass, connect (6)--(7) to the published
   Gaffke validity theorem and state the resulting all-`n` Poisson coverage
   corollary.

Until then, the implemented
[`Stringer--Gaffke safeguard`](GAFFKE-SAFEGUARD.md) is the valid all-sample-
size reporting rule; ordinary Poisson Stringer remains proved only in the
ranges stated in the manuscript.
