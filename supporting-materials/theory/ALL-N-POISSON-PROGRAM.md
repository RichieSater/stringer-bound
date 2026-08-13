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

### What established gamma-tail theory does and does not supply

There is a substantial comparison theory for weighted gamma sums, but its
strongest results stop short of (A).  Bock, Diaconis, Huffer, and Perlman
proved upper- and lower-tail Schur orderings in explicit regions, and
Diaconis and Perlman studied the location of the crossing point.  Yu later
proved their unique-crossing conjecture for iid gamma variables of shape at
least one, which includes the exponentials in (3).  Unique crossing says that
majorization eventually orders the upper tails.  It does not locate every
vector-dependent crossing below the quantile in (3), and tail ordering under
majorization is weaker than the Jensen inequality required for convexity of
`q_alpha`.

The most useful dimension-free tail region presently available is also too
far into the tail.  Specializing the theorem of Roosta-Khorasani and Székely
to unit exponentials, if two weight vectors have common sum `s` and one
majorizes the other, the required upper-tail ordering holds for

\[
 x>\frac32s.
\]

For `m` equal weights, however, the ratio of the fixed-level upper quantile
to `s` is the corresponding `Gamma(m,1)` quantile divided by `m`, which tends
to one as `m` grows.  Thus even the 90%, 95%, and 99% quantiles eventually
fall below the `3s/2` region.  The theorem cannot yield a dimension-free
conventional-level result.

Geometrically, (A) is a one-sided quantile-surface tangency problem for the
product exponential measure: positive homogeneity is automatic, and the
missing property is subadditivity.  Floating-body theory gives a closely
related tangency theorem for symmetric log-concave measures.  The product
exponential law is log-concave but not symmetric, so that theorem does not
apply.  These comparisons explain why route A remains plausible without
turning the existing crossing or floating-body results into a proof.

## 4. A constrained divided-difference target for route B

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
`1/n`-concave. A generic one-dimensional `s`-concave localization theorem
reduces a linear extremal problem to densities proportional to

\[
(A+By)^{n-1}\mathbf1_{[\ell,r]}(y).                              \tag{15}
\]

For (15), integrating (14) by parts `n` times leaves only endpoint terms.
That enlarged class is nevertheless too broad, even after retaining the
necessary mean constraint. At `n=10` the explicit density proportional to

\[
 \left(\frac34+\frac{75y}{388}\right)^9
 \mathbf 1_{[0,97/75]}(y)
\]

has `E Y<10/11`, `Pr(Y>1)>7/15`, and `Pr(SY>10)<7/15` for an
independent `S~Gamma(11,1)`. All three comparisons are certified exactly.
Thus unrestricted `s`-affine endpoint checks cannot prove (B); the
B-spline structure must be retained.

The actual Dirichlet-average structure yields a sharper target. After
scaling a threshold to `n`, write the simplex-projection knots as
`y_0,...,y_n`. At conventional levels they satisfy `sum_i y_i<=n`, and the
Hermite--Genocchi formula turns (13)--(14) into

\[
 [y_0,\ldots,y_n]H_n\ge0,
 \qquad
 H_n(y)=y^ne^{-n/y}-(y-1)_+^n.                         \tag{16}
\]

There is also a general boundary reduction.  The radial deformation

\[
 y_i(t)=1+t(y_i-1)
\]

leaves the simplex-cap event unchanged and makes the poissonized tail
nonincreasing until a knot reaches zero.  Consequently, any counterexample
would have a counterexample on a coordinate face.  For a profile with at
most two distinct knot values, that endpoint is an equal-positive-block
profile, where the result reduces to the Anderson--Samuels
binomial--Poisson comparison.  This proves the constrained inequality for
every two-level profile throughout the sum-constrained domain.  In
particular, it includes every active-boundary profile

\[
 y=(\underbrace{a,\ldots,a}_{k},
    \underbrace{b,\ldots,b}_{n+1-k}),
 \qquad ka+(n+1-k)b=n.
\]

The radial derivative is the density derivative at the threshold.  A
Laplace-transform convolution identity shows that the mode of every
weighted exponential sum is no larger than its mean, giving the required
sign.  Coordinate-face profiles with three or more distinct positive knots
remain open.  See
[`DIRICHLET-POISSONIZATION.md`](DIRICHLET-POISSONIZATION.md) for the proof,
the full reduction, and the exact localization obstruction.

## 5. Shortcuts that do not work

The valid SymPol upper interval is easier to compute than the Gaffke
quantile, but Poisson Stringer does not dominate it for all `n`, even at
conventional levels.  For example, at `alpha=0.05`, `n=15`, and a sample
containing one taint equal to one and fourteen zeros, the numerical values
are

\[
U_{\rm SymPol}\approx0.3187079309,
\qquad
\operatorname{SB}_{\rm P}\approx0.3162576346.                   \tag{17}
\]

Thus SymPol cannot supply the desired pointwise comparison.  A
mean-only Bernoulli-KL confidence bound also exceeds Poisson Stringer on
some samples at 95% confidence.  Generic Chernoff and simultaneous-band
bounds are still wider.  These checks prevent replacing (A)--(B) with a
seemingly simpler but false domination claim.

## 6. Next proof tasks

1. Prove the high-tail monotonicity in (11), or find the first exact
   counterexample to (A) below `4e^-3`.
2. Prove the constrained divided-difference inequality (16), or find an
   actual B-spline counterexample satisfying the sum-of-knots constraint.
3. Use the zero-knot reduction recursively, or prove that a negative minimum
   on a coordinate face must have at most two distinct knot values.
4. Only after both analytic steps pass, connect (6)--(7) to the published
   Gaffke validity theorem and state the resulting all-`n` Poisson coverage
   corollary.

Until then, the implemented
[`Stringer--Gaffke safeguard`](GAFFKE-SAFEGUARD.md) is the valid all-sample-
size reporting rule; ordinary Poisson Stringer remains proved only in the
ranges stated in the manuscript.

## References

M. E. Bock, P. Diaconis, F. W. Huffer, and M. D. Perlman, “Inequalities for
linear combinations of gamma random variables,” *Canadian Journal of
Statistics* 15 (1987), 387--395,
<https://doi.org/10.2307/3315257>.

S. G. Bobkov, “Convex bodies and norms associated to convex measures,”
*Probability Theory and Related Fields* 147 (2010), 303--332,
<https://doi.org/10.1007/s00440-009-0209-7>.

P. Diaconis and M. D. Perlman, “Bounds for tail probabilities of weighted
sums of independent gamma random variables,” in *Topics in Statistical
Dependence*, IMS Lecture Notes--Monograph Series 16 (1990), 147--166,
<https://doi.org/10.1214/lnms/1215457557>.

F. Roosta-Khorasani and G. J. Székely, “Schur properties of convolutions of
gamma random variables,” *Metrika* 78 (2015), 997--1014,
<https://doi.org/10.1007/s00184-015-0537-9>.

Y. Yu, “On the unique crossing conjecture of Diaconis and Perlman on
convolutions of gamma random variables,” *The Annals of Applied Probability*
27 (2017), 3893--3910, <https://doi.org/10.1214/17-AAP1304>.
