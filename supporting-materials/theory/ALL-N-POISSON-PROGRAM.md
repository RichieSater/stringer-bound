# Program for an all-`n` Poisson-Stringer theorem

> **Status.** This note separates proved reductions from two analytic targets.
> Exponential-quantile convexity is now proved sharply in two coordinates but
> remains open in three or more; Dirichlet Poissonization remains open in
> general, although it is proved for every profile with at most two nonzero
> coefficients, for dimension-free convex-core, far-cap, and middle-knot
> budget regions with three nonzero coefficients, and completely for `n=2`
> and `n=3`.  This is a research roadmap,
> not a coverage claim.  The paper
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

The two-coordinate case of (A) is proved sharply below.  Inequality (A) in
three or more coordinates and inequality (B) in general are strongly
supported at conventional levels but remain unproved.

## 3. The exact two-coordinate boundary for route A

The following exact theorem settles the first nontrivial dimension of (A).

\[
\alpha\le \alpha_\star:=4e^{-3}\approx0.1991482735.              \tag{8}
\]

More precisely, for independent unit exponentials `E_1,E_2`, the map

\[
(a,b)\longmapsto Q_{1-\alpha}(aE_1+bE_2)
\]

is convex on the nonnegative quadrant if and only if (8) holds.  Thus route
A is proved for two coordinates at 90%, 95%, and 99% confidence, and the
constant cannot be increased in any dimension-free version of (A).

Here is the local calculation that gives necessity.  Let

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

Sufficiency is global.  The exact hypoexponential tail and the conditional-
variance identity below reduce convexity on `0<t<1` to a one-variable tail
comparison.  The logarithmic comparison margin has a derivative whose sign
reduces to

\[
K(z)=2z^2e^z+ze^{2z}-z-4e^{2z}+8e^z-4>0.
\]

The coefficients of `K` vanish through degree five and are positive from
degree six onward.  Symmetry and the perspective construction then give
convexity on the full quadrant.  The complete proof, including the exact
factorization and all endpoint arguments, is in
[`TWO-EXPONENTIAL-QUANTILE.md`](TWO-EXPONENTIAL-QUANTILE.md); the identities
are regenerated by `make all-n-reduction-check`.  No floating-point sign
decision enters that theorem.

### The equal-weight Hessian in every dimension

The local calculation also has a closed form for every `m>=2`.  Let `x_m`
be the dimensionless gamma quantile

\[
 \Pr\{\Gamma(m,1)>x_m\}=\alpha.
\]

At the equal-weight vector `c 1`, conditionally on
`c\sum_iE_i=cx_m`, the normalized exponential vector is uniform Dirichlet.
For a direction `z=(z_1,...,z_m)`, its conditional variance is

\[
 \operatorname{Var}\!\left(\left.\sum_i z_iE_i\,\right|
 c\sum_iE_i=cx_m\right)
 =\frac{x_m^2}{m^2(m+1)}
 \left\{m\sum_i z_i^2-\left(\sum_i z_i\right)^2\right\}.
\]

The density of `c\sum_iE_i` is proportional to
`x^(m-1)e^(-x/c)`.  Substitution in (11) therefore gives the exact Hessian
quadratic form

\[
 z^{\mathsf T}\nabla^2q_\alpha(c\mathbf1)z
 =\frac{x_m(x_m-m-1)}{c\,m^2(m+1)}
 \left\{m\sum_i z_i^2-\left(\sum_i z_i\right)^2\right\}.       \tag{9a}
\]

Thus the Hessian at equal weights is positive semidefinite exactly when

\[
 \alpha\le\alpha_m
 :=e^{-(m+1)}\sum_{j=0}^{m-1}\frac{(m+1)^j}{j!}.                \tag{9b}
\]

These thresholds are strictly increasing for `m>=2`.  To see this, put
`n=m+1`, let `F_lambda(k)=Pr{Pois(lambda)<=k}`, and let `p_lambda(k)` be the
corresponding point mass.  Then `partial_lambda F_lambda(k)=-p_lambda(k)` and

\[
 \alpha_{m+1}-\alpha_m
 =p_n(n-1)-\int_n^{n+1}p_\lambda(n-1)\,d\lambda>0,
\]

because `p_lambda(n-1)` is strictly decreasing for `lambda>=n`.  Hence the
smallest equal-weight threshold is
`alpha_2=4e^-3`.  In particular, no equal-weight Hessian supplies a stronger
obstruction than the two-coordinate theorem.  This is only a pointwise
Hessian result; the unresolved difficulty is still to control all unequal
weight vectors in three or more coordinates.

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

Thus higher-dimensional (A) is equivalent to showing that the
density-weighted conditional variance in (11) is nonincreasing at the
relevant high quantile for every
coefficient vector and every direction.  At equal weights it is proportional
to `x^(m+1)e^(-x/a)`; for `m=2`, its turning point gives the threshold in
(9).  A successful proof must
show that unequal weights cannot create a worse high-tail turning point.  The
first unresolved dimension is now three.

For that first open dimension, the conditional law can be normalized to an
exponentially tilted uniform law on a triangle.  The resulting two-by-two
curvature matrix and hypoexponential tail reduce global convexity at
`alpha=4e^-3` to one explicit inequality in two nonnegative gap variables;
see
[`THREE-EXPONENTIAL-QUANTILE.md`](THREE-EXPONENTIAL-QUANTILE.md).  That
reduction is exact.  An integration-by-parts identity in the same note
expresses the curvature threshold as the largest generalized eigenvalue of a
boundary Gram matrix against the bulk covariance matrix; hence an upper
threshold can be certified by a two-by-two positive-semidefinite test without
third central moments.  The same note proves the inequality on both
coordinate axes, corresponding to coefficient vectors with two equal maximal
weights, and proves that the logarithmic margin has a strictly positive
fixed-sum inward derivative at every finite point of those axes. It also
proves the diagonal, corresponding to two equal smaller weights, and shows
that every positive finite diagonal point is a strict local maximum of the
margin under fixed-total-gap asymmetry.
The complete infinite-gap boundary reduces to the proved two-exponential
case, and a joint expansion proves positivity in a punctured neighborhood of
its only zero-margin corner.  An analytic boundary-versus-bulk comparison
also proves every point at which both gaps are at least `13`, while a direct
tail comparison, preceded by a fixed-total-gap minimization, proves the
triangle in which the gaps sum to at most `8/9`.  The remaining finite
off-symmetry region has total gap above `8/9` and smaller gap below `13`.

Random Jensen tests conducted during exploration found no violation of (A)
at or below (8) in dimensions through 12 and found the failures predicted by
(9) just above (8).  These are heuristic checks only and are not used in the
two-coordinate proof.

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
apply.  These comparisons explain why higher-dimensional route A remains
plausible without
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

There is a second dimension-free family.  If all but two knots vanish and
the remaining values are `0<=a<=b` with `a+b<=n`, then

\[
 [\underbrace{0,\ldots,0}_{n-1},a,b]H_n
 =[a,b]\frac{H_n(u)}{u^{n-1}}\ge0.                 \tag{16a}
\]

The last sign follows from a one-variable monotonicity argument that uses the
sum constraint to cover all possible positions of `a` and `b`.  Combining
this theorem with the radial lemma proves the complete `n=2` threshold
comparison: every three-knot profile moves to a boundary profile with at
most two nonzero knots.  At the tail levels in this program, this proves (B)
in three coordinates; the general dimension remains open.

The next complete dimension is now proved as well.  For `n=3`, the radial
lemma moves every four-knot profile to an ordered boundary profile
`(0,a,b,c)` with `a+b+c<=3`.  Multiplication by the zero knot reduces (16)
to

\[
 [a,b,c]f_3\ge0,
 \qquad
 f_3(u)=u^2\{e^{-3/u}-(1-u^{-1})_+^3\}.             \tag{16b}
\]

Although `f_3` is not globally convex on `[0,3]`, an exact secant argument
proves (16b) on every ordered triple satisfying the sum constraint.  It uses
strict convexity through `7/5`, the one-turn shape of `f_3'`, and the split
`b=9/10`; every exponential comparison is reduced to a rational Taylor
bound.  Thus (B) is complete in four coordinates too.  This does not prove
route (A) in four coordinates and therefore is not, by itself, a Stringer
coverage theorem.

There is also a dimension-free three-positive region.  With `n-2` zero knots
and ordered positive knots `a<=b<=c`, multiplication reduces (16) to the
second divided difference of

\[
 g_n(u)=u^2\{e^{-n/u}-(1-u^{-1})_+^n\}.
\]

Its curvature has the exact probability representation

\[
 \frac{g_n''(u)}2
 =\Pr\{\operatorname{Pois}(n/u)\le2\}
  -\Pr\{\operatorname{Bin}(n,1/u)\le2\}.           \tag{16c}
\]

The Anderson--Samuels inequality therefore makes `g_n` convex through
`n^2/{2(n+1)}`.  The comparison is proved whenever
`c<=n^2/{2(n+1)}`.  Separate analytic arguments also prove the far cap
`c>=n-1` and the full middle-knot region `b<=n/3`.  The part of the
three-positive face outside those regions remains open in general dimensions.

The radial derivative is the density derivative at the threshold.  A
Laplace-transform convolution identity shows that the mode of every
weighted exponential sum is no larger than its mean, giving the required
sign.  In dimensions `n>=4`, coordinate-face profiles with three or more
nonzero knots and at least three distinct coefficient values remain open
outside the proved convex-core, far-cap, and middle-knot regions.  See
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

1. Prove the explicit tilted-simplex inequality at finite off-symmetry
   points `0<z<w<infinity` in
   [`THREE-EXPONENTIAL-QUANTILE.md`](THREE-EXPONENTIAL-QUANTILE.md), or find
   an exact three-coordinate counterexample to (A) below `4e^-3`; the
   repeated-maximum boundary `zw=0`, equal-smaller diagonal `z=w`, and
   infinite-gap boundary are already proved, as is a punctured neighborhood
   of the sharp corner.  The fixed-total-gap derivative points strictly
   inward from each finite axis, and the corresponding second derivative is
   strictly negative at each positive diagonal point.  The region in which
   both gaps are at least `13` and the triangle in which their sum is at most
   `8/9` are also proved.  Then seek a dimension-free argument for (11).
2. Prove the constrained divided-difference inequality (16), or find an
   actual B-spline counterexample satisfying the sum-of-knots constraint;
   the complete `n=2` and `n=3` cases and every two-level profile are already
   proved.  The three-positive face is also proved whenever its largest knot
   is at most `n^2/{2(n+1)}`.  The first open simplex dimension is `n=4`.
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
