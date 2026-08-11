# The constrained Dirichlet--Poissonization inequality

> **Status.** This note gives an exact reduction, proves a nontrivial family
> of profiles, and rules out an over-broad localization shortcut. The central
> divided-difference inequality remains open. Nothing here is yet an
> all-sample-size coverage theorem.

This is the second analytic component of the dimension-free Poisson
Stringer program. It isolates the part of that program that compares a
uniform-simplex average with its independent-gamma poissonization.

## 1. Threshold form

Let $E_0,\ldots,E_n$ be independent unit exponentials, put

\[
 S=\sum_{i=0}^n E_i,
 \qquad D_i=E_i/S,
\]

and let $y_0,\ldots,y_n\ge0$. Then $D=(D_0,\ldots,D_n)$ is uniform on
the $n$-simplex and is independent of $S\sim\Gamma(n+1,1)$. Define

\[
 T_y=\sum_i y_iD_i,
 \qquad Z_y=\sum_i y_iE_i=S T_y.
\]

After scaling an arbitrary comparison threshold to $n$, the desired
tail comparison becomes

\[
 \Pr\{Z_y>n\}\ \ge\ \Pr\{T_y>1\}.                 \tag{1}
\]

The relevant coefficient vectors satisfy

\[
 \sum_{i=0}^n y_i\le n.                              \tag{2}
\]

Indeed, before scaling, the threshold is an upper quantile of $Z_y$. At
tail probabilities below $e^{-1}$, that quantile lies above the mean: a
weighted sum of exponentials has a log-concave density, and the
one-dimensional Grünbaum inequality gives
$\Pr\{Z_y\ge E Z_y\}\ge e^{-1}$. Thus the scaling produces (2) at 90%,
95%, and 99% confidence.

## 2. An exact divided-difference target

For $u>0$, set

\[
 r_n(u)=\Pr\{\Gamma(n+1,1)>n/u\}
 =e^{-n/u}\sum_{j=0}^n\frac{(n/u)^j}{j!}
\]

and define $H_n(0)=0$ and

\[
 H_n(u)=u^ne^{-n/u}-(u-1)_+^n.                       \tag{3}
\]

Direct differentiation gives

\[
 \frac{H_n^{(n)}(u)}{n!}
 =r_n(u)-\mathbf 1_{\{u>1\}}                         \tag{4}
\]

away from $u=1$, with the identity interpreted through the continuous
$(n-1)$st derivative there. The Hermite--Genocchi formula now yields

\[
 \begin{aligned}
 \Pr\{Z_y>n\}-\Pr\{T_y>1\}
 &=E\!\left[r_n(T_y)-\mathbf1_{\{T_y>1\}}\right]\\
 &=[y_0,\ldots,y_n]H_n .                              \tag{5}
 \end{aligned}
\]

Repeated coefficients are understood as confluent divided differences.
Consequently, a sufficient dimension-free lemma is

\[
 \boxed{
 y_i\ge0,\quad \sum_i y_i\le n
 \quad\Longrightarrow\quad
 [y_0,\ldots,y_n]H_n\ge0.}                            \tag{6}
\]

Equation (6), rather than unrestricted $s$-concavity, is the current
analytic target. Its sum-of-knots hypothesis retains information that is
lost in a generic localization relaxation.

## 3. A proved family: equal positive blocks

The target (6) holds when exactly $k$ coordinates have a common positive
value $c$, and all remaining coordinates are zero. The sum condition is
$kc\le n$. If $c\le1$, the right side of (1) is zero. Otherwise, put
$\lambda=n/c$, so that $k\le\lambda<n$. Then

\[
 Z_y\ \stackrel d=\ c\,\Gamma(k,1),
 \qquad
 T_y\ \stackrel d=\ c\,\operatorname{Beta}(k,n+1-k).
\]

Therefore

\[
 \begin{aligned}
 \Pr\{Z_y>n\}
 &=\Pr\{\operatorname{Pois}(\lambda)\le k-1\},\\
 \Pr\{T_y>1\}
 &=\Pr\{\operatorname{Bin}(n,\lambda/n)\le k-1\}.
 \end{aligned}                                       \tag{7}
\]

For $\lambda>k$, the desired ordering is the Anderson--Samuels
binomial--Poisson inequality because $k-1<\lambda-1$. The endpoint
$\lambda=k$ follows by continuity. This proves (6) for every equal-block
profile, including the sparse boundary vector. At that sparse boundary
$(k=1,c=n)$, the positive gap is explicitly

\[
 e^{-1}-\left(1-\frac1n\right)^n>0.                  \tag{8}
\]

These are the equal-block vertices of the fixed-sum ordered coefficient
polytope. A complete proof must still control mixtures of those vertices.

## 4. Why generic $s$-concave localization is too broad

A one-dimensional localization theorem for $1/n$-concave probability
laws reduces a linear extremal problem to point masses or densities of the
form

\[
 (A+By)^{n-1}\mathbf1_{[\ell,r]}(y).                 \tag{9}
\]

That reduction is valid for the enlarged class, but proving the desired
tail implication for every law in (9) is impossible, even if the necessary
mean constraint $E Y\le n/(n+1)$ is retained. Here is an exact
counterexample to that proposed intermediate statement.

Take $n=10$ and let $Y$ have density proportional to

\[
 \left(\frac34+\frac{75y}{388}\right)^9,
 \qquad 0\le y\le\frac{97}{75}.                       \tag{10}
\]

The normalizing integral is $95984119/196608000$. Direct integration gives

\[
 E Y=\frac{729166363}{816359775}<\frac{10}{11}.
\]

Thus this law satisfies the mean restriction implied by (2). Its upper-tail
probability obeys

\[
 \Pr\{Y>1\}
 =\frac{3108309643939756140704768}
 {6633646218308706152889893}
 \approx0.4685672919>\frac7{15}.                       \tag{11}
\]

For an independent $S\sim\Gamma(11,1)$, repeated integration by parts
using (3)--(4) gives

\[
 \Pr\{SY>10\}
 =\frac{5573507995079350591862317513}
 {5265884111440931688376513}e^{-750/97}
 \approx0.4642055262<\frac7{15}.                       \tag{12}
\]

The mean and both inequalities against $7/15$ are certified with rational
arithmetic;
the exponential in (12) is enclosed by alternating-series bounds after
power-of-two range reduction. Thus the generic $s$-affine class contains
a genuine obstruction even after the necessary mean constraint is imposed.

This does **not** refute (1) for Dirichlet averages. A projection of the
uniform simplex has a B-spline density whose knots are the actual
coefficients $y_i$; not every density in (9) is such a projection with
the constraint (2). The example shows precisely which structural
information a successful localization argument must preserve.

## 5. Remaining proof problem

The cleanest current target is (6). Three plausible routes remain:

1. prove that the constrained divided difference in (6) has no interior
   negative minimum, reducing its extrema to equal-block boundary profiles;
2. exploit total positivity of the B-spline kernel together with the
   single-crossing theory for weighted gamma convolutions; or
3. prove the equivalent uniform-simplex cap inequality directly by a
   knot-merging or majorization argument that preserves $\sum_i y_i\le n$.

The exact certificate and regression checks are generated by
[`dirichlet_poissonization.py`](../computations/python/dirichlet_poissonization.py):

```sh
make dirichlet-poissonization-check
```

The committed output is
[`dirichlet-poissonization-certificate.json`](../computations/certificates/dirichlet-poissonization-certificate.json).
It certifies the obstruction and checks rational equal-block instances. The
all-parameter equal-block result rests on the written Anderson--Samuels
argument, and the general inequality remains explicitly open.

## References

T. W. Anderson and S. M. Samuels, “Some inequalities among binomial and
Poisson probabilities,” *Proceedings of the Fifth Berkeley Symposium on
Mathematical Statistics and Probability*, vol. 1, 1967, pp. 1--12.

M. Fradelizi and O. Guédon, “The extreme points of subsets of
$s$-concave probabilities and a geometric localization theorem,”
*Discrete & Computational Geometry* 31 (2004), 327--335,
<https://doi.org/10.1007/s00454-003-2868-y>.

M. Meyer, F. Nazarov, D. Ryabogin, and V. Yaskin, “Generalized
Grünbaum inequality,” *Bulletin of the London Mathematical Society* 50
(2018), 745--752, <https://doi.org/10.1112/blms.12175>.
