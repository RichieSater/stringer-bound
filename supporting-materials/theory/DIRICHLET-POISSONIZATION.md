# The constrained Dirichlet--Poissonization inequality

> **Status.** This note gives an exact reduction, pushes any possible
> counterexample to a zero-knot boundary, proves every profile having at most
> two distinct coefficient values, proves the complete comparison for
> `n=2`, and rules out an over-broad localization shortcut. The central
> divided-difference inequality remains open for general `n`. Nothing here is
> yet an all-sample-size coverage theorem.

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

## 4. A radial boundary reduction and every two-level profile

The next lemma applies to arbitrary coefficient vectors, not only to two
levels.  It reduces any possible interior counterexample to the boundary of
the nonnegative orthant.

> **Radial boundary lemma.**  Suppose $y_i>0$ and
> $\sum_{i=0}^n y_i\le n$.  Put $c_i=y_i-1$ and
>
> \[
>  y_i(t)=1+tc_i,
>  \qquad
>  1\le t\le t_*:=\min_{c_i<0}\frac{-1}{c_i}.        \tag{9}
> \]
>
> Then
>
> \[
> \Pr\{T_{y(t)}>1\}
> =\Pr\!\left\{\sum_i c_iD_i>0\right\}              \tag{10}
> \]
>
> is independent of $t$, whereas
> $t\mapsto\Pr\{Z_{y(t)}>n\}$ is nonincreasing.  At $t=t_*$,
> at least one coefficient is zero.

At least one $c_i$ is negative because $\sum_i y_i<n+1$, so $t_*$ is
well-defined.  Equation (10) follows immediately from

\[
 T_{y(t)}=1+t\sum_i c_iD_i.
\]

It remains to prove the monotonicity of the poissonized tail.  We first use
an elementary mode bound.  For arbitrary positive weights
$a_0,\ldots,a_n$, let $f$ be the density of $X=\sum_i a_iE_i$ and put
$\mu=EX=\sum_i a_i$.  The Laplace transform

\[
 L(s)=\prod_i(1+a_is)^{-1}
\]

satisfies

\[
 -L'(s)=L(s)\sum_i\frac{a_i}{1+a_is}.
\]

Laplace inversion therefore gives

\[
 x f(x)=\int_0^x f(x-u)k(u)\,du,
 \qquad
 k(u)=\sum_i e^{-u/a_i}.                              \tag{11}
\]

The density $f$ is log-concave because it is a convolution of log-concave
exponential densities.  If $r$ is a mode, (11) gives

\[
 r f(r)
 \le f(r)\int_0^\infty k(u)\,du
 =\mu f(r).
\]

Thus every such weighted exponential sum has a mode no larger than its
mean.

Now let $X_t=Z_{y(t)}$, let $f_t$ be its density, and set
$A(t)=\Pr\{X_t>n\}$.  Its mean is

\[
 \mu_t=n+1+t\left(\sum_i y_i-n-1\right)\le n.
\]

The mode bound and log-concavity imply $f_t'(n)\le0$.  To connect this
density derivative to $A'(t)$, write $S=\sum_iE_i$.  Differentiating the
tail by the coarea formula on the slice $X_t=n$ gives

\[
 A'(t)
 =f_t(n)E\!\left[\left.\sum_i c_iE_i\,\right|X_t=n\right]
 =\frac{f_t(n)}t\{n-E(S\mid X_t=n)\}.                \tag{12}
\]

The change of variables $e_i=xu_i$ in the density integral for $X_t$
shows that

\[
 x\frac{f_t'(x)}{f_t(x)}
 =n-E(S\mid X_t=x).
\]

Substitution into (12) yields

\[
 A'(t)=\frac nt f_t'(n)\le0.                          \tag{13}
\]

Continuity supplies the endpoint $t=t_*$.  This proves the lemma.
Consequently, if (6) had a counterexample with every knot positive, moving
along (9) would produce a counterexample with at least one zero knot.

The radial lemma also completes the two-level case.

> **Two-level corollary.**  Put $m=n+1$.  Let $k,l\ge1$ with $k+l=m$, and
> suppose the coefficient vector consists of $k$ copies of $a$ and $l$
> copies of $b$, where
>
> \[
>  0\le a\le b,
>  \qquad ka+lb\le n.                                \tag{14}
> \]
>
> Then $\Pr\{Z_y>n\}\ge\Pr\{T_y>1\}$.

If $b\le1$, the right side is zero.  If $a=0$, Section 3 applies directly.
Otherwise $0<a<1<b$, and the radial lemma moves the profile to the endpoint
$t_*=1/(1-a)$.  There the low coefficient is zero and the high coefficient
is

\[
 b_*=1+\frac{b-1}{1-a}=\frac{b-a}{1-a}.              \tag{15}
\]

The cap probability is unchanged, the poissonized tail can only decrease,
and the endpoint still satisfies $lb_*\le n$.  The equal-block result of
Section 3 proves the comparison at that endpoint and hence at the original
profile.


## 5. The complete comparison for `n=2`

The radial reduction closes the first nontrivial simplex dimension.

> **Three-knot theorem.**  Let `n=2` and let
> `y_0,y_1,y_2>=0` satisfy `y_0+y_1+y_2<=2`.  Then
>
> \[
>  \Pr\{Z_y>2\}\ge\Pr\{T_y>1\}.                   \tag{15a}
> \]
>
> Equivalently, the divided-difference target (6) holds for every three-knot
> profile.

By the radial boundary lemma and continuity, it is enough to consider an
ordered boundary vector `(0,a,b)` with `0<a<b` and `a+b<=2`.  Put

\[
 f(u)=\frac{H_2(u)}u
 =\begin{cases}
 u e^{-2/u},&0<u\le1,\\
 u e^{-2/u}-u+2-u^{-1},&u\ge1.
 \end{cases}                                       \tag{15b}
\]

Since `H_2(0)=0`, elementary divided-difference algebra gives

\[
 [0,a,b]H_2=\frac{f(b)-f(a)}{b-a}.                 \tag{15c}
\]

The function `f` is strictly increasing on `(0,1]`.  On `[1,2]`, its
derivative can be written as

\[
 f'(u)=g(1/u),
 \qquad
 g(x)=e^{-2x}(1+2x)-1+x^2,                         \tag{15d}
\]

and

\[
 g'(x)=2x\{1-2e^{-2x}\}>0,
 \qquad \frac12\le x\le1,                         \tag{15e}
\]

because `e>2`.  Thus `f'` is nonincreasing on `[1,2]`, so the minimum of the
concave function `f` on that interval occurs at an endpoint.  Those endpoint
values satisfy

\[
 f(1)=e^{-2},
 \qquad
 f(2)=2e^{-1}-\frac12>e^{-2}.                      \tag{15f}
\]

For the last inequality, multiply by `2e^2`: its numerator is
`4e-e^2-2=2-(e-2)^2>0`, using `2<e<3`.

If `b<=1`, monotonicity on `(0,1]` proves `f(b)>=f(a)`.  Otherwise
`a<=1<b<=2`, because `a+b<=2`; equations (15b)--(15f) give
`f(a)<=f(1)<=f(b)`.  Equation (15c) proves the boundary comparison, and the
radial lemma carries it back to every strictly positive three-knot profile.
Degenerate and repeated-knot cases follow by continuity.  This proves
(15a).


## 6. Why generic $s$-concave localization is too broad

A one-dimensional localization theorem for $1/n$-concave probability
laws reduces a linear extremal problem to point masses or densities of the
form

\[
 (A+By)^{n-1}\mathbf1_{[\ell,r]}(y).                 \tag{16}
\]

That reduction is valid for the enlarged class, but proving the desired
tail implication for every law in (16) is impossible, even if the necessary
mean constraint $E Y\le n/(n+1)$ is retained. Here is an exact
counterexample to that proposed intermediate statement.

Take $n=10$ and let $Y$ have density proportional to

\[
 \left(\frac34+\frac{75y}{388}\right)^9,
 \qquad 0\le y\le\frac{97}{75}.                       \tag{17}
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
 \approx0.4685672919>\frac7{15}.                       \tag{18}
\]

For an independent $S\sim\Gamma(11,1)$, repeated integration by parts
using (3)--(4) gives

\[
 \Pr\{SY>10\}
 =\frac{5573507995079350591862317513}
 {5265884111440931688376513}e^{-750/97}
 \approx0.4642055262<\frac7{15}.                       \tag{19}
\]

The mean and both inequalities against $7/15$ are certified with rational
arithmetic;
the exponential in (19) is enclosed by alternating-series bounds after
power-of-two range reduction. Thus the generic $s$-affine class contains
a genuine obstruction even after the necessary mean constraint is imposed.

This does **not** refute (1) for Dirichlet averages. A projection of the
uniform simplex has a B-spline density whose knots are the actual
coefficients $y_i$; not every density in (16) is such a projection with
the constraint (2). The example shows precisely which structural
information a successful localization argument must preserve.

## 7. Remaining proof problem

The cleanest current target is (6).  Section 4 shows that it is enough to
work on coordinate faces, Sections 3--4 establish it for every profile
having at most two distinct coefficient values, and Section 5 proves the
complete target for `n=2`.  Three plausible routes remain for general `n`:

1. exploit the zero-knot reduction recursively, or prove that a negative
   minimum on a coordinate face must have at most two distinct knot values;
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
It certifies the obstruction, checks rational equal-block instances, and
regression-checks the exact algebra used in the two-level lemma.
The all-parameter statements rest on the written arguments; the general
inequality remains explicitly open.

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
