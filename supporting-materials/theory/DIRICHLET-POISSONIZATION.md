# The constrained Dirichlet--Poissonization inequality

> **Status.** This note gives an exact reduction, pushes any possible
> counterexample to a zero-knot boundary, proves every profile having at most
> two distinct coefficient values or at most two nonzero coefficients, proves
> dimension-free convex-core and far-cap regions and the full middle-knot
> region `b<=n/3` for profiles with three nonzero coefficients,
> proves the complete comparison for `n=2` and `n=3`, and rules out an
> over-broad
> localization shortcut. The central
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


## 5. Two arbitrary positive knots in every dimension

The equal-block theorem can be extended to unequal positive coefficients
when all other coefficients vanish.

> **Two-positive-knot theorem.**  Let `n>=2`, and suppose the coefficient
> vector consists of `n-1` zeros and two values `0<=a<=b` satisfying
> `a+b<=n`.  Then
>
> \[
>  \Pr\{Z_y>n\}\ge\Pr\{T_y>1\}.                   \tag{15a}
> \]

The function `H_n` and all of its derivatives vanish at zero from the right.
The standard multiplication rule for divided differences therefore gives

\[
 [\underbrace{0,\ldots,0}_{n-1},a,b]H_n
 =[a,b]f_n,
 \qquad
 f_n(u)=\frac{H_n(u)}{u^{n-1}}
 =u\{e^{-n/u}-(1-u^{-1})_+^n\}.                    \tag{15b}
\]

It remains to prove `f_n(b)>=f_n(a)`.  The function is strictly increasing
on `(0,1]`.  For `1<u<=n`, differentiation gives

\[
 f_n'(u)=e^{-n/u}(1+n/u)
 -(1-u^{-1})^{n-1}\{1+(n-1)/u\}.                  \tag{15c}
\]

With `x=1/u`, the sign of (15c) is the sign of `R_n(x)-1`, where

\[
 R_n(x)=\frac{e^{-nx}(1+nx)}
 {(1-x)^{n-1}\{1+(n-1)x\}}.
\]

Direct differentiation yields

\[
 \frac{d}{dx}\log R_n(x)
 =\frac{-nx\{n(n-1)x^2+nx-1\}}
 {(x-1)(1+nx)\{1+(n-1)x\}}>0,
 \qquad \frac1n\le x<1.                           \tag{15d}
\]

Indeed, the polynomial in braces is increasing and equals `(n-1)/n>0` at
`x=1/n`.  Hence the sign of `f_n'` can change at most once on `[1,n]`, and
only from positive to negative as `u` increases.  Thus the minimum of `f_n`
on `[1,n]` occurs at an endpoint.

The endpoint at `n` is strictly larger.  The power-series identities for the
two logarithms and the termwise inequality `j2^j>=j+1` give

\[
 \left(1-\frac1n\right)^n
 <e^{-1}\left(1-\frac1{2n}\right).                 \tag{15e}
\]

Consequently,

\[
 f_n(n)>\frac1{2e}>e^{-n}=f_n(1),                  \tag{15f}
\]

where the second inequality follows from `e^{n-1}>2`.

One additional monotonicity detail handles the case in which both positive
knots exceed one.  For `n>=3`, put `r=1/(n-1)`.  At `x=r`,

\[
 \begin{aligned}
 \log R_n(r)
 &=-1-r+\log(1+r/2)-r^{-1}\log(1-r)\\
 &\ge -1-r+\left(\frac r2-\frac{r^2}{8}\right)
       +1+\frac r2+\frac{r^2}{3}
 =\frac{5r^2}{24}>0.                               \tag{15g}
 \end{aligned}
\]

The two displayed logarithmic bounds follow from their alternating and
positive power series.  By (15d), `R_n(x)>1` for `x>=1/(n-1)`, so `f_n` is
strictly increasing on `[1,n-1]`.

There are now three cases.  If `b<=1`, monotonicity on `(0,1]` applies.  If
`a<=1<b`, then (15f) and the endpoint-minimum argument give
`f_n(a)<=f_n(1)<=f_n(b)`.  Finally, if `1<a<=b`, the sum constraint gives
`b<=n-a<n-1`, so (15g) applies.  Thus `[a,b]f_n>=0` in every case.
Confluent and zero endpoints follow by continuity, proving (15a).

## 6. The complete comparison for `n=2`

The radial reduction now closes the first nontrivial simplex dimension.

> **Three-knot corollary.**  Let `n=2` and let
> `y_0,y_1,y_2>=0` satisfy `y_0+y_1+y_2<=2`.  Then
>
> \[
>  \Pr\{Z_y>2\}\ge\Pr\{T_y>1\}.                   \tag{15h}
> \]

If the profile is strictly positive, the radial boundary lemma moves it to a
profile with one zero coefficient, without changing the simplex tail and
without increasing the poissonized tail.  The endpoint is covered by the
two-positive-knot theorem.  Profiles already on the boundary and all
degenerate cases follow directly or by continuity.  This proves (15h).


## 7. The complete comparison for `n=3`

The next dimension also admits a closed analytic proof.

> **Four-knot corollary.**  Let `n=3` and let
> `y_0,y_1,y_2,y_3>=0` satisfy
> `y_0+y_1+y_2+y_3<=3`.  Then
>
> \[
>  \Pr\{Z_y>3\}\ge\Pr\{T_y>1\}.                    \tag{15i}
> \]

By the radial boundary lemma, it is enough to consider an ordered endpoint
`(0,a,b,c)` with

\[
 0\le a\le b\le c,
 \qquad a+b+c\le3.                                  \tag{15j}
\]

The divided-difference multiplication rule gives

\[
 [0,a,b,c]H_3=[a,b,c]f,
 \qquad
 f(u)=\frac{H_3(u)}u
 =u^2\{e^{-3/u}-(1-u^{-1})_+^3\},                  \tag{15k}
\]

where `f(0)=0`.  The elementary inequality
`1-u^{-1}<e^{-1/u}` also gives `f(u)>=0`.  We shall prove that
`[a,b,c]f>=0` under (15j),
although `f` is not convex on the whole interval `[0,3]`.

The required derivative structure is elementary.  For `u>1`,

\[
 \begin{aligned}
 f'(u)&=(2u+3)e^{-3/u}-(2u-3+u^{-2}),\\
 f''(u)&=(2+6u^{-1}+9u^{-2})e^{-3/u}-2+2u^{-3},\\
 f'''(u)&=u^{-4}\{27e^{-3/u}-6\}.                  \tag{15l}
 \end{aligned}
\]

On `(0,1]`, the subtraction term in (15k) vanishes and `f''>0`.
On `[1,7/5]`, equation (15l) gives `f'''<0`, while

\[
 f''(7/5)
 =\frac{533}{49}e^{-15/7}-\frac{436}{343}
 >\frac{564}{116963}>0.                             \tag{15m}
\]

Indeed, `9/2<e^(15/7)<341/40`; the upper bound gives
`e^(-15/7)>40/341`, which yields (15m).  Thus `f` is strictly convex on
`[0,7/5]`.

The last expression in (15l) changes sign at most once.  Hence `f''`
first decreases and then increases on `[1,3]`.  Since

\[
 f''(1)=17e^{-3}>0,
 \qquad
 f''(3)=5e^{-1}-\frac{52}{27}<0,                    \tag{15n}
\]

`f''` has exactly one zero on `[0,3]`.  In particular, `f'` first
increases and then decreases, so its minimum on every subinterval is attained
at an endpoint.  The inequality in (15n) follows already from `e>8/3`.

We use the following exact endpoint bounds:

\[
 \frac16<f'(9/10)<\frac7{40}<f'(3)<\frac15<f'(1),  \tag{15o}
\]

and

\[
 f'(3/2)>\frac13.                                   \tag{15p}
\]

For completeness, these follow from

\[
 \frac{192}{7}<e^{10/3}<\frac{144}{5},
 \qquad
 \frac{405}{149}<e<\frac{68}{25}<\frac{11}{4}.     \tag{15q}
\]

For example,
`f'(9/10)=(24/5)e^(-10/3)`,
`f'(3)=9/e-28/9`, and
`f'(3/2)=6/e^2-4/9`.  All bounds in (15m)--(15q) follow from finite
positive Taylor sums; for an upper bound, the remaining terms are dominated
by a geometric series.  Exact rational bounds used by the regression check
are recorded in the accompanying certificate.

We can now compare the two secant slopes.  If `c<=7/5`, strict convexity
already gives

\[
 [a,b]f\le[b,c]f.                                   \tag{15r}
\]

Suppose `c>=7/5`.  If `b<=9/10`, convexity on `[0,b]` gives
`[a,b]f<=f'(b)`.  Here `f'(b)<=f'(9/10)<f'(3)`.  Since
`c>=7/5>1`, the endpoint-minimum property on `[1,3]` and
`f'(1)>f'(3)` give `f'(c)>=f'(3)`.  Unimodality then shows that every
derivative on `[b,c]` is at least `f'(b)`.  Hence (15r) follows again.

It remains to consider `b>=9/10`.  The sum constraint also gives
`b<=3/2`, and

\[
 a\le3-b-c\le\frac85-b,
 \qquad b-a\ge2b-\frac85.                           \tag{15s}
\]

Define

\[
 Q(u)=\frac u3-\frac4{15}-f(u),
 \qquad \frac9{10}\le u\le\frac32.
\]

By (15o)--(15p) and the unimodality of `f'`, the derivative
`Q'(u)=1/3-f'(u)` changes sign at most once, from positive to negative.
Thus the minimum of `Q` occurs at an endpoint.  The same exact bounds give

\[
 Q(9/10)>\frac1{30}-\frac{189}{6400}>0,
 \qquad
 Q(3/2)>\frac{19}{60}-\frac{81}{256}
 =\frac1{3840}>0.                                   \tag{15t}
\]

Consequently, using `f(a)>=0` and (15s),

\[
 [a,b]f
 \le\frac{f(b)}{b-a}
 <\frac{2b-8/5}{6(b-a)}
 \le\frac16.                                       \tag{15u}
\]

On the other hand, (15n)--(15o) and unimodality show that
`f'(u)>1/6` throughout `[9/10,3]`.  Therefore `[b,c]f>1/6`, proving
(15r) in the final case.

For distinct knots, (15r) is exactly `[a,b,c]f>=0`.  Repeated and zero
knots follow by continuity.  Together with (15k) and the radial boundary
lemma, this proves (15i).


## 8. A dimension-free three-positive convex core

The same multiplication step gives a nontrivial all-`n` region with three
arbitrary positive knots.

> **Three-positive convex-core theorem.**  Let `n>=3`, and suppose the
> coefficient vector consists of `n-2` zeros and three values
> `0<=a<=b<=c` satisfying
>
> \[
>  a+b+c\le n,
>  \qquad
>  c\le c_n:=\frac{n^2}{2(n+1)}.                    \tag{15v}
> \]
>
> Then `Pr{Z_y>n}>=Pr{T_y>1}`.

As before, the flatness of `H_n` at zero and the divided-difference
multiplication rule give

\[
 [\underbrace{0,\ldots,0}_{n-2},a,b,c]H_n
 =[a,b,c]g_n,
 \qquad
 g_n(u)=u^2\{e^{-n/u}-(1-u^{-1})_+^n\}.            \tag{15w}
\]

For `u>1`, put `lambda=n/u`.  Direct differentiation yields the exact
probabilistic identity

\[
 \begin{aligned}
 \frac{g_n''(u)}2
 &=e^{-\lambda}\left(1+\lambda+\frac{\lambda^2}{2}\right)\\
 &\quad-\left\{
 (1-u^{-1})^n
 +\frac n u(1-u^{-1})^{n-1}
 +\binom{n}{2}u^{-2}(1-u^{-1})^{n-2}
 \right\}\\
 &=\Pr\{\operatorname{Pois}(\lambda)\le2\}
   -\Pr\{\operatorname{Bin}(n,\lambda/n)\le2\}.   \tag{15x}
 \end{aligned}
\]

The Anderson--Samuels comparison used in Section 3 makes the last difference
strictly positive whenever

\[
 2<\lambda-\frac{\lambda}{n+1},
 \quad\text{equivalently}\quad
 u<\frac{n^2}{2(n+1)}.
\]

At equality the weak inequality follows by continuity.  On `(0,1]`, the
binomial term in (15w) vanishes and direct differentiation gives
`g_n''(u)>0`.  Hence `g_n` is convex throughout `[0,c_n]`.  Under (15v),
all three knots lie in that interval, so `[a,b,c]g_n>=0`.  Equations
(5) and (15w) prove the theorem.  For `n=3`, Section 7 is strictly stronger:
it removes the restriction `c<=9/8`.


## 9. A dimension-free three-positive far cap

The opposite end of the three-positive face can also be settled in every
dimension.

> **Three-positive far-cap theorem.**  Let `n>=4`, and suppose the
> coefficient vector consists of `n-2` zeros and three values
> `0<=a<=b<=c` satisfying
>
> \[
>  a+b+c\le n,
>  \qquad
>  c\ge n-1.                                       \tag{15y}
> \]
>
> Then `Pr{Z_y>n}>=Pr{T_y>1}`.

We continue to use the function `g_n` from (15w).  We first record a
monotonicity fact.  On `(0,1]`, direct differentiation shows that `g_n` is
strictly increasing and strictly convex.  For `1<u<=n`, put
`lambda=n/u`.  Then

\[
 \frac{g_n'(u)}u
 =e^{-\lambda}(2+\lambda)
  -\left(1-\frac\lambda n\right)^{n-1}
   \left(2+\lambda-\frac{2\lambda}n\right).        \tag{15z}
\]

The ratio of the first term in (15z) to the second is

\[
 R_n(\lambda)=
 \frac{e^{-\lambda}(2+\lambda)}
 {(1-\lambda/n)^{n-1}(2+\lambda-2\lambda/n)}.
\]

For `1<=lambda<n`, direct differentiation gives

\[
 \frac{d}{d\lambda}\log R_n(\lambda)
 =\frac{-\lambda^2\{\lambda n-2\lambda+3n-2\}}
 {(\lambda+2)(\lambda-n)(\lambda n-2\lambda+2n)}>0. \tag{15aa}
\]

Moreover,

\[
 \begin{aligned}
 \log R_n(1)
 &=-1-(n-1)\log(1-1/n)-\log(1-2/(3n))\\
 &> -1+(n-1)\left(\frac1n+\frac1{2n^2}\right)
       +\frac2{3n}
 =\frac{n-3}{6n^2}>0.                              \tag{15ab}
 \end{aligned}
\]

Here we used the positive power series for `-log(1-x)`.  Equations
(15z)--(15ab) prove that `g_n` is strictly increasing on `[0,n]` for
`n>=4`.

The remaining ingredient is the endpoint estimate

\[
 g_n(n-1)>(n+1)^2e^{-n},\qquad n\ge4.              \tag{15ac}
\]

For `n=4`, the finite Taylor bounds `e<11/4` and `e^4>54` give

\[
 g_4(3)=9e^{-4/3}-\frac{16}{9}
 >\frac{17}{36}>\frac{25}{54}>25e^{-4}.            \tag{15ad}
\]

For `n>=5`, set

\[
 x=\frac n{n-1},\qquad y=\frac n{2(n-1)^2}.
\]

The first two nonzero terms of the logarithmic series and the inequality
`1-e^{-y}>ye^{-y}` yield

\[
 \begin{aligned}
 g_n(n-1)
 &=(n-1)^2\left\{e^{-x}
       -\left(1-\frac1{n-1}\right)^n\right\}\\
 &>(n-1)^2e^{-x}(1-e^{-y})
 >\frac n2e^{-x-y}.                                \tag{15ae}
 \end{aligned}
\]

Now `x+y=1+3/{2(n-1)}+1/{2(n-1)^2}<=45/32`.  The inequality

\[
 e^{\,n-45/32}>\frac{2(n+1)^2}{n},\qquad n\ge5,   \tag{15af}
\]

follows by induction.  At `n=5`, its left side exceeds
`e^3>131/8>72/5`; thereafter the left side is multiplied by `e>2`,
whereas the ratio of consecutive right sides is

\[
 \frac{n(n+2)^2}{(n+1)^3}<2.
\]

Combining (15ae)--(15af) proves (15ac).

We now compare the two secant slopes.  Condition (15y) implies `a+b<=1`,
so `b<=1`.  Convexity on `[0,1]`, monotonicity on `[0,n]`, and (15ac) give

\[
 \begin{aligned}
 [a,b]g_n
 &\le g_n'(1)=(n+2)e^{-n},\\
 [b,c]g_n
 &\ge\frac{g_n(n-1)-g_n(1)}n
 >(n+2)e^{-n}.                                     \tag{15ag}
 \end{aligned}
\]

The second line uses `c-b<=n`, `g_n(b)<=g_n(1)=e^{-n}`, and
`g_n(c)>=g_n(n-1)`.  Thus `[a,b,c]g_n>0`, with confluent endpoints supplied
by continuity.  Equation (15w) proves the theorem.  Together, Sections 8
and 9 settle both the convex core `c<=n^2/{2(n+1)}` and the far cap
`c>=n-1` of every three-positive face.  Section 10 below covers an additional
part of the intervening region.


## 10. A dimension-free middle-knot region

There is also a dimension-free region that can pass through the nonconvex
part of `g_n`.  It is controlled by the middle knot rather than by the largest
knot alone.

> **Three-positive middle-knot theorem.**  Let `n>=4`, and suppose the
> coefficient vector consists of `n-2` zeros and three values
> `0<=a<=b<=c` satisfying
>
> \[
>  a+b+c\le n,
>  \qquad
>  b\le\frac n3.                                    \tag{15ah}
> \]
>
> Then `Pr{Z_y>n}>=Pr{T_y>1}`.

The proof uses a one-crossing property of the curvature in (15x).  Write

\[
 D_n(\lambda)=
 \Pr\{\operatorname{Pois}(\lambda)\le2\}
 -\Pr\{\operatorname{Bin}(n,\lambda/n)\le2\}.
\]

For `0<lambda<n`, factor the binomial probability in (15x), and let
`S_n(lambda)` be the ratio of the Poisson probability to the binomial
probability.  Direct differentiation gives

\[
 \frac{d}{d\lambda}\log S_n(\lambda)
 =\frac{-\lambda^2P_n(\lambda)}
 {(\lambda-n)(\lambda^2+2\lambda+2)Q_n(\lambda)}, \tag{15ai}
\]

where

\[
 \begin{aligned}
 P_n(\lambda)
 &=(n-1)(n-2)\lambda^3+2n(n-2)\lambda^2+4n\lambda
   -2n(3n-2),\\
 Q_n(\lambda)
 &=(n-1)(n-2)\lambda^2+2n(n-2)\lambda+2n^2>0.
 \end{aligned}
\]

The polynomial `P_n` is strictly increasing on the positive half-line.  Thus
`S_n` first decreases and then increases.  Since `S_n(0)=1` and
`S_n(lambda)` tends to infinity as `lambda` tends to `n`, the sign of `D_n`
can change at most once, from negative to positive.  Anderson--Samuels gives
`D_n(lambda)<0` for `0<lambda<2` and `D_n(lambda)>0` for
`3<=lambda<n`.  Hence `g_n''` changes sign at most once on `[0,n]`, from
positive to negative as its argument increases.  In particular, `g_n'` has
at most one maximum.

We next need four exact endpoint comparisons.  Define

\[
 h_n(\lambda)=g_n'(n/\lambda)
 =\frac n\lambda\left\{
 e^{-\lambda}(\lambda+2)
 -(1-\lambda/n)^{n-1}(\lambda+2-2\lambda/n)
 \right\}.                                        \tag{15aj}
\]

For `z=1/n`, put

\[
 d_\lambda(z)=
 \log\frac{e^{-\lambda}(\lambda+2)}
 {(1-\lambda z)^{1/z-1}(\lambda+2-2\lambda z)}.
\]

The positive difference in braces in (15aj) can then be written as its first
term times `1-e^{-d_lambda(z)}`.  Expanding only logarithms gives the exact
series

\[
 d_\lambda(z)=\sum_{r\ge1}\left\{
 \frac{\lambda^{r+1}}{r+1}-\frac{\lambda^r}{r}
 +\frac1r\left(\frac{2\lambda}{\lambda+2}\right)^r
 \right\}z^r.                                     \tag{15ak}
\]

The following bounds hold for every integer `n>=4`:

\[
 \begin{aligned}
 d_1(z)&>-\log(1-z/6),&
 d_2(z)&>-\log(1-z),\\
 d_3(z)&<-\log(1-16z/5),&
 d_4(z)&<-\log(1-16z/3)\quad(n\ge6).              \tag{15al}
 \end{aligned}
\]

Here are exact coefficient checks for completeness.  For the first bound,
the coefficient difference from `-log(1-z/6)` is zero at order one, equals
`1/24` and `1/72` at orders two and three, and its remaining tail is bounded
below by

\[
 -\frac{z^4}{20(1-z)}
 -\frac{(z/6)^4}{4(1-z/6)}.
\]

This is smaller in magnitude than `z^2/24` when `0<z<=1/4`.  For the second
bound, every coefficient after the first strictly exceeds the corresponding
coefficient of `-log(1-z)`.

For the third bound, let `c_r` denote the coefficient of
`-log(1-16z/5)-d_3(z)`.  Its first eight values are

\[
 \frac12,-\frac1{10},-\frac{271}{300},-\frac{1327}{500},
 -\frac{7861}{1250},-\frac{1636591}{131250},
 -\frac{15185543}{875000},\frac{4360661}{625000}.
\]

The polynomial through order seven is positive on `(0,1/4]`: after division
by `z` it is decreasing, and at `z=1/4` its value before division is
`516239213/6144000000`.  For `r>=8`,

\[
 \frac{r c_r}{3^r}
 =\left(\frac{16}{15}\right)^r-2+\frac3{r+1}
  -\left(\frac25\right)^r>0.                       \tag{15am}
\]

Positivity holds directly at `r=8`; its forward difference is positive
thereafter.  Finally, for the fourth bound, the coefficient difference is
zero at orders one and two, while for `r>=3` its sign is the sign of

\[
 \left(\frac43\right)^r-\frac{3r-1}{r+1}-3^{-r}>0.
\]

The cases `n=4,5` will not require this last logarithmic bound.

Substitution into (15aj) now gives

\[
 h_n(1)>\frac1{2e}>\frac8{e^4}>h_n(4),
 \qquad
 h_n(2)>\frac2{e^2}>\frac{16}{3e^3}>h_n(3).       \tag{15an}
\]

For `h_n(4)` with `n=4,5`, simply omit the nonnegative binomial term in
(15aj) to get `h_n(4)<=3ne^{-4}/2<8e^{-4}`.  The middle inequalities in
(15an) use the finite Taylor bounds `e^3>131/8>16` and `e>8/3`.

We can now prove the derivative reflection

\[
 g_n'(n-2b)\ge g_n'(b),\qquad 0\le b\le n/3.      \tag{15ao}
\]

First suppose `b>=1`, and put `r=n/b>=3`.  The argument on the left of
(15ao) is `n-2b`, so its `lambda` coordinate is
`s=r/(r-2)`.  If `3<=r<=4`, then `2<=s<=3`.  The one-maximum property and
(15an) give `h_n(s)>=h_n(3)>=h_n(r)`.  If `r>=4`, then `1<=s<=2`, and
`h_n(s)>=h_n(1)>h_n(4)>=h_n(r)`.  This proves (15ao) for `b>=1`.

If `0<=b<=1`, then `n-2b` has `lambda` coordinate in `[1,2]`, so
`g_n'(n-2b)>=h_n(1)>1/(2e)`.  Convexity of `u^2e^{-n/u}` on `[0,1]` gives
`g_n'(b)<=g_n'(1)=(n+2)e^{-n}<1/(2e)`.  The last inequality starts from
`12<e^3` at `n=4` and propagates by induction.  Thus (15ao) holds in every
case.

The reflection already proves the theorem in the subregion `2b+c<=n`.  To
cover the rest of (15ah), we first strengthen the curvature estimate to

\[
 (n-2b)g_n''(b)\ge g_n'(b),
 \qquad 0\le b\le n/3.                             \tag{15ap}
\]

For `0<b<=1`, put `lambda=n/b>=n`.  The left side of (15ap) minus the
right side is

\[
 \frac n\lambda e^{-\lambda}
 (\lambda^3-3\lambda-6)>0.
\]

At `b=0`, (15ap) follows by continuity.

For `b>1`, let `A_n(lambda)` denote the difference in braces in (15aj).
By (15x), inequality (15ap) is equivalent to

\[
 E_n(\lambda):=2(\lambda-2)D_n(\lambda)-A_n(\lambda)\ge0,
 \qquad 3\le\lambda<n.                             \tag{15aq}
\]

The Poisson part of `E_n` is
`(lambda^3-3lambda-6)e^{-lambda}`.  Its binomial part is

\[
 \left(1-\frac\lambda n\right)^{n-2}
 \frac{R_n(\lambda)}{n^2},
\]

where

\[
 \begin{aligned}
 R_n(\lambda)={}&(n-1)(n-2)\lambda^3
 +3(n-2)\lambda^2\\
 &-3n(n-4)\lambda-6n^2.
 \end{aligned}
\]

This polynomial is positive for `3<=lambda<=n`.  The logarithmic derivative
of the ratio of the Poisson part to the binomial part is

\[
 \frac{-\lambda^2W_n(\lambda)}
 {(\lambda-n)(\lambda^3-3\lambda-6)R_n(\lambda)}, \tag{15ar}
\]

and `W_n(lambda)>0` on the same interval.  Both positivity claims are
transparent after putting `x=lambda-3` and `m=n-lambda`: the expansions of
`R_n` and `W_n` have only positive coefficients.  For reference,

\[
 \begin{aligned}
 R_n={}&x^5+2x^4m+12x^4+x^3m^2+21x^3m+56x^3\\
 &+9x^2m^2+78x^2m+126x^2+24xm^2+117xm+135x\\
 &+12m^2+54m+54.
 \end{aligned}
\]

Also,

\[
 \begin{aligned}
 W_n={}&x^7+2x^6m+16x^6+x^5m^2+29x^5m+104x^5\\
 &+13x^4m^2+162x^4m+350x^4+60x^3m^2+432x^3m+651x^3\\
 &+111x^2m^2+573x^2m+666x^2+84xm^2+372xm+348x\\
 &+36m^2+108m+72.
 \end{aligned}
\]

At `lambda=3`, the Poisson-to-binomial ratio whose logarithmic derivative is
displayed in (15ar) is at least one because

\[
 \left(1-\frac3{2n}\right)
 \left(1-\frac3n\right)^{n-2}<e^{-3}.             \tag{15as}
\]

Indeed, with `t=3/n`, the logarithm of the left side plus three is

\[
 -\sum_{r\ge2}\left\{
 \frac{r-2}{r(r+1)}+\frac1{r2^r}
 \right\}t^r<0.
\]

On `3<=lambda<n`, both the numerator and denominator in (15ar) are
negative.  The ratio is therefore increasing, so (15as) proves (15aq), and
hence (15ap).

One more endpoint comparison will be useful.  The series (15ak) gives

\[
 \begin{aligned}
 d_{3/2}(z)
 &>-\log\left(1-\frac{27}{56}z-\frac{225}{896}z^2\right),\\
 d_3(z)
 &<-\log\left(1-\frac{27}{10}z-\frac74z^2\right).
                                                               \tag{15at}
 \end{aligned}
\]

For the first inequality, the two series agree through order two.  If
`p,q` are defined by

\[
 (1-pz)(1-qz)=1-\frac{27}{56}z-\frac{225}{896}z^2,
\]

then `0<p<6/7` and `|q|<6/7`; these bounds follow from
`p=(27+sqrt(3879))/112` and `3879<69^2`.  Moreover,

\[
 r[z^r]d_{3/2}(z)
 =\left(\frac32\right)^r\frac{r-2}{2(r+1)}
  +\left(\frac67\right)^r.
\]

For odd `r>=3`, this exceeds `p^r+q^r` because `q<0`.  For even `r>=4`,
use

\[
 \left(\frac32\right)^r\frac{r-2}{2(r+1)}
 >\left(\frac67\right)^r,
\]

which holds at `r=4` and propagates because `(7/4)^r(r-2)/(2(r+1))`
is increasing.

For the second inequality in (15at), define `c_r` as the coefficient of the
right side minus `d_3`.  The coefficients `c_2,c_6,\ldots,c_{10}` are
positive, while only `c_3,c_4,c_5` are negative among the first ten orders.
The polynomial

\[
 \sum_{r=2}^9c_r z^{r-2}
\]

has positive degree-eight Bernstein coefficients on `[0,1/4]`; the smallest
is `264755763361/68812800000000`.  For all `r>=11`, factor

\[
 1-\frac{27}{10}z-\frac74z^2=(1-pz)(1-qz).
\]

Here `p=(27+sqrt(1429))/20`, so `37^2<1429<39^2` gives
`16/5<p<33/10` and `|q|<3/5`.  Both logarithmic series therefore converge
on `[0,1/4]`.  Direct coefficient extraction gives

\[
 rc_r=p^r+q^r-3^r\frac{2r-1}{r+1}-\left(\frac65\right)^r.
\]

The elementary inequality

\[
 \left(\frac{16}{5}\right)^r-\left(\frac35\right)^r
 >2\,3^r+\left(\frac65\right)^r
\]

holds at `r=11`; subtracting three times its left-minus-right difference at
order `r` from that at order `r+1` leaves

\[
 \frac15\left(\frac{16}{5}\right)^r
 +\frac{12}{5}\left(\frac35\right)^r
 +\frac95\left(\frac65\right)^r>0.
\]

It therefore holds at every later order.  This proves `c_r>0` for the
remaining orders and completes (15at).

Substitution in (15aj) yields

\[
 \begin{aligned}
 h_n(3/2)&>e^{-3/2}\left(\frac98+\frac{75}{128}z\right),\\
 h_n(3)&<e^{-3}\left(\frac92+\frac{35}{12}z\right).
 \end{aligned}                                      \tag{15au}
\]

The ratio of the second parenthesis to the first is increasing in `z` and is
at most `8032/1953` for `0<=z<=1/4`.  Its square is smaller than `92/5`,
whereas the degree-five Taylor sum gives `e^3>92/5`.  Consequently,

\[
 h_n(3/2)>h_n(3).                                  \tag{15av}
\]

We can now finish the proof.  For `0<=b<=n/3`, define the endpoint tangent
remainder

\[
 K_n(b)=g_n(n-b)-g_n(b)-(n-2b)g_n'(b).
\]

Differentiation, (15ap), and the strict positivity of `g_n'` proved in
Section 9 give

\[
 K_n'(b)=g_n'(b)-g_n'(n-b)-(n-2b)g_n''(b)<0.
\]

At `b=n/3`, the two endpoint derivatives are `h_n(3)` and `h_n(3/2)`.
The one-maximum property and (15av) therefore show that
`g_n'(u)>=g_n'(n/3)` throughout `[n/3,2n/3]`, so `K_n(n/3)>0`.
Hence `K_n(b)>0` for every `0<=b<=n/3`.

For fixed `b`, put

\[
 F_b(v)=g_n(v)-g_n(b)-(v-b)g_n'(b).
\]

Because `g_n'` has at most one maximum, `F_b'` is either nonpositive
throughout `[b,n-b]` or is first nonnegative and then nonpositive.  Thus the
minimum of `F_b` occurs at an endpoint.  We have `F_b(b)=0` and
`F_b(n-b)=K_n(b)>0`, so `F_b(v)>=0` throughout that interval.

Finally, (15ah) and the sum constraint give `c<=n-a-b<=n-b`.  The curvature
is positive on `[0,b]`, so

\[
 [a,b]g_n\le g_n'(b).
\]

The nonnegativity of `F_b(c)` gives

\[
 [b,c]g_n\ge g_n'(b)\ge[a,b]g_n.
\]

Hence `[a,b,c]g_n>=0`; confluent cases follow by continuity.  Equation (15w)
proves the middle-knot theorem.


## 11. Why generic $s$-concave localization is too broad

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

## 12. Remaining proof problem

The cleanest current target is (6).  Section 4 shows that it is enough to
work on coordinate faces, Sections 3--5 establish it for every profile
having at most two distinct coefficient values or at most two nonzero
coefficients, and Sections 6--7 prove the complete target for `n=2` and
`n=3`.  Section 8 also proves every three-positive profile whose largest
knot is at most `n^2/{2(n+1)}`, Section 9 proves every such profile whose
largest knot is at least `n-1`, and Section 10 proves the full middle-knot
region `b<=n/3`.  Three-positive profiles outside those regions and all
profiles with four or more nonzero knots remain open.  Three plausible routes
remain for general `n`:

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
regression-checks the exact algebra used in the two-level, two-positive-knot,
three-positive convex-core, far-cap, and middle-knot regions, and complete
`n=2` and `n=3` proofs.
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
