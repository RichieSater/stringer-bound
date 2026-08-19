# The constrained Dirichlet--Poissonization inequality

> **Status.** This note gives an exact reduction, pushes any possible
> counterexample to a zero-knot boundary, proves every profile having at most
> two distinct coefficient values or at most two nonzero coefficients, proves
> every profile with three nonzero coefficients (including analytic
> convex-core, far-cap, and middle-knot subregions), proves a sparse
> convex-core region on every coordinate face (including a four-positive
> region in every dimension `n>=4`), proves the four-positive far cap
> `d>=n-1`,
> proves the complete comparison for `n=2`, `n=3`, `n=4`, and `n=5`, and
> rules out an over-broad localization shortcut. The central
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


## 8. A sparse convex-core theorem

The multiplication step and the Anderson--Samuels comparison extend to
every derivative order.  This gives a dimension-free core on every
coordinate face.

> **Sparse convex-core theorem.**  Let `n>=2` and `1<=k<=n-1`.  Suppose the
> coefficient vector consists of `n-k` zeros and `k+1` further values
> `0<=x_0<=...<=x_k` satisfying
>
> \[
>  \sum_{j=0}^k x_j\le n,
>  \qquad
>  x_k\le c_{n,k}:=\frac{n^2}{k(n+1)}.
> \]
>
> Then `Pr{Z_y>n}>=Pr{T_y>1}`.

The flatness of `H_n` at zero and `n-k` applications of the
divided-difference multiplication rule give

\[
 [\underbrace{0,\ldots,0}_{n-k},x_0,\ldots,x_k]H_n
 =[x_0,\ldots,x_k]F_{n,k},
 \qquad
 F_{n,k}(u)=u^k\{e^{-n/u}-(1-u^{-1})_+^n\}.
\]

For `u>1`, put `lambda=n/u`.  Direct induction in `k` gives

\[
 \begin{aligned}
 \frac{F_{n,k}^{(k)}(u)}{k!}
 &=e^{-\lambda}\sum_{j=0}^k\frac{\lambda^j}{j!}
   -\sum_{j=0}^k\binom nj u^{-j}(1-u^{-1})^{n-j}\\
 &=\Pr\{\operatorname{Pois}(\lambda)\le k\}
   -\Pr\{\operatorname{Bin}(n,\lambda/n)\le k\}.
 \end{aligned}
\]

For completeness, the induction uses the following elementary recurrence.
If

\[
 A_j(u)=\frac1{j!}\frac{d^j}{du^j}\{u^jh(u)\},
\]

then the product rule gives

\[
 A_{j+1}(u)=A_j(u)+\frac{u}{j+1}A_j'(u).
\]

Starting with `h(u)=e^{-n/u}` or `h(u)=(1-u^{-1})^n`, the corresponding
Poisson or binomial lower-tail CDF satisfies this recurrence because

\[
 -\lambda\frac d{d\lambda}
   \Pr\{\operatorname{Pois}(\lambda)\le j\}
 =(j+1)\Pr\{\operatorname{Pois}(\lambda)=j+1\}
\]

and, with `p=1/u`,

\[
 -p\frac d{dp}\Pr\{\operatorname{Bin}(n,p)\le j\}
 =(j+1)\Pr\{\operatorname{Bin}(n,p)=j+1\}.
\]

On `(0,1]`, the binomial term in `F_{n,k}` vanishes, and the same derivative
identity makes `F_{n,k}^{(k)}` positive.  The binomial component switched on
for `u>1` is `(u-1)^n/u^{n-k}`.  Because `k<=n-1`, its derivatives through
order `k` vanish at one, so `F_{n,k}^{(k)}` is continuous there.  Notice also
that `c_{n,k}>1`, since `k<=n-1`.

For `1<u<c_{n,k}`, the Anderson--Samuels inequality applies because

\[
 k<\lambda-\frac{\lambda}{n+1}.
\]

The endpoint `u=c_{n,k}` follows by continuity.  Consequently,
`F_{n,k}^{(k)}>=0` throughout `[0,c_{n,k}]`.  The Hermite--Genocchi formula
now gives `[x_0,...,x_k]F_{n,k}>=0`, including repeated knots by continuity,
and (5) proves the theorem.

The first consequence on a previously untreated coordinate face is worth
recording explicitly.

> **Four-positive convex-core corollary.**  Let `n>=4`, and suppose the
> coefficient vector consists of `n-3` zeros and four ordered nonnegative
> values `a<=b<=c<=d` whose sum is at most `n`.  If
>
> \[
>  d\le\frac{n^2}{3(n+1)},
> \]
>
> then `Pr{Z_y>n}>=Pr{T_y>1}`.  For the first open dimension `n=4`, the
> condition is `d<=16/15`.

This corollary is a proper region of the four-positive face, not a proof of
that complete face.  Taking `k=2` in the sparse theorem gives the following
three-positive subcase, for which it is useful to display the curvature
identity separately.

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

The same endpoint machinery also proves a far cap on the next coordinate
face.

> **Four-positive far-cap theorem.**  Let `n>=4`, and suppose the coefficient
> vector consists of `n-3` zeros and four values
> `0<=a<=b<=c<=d` satisfying
>
> \[
>  a+b+c+d\le n,
>  \qquad
>  d\ge n-1.
> \]
>
> Then `Pr{Z_y>n}>=Pr{T_y>1}`.

Write `F_{n,3}=u g_n` and put

\[
 Q(x,y,z)=[x,y,z]g_n.
\]

For distinct knots, two applications of the multiplication rule give

\[
 [a,b,c,d]F_{n,3}
 =\frac{dQ(b,c,d)-aQ(a,b,c)}{d-a}.                 \tag{15ag1}
\]

We first establish the endpoint estimate needed to compare the two terms.
Set

\[
 q_n=\frac{g_n''(1)}2
 =e^{-n}\left(1+n+\frac{n^2}{2}\right)
\]

and

\[
 \begin{aligned}
 T_n
 &=e^{-n}\left\{\frac1{n-1}+(n+2)
        +\frac{1+n+n^2/2}{3}\right\}\\
 &=e^{-n}\left\{\frac{n^2}{6}+\frac{4n}{3}
        +\frac73+\frac1{n-1}\right\}.
 \end{aligned}
\]

Then

\[
 \min\{f_n(n-1),f_n(n)\}>T_n,
 \qquad n\ge4.                                    \tag{15ag2}
\]

For `n=4`, the coefficient in braces is `32/3`.  Positive Taylor sums and a
geometric bound for the positive tail give

\[
 e^3>\frac{89641}{4480},\qquad
 e^{8/3}>\frac{12045015679}{837019575},\qquad
 e^4<\frac{553360529}{10135125}.
\]

Substitution in the two endpoint inequalities leaves the exact positive
margins

\[
 \begin{aligned}
 9e^{8/3}-\frac{16}{9}e^4-32
 &>\frac{29899236229}{66496555125},\\
 4e^3-\frac{81}{64}e^4-\frac{32}{3}
 &>\frac{6461863}{24024000}.
 \end{aligned}
\]

These are respectively `f_4(3)>T_4` and `f_4(4)>T_4` after multiplication
by the relevant positive power of `e`.

For `n>=5`, write

\[
 P_n=\frac{n^2}{6}+\frac{4n}{3}+\frac73+\frac1{n-1}.
\]

The exact identity

\[
 2P_n-P_{n+1}
 =\frac{(n+1)(n^3+4n^2-5n+6)}{6n(n-1)}>0
\]

and `e>2` show that `T_n=P_ne^{-n}` decreases.  The degree-eight Taylor
sum gives `e^5>1115309/8064>805/6`, and hence `T_5<1/10`.  On the other
hand, (15ae) gives

\[
 f_n(n-1)>\frac12e^{-45/32}>T_5;
\]

the last comparison follows from the degree-five bound

\[
 e^{115/32}>\frac{24748696103}{805306368}>\frac{161}{6}.
\]

Finally, (15f) gives `f_n(n)>1/(2e)>1/6>1/10`.  This proves (15ag2).

We now return to (15ag1).  The hypothesis implies
`a+b+c<=1`, hence `c<=1` and `a<=1/3`.  On `[0,1]`, the quantity
`g_n''(u)/2=Pr{Pois(n/u)<=2}` increases with `u`.  The
Hermite--Genocchi formula therefore gives

\[
 0\le Q(a,b,c)\le q_n.                             \tag{15ag3}
\]

Also, convexity on `[0,1]` gives

\[
 [b,c]g_n\le g_n'(1)=(n+2)e^{-n}.                 \tag{15ag4}
\]

Put `R_n(d)=\{g_n(d)-g_n(1)\}/d`.  The one-turn derivative argument in
Section 5 says that the minimum of `f_n` on `[n-1,n]` occurs at an endpoint.
Since `g_n=uf_n`, (15ag2) yields

\[
 \begin{aligned}
 R_n(d)
 &=f_n(d)-\frac{e^{-n}}d\\
 &>T_n-\frac{e^{-n}}{n-1}
  =(n+2)e^{-n}+\frac{q_n}{3}.                     \tag{15ag5}
 \end{aligned}
\]

Monotonicity of `g_n` from (15z)--(15ab), together with `c<=1`, now gives

\[
 [c,d]g_n
 =\frac{g_n(d)-g_n(c)}{d-c}
 \ge R_n(d).
\]

Combining this inequality with (15ag4)--(15ag5),

\[
 \begin{aligned}
 dQ(b,c,d)
 &=\frac d{d-b}\{[c,d]g_n-[b,c]g_n\}
 >\frac{q_n}{3}\\
 &\ge a q_n\ge aQ(a,b,c).
 \end{aligned}
\]

Thus (15ag1) is positive.  Repeated and zero knots follow by continuity,
proving the theorem.  This settles the region `d>=n-1` of every
four-positive face; it does not settle the intervening region between that
far cap and the sparse convex core in Section 8.


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


## 11. The complete three-positive face

The preceding regions can be joined by a derivative estimate that closes the
entire face.

> **Three-positive face theorem.**  Let `n>=4`, and suppose the coefficient
> vector consists of `n-2` zeros and three values `0<=a<=b<=c` satisfying
> `a+b+c<=n`.  Then
>
> \[
>  \Pr\{Z_y>n\}\ge\Pr\{T_y>1\}.
> \]

Put `f_n(u)=g_n(u)/u`, with `f_n(0)=0`, as in Section 5.  The divided-
difference multiplication rule gives

\[
 [a,b,c]g_n
 =\frac{c[b,c]f_n-a[a,b]f_n}{c-a}.               \tag{15aw}
\]

Confluent cases follow by continuity.  Sections 8--10 already cover
`c<=n^2/{2(n+1)}`, `c>=n-1`, and `b<=n/3`.  It remains to consider

\[
 \frac n3<b\le c,\qquad
 \frac{n^2}{2(n+1)}<c<n-1,
 \qquad a+b+c\le n.                               \tag{15ax}
\]

In particular,

\[
 c<\frac{2n}{3},
 \qquad
 a\le\frac{2n}{3}-c.                             \tag{15ay}
\]

We need two bounds for `phi_n=f_n'`.  With `lambda=n/u`, formula (15c)
becomes

\[
 \phi_n(n/\lambda)
 =e^{-\lambda}(1+\lambda)
  -\left(1-\frac\lambda n\right)^{n-1}
   \left(1+\lambda-\frac\lambda n\right).        \tag{15az}
\]

First, `phi_n` has only one maximum on `[0,n-1]`.  For `u>1`, direct
differentiation gives

\[
 f_n''(u)=\frac{\lambda^2}{u}
 \left\{e^{-\lambda}
 -(1-1/n)(1-\lambda/n)^{n-2}\right\}.             \tag{15ba}
\]

The logarithmic derivative, with respect to `lambda`, of the ratio of the
two terms in braces is

\[
 \frac{\lambda-2}{n-\lambda}.                    \tag{15bb}
\]

At `lambda=1` the ratio is less than one because
`(1-1/n)^{n-1}>e^{-1}`, whereas it tends to infinity as `lambda` tends to
`n`.  Thus `f_n''` changes sign only once as `u` increases, from positive to
negative.  On `(0,1]`, the same conclusion follows directly from
`f_n(u)=u e^{-n/u}`.

For `n>=5`, the maximum occurs before `n/3`.  Indeed,

\[
 (1-1/n)(1-3/n)^{n-2}>e^{-3}.                    \tag{15bc}
\]

For `5<=n<=12`, this follows from the degree-seven lower Taylor bound
`e^3>5557/280`; the smallest resulting rational margin is `353/21875`.
For `n>=13`, put `z=1/n`.  The elementary bound
`log(1-x)>=-x-x^2/{2(1-x)}` gives

\[
 \begin{aligned}
 &\log(1-z)+(z^{-1}-2)\log(1-3z)+3\\
 &\quad\ge
 \frac{z(15z^2-14z+1)}{2(z-1)(3z-1)}>0,
 \end{aligned}
\]

where the numerator at `z=1/13` is `2/169` and is larger for smaller
positive `z`.  Hence `phi_n` is decreasing on `[n/3,n-1]` when `n>=5`.

For `n=4`, its maximum lies before `8/5`, and
`phi_4(4/3)>phi_4(8/5)`.  The first assertion follows from
`e^{5/2}>256/27`.  After multiplication by `e^3`, the second follows from

\[
 4-\frac72e^{1/2}+\frac{413}{4096}e^3>0;
\]

the Taylor bounds used in the certificate leave the rational margin
`1297/15360`.  The one-maximum property now gives, under (15ax),

\[
 \phi_n(b)\ge\phi_n(c),
 \qquad
 [b,c]f_n\ge\phi_n(c).                            \tag{15bd}
\]

The second ingredient is the pair of uniform estimates

\[
 \phi_n(u)<\frac{3}{5n}\qquad(0<u\le n-1),      \tag{15be}
\]

and

\[
 n\phi_n(n/\lambda)\ge\frac{2\lambda-3}{5},
 \qquad
 \frac32\le\lambda\le2+\frac2n.                 \tag{15bf}
\]

Here are complete certification details.  For (15be), the finite prefix
`4<=n<=21` is covered by a 128-bit Arb interval proof.  It bisects
`1<=lambda<=n` until the whole image of (15az) lies strictly below the exact
rational `3/(5n)`.  Regeneration produces 2,006 certified terminal intervals;
the certificate records their counts by `n`, and the maximum bisection depth
is ten.  On `0<u<=1`, the
bound follows from `(1+lambda)e^{-lambda}<3/(5n)` at `lambda>=n`, starting
from `e^4>54`.

For the analytic tail `n>=22`, put `p=lambda/n` and let `d` be the logarithm
of the first term in (15az) divided by the second.  If `phi_n>0`, then
`d>0`.  The positive logarithmic series and `log(1+x)<=x` give

\[
 d\le
 \frac{(n-1)p^2(\lambda-1+p)}
 {2(1-p)(1+\lambda-p)}.
\]

Since `1-e^{-d}<=d`, it follows that

\[
 n\phi_n(n/\lambda)\le e^{-\lambda}U_n(\lambda),
\]

where

\[
 U_n(\lambda)=
 \frac{(1-1/n)(1+\lambda)\lambda^2
       (\lambda-1+\lambda/n)}
 {2(1-\lambda/n)(1+\lambda-\lambda/n)}.           \tag{15bg}
\]

If `p>=1/4`, the Poisson term alone suffices:
`n(1+n/4)e^{-n/4}<3/5` at `n=22`, and the left side decreases thereafter.
If `p<=1/4` and `lambda>=9/2`, then

\[
 U_n(\lambda)\le
 \frac23(1+\lambda)\lambda^2
 \frac{\lambda-3/4}{\lambda+3/4}.
\]

The right side after multiplication by `e^{-lambda}` is decreasing from
`lambda=9/2`; the exact bounds
`e^{9/2}>202948427/2293760>2475/28` put it below `3/5` there.
Finally, if `p<=1/4` and `lambda<=9/2`, the rational factor in (15bg) is
increasing in `1/n`, so `n=22` is worst.  After replacing `e^lambda` by its
degree-seven Taylor polynomial, the remaining cleared polynomial has
positive Bernstein coefficients on

\[
 [1,15/8],\ [15/8,11/4],\ [11/4,29/8],\ [29/8,9/2].
\]

The smallest is `6341019244847/25165824`.  This proves (15be).

For (15bf), the exact logarithmic series is

\[
 d=\sum_{r\ge1}\left\{
 \frac{\lambda^{r+1}}{r+1}-\frac{\lambda^r}{r}
 +\frac1r\left(\frac{\lambda}{\lambda+1}\right)^r
 \right\}n^{-r}.                                  \tag{15bh}
\]

Every coefficient is positive for `lambda>=3/2`.  If `n>=9`, retain only

\[
 c_1=\frac{\lambda^2(\lambda-1)}{2(\lambda+1)}.
\]

Using `1-e^{-d}>=d/(1+d)` reduces (15bf) to a decreasing function on
`[3/2,20/9]`.  The derivative claim is certified by four exact Bernstein
blocks with smallest coefficient `26185/189`.  At the right endpoint it
reduces to

\[
 e^{20/9}<\frac{2871000}{303433}.
\]

This follows from `e<11/4`, `e^{2/9}<5/4`, and
`605/64<2871000/303433`.  For `4<=n<=8`, retaining the first three terms in
(15bh) gives the same monotonicity.  Exact degree-five Taylor--Bernstein
checks and endpoint exponential enclosures are recorded in the certificate.
This proves (15bf).

We can now finish.  Put `lambda=n/c`.  Conditions (15ax)--(15ay) give the
range in (15bf).  Equations (15bd)--(15bf) yield

\[
 \begin{aligned}
 c[b,c]f_n
 &\ge c\phi_n(c)
 \ge\frac{2\lambda-3}{5\lambda},\\
 a[a,b]f_n
 &\le \frac{3a}{5n}
 \le\frac{2\lambda-3}{5\lambda}.
 \end{aligned}
\]

Thus the numerator in (15aw) is nonnegative.  Together with Sections 8--10,
this proves the three-positive face theorem.


## 12. The complete comparison for `n=4`

The next simplex dimension can also be closed.  The proof uses
exact scalar inequalities and two finite real-ball subdivisions; it is not a
floating-point search.

> **Five-knot theorem.**  Let `n=4` and let
> `y_0,...,y_4>=0` satisfy `sum_i y_i<=4`.  Then
>
> \[
>  \Pr\{Z_y>4\}\ge\Pr\{T_y>1\}.                  \tag{15bi}
> \]

By the radial boundary lemma, it is enough to treat
`(0,a,b,c,d)` with

\[
 0\le a\le b\le c\le d,
 \qquad a+b+c+d\le4.
\]

Put

\[
 F(u)=u^3\{e^{-4/u}-(1-u^{-1})_+^4\}.
\]

Use the continuous value `F(0)=0`.

The multiplication rule gives

\[
 [0,a,b,c,d]H_4=[a,b,c,d]F.                       \tag{15bj}
\]

If `W=(W_0,...,W_3)` is uniform on the three-simplex and
`U=aW_0+bW_1+cW_2+dW_3`, the Hermite--Genocchi formula gives

\[
 [a,b,c,d]F=E q(U),
 \qquad E U=\frac{a+b+c+d}{4}\le1,                \tag{15bk}
\]

where, with `lambda=4/u`,

\[
 q(u)=\frac{F'''(u)}6
 =\begin{cases}
 e^{-\lambda}(1+\lambda+\lambda^2/2+\lambda^3/6),&0<u\le1,\\
 e^{-\lambda}(1+\lambda+\lambda^2/2+\lambda^3/6)-1+u^{-4},&u>1.
\end{cases}                                      \tag{15bl}
\]

The two pieces make `F` a `C^3` function, so repeated knots are obtained by
continuity from the strict-knot formulas below.

We first remove two analytic regions.  For `u>1`, differentiation gives

\[
 u q'(u)=4\{\Pr(\operatorname{Pois}(\lambda)=4)
              -\Pr(\operatorname{Bin}(4,1/u)=4)\},
\]

and the ratio of the two masses is `(32/3)e^{-lambda}`.  Hence `q` decreases
on `[1,3/2]`.  At `u=19/15`, the inequality `q(u)>0` is equivalent to

\[
 e^{60/19}<\frac{110333}{4688}.
\]

The degree-seven Taylor upper bound is

\[
 e^{60/19}<\frac{778209455563}{33073254343},
\]

leaving the exact margin
`825443746875/155047416359984`.  Thus `q>=0` through `19/15`, and (15bk)
proves every profile with `d<=19/15`.

The second scalar bound is

\[
 q(u)\ge\frac{1-u}{6},
 \qquad \frac9{16}\le u\le\frac{37}{16}.         \tag{15bm}
\]

On `[9/16,1]`, the left side increases and the right side decreases, so the
claim follows from one rational Taylor upper bound at `9/16`.  For `u>1`,
putting `lambda=4/u` reduces (15bm) to

\[
 1+\lambda+\frac{\lambda^2}{2}+\frac{\lambda^3}{6}
 \ge e^\lambda\left(
 \frac76-\frac{2}{3\lambda}-\frac{\lambda^4}{256}
 \right).                                         \tag{15bn}
\]

The factor in parentheses is nonnegative for `64/37<=lambda<=4`.
Replacing `e^lambda` by its degree-six Taylor sum plus the geometric tail
majorant produces a rational polynomial after multiplication by the positive
denominator `3870720 lambda(8-lambda)`.  Its exact Bernstein coefficients on

\[
 [64/37,2], [2,5/2], [5/2,3], [3,7/2], [7/2,4]
\]

are positive; the smallest is `11480818817/25344`.  This proves (15bm).
If `a>=9/16`, the sum constraint gives `d<=37/16`, so (15bk)--(15bm) yield

\[
 [a,b,c,d]F\ge\frac{1-EU}{6}\ge0.                \tag{15bo}
\]

We next record the quantitative comparison used to prune another part of
the face.  Write

\[
 g(u)=u^2\{e^{-4/u}-(1-u^{-1})_+^4\},\quad
 f(u)=u\{e^{-4/u}-(1-u^{-1})_+^4\},\quad
 h(u)=e^{-4/u}-(1-u^{-1})_+^4.
\]

Each function is assigned its continuous value zero at `u=0`.

Exact Taylor--Bernstein bounds give

\[
 \frac{g''(u)}2\le\frac u4,
 \qquad u\ge0.                                    \tag{15bp}
\]

For `u<=1`, this follows because
`lambda Pr{Pois(lambda)<=2}` decreases for `lambda>=4` and is below one at
four.  For `u>1`, the multiplier
`1/lambda+1-lambda^3/16+3lambda^4/256` is positive on `0<lambda<=4`.
Multiplication by `e^lambda` and substitution of the degree-eight Taylor
lower bound therefore leave a polynomial whose exact Bernstein coefficients
on the eight half-unit subintervals of `[0,4]` are positive; the smallest is
`90499/640640`.

The first interval stage certifies the following auxiliary inequality:

\[
 [b,c,d]f\ge0
 \quad\text{whenever}\quad
 0\le b\le c\le d,\qquad b+c+d\le4.             \tag{15bq}
\]

Here `r=f''/2` satisfies `|r'|<3/5`.  On `u<=1`, this follows by maximizing
`e^{-lambda}lambda^4(lambda-3)/32` at `lambda=6`.  On `u>1`, the exact
formula is

\[
 r'(u)=\frac{15\lambda^6}{2048}-\frac{3\lambda^5}{64}
       +\frac{9\lambda^4}{128}
       +e^{-\lambda}\left(\frac{\lambda^5}{32}
                          -\frac{3\lambda^4}{32}\right).
\]

On `0<=lambda<=4`, 438 directed 160-bit Arb intervals enclose this expression
strictly inside `(-3/5,3/5)`, to maximum bisection depth 12.  The region
`d<=4/3` is analytic: for `u<=1`, `r>0`, while for `1<u<=4/3` the ratio of
the positive and negative mass terms is

\[
 \frac{4e^{-\lambda}}{3(1-\lambda/4)^2}.
\]

Its logarithmic derivative is `(lambda-2)/(4-lambda)`, and its value at
`lambda=3` exceeds one because `e^3<64/3`.  Hence `f''>=0` throughout the
core.
The function `f` is `C^2`, so confluent divided differences again follow by
continuity.

For the remainder, set

\[
 (b,c,d)=(b,b+s,b+s+t),qquad 3b+2s+t\le4.
\]

Since a partial derivative of `E r(bW_0+cW_1+dW_2)` with respect to any
knot has absolute value at most `(3/5)E W_i=1/5`, a parameter box has the
certified Lipschitz error

\[
 \frac15(3\Delta b+2\Delta s+\Delta t).           \tag{15br}
\]

Here each `Delta` is the larger distance from the chosen evaluation point to
the two endpoints in that coordinate; this remains valid if the radial
scaling moves the evaluation point away from the coordinate midpoint.

The initial box is `[0,4/3] x [0,2] x [0,4]`.  Before testing a box, each
upper endpoint is tightened using the other lower endpoints and
`3b+2s+t<=4`.  The center is the coordinate midpoint, scaled radially to the
budget boundary if necessary.  An unresolved box is bisected in the first
coordinate maximizing `w_i(U_i-L_i)`, with weights `(3,2,1)`; the children
are inserted lower then upper on the last-in, first-out stack.  These rules
and the pinned Arb precision make the transcript deterministic.

The deterministic longest-weighted-side subdivision terminates after
24,479 branch calls.  It has 53 analytic-core terminal boxes and 12,187
Lipschitz-certified terminal boxes, with maximum depth 20.  Every center is
evaluated by directed 160-bit Arb arithmetic.  The SHA-256 digest of the
ordered terminal transcript is

```text
baf76e5da205718ac2f7e7037bde03a4d02e5cacfa9351c14254c24f1ca31dfe
```

Thus (15bq) is an exhaustive interval proof, not evidence from sampled
points.

We use (15bq) to obtain a simple lower bound.  Two knot-insertion identities
give

\[
 [b,c,d]g-h(d)
 =b[b,d]h+c[b,c,d]f.                              \tag{15bs}
\]

The sum constraint implies `2b+d<=4`.  The derivative is

\[
 h'(u)=\frac4{u^2}\{e^{-4/u}-(1-u^{-1})^3\},
\]

and the logarithmic derivative of the ratio inside braces is
`(lambda-1)/(4-lambda)`.  Since `e^2<8`, the function `h` increases through
two and thereafter has at most one turn, necessarily a maximum, on `[2,4]`.
If `d<=2`, this gives `h(d)>=h(b)`
directly.  If `d>2`, then `b<1`.  Moreover, rational Taylor bounds give

\[
 e^3\left(1-\frac{81e}{256}\right)-1
 >\frac{213}{2048}>0,
\]

which is exactly `h(4)>h(1)` after multiplication by `e^4`.  Thus the
one-turn shape gives
`h(d)>=min\{h(2),h(4)\}>h(1)>=h(b)`.  Therefore both terms on the right of
(15bs) are nonnegative, and

\[
 [b,c,d]g\ge h(d).                                \tag{15bt}
\]

On the other hand, (15bp) and Hermite--Genocchi give

\[
 [a,b,c]g\le\frac{a+b+c}{12}\le\frac{4-d}{12}.
                                                               \tag{15bu}
\]

The scalar test below uses the value of `f` at a lower bound for `d`.  This
is legitimate in the branch where that lower bound is below three.  Indeed,
the one-maximum analysis of `f'` in Section 11, together with the first
bound below, shows that `f` increases through three.  On `[3,4]`, equations
(15ba)--(15bb) and `e>64/27` give `f''<0`, so `f'` can cross zero at most
once and `f` has no interior minimum.  Rational Taylor bounds give

\[
 e^{4/3}<\frac{59}{15}<\frac{63}{16},
 \qquad
 e>\frac52>\frac{64}{27},
 \qquad
 f(4)-f(3)>
 \frac{4159877}{1167880896}>0.
\]

Consequently `f(d)>=f(d_0)` whenever `1<d_0<3` and `d_0<=d<=4`.
Using the multiplication identity (15ag1), equations (15bt)--(15bu) now
prove the target whenever

\[
 f(d)\ge\frac{a(4-d)}{12}.                         \tag{15bv}
\]

It remains to certify the boxes not covered by `d<=19/15`, `d>=3`,
`a>=9/16`, or (15bv).  Parameterize them by

\[
 (a,b,c,d)=(a,a+r,a+r+s,a+r+s+t),qquad
 4a+3r+2s+t\le4.                                  \tag{15bw}
\]

For `u>1`, the mass-difference formula gives `|q'(u)|<4` directly.  For
`u<=1`, it reduces to
`e^{-lambda}lambda^5/24`, whose maximum on `lambda>=4` occurs at five and is
also below four.  Thus `q` is globally four-Lipschitz.  Coupling the
Hermite--Genocchi expectations shows that changing one knot by `delta`
changes (15bk) by at most `|delta|`.  A parameter box therefore has
Lipschitz error

\[
 4\Delta a+3\Delta r+2\Delta s+\Delta t.          \tag{15bx}
\]

The main initial box is
`[0,1] x [0,4/3] x [0,2] x [0,4]`.  The same tightening, center-scaling,
stack, and tie-breaking rules apply, now with weights `(4,3,2,1)`.  This
deterministic subdivision closes after 77,401 branch calls.
Its terminal boxes comprise 6,190 convex-core boxes, 23 central boxes, 7,599
boxes certified by (15bv), and 24,889 boxes certified directly by (15bx).
The maximum depth is 28, and the ordered terminal-transcript digest is

```text
66dfba78573895b27c4ba6dd7616b68ed96a377f0c767107d30a5f08de8fbf94
```

At every directly certified box, the Arb lower bound at the rational center
strictly exceeds the exact Lipschitz radius.  The source regenerates the
partition, all scalar checks, both digests, and every count.  Thus
`[a,b,c,d]F>=0` on the complete ordered face.  Equations (15bj) and (5),
followed by the radial boundary lemma, prove (15bi).

This theorem completes the Dirichlet--Poissonization comparison for `n=4`.
It does not prove the separate exponential-quantile convexity target and is
not, by itself, a general-`n` Stringer coverage theorem.


## 13. The complete four-positive face for `n=5`

The same expectation-and-subdivision method closes the next complete
coordinate face.

> **`n=5` four-positive-face theorem.**  Suppose the coefficient profile has
> two zero coordinates and four nonnegative coordinates.  If its coordinates
> sum to at most five, then
>
> \[
>  \Pr\{Z_y>5\}\ge\Pr\{T_y>1\}.                  \tag{15by}
> \]

After ordering, write the positive face as `(0,0,a,b,c,d)`, where

\[
 0\le a\le b\le c\le d,
 \qquad a+b+c+d\le5.
\]

Put, with its continuous value at zero,

\[
 F(u)=u^3\{e^{-5/u}-(1-u^{-1})_+^5\}.
\]

Repeated multiplication by the two zero knots reduces the target to
`[a,b,c,d]F>=0`.  If `U` is the uniform-simplex convex combination of these
four knots, then

\[
 [a,b,c,d]F=E q(U),
 \qquad EU=\frac{a+b+c+d}{4}\le\frac54,          \tag{15bz}
\]

where, for `lambda=5/u`,

\[
 q(u)=\frac{F'''(u)}6
 =\begin{cases}
 e^{-\lambda}(1+\lambda+\lambda^2/2+\lambda^3/6),&0<u\le1,\\
 e^{-\lambda}(1+\lambda+\lambda^2/2+\lambda^3/6)
 -1+5u^{-4}-4u^{-5},&u>1.
 \end{cases}                                      \tag{15ca}
\]

The two pieces make `F` a `C^4` function on `[0,5]`, so `q` is continuously
differentiable and confluent divided differences follow by continuity.

The sparse-core theorem proves the face when `d<=25/18`, and the
four-positive far-cap theorem proves it when `d>=4`.  A third analytic
region follows from the affine minorant

\[
 q(u)\ge\frac2{25}\left(\frac54-u\right),
 \qquad \frac{13}{20}\le u\le\frac{61}{20}.     \tag{15cb}
\]

The two pieces of (15ca) are evaluated with directed 160-bit Arb arithmetic.
Dyadic bisection proves (15cb) with seven terminal intervals and maximum
depth five on `[13/20,1]`, and with 256 terminal intervals and maximum depth
ten on `[1,61/20]`.  These are strict interval enclosures of the complete
pieces, not point samples.  Their ordered transcript, combined with the
derivative transcript below, has SHA-256 digest

```text
b5cc3d8cbda175722021249a795d275b0f29adbb6b76c658dada7a3903b2db2e
```

If `a>=13/20`, the sum constraint gives `d<=5-3a<=61/20`.  Equations
(15bz)--(15cb) therefore yield

\[
 [a,b,c,d]F
 \ge\frac2{25}\left(\frac54-EU\right)\ge0.       \tag{15cc}
\]

For the remaining boxes, we use the global bound

\[
 |q'(u)|<1,
 \qquad 0\le u\le5.                              \tag{15cd}
\]

On `u<=1`, the derivative is
`e^{-lambda}lambda^5/30`.  It decreases for `lambda>=5`, and its value at
five is below one because the degree-six Taylor lower bound gives
`e^5>16289/144>625/6`.  On `u>1`, direct differentiation gives

\[
 q'(u)=\frac{e^{-\lambda}\lambda^5}{30}
       -\frac{4\lambda^5(5-\lambda)}{3125}.
\]

Seventeen directed intervals certify that this expression lies in `(-1,1)`
for `1<=lambda<=5`, with maximum bisection depth five.  This completes
(15cd), including the matching value at one.

Parameterize the ordered face by

\[
 (a,b,c,d)=(a,a+r,a+r+s,a+r+s+t),
 \qquad 4a+3r+2s+t\le5.                          \tag{15ce}
\]

Because the Hermite--Genocchi simplex now has four vertices, (15cd) makes
the expectation in (15bz) `1/4`-Lipschitz in each knot.  A parameter box
therefore has error at most

\[
 \frac14(4\Delta a+3\Delta r+2\Delta s+\Delta t). \tag{15cf}
\]

The initial box is
`[0,5/4] x [0,5/3] x [0,5/2] x [0,5]`.  Upper endpoints are tightened from
the other lower endpoints and the budget in (15ce).  The midpoint is scaled
radially to the budget boundary when necessary.  As in Section 12, each
`Delta` is the maximum distance from this evaluation point to the two box
endpoints.  Unresolved boxes are bisected in the first coordinate maximizing
`w_i(U_i-L_i)`, for weights `(4,3,2,1)`, and inserted lower then upper on the
last-in, first-out stack.

The Lipschitz conclusion is applied only to budget-feasible points in each
box.  Both a feasible target and the radially scaled evaluation point satisfy
(15ce), so their connecting segment is feasible and every knot on it lies in
`[0,5]`, the certified domain of (15cd).  Using the full box endpoints in
(15cf) merely enlarges the distance bound.

The subdivision terminates after 81,703 branch calls.  Its terminal boxes
comprise 1,719 sparse-core boxes, 350 far-cap boxes, 82 central-minorant
boxes, and 38,701 boxes certified directly by (15cf).  The maximum depth is
28.  Every direct center value is enclosed by directed 160-bit Arb
arithmetic and strictly exceeds its exact rational Lipschitz radius.  The
ordered terminal-transcript digest is

```text
f4bae2dcb30244286c97891596caffa5c8b1580f41aa80c600fba63c5af7d5a3
```

The source regenerates the scalar enclosures, subdivision, counts, depths,
and both digests.  Hence `[a,b,c,d]F>=0` throughout the ordered face, with
confluent cases supplied by continuity, proving (15by).

Together with the lower-dimensional face theorems, this proves the
comparison for every `n=5` profile having at most four nonzero coordinates.
By itself it does not close the five-positive boundary; Section 14 uses this
face theorem recursively to complete the `n=5` comparison.


## 14. The complete comparison for `n=5`

The four-positive-face theorem supplies both a boundary result and the
recursive estimate needed to close the remaining face.

> **Complete `n=5` theorem.**  For every nonnegative six-coordinate vector
> `y` with `sum_i y_i<=5`, one has
>
> \[
>  \Pr\{Z_y>5\}\ge\Pr\{T_y>1\}.                 \tag{15cg}
> \]

By the radial boundary lemma, it is enough to consider an ordered boundary
profile `(0,a,b,c,d,e)` with

\[
 0\le a\le b\le c\le d\le e,
 \qquad a+b+c+d+e\le5.
\]

Put, with its continuous value at zero,

\[
 F(u)=u^4\{e^{-5/u}-(1-u^{-1})_+^5\}.
\]

Multiplication by the zero knot reduces the target to

\[
 [a,b,c,d,e]F\ge0.                               \tag{15ch}
\]

For a uniform-simplex convex combination `U` of these five knots,
Hermite--Genocchi gives `[a,b,c,d,e]F=E q(U)` and `EU<=1`, where, with
`lambda=5/u`,

\[
 q(u)=\frac{F^{(4)}(u)}{24}
 =\begin{cases}
 e^{-\lambda}\displaystyle\sum_{j=0}^4\frac{\lambda^j}{j!},
       &0<u\le1,\\[2mm]
 e^{-\lambda}\displaystyle\sum_{j=0}^4\frac{\lambda^j}{j!}
       -1+u^{-5},&u>1.
 \end{cases}                                      \tag{15ci}
\]

### Scalar regions

The integrand is nonnegative through `u=19/16`.  On `[1,19/16]` it is
decreasing: `q'(u)<0` is equivalent to `e^lambda>625/24`, and the
degree-four Taylor lower bound at `lambda=80/19` is
`5162241/130321>625/24`.  At the right endpoint, positivity is equivalent to

\[
 e^{80/19}<\frac{32694193}{475841}.
\]

The degree-six geometric-tail Taylor upper bound is
`1819400743567/26675014527`, leaving the exact positive margin

\[
 \frac{6372604003876864}{12693065587542207}.
\]

A second scalar region follows from the directed affine minorant

\[
 q(u)\ge\frac{13}{50}(1-u),
 \qquad \frac58\le u\le\frac52.                 \tag{15cj}
\]

At 160-bit Arb precision, complete dyadic interval evaluation proves
(15cj) with nine terminal intervals and maximum depth seven on `[5/8,1]`,
and 520 intervals and maximum depth 14 on `[1,5/2]`.  The ordered scalar
transcript has SHA-256 digest

```text
bdf4f9223239168671e7eb3517d210b90f270b6394c668e6b21703065234fe97
```

Thus (15ch) holds if `e<=19/16`.  It also holds if `a>=5/8`, because the sum
constraint then gives `e<=5-4a<=5/2`, and (15cj) yields
`E q(U)>=(13/50)(1-EU)>=0`.

### Recursive pruning from the four-positive face

Let `f_3(u)=F(u)/u`, and define

\[
 A=[a,b,c,d]f_3,
 \qquad B=[b,c,d,e]f_3.
\]

The divided-difference multiplication identity gives

\[
 [a,b,c,d,e]F=\frac{eB-aA}{e-a}.                 \tag{15ck}
\]

Section 13 proves `A,B>=0`.  Its Hermite--Genocchi integrand is at most one:
for `u<=1` it is a Poisson CDF, while for `u>1` its correction satisfies

\[
 5u^{-4}-4u^{-5}\le1
\]

because its derivative is `20(1-u)/u^6` and its value at one is one.
Consequently,

\[
 0\le A\le1.                                     \tag{15cl}
\]

Section 13 also proves that this integrand has derivative of absolute value
less than one on `[0,5]`.  For a cumulative-parameter box with coordinate
radii `(Delta a,Delta r,Delta s,Delta t,Delta v)`, directed center values
`A_0,B_0` therefore give

\[
\begin{aligned}
 A&\le A_0+\Delta a+\tfrac34\Delta r
                +\tfrac12\Delta s+\tfrac14\Delta t,\\
 B&\ge B_0-\Delta a-\Delta r-\tfrac34\Delta s
                -\tfrac12\Delta t-\tfrac14\Delta v.
\end{aligned}
\]

The upper bound for `A` is intersected with (15cl).  If `e_-` and `a_+` are
the box lower bound for `e` and upper bound for `a`, respectively, the box
is certified whenever the resulting bounds give `e_- B_->a_+ A_+`.
Equation (15ck) then proves strict positivity throughout the feasible part
of the box.  This recursive test closes most of the residual region.

### Directed Lipschitz bound for the remaining boxes

The direct certificate uses two derivative bounds for (15ci):

\[
 \operatorname{Lip}_{[0,5]}(q)<\frac{33}{8},
 \qquad
 |q'(u)|<1\quad\left(0<u<1\ \hbox{or}\ 6/5\le u\le5\right),
                                                        \tag{15cm}
\]

where the first statement gives a global Lipschitz bound by continuity at
the knot `u=1`.  On `u<=1`,
`q'(u)=e^{-lambda}lambda^6/120`; its maximum occurs at `lambda=6`, and the
degree-eleven lower bound gives `e^6>1944/5`.  On `u>1`,

\[
 q'(u)=e^{-\lambda}\frac{\lambda^6}{120}
       -\frac{\lambda^6}{3125}.
\]

As a function of `lambda`, its derivative has the sign of
`e^{-lambda}(6-lambda)/120-6/3125`, a strictly decreasing function on
`[1,5]`.  Its minimum on each subinterval is therefore at an endpoint.
The exact bounds

\[
 e^5<\frac{129437747}{870912}<\frac{3125}{21}
\]

and

\[
 e^{25/6}<\frac{55091469462637}{853298675712}
          <\frac{48828125}{755256}
\]

prove the lower sides of the two inequalities in (15cm); the upper sides
follow from the same one-turn argument and the elementary endpoint bounds.

Let `W=(W_0,...,W_4)` be uniform Dirichlet and put
`U=sum_i x_iW_i`, where `(x_0,...,x_4)=(a,b,c,d,e)`.  Along a segment inside
a parameter box, (15cm) gives

\[
 \left|\frac{\partial}{\partial x_i}E q(U)\right|
 \le \frac15+\frac{25}{8}
       E\{W_i\mathbf1_{\{1<U<6/5\}}\}.           \tag{15cn}
\]

The event moments are bounded without sampling.  Suppose the first `m`
knot upper bounds are at most `D<1`, all remaining knot upper bounds are at
most `E`, and let `S` be the aggregate Dirichlet weight on the `h=5-m` high
knots.  Then `U>1` implies

\[
 S>\theta=\frac{1-D}{E-D}.                       \tag{15co}
\]

Here `S` has the `Beta(h,5-h)` law.  With `R=1-theta`, the individual
low- and high-coordinate tail moments are as follows:

\[
\begin{array}{c|c|c}
 h & E[W_{\rm low}\mathbf1_{\{S>\theta\}}]
   & E[W_{\rm high}\mathbf1_{\{S>\theta\}}]\\ \hline
1 & R^5/5 & R^4(5-4R)/5\\
2 & R^4-4R^5/5 & 2R^3-3R^4+6R^5/5\\
3 & 2R^3-3R^4+6R^5/5 & 2R^2-4R^3+3R^4-4R^5/5\\
4 & 2R^2-4R^3+3R^4-4R^5/5
  & R-2R^2+2R^3-R^4+R^5/5.
\end{array}                                      \tag{15cp}
\]

The source derives all eight formulas symbolically.  Bounds from the other
side are obtained similarly: if the last `h` knot lower bounds are at least
`L>6/5`, then `U<6/5` implies `S<(6/5)/L`; subtracting the corresponding
moments in (15cp) from `E W_i=1/5` gives another valid bound.  For each knot,
the certificate takes the smallest directed bound supplied by every
admissible split and substitutes it in (15cn).

Finally parameterize the ordered face by

\[
 (a,b,c,d,e)=(a,a+r,a+r+s,a+r+s+t,a+r+s+t+v),
 \qquad 5a+4r+3s+2t+v\le5.                      \tag{15cq}
\]

If the five knot-coordinate bounds from (15cn) are `L_0,...,L_4`, the
cumulative parameter coefficients are

\[
 L_0+\cdots+L_4,\quad L_1+\cdots+L_4,\quad
 L_2+L_3+L_4,\quad L_3+L_4,\quad L_4.           \tag{15cr}
\]

The initial box is
`[0,1] x [0,5/4] x [0,5/3] x [0,5/2] x [0,5]`.  Coordinate upper endpoints
are tightened using the other lower endpoints and the budget in (15cq).
If the rational midpoint exceeds the budget, it is moved from the lower box
corner along their connecting segment to the budget boundary.  The resulting
evaluation point is both feasible and inside the box.  Whole boxes are first
pruned by the convex core, the affine minorant, and the recursive test.  A
remaining box is certified when the
160-bit Arb lower bound for its center value exceeds the directed radius
from (15cr).  Otherwise, the first parameter maximizing
`(5,4,3,2,1)_i` times its width is bisected.  Lower and upper children are
inserted in that order on a last-in, first-out stack.

The subdivision terminates after 1,669,573 branch calls and maximum depth
38.  Its 834,787 terminal boxes comprise 43,171 convex-core boxes, six
affine-minorant boxes, 105,885 boxes closed recursively with `A<=1`,
558,019 boxes closed with the local bound for `A`, and 127,706 direct
Lipschitz boxes.  No infeasible terminal box occurs after coordinatewise
tightening.  The ordered terminal-transcript digest is

```text
c07fa042526bb48d7a79cf67c80297fcbcd59481ac8b3560645bafd80e23a5f2
```

Every numerical inequality in this proof is a directed Arb enclosure at
160-bit precision; every scalar Taylor comparison and all subdivision
endpoints are exact rationals.  The source regenerates the scalar partition,
the beta-moment identities, the complete branch partition, all counts and
depths, and both transcript digests.  Continuity supplies repeated-knot
cases.  This proves (15ch), and the radial boundary lemma proves (15cg).

Thus the Dirichlet--Poissonization comparison is complete through `n=5`.
Section 15 closes the next dimension. The `n=5` result does not by itself
prove the separate exponential-quantile convexity target or a
general-sample-size Stringer coverage theorem.


## 15. The complete comparison for `n=6`

A nested use of the preceding face argument closes the next complete
simplex dimension.

> **Complete `n=6` theorem.** For every nonnegative seven-coordinate vector
> `y` with `sum_i y_i<=6`, one has
>
> \[
>  \Pr\{Z_y>6\}\ge\Pr\{T_y>1\}.                 \tag{15cs}
> \]

The radial boundary lemma reduces the theorem to profiles having at least
one zero coordinate. Lower-dimensional faces are already covered by the
results above, so three ordered faces remain. They are certified in order,
because each face supplies the recursive estimate for the next one.

### Four positive knots

For `(0,0,0,a,b,c,d)`, three zero-knot multiplications reduce the target to
`[a,b,c,d]F_3>=0`, where

\[
 F_3(u)=u^3\{e^{-6/u}-(1-u^{-1})_+^6\}.
\]

Hermite--Genocchi represents this divided difference as `E q_3(U)`, with
`EU<=(a+b+c+d)/4<=3/2`, and

\[
 q_3(u)=\begin{cases}
 e^{-\lambda}\displaystyle\sum_{j=0}^3\lambda^j/j!,&0<u\le1,\\
 e^{-\lambda}\displaystyle\sum_{j=0}^3\lambda^j/j!
 -1+15u^{-4}-24u^{-5}+10u^{-6},&u>1,
 \end{cases}\qquad \lambda=6/u.                 \tag{15ct}
\]

Directed 160-bit Arb subdivisions certify

\[
\begin{aligned}
 q_3(u)&\ge0 &&(0\le u\le31/16),\\
 q_3(u)&\ge\frac1{20}(3/2-u)
   &&(3/4\le u\le15/4),\\
 |q_3'(u)|&<11/20 &&(0\le u\le6),\\
 q_3(u)&<41/200 &&(0\le u\le6).
\end{aligned}                                      \tag{15cu}
\]

The lower exponential pieces are completed analytically from exact Taylor
bounds; the certificate evaluates every compact upper piece as a directed
interval, not as a grid. The scalar nonnegativity bound and the all-`n` far
cap cover `d<=31/16` and `d>=5`, respectively. (The dimension-free sparse
core separately reaches `d<=12/7`.) If `a>=3/4`, the budget implies
`d<=15/4`, and the affine minorant in (15cu) applies.

On the residual region, use cumulative parameters

\[
 (a,b,c,d)=(a,a+r,a+r+s,a+r+s+t),\qquad
 4a+3r+2s+t\le6.
\]

The derivative bound in (15cu) gives the exact cumulative-parameter radius
`(11/80)(4 Delta a+3 Delta r+2 Delta s+Delta t)`. Coordinatewise budget
tightening, a budget-clipped midpoint, and deterministic longest-weighted-
edge bisection close the complete box. The certificate records 36,203 branch
calls and maximum depth 24. Its 18,102 terminal boxes comprise 334 scalar-
core boxes, 152 far-cap boxes, one central-minorant box, and 17,615 boxes
closed by the directed Lipschitz test. The terminal-transcript digest is

```text
796d71cff3aebbbaf38f6b63814f8a25d03eaa2a0baf701d9b62e23689ec888b
```

### Five positive knots

For `(0,0,a,b,c,d,e)`, put

\[
 F_4(u)=u^4\{e^{-6/u}-(1-u^{-1})_+^6\}.
\]

The target is `[a,b,c,d,e]F_4>=0`, with scalar integrand

\[
 q_4(u)=\begin{cases}
 e^{-\lambda}\displaystyle\sum_{j=0}^4\lambda^j/j!,&0<u\le1,\\
 e^{-\lambda}\displaystyle\sum_{j=0}^4\lambda^j/j!
 -1+6u^{-5}-5u^{-6},&u>1.
 \end{cases}                                      \tag{15cv}
\]

The directed scalar certificate proves

\[
\begin{aligned}
 q_4(u)&\ge0 &&(0\le u\le23/16),\\
 q_4(u)&\ge\frac{13}{100}(6/5-u)
   &&(7/10\le u\le16/5),\\
 |q_4'(u)|&<24/25 &&(0\le u\le6),\\
 q_4(u)&<149/500 &&(0\le u\le6).
\end{aligned}                                      \tag{15cw}
\]

Let `f_3=F_4/u` and set

\[
 A=[a,b,c,d]f_3,\qquad B=[b,c,d,e]f_3.
\]

The four-positive theorem gives `A,B>=0`, and knot insertion gives

\[
 [a,b,c,d,e]F_4=\frac{eB-aA}{e-a}.                \tag{15cx}
\]

The pointwise bound in (15cu) yields `A<=41/200`; its derivative bound gives
directed local radii for `A` and `B`. The recursive test `e_-B_->a_+A_+`
closes almost the entire residual region, and the global derivative bound in
(15cw) closes the remainder directly. With cumulative weights `(5,4,3,2,1)`,
the deterministic subdivision terminates after 725,061 calls and maximum
depth 37. Its 362,531 terminal boxes comprise 21,295 scalar-core boxes, 164
central-minorant boxes, 189,416 global-`A` recursive boxes, 144,720 local-`A`
recursive boxes, and 6,936 direct Lipschitz boxes. The terminal-transcript
digest is

```text
d3442eed3224f9c4dc1fb65624c1d375aa5d223670830948f7a13009ac9039b8
```

### Six positive knots

It remains to consider `(0,a,b,c,d,e,f)`. Put

\[
 F_5(u)=u^5\{e^{-6/u}-(1-u^{-1})_+^6\}.
\]

The target is `[a,b,c,d,e,f]F_5>=0`. Its Hermite--Genocchi integrand is

\[
 q_5(u)=\begin{cases}
 e^{-\lambda}\displaystyle\sum_{j=0}^5\lambda^j/j!,&0<u\le1,\\
 e^{-\lambda}\displaystyle\sum_{j=0}^5\lambda^j/j!-1+u^{-6},&u>1.
 \end{cases}                                      \tag{15cy}
\]

Directed interval evaluation proves

\[
 q_5(u)\ge0\quad(0\le u\le37/32),\qquad
 q_5(u)\ge\frac38(1-u)\quad(27/40\le u\le21/8). \tag{15cz}
\]

These bounds close the largest-knot core and, using `EU<=1`, every box with
`a>=27/40`.

For the residual region, let `f_4=F_5/u` and define

\[
 A=[a,b,c,d,e]f_4,\qquad B=[b,c,d,e,f]f_4.
\]

The five-positive theorem gives `A,B>=0`, (15cw) gives `A<=149/500`, and

\[
 [a,b,c,d,e,f]F_5=\frac{fB-aA}{f-a}.              \tag{15da}
\]

A sharper localized radius makes this recursion effective. In addition to
the global bound in (15cw), exact scalar estimates give the nested envelope

\[
 |q_4'(u)|<\frac1{20}
 +\frac15\mathbf1_{\{2/5\le u\le19/10\}}
 +\frac3{10}\mathbf1_{\{11/20\le u\le8/5\}}
 +\frac{13}{50}\mathbf1_{\{7/10\le u\le1\}}
 +\frac{41}{100}\mathbf1_{\{1\le u\le3/2\}}.    \tag{15db}
\]

Thus the successive bounds are `1/20`, `1/4`, and `11/20` outside the
nested intervals, `81/100` on `[7/10,1]`, and `24/25` on `[1,3/2]` (with an
irrelevant larger value at the shared endpoint in the displayed envelope).
On `u<=1`, one has
`q_4'(u)=e^{-lambda}lambda^6/144`; its derivative in `lambda` is negative
for `lambda>=6`. Exact Taylor lower bounds at
`lambda=15,120/11,60/7,6` prove the four successive lower-piece bounds. On
`u>1`, the certificate evaluates the exact derivative over the three outer
`lambda` intervals. The upper and lower signs require respectively 102 and
13 terminal intervals for the `1/20` bound, three and nine for the `1/4`
bound, and one and two for the `11/20` bound; their maximum depths are nine,
six, two, five, zero, and one. The global directed bound in (15cw) supplies
the upper inner interval.

For five Dirichlet coordinates,

\[
\begin{aligned}
 \left|\frac{\partial}{\partial x_i}E q_4(U)\right|
 \le{}&\frac1{100}
   +\frac15E\{W_i\mathbf1_{\{2/5\le U\le19/10\}}\}\\
 &+\frac3{10}E\{W_i\mathbf1_{\{11/20\le U\le8/5\}}\}
   +\frac{13}{50}E\{W_i\mathbf1_{\{7/10\le U\le1\}}\}\\
 &+\frac{41}{100}E\{W_i\mathbf1_{\{1\le U\le3/2\}}\}.
                                                               \tag{15dc}
\end{aligned}
\]

Every event moment in (15dc) is bounded by the exact `Beta(h,5-h)` tail
formulas already used in Section 14. The source rederives those formulas
symbolically and takes the smallest directed bound supplied by every
admissible low/high knot split. Support-disjoint boxes receive moment zero.
The resulting local upper bound for `A` is intersected with `149/500`, and
the lower bound for `B` is substituted in (15da).

With cumulative parameters and weights `(6,5,4,3,2,1)`, the exact initial
box is first divided into 1,024 fixed rational base boxes. Along the
all-zero corner, the first and most expensive box is refined twice by the
same exact budget-aware 64-way rule. This yields 1,150 fixed roots for
parallel regeneration. Within each root, an unresolved box is split in the
coordinate with the largest directed
contribution to the current lower-bound radius. This recursive test, together
with the two scalar regions in (15cz), certifies the complete six-positive
boundary after 46,317,164 branch calls and maximum depth 45. Its 23,159,157
terminal boxes comprise 809,393 scalar-core boxes, 348 central-minorant
boxes, 10,442,544 boxes closed by the global cap for `A`, and 11,906,872
boxes closed by its local directed bound. No infeasible terminal box occurs.
The fixed-root partition digest is

```text
2c8df24aedb29c4b796687d8f79ebeed0bc85b6c6784682da480deef727fbf36
```

The ordered aggregate terminal-transcript digest is

```text
0c6ff61e6a587c8e00975053ad00848c008ee530c6bfcab961175397ee801431
```

The 21 directed scalar partitions used across the three faces contain 9,012
terminal intervals and have maximum depth 18. Their aggregate transcript
digest is

```text
45911d92da792ce6dd0ce418ed617c0530f033ac98f80a4a5d0ece0c8d909db1
```

Every numerical sign in the three subdivisions is a directed 160-bit Arb
inclusion; all box endpoints, budgets, derivative constants, and Taylor
comparisons are exact rationals. The source regenerates the scalar interval
partition, symbolic beta moments, all terminal counts and depths, and every
transcript digest. Continuity supplies repeated knots. The radial boundary
lemma then proves (15cs).

Thus the Dirichlet--Poissonization comparison is complete through `n=6`.
The first open complete simplex dimension is `n=7`. This theorem does not
prove the separate exponential-quantile convexity target and is not, by
itself, a general-sample-size Stringer coverage theorem.


## 16. Why generic $s$-concave localization is too broad

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

## 17. Remaining proof problem

The cleanest current target is (6).  Section 4 shows that it is enough to
work on coordinate faces, Sections 3--5 establish it for every profile
having at most two distinct coefficient values or at most two nonzero
coefficients, and Sections 6--7 prove the complete target for `n=2` and
`n=3`, while Section 12 proves the complete target for `n=4`, Section 14
proves it for `n=5`, and Section 15 proves it for `n=6`.  Sections 8--11
prove every profile having exactly three nonzero coefficients, and Sections
13 and 15 prove the additional coordinate faces used in the two recursive
complete theorems.
Section 8 also proves a convex core on every coordinate face, including the
region `d<=n^2/{3(n+1)}` for four nonzero knots, and Section 9 proves the
opposite four-positive far cap `d>=n-1`.  In dimensions `n>=7`, profiles
with four or more nonzero knots outside the proved regions remain open.  The
first open complete simplex dimension is now `n=7`.  Three plausible routes
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
It certifies the obstruction, checks rational equal-block instances,
regression-checks the exact algebra used in the two-level, two-positive-knot,
arbitrary-order sparse-core, four-positive far-cap, three-positive subregion,
and complete `n=2`, `n=3`, `n=4`, `n=5`, and `n=6` proofs.  It also supplies the
directed interval components of the complete three-positive-face, `n=4`,
both `n=5` face proofs, and all three nested `n=6` face proofs, including
deterministic terminal-box counts and transcript digests.  The general
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
