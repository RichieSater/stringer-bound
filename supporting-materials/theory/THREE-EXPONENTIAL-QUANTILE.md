# Exact reduction of the three-exponential convexity problem

> **Status.** This note reduces three-coordinate exponential-quantile
> convexity to one explicit two-variable inequality, gives a boundary-trace
> variational formula for its curvature threshold, and proves the inequality
> on the repeated-maximum boundary and the equal-smaller symmetry line.  It
> also proves strict fixed-sum inward transversality at the former and a
> strictly negative fixed-sum transverse second derivative at the latter,
> reduces the complete
> infinite-gap boundary to the two-exponential theorem, and proves positivity
> in a punctured neighborhood of the only zero-margin corner.  It also proves
> the full regions in which both gaps are at most `2/5` or both gaps are at
> least `13`.  The remaining finite off-symmetry region with larger gap above
> `2/5` and smaller gap below `13` is open.  Nothing in this
> note is a general-`n` Stringer coverage theorem.

Let `E_1,E_2,E_3` be independent unit exponentials and put

\[
 q_\alpha(a,b,c)
 =Q_{1-\alpha}(aE_1+bE_2+cE_3),
 \qquad a,b,c\ge0.                                      \tag{1}
\]

The sharp two-coordinate theorem shows that convexity of (1) in any
dimension can hold only when `alpha<=4e^-3`.  The purpose of this note is to
turn sufficiency at that same threshold in dimension three into a concrete
two-variable sign problem.

## 1. Normalization and the tilted simplex

Positive homogeneity and permutation symmetry let us work locally in a cone
where the third coefficient is largest and normalize it to one.  Write

\[
 X=aE_1+bE_2+E_3,
 \qquad 0<a,b\le1,
 \qquad x=q_\alpha(a,b,1),                            \tag{2}
\]

and introduce the nonnegative gap coordinates

\[
 z=\frac{x(1-a)}a,
 \qquad
 w=\frac{x(1-b)}b.
                                                               \tag{3}
\]

Equivalently, `a=x/(x+z)` and `b=x/(x+w)`.
If `g(a,b)=q_alpha(a,b,1)`, then
`q_alpha(a,b,c)=c g(a/c,b/c)` for `c>0`.  Thus the usual perspective
criterion makes positive semidefiniteness of the two-by-two Hessian of `g`
equivalent to positive semidefiniteness of the full Hessian transverse to
the radial direction.  Choosing a largest coordinate and permuting covers
every point in the positive orthant.

Let

\[
 \Delta=\{(u,v):u\ge0,\ v\ge0,\ u+v\le1\}
\]

and define the partition function

\[
 Z(z,w)=\int_\Delta e^{-zu-wv}\,du\,dv.               \tag{4}
\]

If

\[
 B(r)=\frac{1-e^{-r}}r,
\]

with its removable value `B(0)=1`, then direct integration gives

\[
 Z(z,w)=\frac{B(w)-B(z)}{z-w},                        \tag{5}
\]

where equal arguments are interpreted confluently.

Conditionally on `X=x`, the vector

\[
 R=(U,V)=\left(\frac{aE_1}{x},\frac{bE_2}{x}\right)
\]

has density `Z(z,w)^(-1)e^(-zu-wv)` on `Delta`.  This follows by eliminating
`E_3=x-aE_1-bE_2` and applying the change of variables
`E_1=xu/a`, `E_2=xv/b`.  In particular, if `C` is the covariance matrix of
`R` under this tilted law and

\[
 K(z,w)=Z(z,w)C(z,w)
       =Z(z,w)\nabla^2\log Z(z,w),                       \tag{6}
\]

then `K` is positive definite for finite `z,w`.

## 2. An exact curvature matrix

Let `W=diag(a,b)` and let

\[
 \mathcal D=z\frac{\partial}{\partial z}
             +w\frac{\partial}{\partial w}.           \tag{7}
\]

The density of `X` at `x` and the density-weighted conditional covariance of
`(E_1,E_2)` are

\[
 f_X(x)=\frac{e^{-x}x^2}{ab}Z(z,w),                    \tag{8}
\]

and

\[
 f_X(x)\operatorname{Cov}\{(E_1,E_2)\mid X=x\}
 =\frac{e^{-x}x^4}{ab}W^{-1}K(z,w)W^{-1}.              \tag{9}
\]

When `a,b` are held fixed and `x` varies, (3) gives
`x\,d/dx=\mathcal D` on functions of `(z,w)`.  Applying the
implicit-quantile curvature identity componentwise to (9) therefore yields

\[
 \boxed{
 \nabla^2 q_\alpha(a,b,1)
 =\frac{x}{Z(z,w)}W^{-1}
 \{(x-4)K(z,w)-\mathcal D K(z,w)\}W^{-1}.}             \tag{10}
\]

Thus the full two-by-two Hessian question has no remaining implicit
differentiation.

There is also a useful probabilistic form.  Put

\[
 T=zU+wV,
 \qquad \mu=E_{z,w}R.
\]

Differentiating the raw tilted moments in (6) gives

\[
 \mathcal D K
 =-Z\,E_{z,w}\!\left[T(R-\mu)(R-\mu)^{\mathsf T}\right]. \tag{11}
\]

For a nonzero vector `r`, set

\[
 \mathcal R_{z,w}(r)
 =\frac{E_{z,w}[T\{r^{\mathsf T}(R-\mu)\}^2]}
        {E_{z,w}[\{r^{\mathsf T}(R-\mu)\}^2]},
\qquad
 \rho(z,w)=\inf_{r\ne0}\mathcal R_{z,w}(r).           \tag{12}
\]

Equivalently, `rho` is the smaller generalized eigenvalue of the matrix pair
`(-mathcal D K,K)`.  Equations (10)--(12) show that the Hessian in (10) is
positive semidefinite exactly when

\[
 x\ge h(z,w):=4-\rho(z,w).                             \tag{13}
\]

If `h(z,w)<=0`, this condition is automatic because `x>0`.

There is a second, boundary-only representation of the same threshold that
removes the radial derivative and the third central moments.  For a vector
`c in R^2`, put

\[
 \ell_c(u,v)=c^{\mathsf T}\{(u,v)-\mu\}
\]

and define the boundary Gram matrix

\[
 L(z,w)=\int_0^1
 \left\{\binom{u}{1-u}-\mu\right\}
 \left\{\binom{u}{1-u}-\mu\right\}^{\!\mathsf T}
 e^{-zu-w(1-u)}\,du.                                  \tag{13a}
\]

**Trace identity.**  The curvature threshold is the largest generalized
eigenvalue of the boundary and bulk Gram matrices:

\[
 \boxed{
 h(z,w)=\sup_{c\ne0}
 \frac{c^{\mathsf T}L(z,w)c}{c^{\mathsf T}K(z,w)c}
 =\lambda_{\max}\{L(z,w),K(z,w)\}.}                  \tag{13b}
\]

To prove this, write `T=zu+wv` and
`d_c=c^T mu`.  Direct differentiation gives

\[
 \nabla\!\cdot\!\left[(u,v)\ell_c(u,v)^2e^{-T}\right]
 =\{(4-T)\ell_c^2+2d_c\ell_c\}e^{-T}.                \tag{13c}
\]

The integral of the last linear term vanishes by the definition of `mu`.
The radial flux vanishes on the two coordinate edges of `Delta`; on
`u+v=1`, the product of the outward-normal component and arclength is `du`.
The divergence theorem therefore gives

\[
 \int_0^1\ell_c(u,1-u)^2e^{-zu-w(1-u)}\,du
 =\int_\Delta(4-T)\ell_c(u,v)^2e^{-zu-wv}\,du\,dv.    \tag{13d}
\]

The two quadratic forms on the right are `4c^T Kc` and the unnormalized
numerator in (12).  Taking the supremum proves (13b).  Equivalently,

\[
 L(z,w)=4K(z,w)+\mathcal D K(z,w).                    \tag{13e}
\]

This representation gives a particularly simple route for a rigorous
interior certificate.  To prove `h(z,w)<=x`, it is enough to prove the
two-by-two matrix inequality `xK(z,w)-L(z,w)\succeq0`.  If the tail in
(15) at the same `x` is at least `4e^{-3}`, monotonicity then proves (16).
No such full interior certificate is asserted here.

## 3. The remaining two-variable inequality

The tail in the same coordinates is especially simple.  Define

\[
 A(z,w)=\frac{zB(w)-wB(z)}{z-w}
       =B(z)+zZ(z,w)=B(w)+wZ(z,w).                    \tag{14}
\]

The exact hypoexponential tail, after substituting (3), is

\[
 S_{z,w}(x)
 :=\Pr\!\left\{
 \frac{x}{x+z}E_1+\frac{x}{x+w}E_2+E_3>x
 \right\}
 =e^{-x}\{1+A(z,w)x+Z(z,w)x^2\}.                     \tag{15}
\]

For fixed `z,w`, this is strictly decreasing in positive `x`.  A pathwise
proof avoids any differentiation: after division by `x`, the event in (15)
is

\[
 \frac{E_1}{x+z}+\frac{E_2}{x+w}+\frac{E_3}{x}>1,
\]

and its left side decreases strictly with `x` almost surely.

It follows that the three-coordinate theorem at the sharp candidate
threshold is equivalent to the following explicit inequality:

\[
 \boxed{
 e^{-h(z,w)}
 \{1+A(z,w)h(z,w)+Z(z,w)h(z,w)^2\}
 \ge4e^{-3}}                                         \tag{16}
\]

for every `z,w>=0` with `h(z,w)>0`, with all confluent and infinite-boundary
values taken by continuity.

It is convenient to denote the corresponding logarithmic margin by

\[
 M(z,w)=3-h(z,w)
 +\log\{1+A(z,w)h(z,w)+Z(z,w)h(z,w)^2\}-\log4.       \tag{16a}
\]

Thus (16) is equivalent to `M(z,w)>=0`.

Indeed, for the actual quantile point in (2)--(3), (15) gives
`S_{z,w}(x)=alpha`.  If `alpha<=4e^-3` and (16) holds, strict monotonicity in
`x` implies `x>=h(z,w)`.  Equation (10) then makes the Hessian positive
semidefinite.  Permutation, homogeneity, and continuity cover all cones and
their boundaries.  Conversely, the two-coordinate face already rules out
every `alpha>4e^-3`.

At the equal-weight point `z=w=0`, equations (12)--(16) give

\[
 h(0,0)=4,
 \qquad
 S_{0,0}(4)=13e^{-4}>4e^{-3}.                         \tag{17}
\]

At the sharp coordinate-face corner, one gap tends to zero and the other to
infinity; the limit is the equal-weight two-exponential problem,

\[
 h\longrightarrow3,
 \qquad
 S(h)\longrightarrow4e^{-3}.                         \tag{18}
\]

These endpoint calculations explain both the candidate constant and where
equality must occur.  Sections 4--6 prove the coordinate axes, the diagonal,
and the infinite-gap boundary; the finite off-symmetry interior remains.
Because `Z`, `A`, and every entry of `K` are elementary divided differences
of the exponential, (16) is an explicit two-variable analytic inequality
rather than an optimization over distributions.

There is also a direct proof on a nontrivial neighborhood of the equal-weight
point.

**Small-gap theorem.** Inequality (16) holds strictly whenever

\[
 0\le z,w\le\frac25.                                \tag{18a}
\]

Indeed, `T=zU+wV>=0` in (12), so `rho(z,w)>=0` and hence `h(z,w)<=4`.
The case `h<=0` is automatic.  If `0<h<=4`, monotonicity of (15) gives

\[
 S_{z,w}(h)\ge S_{z,w}(4)
 =e^{-4}\{1+4A(z,w)+16Z(z,w)\}.                    \tag{18b}
\]

Both coefficients in braces are coordinatewise nonincreasing.  This follows
immediately for `Z` from its integral representation.  For `A`, define the
boundary integral

\[
 J(z,w)=\int_0^1e^{-zu-w(1-u)}\,du.
\]

The identity `A(z,w)=B(z)+B(w)-J(z,w)` gives

\[
 \frac{\partial A}{\partial z}
 =-\int_0^1u e^{-zu}\,du
   +\int_0^1u e^{-zu-w(1-u)}\,du\le0,               \tag{18c}
\]

and the derivative in `w` is handled symmetrically.  Therefore, with
`r=2/5` and `t=e^{-2/5}`, condition (18a) implies

\[
 \begin{aligned}
 A(z,w)&\ge A(r,r)=5-6t,\\
 Z(z,w)&\ge Z(r,r)=\frac{25}{4}-\frac{35}{4}t,\\
 1+4A(z,w)+16Z(z,w)&\ge121-164t.                   \tag{18d}
 \end{aligned}
\]

All constants can be compared rationally.  The positive exponential series
and its tail bound give

\[
 e^{2/5}>\sum_{k=0}^4\frac{(2/5)^k}{k!}
 =\frac{2797}{1875}>\frac{1000}{671},
 \qquad
 e<\sum_{k=0}^6\frac1{k!}+\frac1{4410}
 =\frac{31967}{11760}<\frac{2719}{1000}.            \tag{18e}
\]

For the second inequality, use
`k!>=7! 8^(k-7)` for `k>=7` and sum the resulting geometric tail.
Consequently,

\[
 121-164e^{-2/5}
 >121-164\frac{671}{1000}
 =\frac{2739}{250}
 >\frac{2719}{250}>4e.                              \tag{18f}
\]

Equations (18b)--(18f) yield `S_{z,w}(h)>4e^{-3}`, proving the theorem.

## 4. A proved boundary: two equal maximal weights

**Boundary theorem.** Inequality (16) holds whenever `z=0` or `w=0`, and it
is strict at every finite point on those axes.  Consequently, for
`alpha<=4e^-3`, the Hessian criterion (10) is positive semidefinite at every
positive coefficient vector having two equal maximal coordinates.

By symmetry, it is enough to put `z=0` and `w=s>0`.  The corresponding
coefficient vector has two equal maximal weights.  Under the tilted-simplex
law, the marginal density of `V` is proportional to

\[
 (1-v)e^{-sv},\qquad 0\le v\le1,
\]

and, conditionally on `V=v`, the coordinate `U` is uniform on `[0,1-v]`.
The symmetry interchanging `U` and `1-U-V` splits the two transverse
curvature directions into two eigenspaces.

More explicitly, write `R_3=1-U-V`.  The two centered contrasts are
`U-R_3=2U+V-1` and `V-EV`, corresponding to coordinate directions `(2,1)`
and `(0,1)`.  Because `E[U-R_3\mid V]=0`, their cross term vanishes for both
the covariance form in the denominator of (12) and its `T=sV`-weighted
counterpart in the numerator.  They are therefore the two generalized
eigendirections.  Moreover,

\[
 \operatorname{Var}(U-R_3\mid V=v)=\frac{(1-v)^2}{3}.
\]

For the antisymmetric direction, put

\[
 I_j(s)=\int_0^1v^j(1-v)^3e^{-sv}\,dv.
\]

Equation (12) gives

\[
 \mathcal R_-(s)=s\frac{I_1(s)}{I_0(s)}.              \tag{19}
\]

For the other transverse direction, let expectation refer to the normalized
density proportional to `(1-v)e^(-sv)`.  Its generalized eigenvalue is

\[
 \mathcal R_+(s)
 =s\frac{E[V\{V-EV\}^2]}{\operatorname{Var}(V)}.      \tag{20}
\]

The first eigenvalue is always smaller.  Define

\[
 \begin{aligned}
 N(s)&=s^2e^s-2se^s+2e^s-2,\\
 D(s)&=s^3e^s-3s^2e^s+6se^s-6e^s+6,\\
 L(s)&=s^4e^s+2s^2e^{2s}+8s^2e^s+2s^2
       -12se^{2s}+12s\\
     &\hspace{4.5em}+12e^{2s}-24e^s+12.
 \end{aligned}                                       \tag{21}
\]

Also put

\[
 J(s)=\int_0^1(1-v)e^{-sv}\,dv,
 \qquad G(s)=s^2e^sJ(s),
 \qquad H(s)=s^2G(s)^2\operatorname{Var}(V).
\]

All three quantities are positive.  A direct common-denominator calculation
in (19)--(20) gives

\[
 \mathcal R_+(s)-\mathcal R_-(s)
 =\frac{s^2N(s)L(s)e^s}{G(s)D(s)H(s)}.               \tag{22}
\]

Indeed, `D(s)=s^4e^sI_0(s)>0`.  Moreover, `N(0)=0` and
`N'(s)=s^2e^s>0`.  The power-series coefficients of `L` vanish through
degree seven.  For every integer `j>=4`, its coefficient of `s^j`,
multiplied by `j!`, is

\[
 j(j-1)(j-2)(j-3)+8j(j-1)-24
 +2^{j-1}(j^2-13j+24).                               \tag{23}
\]

This is `56`, `504`, and `2664` at `j=8,9,10`, respectively.  For `j>=11`,
both the first polynomial part and `j^2-13j+24` are positive.  Thus
`L(s)>0`, proving
`mathcal R_+(s)>mathcal R_-(s)` for `s>0`.

Consequently, (13) and direct evaluation of (19) give

\[
 h(0,s)=4-\mathcal R_-(s)=\frac{3sN(s)}{D(s)}.         \tag{24}
\]

On this axis, `A(0,s)=1` and `Z(0,s)=\{1-B(s)\}/s`.  Define the logarithmic
margin

\[
 M_0(s)=3-h(0,s)
 +\log\{1+h(0,s)+Z(0,s)h(0,s)^2\}-\log4.             \tag{25}
\]

Exact differentiation simplifies to

\[
 M_0'(s)
 =-\frac{27s^2N(s)^2L(s)}
 {D(s)^4\{1+h(0,s)+Z(0,s)h(0,s)^2\}}<0.              \tag{26}
\]

Every factor in the denominator is positive.  Finally,

\[
 h(0,s)\longrightarrow3,
 \qquad Z(0,s)\longrightarrow0,
 \qquad M_0(s)\longrightarrow0
 \quad(s\longrightarrow\infty).                     \tag{27}
\]

Hence `M_0(s)>0` for every finite `s`, which proves (16) on `z=0` and, by
symmetry, on `w=0`.  Equality is approached only at the coordinate-face
corner in (18).

The same factorization also controls the direction transverse to the axis.
For fixed total gap `s`, set

\[
 \Phi_s(\varepsilon)=M(\varepsilon,s-\varepsilon),
 \qquad 0\le\varepsilon\le s/2.
\]

At `epsilon=0`, the smaller generalized eigenvalue is simple by (22),
and its eigendirection is the contrast `X=2U+V`.  Under the axis law,
`E[X]=1`.  Put `Y=U-V`; along this path the tilted density is proportional
to

\[
 e^{-sV-\varepsilon Y}.
\]

Differentiating the centered Rayleigh quotient for `X` at zero, together
with the path derivatives `A'(0)=-1/2+Z(0,s)` and
`Z'(0)=-\int_\Delta Y e^{-sV}\,du\,dv`, gives

\[
 \boxed{
 \Phi_s'(0)
 =\frac{81s^2N(s)^2L(s)}
 {2D(s)^4\{1+h(0,s)+Z(0,s)h(0,s)^2\}}
 =-\frac32M_0'(s)>0.}                                \tag{27a}
\]

Every factor in the first denominator is positive, and the final equality
uses (26).  Thus each finite repeated-maximum point is not only positive but
a strict one-sided local minimum of the margin under fixed-sum splitting of
the two gaps.  The remaining open domain can therefore be restricted to
`z,w>0`.

## 5. A proved symmetry line: two equal smaller weights

**Symmetry-line theorem.** Inequality (16) also holds whenever `z=w`, and it
is strict at every finite point on that ray.  Consequently, for
`alpha<=4e^-3`, the Hessian criterion (10) is positive semidefinite at every
positive coefficient vector having two equal coordinates no larger than the
third.

Put `z=w=s>0` and let `Y=U+V`.  Under the tilted-simplex law, `Y` has density
proportional to `ye^(-sy)` on `[0,1]`; conditionally on `Y=y`, the coordinate
`U` is uniform on `[0,y]`.  If

\[
 J_j(s)=\int_0^1y^je^{-sy}\,dy,
\]

the antisymmetric contrast `U-V` and the symmetric contrast `Y-EY` are the
two generalized eigendirections in (12).  Their eigenvalues are

\[
 \mathcal R_{\rm a}(s)=s\frac{J_4(s)}{J_3(s)},
 \qquad
 \mathcal R_{\rm s}(s)
 =s\frac{E[Y\{Y-EY\}^2]}{\operatorname{Var}(Y)},       \tag{28}
\]

where the second expectation uses the normalized density proportional to
`ye^(-sy)`.  Clearing the positive moment denominators in (28) shows that
the numerator of `mathcal R_a-mathcal R_s` is

\[
 s^2\{2e^s-2-2s-s^2\}L(s).                           \tag{29}
\]

The first factor is positive by the exponential series, and `L(s)>0` was
proved in Section 4.  Hence `mathcal R_s<mathcal R_a`, so the smaller
generalized eigenvalue is `mathcal R_s`.

For the monotonicity calculation, define

\[
 \begin{aligned}
 C(s)&=s^2+2s-2e^s+2,\\
 P(s)&=se^s+s-2e^s+2,\\
 Q(s)&=s^3e^s-s^2e^s-s^2+4se^s-4s
       -2e^{2s}+4e^s-2.
 \end{aligned}                                       \tag{30}
\]

Here `C(s)<0`, while `P(s)>0` because its coefficient of `s^j` is
`(j-2)/j!` for every `j>=3`.  Direct simplification of (28) gives

\[
 h(s,s)=\frac{s^2P(s)^2}{\{s-e^s+1\}Q(s)}.           \tag{31}
\]

The variance in (28) is

\[
 \operatorname{Var}(Y)
 =-\frac{Q(s)}{s^2\{s-e^s+1\}^2},
\]

so `Q(s)<0`.  Differentiating (31) yields

\[
 \frac{d}{ds}h(s,s)=
 \frac{sC(s)P(s)R(s)}
 {\{s-e^s+1\}^2Q(s)^2},                              \tag{32}
\]

where

\[
 \begin{aligned}
 R(s)={}&s^4e^{2s}-s^4e^s-2s^3e^{2s}-2s^3e^s
 +s^2e^{3s}+9s^2e^{2s}-9s^2e^s-s^2\\
 &-6se^{3s}+6se^{2s}+6se^s-6s
 +4e^{3s}-12e^{2s}+12e^s-4.
 \end{aligned}                                       \tag{33}
\]

The coefficients of `R` vanish through degree eight.  For `j>=4`, the
coefficient of `s^j`, multiplied by `j!`, is

\[
 \begin{aligned}
 r_j={}&\frac{2^j}{16}
 (j^4-10j^3+59j^2-2j-192)
 +\frac{3^j}{9}(j^2-19j+36)\\
 &-j^4+4j^3-14j^2+17j+12.                            \tag{34}
 \end{aligned}
\]

Direct substitution gives `r_j>0` for `9<=j<=16`.  For `j>=17`, the first
polynomial in (34) exceeds `j^4/2`, the second is positive, and the final
polynomial exceeds `-j^4`.  The first term in (34) alone then exceeds
`j^4`.  Thus `R(s)>0`; equation (32) proves that `h(s,s)` is strictly
decreasing.

On the diagonal,

\[
 Z(s,s)=\frac{1-(s+1)e^{-s}}{s^2},
 \qquad
 A(s,s)=\frac{2(1-e^{-s})}{s}-e^{-s}.                \tag{35}
\]

Both functions are strictly decreasing.  For `Z` this follows from its
integral representation `Z=J_1`; for `A=B+sZ`, differentiation under the
integral gives `A'(s)=-sJ_2(s)<0`.  Define

\[
 M_{\rm d}(s)=3-h(s,s)
 +\log\{1+A(s,s)h(s,s)+Z(s,s)h(s,s)^2\}-\log4.       \tag{36}
\]

To locate its minimum, put

\[
 \begin{aligned}
 F(s)={}&(-2s^2+12s-4)e^{4s}
 +(-s^4-4s^3-20s^2-32s+16)e^{3s}\\
 &+(2s^6-2s^5+18s^4+24s^3+48s^2+24s-24)e^{2s}\\
 &+(-s^6-6s^5-21s^4-20s^3-28s^2+16)e^s
 +2s^2-4s-4.
 \end{aligned}                                       \tag{37}
\]

Exact differentiation and cancellation give

\[
 M_{\rm d}'(s)=
 \frac{C(s)^2P(s)^2\{s^2e^s-(e^s-1)^2\}F(s)}
 {\{s-e^s+1\}^2Q(s)^4e^s
  \{1+A(s,s)h(s,s)+Z(s,s)h(s,s)^2\}}.               \tag{38}
\]

Every denominator factor is positive, while
`s^2e^s-(e^s-1)^2<0` because `2sinh(s/2)>s`.  Therefore the sign of
`M_d'` is the opposite of the sign of `F`.

Write `F(s)=sum_(j>=0) f_j s^j/j!`.  Its coefficients vanish through
`j=11`, and

\[
 f_{12}=18480,qquad f_{13}=240240,qquad f_{14}=480480.       \tag{39}
\]

For every `j>=4`, direct expansion gives

\[
 \begin{aligned}
 f_j={}&\frac{4^j}{8}(-j^2+25j-32)\\
 &+\frac{2^j}{32}
 (j^6-17j^5+141j^4-415j^3+866j^2-192j-768)\\
 &-\frac{3^j}{81}(j^4+6j^3+155j^2+702j-1296)\\
 &-j^6+9j^5-46j^4+121j^3-173j^2+90j+16.
 \end{aligned}                                       \tag{40}
\]

Substitution gives `f_j<0` for `15<=j<=23`.  For `j>=24`, the negative
`4^j` term in (40) dominates the only positive exponential term: indeed,
the second polynomial is less than `j^6`, and

\[
 4\,2^j(j^2-25j+32)>j^6
\]

holds at `j=24` and then inductively.  The `3^j` and final polynomial terms
are also negative.  For the induction, `j^2-25j+32` is positive and
increasing from its value `8` at `j=24`, while `((j+1)/j)^6<2`.  Thus every
coefficient from degree 15 onward is negative.  The two auxiliary polynomial
comparisons just used can also be checked by substituting `j=24+k`; their
differences then have positive coefficients in `k`.  It follows directly
that `F(s)/s^14` is strictly decreasing on
`s>0`: its first two terms are positive multiples of `s^{-2}` and `s^{-1}`,
its third is constant, and all remaining terms are negative multiples of
positive powers of `s`.  Hence `F` has at most one positive zero.

The following short rational bounds isolate that zero and bound the margin.
They use alternating-series enclosures for the two exponentials; the checker
regenerates every inequality:

\[
 \begin{gathered}
 \frac{3678794411}{10^{10}}<e^{-1}
 <\frac{3678794412}{10^{10}},\\
 \frac{3011942119}{10^{10}}<e^{-6/5}
 <\frac{3011942120}{10^{10}},                         \tag{41}\\
 F(1)>0,qquad F(6/5)<0,\\
 h(1,1)<\frac{347}{100},qquad
 h(6/5,6/5)>\frac{84}{25},\\
 A(6/5,6/5)>\frac{43}{50},qquad
 Z(6/5,6/5)>\frac{117}{500}.
 \end{gathered}
\]

Thus the unique minimum of `M_d` occurs at some `s_0 in (1,6/5)`.  At that
point, monotonicity and (41) give

\[
 \begin{aligned}
 1+A h+Zh^2
 &>1+\frac{43}{50}\frac{84}{25}
   +\frac{117}{500}\left(\frac{84}{25}\right)^2\\
 &=\frac{510263}{78125}>\frac{32}{5}.                \tag{42}
 \end{aligned}
\]

Finally, for `x=47/100`, the exponential series and
`j!>=24*5^(j-4)` for `j>=4` give

\[
 e^x
 \le1+x+\frac{x^2}{2}+\frac{x^3}{6}
 +\frac{x^4}{24(1-x/5)}
 =\frac{17395178081}{10872000000}<\frac85.           \tag{43}
\]

Equations (36), (41)--(43) now imply

\[
 M_{\rm d}(s_0)
 >-\frac{47}{100}+\log\frac85>0.                    \tag{44}
\]

This proves (16) on the full diagonal, including `s=0` by (17).  Together
with Section 4, this leaves finite points with `z,w>0` and `z\ne w`.
Section 6 checks the remaining boundary at infinity.

The margin also has a strict transverse sign on this symmetry line.

**Diagonal-transversality theorem.** For every `s>0`, set

\[
 \Psi_s(\varepsilon)=M(s-\varepsilon,s+\varepsilon),
 \qquad |\varepsilon|<s.                              \tag{44a}
\]

Then `Psi_s` is even and

\[
 \boxed{\Psi_s''(0)<0.}                              \tag{44b}
\]

Thus every positive finite diagonal point is a strict local maximum of the
margin under fixed-total-gap asymmetry.  This is the complementary local
sign to the strict inward minimum at the coordinate axes in (27a).

Here is a direct proof using the trace identity.  Put `t=e^{-s}` and

\[
 j_r=j_r(s)=\int_0^1y^re^{-sy}\,dy.
\]

Along (44a), use the coordinates `Y=U+V` and `D=U-V`.  The tilted density is
proportional to `e^{-sY+\varepsilon D}`.  In the `(Y,D)` basis, direct
integration gives the expansions

\[
 \begin{aligned}
 K_\varepsilon
 &=\begin{pmatrix}
 k_0+k_2\varepsilon^2+O(\varepsilon^4)&
 r_1\varepsilon+O(\varepsilon^3)\\
 r_1\varepsilon+O(\varepsilon^3)&r_0+O(\varepsilon^2)
 \end{pmatrix},\\
 L_\varepsilon
 &=\begin{pmatrix}
 \ell_0+\ell_2\varepsilon^2+O(\varepsilon^4)&
 m_1\varepsilon+O(\varepsilon^3)\\
 m_1\varepsilon+O(\varepsilon^3)&m_0+O(\varepsilon^2)
 \end{pmatrix},                                      \tag{44c}
 \end{aligned}
\]

where

\[
 \begin{aligned}
 k_0&=j_3-\frac{j_2^2}{j_1},&
 k_2&=\frac{j_5}{6}-\frac{j_2j_4}{3j_1}
       +\frac{j_2^2j_3}{6j_1^2},\\
 r_0&=\frac{j_3}{3},&
 r_1&=\frac{j_4}{3}-\frac{j_2j_3}{3j_1},\\
 a_0&=1-\frac{j_2}{j_1},&
 a_2&=-\frac{j_4}{6j_1}+\frac{j_2j_3}{6j_1^2},\\
 \ell_0&=t a_0^2,&
 \ell_2&=t\left(\frac{a_0^2}{6}+2a_0a_2\right),\\
 m_0&=\frac t3,&
 m_1&=ta_0\left(\frac13-\frac{j_3}{3j_1}\right).
                                                               \tag{44d}
 \end{aligned}
\]

For example, the partition function is
`j_1+\varepsilon^2j_3/6+O(\varepsilon^4)`, the bulk means are

\[
 E[Y]=\frac{j_2}{j_1}
 +\left(\frac{j_4}{6j_1}-\frac{j_2j_3}{6j_1^2}\right)
 \varepsilon^2+O(\varepsilon^4),
 \qquad
 E[D]=\frac{j_3}{3j_1}\varepsilon+O(\varepsilon^3),
\]

and (44c)--(44d) follow by centering the bulk and boundary second moments.

At `epsilon=0`, Section 5 proves that the largest generalized eigenvalue is
simple and lies in the `Y` direction.  Write

\[
 h_0=\frac{\ell_0}{k_0},
 \qquad d_0=m_0-h_0r_0<0.
\]

The standard two-by-two simple-eigenvalue expansion, or direct expansion of
`det(L_epsilon-h K_epsilon)=0`, gives

\[
 h(s-\varepsilon,s+\varepsilon)
 =h_0+h_2\varepsilon^2+O(\varepsilon^4),
 \qquad
 h_2=\frac{\ell_2-h_0k_2-(m_1-h_0r_1)^2/d_0}{k_0}.   \tag{44e}
\]

The other two functions in the tail comparison satisfy

\[
 \begin{aligned}
 Z(s-\varepsilon,s+\varepsilon)
 &=j_1+\frac{j_3}{6}\varepsilon^2+O(\varepsilon^4),\\
 A(s-\varepsilon,s+\varepsilon)
 &=j_0+sj_1+\left(\frac{j_2}{2}+\frac{sj_3}{6}\right)
   \varepsilon^2+O(\varepsilon^4).                  \tag{44f}
 \end{aligned}
\]

Consequently, if

\[
 T_0(s)=1+A(s,s)h_0+Z(s,s)h_0^2>0,
\]

then `Psi_s(epsilon)=Psi_s(0)+c(s)epsilon^2+O(epsilon^4)`, where substitution
of (44e)--(44f) and the exact moments

\[
 j_r=\frac{r!}{s^{r+1}}
 \left(1-e^{-s}\sum_{q=0}^r\frac{s^q}{q!}\right)
\]

simplifies to

\[
 \boxed{
 c(s)=
 \frac{e^{-s}P(s)^2C(s)\mathcal E(s)}
 {6\{s-e^s+1\}^2Q(s)^4T_0(s)}.}                    \tag{44g}
\]

Here `C`, `P`, and `Q` are the functions in (30), and

\[
 \mathcal E(s)=\sum_{k=0}^7 E_k(s)e^{ks},            \tag{44h}
\]

with

\[
 \begin{aligned}
 E_7={}&8(s^2-5s-3),\\
 E_6={}&-4(s^4-8s^3-62s^2-52s-54),\\
 E_5={}&-2(4s^6+19s^5+107s^4+400s^3+504s^2+348s+396),\\
 E_4={}&9s^8+26s^7+142s^6+544s^5+1370s^4+2200s^3
         +1760s^2+1760s+1560,\\
 E_3={}&-2(s^{10}+5s^9+16s^8+93s^7+282s^6+678s^5
         +1090s^4+1200s^3+1100s^2+1420s+900),\\
 E_2={}&s^{10}+12s^9+71s^8+238s^7+520s^6+968s^5
         +1240s^4+1520s^3+2232s^2+2640s+1224,\\
 E_1={}&-2(3s^7+22s^6+71s^5+163s^4+416s^3+704s^2+644s+228),\\
 E_0={}&2(s^6+12s^5+57s^4+140s^3+184s^2+128s+36).
                                                               \tag{44i}
 \end{aligned}
\]

It remains only to prove that `mathcal E(s)>0`.  The following coefficient
argument does so without numerical approximation.  Write

\[
 \mathcal E(s)=\sum_{n\ge0}e_n\frac{s^n}{n!},
 \qquad
 p_k(n)=\sum_r [s^r]E_k(s)\frac{(n)_r}{k^r},         \tag{44j}
\]

where `(n)_r=n(n-1)\cdots(n-r+1)`.  Since `E_0` has degree six, for
`n>=18`

\[
 e_n=\sum_{k=1}^7 k^n p_k(n).                       \tag{44k}
\]

Exact expansion gives `e_0=\cdots=e_{17}=0` and

\[
 \min_{18\le n\le40}e_n=e_{18}=61751289600>0.       \tag{44l}
\]

For the infinite tail of coefficients, the sign bookkeeping is short.
The shifted polynomials `-p_k(18+m)` have positive coefficients in `m` for
`k=1,3,5`, while `p_k(18+m)` has positive coefficients for `k=2,4`.
Moreover,

\[
 p_7(n)=\frac8{49}(n^2-36n-147)>0\quad(n\ge40),      \tag{44m}
\]

and

\[
 p_6(n)=\frac{-n^4+54n^3+2077n^2+9102n+69984}{324}.
                                                               \tag{44n}
\]

The last polynomial is positive through `n=81`: for `n<=54` this is
immediate; for `55<=n<=79`, group its numerator as
`n^2(-n^2+54n+2077)+9102n+69984` and note that the quadratic factor is at
least its value `102` at `n=79`; finally,
`p_6(80)=64912/27` and `p_6(81)=264`.  The polynomial
`-p_6(82+m)` has positive coefficients, so `p_6(n)<0` for `n>=82`.

For `k=1,3,5`, put `q_k=-p_k` and

\[
 R_k(n)=\frac{q_k(n)k^n}{p_7(n)7^n}.
\]

Each `R_k` decreases for `n>=41`, because the polynomial

\[
 7q_k(n)p_7(n+1)-kq_k(n+1)p_7(n)                  \tag{44o}
\]

has positive coefficients after substituting `n=41+m`.  Direct exact
substitution gives

\[
 R_1(41)+R_3(41)+R_5(41)<\frac12.                  \tag{44p}
\]

For `q_6=-p_6`, the analogous polynomial in (44o) has positive coefficients
after `n=88+m`.  Exact substitution for `82<=n<=88` gives

\[
 0<R_6(n)<10^{-4},
\]

so this bound persists for every `n>=82`.  Thus, from `n=41` onward, the
total magnitude of every negative term in (44k) is less than
`(1/2+10^{-4})p_7(n)7^n`.  The omitted `k=2,4` terms, and the `k=6` term
through `n=81`, are positive.  Together with (44l), this proves
`e_n>0` for every `n>=18`, and hence `mathcal E(s)>0` for every `s>0`.

Section 5 already proved `P(s)>0`, `C(s)<0`, and `Q(s)<0`.  Every other
factor in the denominator of (44g) is positive.  Therefore `c(s)<0`, so
`Psi_s''(0)=2c(s)<0`, proving (44b).  This transverse result is local: it
does not assert monotonicity across the finite off-symmetry interior.

## 6. The infinite-gap boundary

The other boundary of the compactified gap domain also reduces exactly to
the proved two-exponential theorem.  Fix `z>=0` and let `w` tend to infinity.
After the rescaling `Y=wV`, the tilted-simplex law converges to the product
of

\[
 U\ hbox{with density}\ \frac{e^{-zu}}{B(z)},
 \quad 0\le u\le1,
 \qquad\hbox{and}\qquad
 Y\sim\operatorname{Exp}(1).                         \tag{45}
\]

Indeed, the rescaled domain is
`{(u,y):0<=u<=1, 0<=y<=w(1-u)}`, its density is proportional to
`e^(-zu-y)`, and dominated convergence applies to every moment used in
(12).  Generalized eigenvalues are unchanged by the coordinate rescaling
from `(U,V)` to `(U,Y)`.

Let `mu_z=EU` and `sigma_z^2=Var(U)` under (45).  For
`T=zU+Y`, independence makes both limiting covariance matrices diagonal.
Their generalized eigenvalues are

\[
 \rho_U(z)=1+z\frac{E[U\{U-\mu_z\}^2]}{\sigma_z^2},
 \qquad
 \rho_Y(z)=3+z\mu_z.                                \tag{46}
\]

The first expression is `1` plus the tilted-interval Rayleigh quotient in
the two-exponential proof.  Consequently,

\[
 \rho_U(z)=4-h_2(z),                                 \tag{47}
\]

where `h_2` is the function in equation (11) of
[`TWO-EXPONENTIAL-QUANTILE.md`](TWO-EXPONENTIAL-QUANTILE.md).  Also,

\[
 \mu_z=\frac1z-\frac1{e^z-1}
\]

by continuous extension at zero.  The other eigenvalue is no smaller,
because

\[
 \begin{aligned}
 \rho_Y(z)-\rho_U(z)
 &=h_2(z)-\frac{z}{e^z-1}\\
 &=\frac{-z^2e^z\{ze^z+z-2e^z+2\}}
 {(e^z-1)\{z^2e^z-(e^z-1)^2\}}>0                  \tag{48}
 \end{aligned}
\]

for finite `z>0`.  The numerator's braced factor has positive power-series
coefficients from degree three onward, and the last denominator factor is
negative by (10) of the two-exponential note.  Thus

\[
 h(z,w)\longrightarrow h_2(z).                       \tag{49}
\]

Finally, (5) and (14) give

\[
 wZ(z,w)\longrightarrow B(z),
 \qquad A(z,w)\longrightarrow B(z),                 \tag{50}
\]

and therefore

\[
 S_{z,w}\{h(z,w)\}
 \longrightarrow e^{-h_2(z)}\{1+B(z)h_2(z)\}
 \ge4e^{-3}.                                        \tag{51}
\]

The inequality in (51) is exactly the sharp comparison (16) in the
two-exponential proof.  It is strict for `z>0` and becomes equality only at
`z=0`, the sharp corner already identified in (18).

For completeness, the case in which both gaps diverge is uniformly benign.
Set `X=zU` and `Y=wV`.  If `min(z,w)` tends to infinity, the rescaled domain

\[
 \{(x,y):x\ge0,\ y\ge0,\ x/z+y/w\le1\}
\]

exhausts the nonnegative quadrant, and the tilted law converges in all
moments to two independent unit exponentials.  In these coordinates,
write `R_hat=(X,Y)`.  Then

\[
 \operatorname{Cov}(\widehat R)\longrightarrow I_2,
 \qquad E[(X+Y)(\widehat R-E\widehat R)
                 (\widehat R-E\widehat R)^{\mathsf T}]
 \longrightarrow4I_2.                               \tag{52}
\]

Thus `rho(z,w)->4`, `h(z,w)->0`, `A(z,w)->0`, and `Z(z,w)->0`; the left side
of (16) tends to one.

There is also a useful nonasymptotic version of this two-large-gap
conclusion.

**Two-large-gap theorem.** Inequality (16) holds throughout

\[
 \min(z,w)\ge13.                                     \tag{52a}
\]

To prove this, assume by symmetry that `1<=z<=w` and put `delta=w-z`.
The tilted-simplex law is the law of independent exponential variables of
rates `z,w`, conditioned on their sum being at most one.  This decreasing
conditioning event gives

\[
 0\le\mu_1\le\frac1z,
 \qquad 0\le\mu_2\le\frac1w.                        \tag{52b}
\]

For completeness, condition first on the first exponential.  The
conditional probability of the event is a decreasing function of that
coordinate, so its covariance with the coordinate is nonpositive; division
by the event probability proves the first bound.  The second is symmetric.

On the boundary `u+v=1`, set `r=1-u`.  Its density is
`e^{-z}e^{-\delta r}`.  For `c=(c_1,c_2)`, the elementary inequality
`(x+y)^2<=2x^2+2y^2`, (52b), and

\[
 \begin{aligned}
 \int_0^1e^{-\delta r}\,dr
 &\le\frac2{1+\delta},\\
 \int_0^1r^2e^{-\delta r}\,dr
 &\le\frac4{(1+\delta)^3}.                          \tag{52c}
 \end{aligned}
\]

The first bound follows by combining the elementary bounds `1` and
`1/delta`.  For the second, multiply the desired inequality by
`delta^3(1+delta)^3e^delta`.  The resulting exponential polynomial is

\[
 H(\delta)=\delta^5+5\delta^4+11\delta^3+13\delta^2+8\delta+2
 +(2\delta^3-6\delta^2-6\delta-2)e^\delta.
\]

Its coefficients vanish through degree two; after multiplication by `n!`,
the coefficients in degrees `3,4,5` are `22,70,88`, and for `n>=6` they are
`2n^3-12n^2+4n-2=2n^2(n-6)+4n-2>0`.  Thus `H(delta)>0`.

give

\[
 \begin{aligned}
 c^{\mathsf T}Lc
 \le{}&\frac{4e^{-z}(1+1/z)^2}{1+\delta}\,c_1^2\\
 &+e^{-z}\left\{
 \frac{16}{(1+\delta)^3}
 +\frac{8}{w^2(1+\delta)}\right\}c_2^2.            \tag{52d}
 \end{aligned}
\]

These estimates use only `(r+1/w)^2<=2r^2+2/w^2` in addition to (52c).

For a sharper bulk lower bound, extend the integral temporarily to the full
nonnegative quadrant.  Set

\[
 r_1=\frac{c_1}{z},\qquad r_2=\frac{c_2}{w},
 \qquad d=c^{\mathsf T}\mu.
\]

After `x=zu` and `y=wv`, direct integration over the quadrant gives

\[
 \int_{\mathbb R_+^2}\{c^{\mathsf T}((u,v)-\mu)\}^2
 e^{-zu-wv}\,du\,dv
 =\frac{r_1^2+r_2^2+(r_1+r_2-d)^2}{zw}
 \ge\frac{r_1^2+r_2^2}{zw}.                         \tag{52e}
\]

The complement of `Delta` is contained, after this scaling, in
`{x+y>z}`.  By (52b), `|d|<=|r_1|+|r_2|`.  Hence

\[
 (r_1x+r_2y-d)^2
 \le2(r_1^2+r_2^2)(x^2+y^2+2).
\]

Writing `s=x+y`, integrating first along the segment `x+y=s`, and using
`x^2+y^2<=s^2`, the complement contributes at most

\[
 \frac{2(r_1^2+r_2^2)}{zw}
 \int_z^\infty(s^3+2s)e^{-s}\,ds
 =\frac{\theta(z)(r_1^2+r_2^2)}{zw},                \tag{52f}
\]

where

\[
 \theta(z)=2e^{-z}(z^3+3z^2+8z+8).
\]

Subtracting (52f) from (52e) therefore yields

\[
 c^{\mathsf T}Kc
 \ge\frac{1-\theta(z)}{zw}
 \left(\frac{c_1^2}{z^2}+\frac{c_2^2}{w^2}\right). \tag{52g}
\]

Since `w/(1+delta)<=z` for `z>=1`, comparing (52d) and (52g) gives

\[
 h(z,w)=\sup_{c\ne0}\frac{c^{\mathsf T}Lc}{c^{\mathsf T}Kc}
 \le\frac{24z^4e^{-z}}{1-\theta(z)}.                \tag{52h}
\]

Both `theta(z)` and `z^4e^{-z}` decrease for `z>=4`; for the first function,
differentiate its polynomial factor to obtain
`e^z\theta'(z)/2=-z^3-2z<0`.  The exponential series gives
`e>163/60>19/7`.  Since the polynomial in `theta` equals `2816` at `z=13`,
the exact integer inequalities

\[
 163^{13}>100\cdot2816\cdot60^{13},
 \qquad
 49\,163^{13}>750\cdot13^4\,60^{13}                \tag{52i}
\]

imply `theta(13)<1/50` and make the right side of (52h) less than `8/5`
at `z=13`.  It remains below `8/5` thereafter.  Thus `h(z,w)<8/5` under
(52a).  If `h<=0`, the curvature condition is automatic.  Otherwise,

\[
 e^{-h}\{1+Ah+Zh^2\}\ge e^{-8/5}>4e^{-3},           \tag{52j}
\]

The last inequality is equivalent to `e^7>4^5`; it follows from
`e>19/7` and the exact inequality `19^7>4^5 7^7`.  This proves the theorem.

Hence every finite point on the coordinate axes, every point on the
diagonal, the complete infinite-gap boundary, the small-gap square (18a),
and the two-large-gap region (52a) satisfy the target.  What remains is the
finite off-symmetry domain
`0<z<w<infinity` with `z<13` and `w>2/5`, up to interchange of `z` and `w`.

## 7. A positive neighborhood of the sharp corner

The zero margin at `(z,w)=(0,infinity)` cannot conceal a sequence of
violations approaching that corner.  Put `t=1/w` and rescale the second
simplex coordinate to `Y=wV`.  Since

\[
 B(w)=t\{1-e^{-1/t}\},
\]

the exponentially small term is flat at `t=0`.  After dropping only that
flat remainder, the partition function is

\[
 Z_0(z,t)=\frac{t\{B(z)-t\}}{1-zt}.                  \tag{53}
\]

Use the congruence transformation from `(U,V)` to `(U,Y)` in both matrices
of the generalized-eigenvalue problem.  The rescaled matrices extend
smoothly to `(z,t)=(0,0)`, where the limiting covariance and
weighted-covariance matrices are

\[
 C_0=\begin{pmatrix}1/12&0\\0&1\end{pmatrix},
 \qquad
 H_0=\begin{pmatrix}1/12&0\\0&3\end{pmatrix}.
\]

The smaller generalized eigenvalue is therefore simple.  Taylor expansion
of that eigenvalue, using (53), gives

\[
 \begin{aligned}
 \rho(z,1/t)={}&1+\frac z2-3t+\frac{z^2}{60}+3t^2
 -\frac{z^2t}{60}-\frac32zt^2+9t^3\\
 &+O\{(z+t)^4\}.                                    \tag{54}
 \end{aligned}
\]

The flat `e^(-1/t)` contribution and all of its polynomially rescaled
derivatives are absorbed by the remainder.  Substituting `h=4-rho` together
with the corresponding expansions of `A` and `Z` into the logarithmic
margin gives the especially simple expression

\[
 M(z,1/t)
 =\frac{3}{40}z^2-\frac{3}{20}z^2t+\frac92t^3
  +O\{(z+t)^4\}.                                    \tag{55}
\]

This expansion proves a genuine two-variable neighborhood result.  For
`0<=z,t<=delta<=1/4`, the displayed part of (55) is at least

\[
 \frac{3}{80}z^2+\frac92t^3.
\]

If the absolute remainder is at most `C(z+t)^4`, then

\[
 (z+t)^4\le8(z^4+t^4)
 \le8\delta(z^2+t^3).
\]

Choosing `delta>0` sufficiently small makes the remainder less than half of
the displayed positive lower bound.  Hence `M(z,1/t)>0` throughout a
punctured neighborhood of the sharp corner.  In particular, any sequence of
failures of (16) would have to remain a positive distance from every
compactified boundary family already treated; the unresolved task is an
interior sign problem, not a singular-boundary problem.
