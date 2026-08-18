# Exact reduction of the three-exponential convexity problem

> **Status.** This note reduces three-coordinate exponential-quantile
> convexity to one explicit two-variable inequality and proves that
> inequality on the repeated-maximum boundary.  The reduction and boundary
> result are proved; the strict two-variable interior remains open.  Nothing
> in this note is a general-`n` Stringer coverage theorem.

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
equality must occur.  Section 4 rules out every finite unequal pair on the
coordinate axes; the strict interior remains.  Because `Z`, `A`, and every
entry of `K` are elementary divided differences of the exponential, (16) is
an explicit two-variable analytic inequality rather than an optimization
over distributions.

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
corner in (18).  The remaining open domain can therefore be restricted to
`z,w>0`.
