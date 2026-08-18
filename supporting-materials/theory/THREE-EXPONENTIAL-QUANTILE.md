# Exact reduction of the three-exponential convexity problem

> **Status.** This note reduces three-coordinate exponential-quantile
> convexity to one explicit two-variable inequality.  The reduction is
> proved; the final inequality is open.  Nothing in this note is a
> general-`n` Stringer coverage theorem.

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
equality must occur.  What remains is a proof that no finite unequal pair
`(z,w)` makes the left side of (16) smaller.  Because `Z`, `A`, and every
entry of `K` are elementary divided differences of the exponential, (16) is
an explicit two-variable analytic inequality rather than an optimization
over distributions.
