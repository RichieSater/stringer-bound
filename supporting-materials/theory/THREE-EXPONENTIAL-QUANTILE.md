# Exact reduction of the three-exponential convexity problem

> **Status.** This note reduces three-coordinate exponential-quantile
> convexity to one explicit two-variable inequality and proves that
> inequality on the repeated-maximum boundary and the equal-smaller symmetry
> line.  Those results are proved; the off-symmetry strict interior remains
> open.  Nothing in this note is a general-`n` Stringer coverage theorem.

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
equality must occur.  Sections 4 and 5 prove the coordinate axes and the
diagonal; the off-symmetry strict interior remains.  Because `Z`, `A`, and
every entry of `K` are elementary divided differences of the exponential,
(16) is an explicit two-variable analytic inequality rather than an
optimization over distributions.

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
with Section 4, the remaining open domain is
`z,w>0` with `z\ne w`.
