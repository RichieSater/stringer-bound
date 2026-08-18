# The exact convexity threshold for two weighted exponentials

## Status and result

Let `E1,E2` be independent unit exponentials and, for `0<alpha<1`, put

\[
 q_\alpha(a,b)
 =Q_{1-\alpha}(aE_1+bE_2),\qquad a,b\ge0.                 \tag{1}
\]

The following theorem completely resolves the two-coordinate case of the
exponential-quantile convexity problem in
[`ALL-N-POISSON-PROGRAM.md`](ALL-N-POISSON-PROGRAM.md).

> **Theorem.** The map `(a,b) -> q_alpha(a,b)` is convex on the
> nonnegative quadrant if and only if
>
> \[
>  \alpha\le 4e^{-3}.                                  \tag{2}
> \]

Thus the three conventional audit tail probabilities `0.10`, `0.05`, and
`0.01` all lie in the convexity range. The constant in (2) is exact, not a
numerical threshold. This theorem does **not** prove convexity in three or
more coordinates and therefore does not establish an all-sample-size
Stringer coverage theorem.

## 1. Reduction to one variable and necessity

Write

\[
 q(t)=q_\alpha(1,t),\qquad t\ge0.
\]

Positive homogeneity gives

\[
 q_\alpha(a,b)=a q(b/a),\qquad a>0.                    \tag{3}
\]

Consequently, convexity of `q` on the half-line is equivalent, through the
perspective construction and continuous extension to the axes, to convexity
of (1). Symmetry also gives

\[
 q(t)=tq(1/t),qquad t>0.                               \tag{4}
\]

Let `q0=q(1)`. Since `E1+E2` is Erlang with shape two,

\[
 (1+q_0)e^{-q_0}=\alpha.                               \tag{5}
\]

Writing `t=1+r`, the exact tail has the fixed-`x` expansion

\[
 e^{-x}(1+x)+r\frac{x^2e^{-x}}2
 +r^2e^{-x}\left(\frac{x^3}{6}-\frac{x^2}{2}\right)+O(r^3).
\]

Substitution of the implicit quantile path into this expansion gives

\[
 q'(1)=\frac{q_0}{2},
 \qquad
 q''(1)=\frac{q_0(q_0-3)}{12}.                         \tag{6}
\]

The left side of (5) is strictly decreasing for positive `q0`, and it
equals `4e^-3` at `q0=3`. Hence `alpha>4e^-3` makes `q''(1)<0`, which proves
the necessity of (2). The remainder of the note proves sufficiency globally,
not merely near equal weights.

## 2. A curvature criterion

We first record the implicit-quantile identity used below. If
`X_t=X+tY` has a positive smooth density `f_t` at its `(1-alpha)` quantile,
then differentiation of `Pr(X_t<=q(t))=1-alpha` gives

\[
 q'(t)=E(Y\mid X_t=q(t))
\]

and a second differentiation gives

\[
 q''(t)
 =-\frac1{f_t(q(t))}
 \left.\frac{d}{dx}
 \left\{f_t(x)\operatorname{Var}(Y\mid X_t=x)\right\}
 \right|_{x=q(t)}.                                    \tag{7}
\]

For completeness, let
`m_k(t,x)=f_t(x)E(Y^k\mid X_t=x)`.  Differentiating the translated joint
density gives the continuity identities
`\partial_t m_k=-\partial_x m_{k+1}`.  The first derivative of the quantile
equation is therefore `f_t(q)q'-m_1(t,q)=0`.  Differentiating that equality
once more and substituting `q'=m_1/f_t` leaves

\[
 f_t(q)q''
 =-\partial_x\left(m_2-\frac{m_1^2}{f_t}\right)_{x=q},
\]

which is (7).  All densities in the application below are positive and
smooth for `t>0` and positive quantiles, including at `t=1` by continuous
extension.

Fix `0<t<1`, put `delta=1-t`, and, temporarily, condition at
`E1+tE2=x`. The conditional density of `E2` on `[0,x/t]` is proportional
to `exp(-delta y)`. Define

\[
 z=\frac{\delta x}{t},
 \qquad
 A(z)=1-e^{-z}-\frac{z^2}{e^z-1}.                      \tag{8}
\]

The density times the conditional variance in (7) is

\[
 f_t(x)\operatorname{Var}(E_2\mid E_1+tE_2=x)
 =e^{-x}\delta^{-3}A(z).                              \tag{9}
\]

Indeed, after scaling by `delta`, the conditional law is a unit exponential
truncated to `[0,z]`, whose variance is

\[
 1-\frac{z^2e^z}{(e^z-1)^2}.
\]

The factor `A(z)` is positive for `z>0`, because

\[
 A(z)=\frac{(e^z-1)^2-z^2e^z}{e^z(e^z-1)}
 =\frac{4\sinh^2(z/2)-z^2}{e^z-1}>0.                  \tag{10}
\]

Put

\[
 h(z)=\frac{zA'(z)}{A(z)}
 =\frac{z(ze^z-e^z+1)^2}
 {(e^z-1)\{(e^z-1)^2-z^2e^z\}}.                      \tag{11}
\]

Differentiating (9) and using `z=delta x/t` shows that its derivative with
respect to `x` is nonpositive exactly when

\[
 h(z)\le x.                                           \tag{12}
\]

Therefore (7) reduces the desired curvature inequality to

\[
 q(t)\ge h\!\left(\frac{(1-t)q(t)}t\right).           \tag{13}
\]

## 3. The sharp tail comparison

Set

\[
 B(z)=\frac{1-e^{-z}}z.
\]

The exact upper tail of `E1+tE2` is

\[
 \Pr(E_1+tE_2>x)
 =\frac{e^{-x}-te^{-x/t}}{1-t}.                       \tag{14}
\]

When `z=(1-t)x/t`, equivalently `t=x/(x+z)`, (14) becomes

\[
 s_z(x):=e^{-x}\{1+xB(z)\}.                           \tag{15}
\]

For every fixed positive `z`, the function `s_z` is strictly decreasing on
the nonnegative half-line, since `0<B(z)<1` and
`s_z'(x)=-e^{-x}\{1-B(z)+xB(z)\}`. It is therefore enough to prove

\[
 s_z(h(z))\ge4e^{-3}.                                 \tag{16}
\]

Take logarithms and define the comparison margin

\[
 M(z)=3-h(z)+\log\{1+h(z)B(z)\}-\log4.                \tag{17}
\]

The removable endpoint values in (8), (11), and (17) give

\[
 h(0)=3,\qquad M(0)=0.
\]

Direct differentiation and exact factorization yield

\[
 M'(z)=
 \frac{(e^z-1-z)^2(ze^z-e^z+1)^2K(z)e^z}
 {(e^z-1)^2(e^{2z}-1-2ze^z)
  (z^2e^z-e^{2z}+2e^z-1)^2},                         \tag{18}
\]

where

\[
 K(z)=2z^2e^z+ze^{2z}-z-4e^{2z}+8e^z-4.              \tag{19}
\]

Every denominator factor in (18) is positive.  Indeed, (10) shows that
`z^2e^z-(e^z-1)^2<0`, so the last squared factor is nonzero, while

\[
 e^{2z}-1-2ze^z=2e^z(\sinh z-z)>0.                   \tag{20}
\]

Also, `e^z-1-z>0`, and `ze^z-e^z+1>0` because the latter function vanishes
at zero and has derivative `ze^z>0`.  The remaining sign is elementary. The
coefficients of `K` through degree
five vanish. For every integer `m>=2`, the coefficient of `z^m`, multiplied
by `m!`, is

\[
 2m(m-1)+(m-8)2^{m-1}+8.                              \tag{21}
\]

It equals zero for `2<=m<=5`, equals `4` and `28` at `m=6` and `m=7`,
respectively, and is plainly positive for every `m>=8`. Thus `K(z)>0` for
`z>0`, so (18) gives
`M'(z)>0`. Equations (17)--(20) prove (16).

Now suppose `alpha<=4e^-3` and use the actual pair `(q(t),z)` in (15).
Because `s_z(q(t))=alpha`, (16) and strict monotonicity of `s_z` imply
`q(t)>=h(z)`. By (7)--(13), `q''(t)>=0` for `0<t<1`. At `t=1`, equations
(5)--(6) give the same conclusion. Finally, differentiating (4) twice gives

\[
 q''(t)=t^{-3}q''(1/t),
\]

so convexity extends to `t>1`. The perspective relation (3) and continuity
on the axes complete the sufficiency proof.

## 4. Reproduction and remaining boundary

The exact symbolic identities in (6), (11), and (18), together with the
coefficient formula (21), are checked by
[`all_n_poisson_reductions.py`](../computations/python/all_n_poisson_reductions.py):

```sh
make all-n-reduction-check
```

No floating-point sign decision enters the theorem. The next unresolved
case is convexity of the same quantile map in three coordinates. Even a
proof on every two-coordinate face does not by itself imply convexity in
higher dimension; that distinction is retained explicitly in the project
status language.
