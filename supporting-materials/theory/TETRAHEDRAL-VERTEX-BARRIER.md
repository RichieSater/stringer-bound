# The four-coordinate section-centroid barrier

## Status

This note proves, by ordinary reduction plus an exact rational Bernstein
certificate, the first genuinely higher-dimensional case in this program of the
section-centroid vertex barrier.  As a consequence, every cap whose threshold
is a nondecreasing barycentric average of four ordered knots is maximized at a
step vertex.  Applied to the Clopper--Pearson weights, this extends pointwise
Stringer--Gaffke domination at sample size `n=3` from the previously certified
range `0.01 <= alpha <= 0.20` to the analytic range

\[
 0<\alpha\le \alpha_*
 :=\left(\frac{19+\sqrt{21}}{34}\right)^3
   =0.3336852118672717\ldots .
\]

Thus binomial Stringer is conservative at `n=3` for every nominal confidence
at least `66.6314788...%`, and this endpoint is sharp for pointwise
Stringer--Gaffke domination. This does **not** prove the all-sample-size cap
conjecture. The proposed five- and six-coordinate global extensions
remain conjectural: their local bounds do not yet establish the required
active-prefix comparison on repeated-lowest-knot faces. See
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).
That obligation does not affect the four-coordinate proof in this note.

## 1. Section centroids and the global mechanism

Let

\[
 D=(D_0,\ldots,D_{N-1})\sim\operatorname{Dirichlet}(1,\ldots,1),
 \qquad x_0\le\cdots\le x_{N-1},
\]

and put `T=x dot D`.  At a regular level `q`, write

\[
 p=\Pr\{T>q\},\qquad
 \mu_i=\mathbb E(D_i\mid T=q),\qquad
 M_r=\sum_{i<r}\mu_i.                                      \tag{1}
\]

The cap at the `r`-th step vector has value

\[
 V_r(\mu)=I_{M_r}(r,N-r),                                  \tag{2}
\]

where `I` is the regularized incomplete beta function.  The relevant local
statement is

\[
 \mu_0\le\cdots\le\mu_{N-1}
 \quad\Longrightarrow\quad
 p\le\max_{r:\,x_{r-1}<x_r}V_r(\mu).                       \tag{3}
\]

Here is why (3) is a maximum principle rather than an isolated section
inequality.  Let `c_0 <= ... <= c_(N-1)` be positive and sum to one, and set

\[
 \Phi_c(x)=\Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}.
\]

If `f_x` is the density of `T` at `q=c dot x`, differentiation gives

\[
 \frac{\partial\Phi_c}{\partial x_i}
 =f_x(q)(\mu_i-c_i).                                      \tag{4}
\]

On a knot face, collect equal consecutive knots into blocks.  At an interior
stationary point, (4) says that each block's conditional Dirichlet mass equals
the corresponding block sum of the `c_i`.  Symmetry makes the individual
centroid coordinates within that block equal to the average of those `c_i`.
Consecutive averages of a nondecreasing vector are nondecreasing, so (3)
applies.  At every active block boundary it has

\[
 M_r=\sum_{i<r}c_i=:C_r,
\]

and hence the stationary cap is no larger than
`max_r I_(C_r)(r,N-r)`.  If a face has no interior maximum, its maximum lies
on a proper face.  Induction over faces, followed by continuity at singular
sections, proves

\[
 \boxed{\displaystyle
 \Phi_c(x)\le\max_{1\le r<N}I_{C_r}(r,N-r).}                \tag{5}
\]

It remains to establish (3).  This note does so for `N=4`.

## 2. The easy crossing types

Suppose first that only the largest knot lies above `q`.  Define

\[
 a_i=\frac{x_3-q}{x_3-x_i},\qquad i=0,1,2.
\]

The upper cap is a tetrahedron and its section is its base, so

\[
 p=a_0a_1a_2,\qquad
 \mu_i=\frac{a_i}{3}\ (i<3),\qquad
 M_3=\frac{a_0+a_1+a_2}{3}.
\]

AM--GM gives

\[
 p\le M_3^3=I_{M_3}(3,1),                                  \tag{6}
\]

which is one of the active step values.

At the other extreme, suppose only the smallest knot lies below `q`, and put

\[
 b_i=\frac{q-x_0}{x_i-x_0},\qquad i=1,2,3.
\]

Then

\[
 p=1-b_1b_2b_3,\qquad
 \mu_0=1-\frac{b_1+b_2+b_3}{3},\qquad
 \mu_i=\frac{b_i}{3}\ (i>0).                               \tag{7}
\]

The ordered knots make `b_1 >= b_2 >= b_3`, whereas the ordered centroid
requires the reverse inequalities.  Thus all three are a common value `b`,
the upper knots coincide, and

\[
 p=1-b^3=I_{1-b}(1,3)=I_{M_1}(1,3).                         \tag{8}
\]

Only the two-versus-two crossing remains.

## 3. The two-versus-two value barrier

### Theorem 1

Let `N=4`, let `x_0<x_1<q<x_2<x_3`, and suppose the section centroid in (1)
is ordered.  Then

\[
 \boxed{p\le I_{M_2}(2,2)=3M_2^2-2M_2^3.}                  \tag{9}
\]

The same assertion holds by continuity when two lower knots coincide.  When
the two upper knots coincide, (9) follows from the separate boundary
certificate in Section 4.

#### Exact reduction

Translation and positive scaling put the centered knots in the form

\[
 (-1,-u,tv,t),\qquad 0<u,v<1,\quad t>0.                    \tag{10}
\]

For temporary positive variables `A,B,C,D`, the divided-difference formula
for the cap above zero is

\[
 p(A,B,C,D)=
 \frac{C^3}{(C+A)(C+B)(C-D)}
 +\frac{D^3}{(D+A)(D+B)(D-C)}.                             \tag{11}
\]

Its density at zero is obtained by replacing the two cubes by three times
their squares.  Differentiation therefore gives the section centroid as

\[
 \mu=\frac{1}{f(0)}
 \left(-p_A,-p_B,p_C,p_D\right),                            \tag{12}
\]

followed by `(A,B,C,D)=(1,u,tv,t)`.  Exact simplification of the final upper
centroid inequality yields

\[
 \mu_3-\mu_2=\frac{(v-1)R(u,v,t)}{\Delta_0(u,v,t)},
 \qquad \Delta_0>0,                                        \tag{13}
\]

where

\[
\begin{aligned}
 R(u,v,t)={}&t^4v^2(u^2+u+1)
 +2t^3uv(u+1)(v+1)\\
 &+t^2u^2(v^2+4v+1)-u^3.                                  \tag{14}
\end{aligned}
\]

Since `v<1`, orderedness implies `R<=0`.  Dropping positive terms in (14)
also gives

\[
 t^2(v^2+4v+1)\le u<1,                                    \tag{15}
\]

so all three variables lie in the unit cube.

Put `M=mu_0+mu_1`.  The target margin has the exact form

\[
 3M^2-2M^3-p=\frac{t^2P(u,v,t)}{\Delta(u,v,t)},             \tag{16}
\]

with

\[
\begin{aligned}
 \Delta={}&27(t+1)^3(t+u)^3(tv+1)^3(tv+u)^3\\
 &\qquad\times(tuv+tv+uv+u)^3>0,                            \tag{17}
\end{aligned}
\]

and tensor degree `deg(P)=(9,9,11)`.  Expanding the 436-term polynomial `P`
in the paper would obscure the idea.  The next section gives the complete,
smaller proof witness.

## 4. Exact Bernstein-multiplier certificate

Write

\[
 B_{i,d}(z)=\binom di z^i(1-z)^{d-i}.
\]

The checked certificate supplies

\[
 S=\sum_{i=0}^{6}\sum_{j=0}^{7}\sum_{k=0}^{7}
 s_{ijk}B_{i,6}(u)B_{j,7}(v)B_{k,7}(t),                    \tag{18}
\]

where 79 coefficients are positive and the other 369 are zero.  It then
performs the exact degree elevation and verifies

\[
 H:=P+RS
 =\sum_{i=0}^{9}\sum_{j=0}^{9}\sum_{k=0}^{11}
 h_{ijk}B_{i,9}(u)B_{j,9}(v)B_{k,11}(t),                   \tag{19}
\]

where 631 coefficients are positive and the other 569 are zero.  None is
negative.  The smallest positive coefficient of `S` is `561/10640`; that of
`H` is `5/99792`.  Therefore, on `R<=0`,

\[
 P=H-RS\ge0,                                                \tag{20}
\]

which proves (9).  In the open domain the right side is strictly positive.

There is one boundary not obtained by setting `v=1` in (13), because equal
upper knots make `mu_3-mu_2=0` identically.  On this boundary, exact
simplification gives

\[
 \mu_2-\mu_1=\frac{G(u,t)}{\Delta_2(u,t)},\qquad
 \Delta_2>0,                                                \tag{21}
\]

where

\[
 G=-2t^3u-t^3+t^2u^2-4t^2u+3tu^2+3u^2.                    \tag{22}
\]

If `t>=1`, then `t^3>=t`, `t^2>=1`, and

\[
 G\le (u-1)\{(3u+1)t+4u\}\le0.                            \tag{23}
\]

Thus orderedness confines the nondegenerate case to the closed unit square.
After removing the positive factor `t^2(1-u)^2` from (16) at `v=1`, call the
remaining numerator (P_\partial).  The boundary witness is the single
multiplier

\[
 S_\partial=2B_{0,2}(u)B_{2,2}(t)=2(1-u)^2t^2.             \tag{24}
\]

All degree `(4,5)` Bernstein coefficients of

\[
 H_\partial=P_\partial-GS_\partial                         \tag{25}
\]

are nonnegative: 19 are positive and 11 are zero, with smallest positive
coefficient `7/5`.  Since orderedness gives `G>=0`, equations (24)--(25) give
(P_\partial=H_\partial+G S_\partial\ge0).  This completes every limiting
two-versus-two section.

The generator
[`tetrahedral_vertex_barrier.py`](../computations/python/tetrahedral_vertex_barrier.py)
derives (11)--(17) symbolically, converts power coefficients to Bernstein
coefficients by exact binomial identities, multiplies the two Bernstein
expansions exactly, checks every sign with rational arithmetic, and verifies
the three polynomial identities locating the `n=3` weight-order endpoint.
The
machine-readable witness is
[`tetrahedral-vertex-barrier-certificate.json`](../computations/certificates/tetrahedral-vertex-barrier-certificate.json).
The command

```bash
make tetrahedral-vertex-barrier-check
```

regenerates the certificate and compares it byte for byte.  Floating-point
optimization proposed the sparse multiplier but makes no proof decision.

### Corollary 2 (four-coordinate monotone-weight cap theorem)

For every nondecreasing probability vector `c=(c_0,c_1,c_2,c_3)` and every
ordered four-knot vector `x`,

\[
 \Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}
 \le \max_{r=1,2,3}I_{C_r}(r,4-r),
 \qquad C_r=\sum_{i<r}c_i.                                  \tag{26}
\]

Indeed, (6), (8), and Theorem 1 prove the section barrier (3) in every
crossing stratum.  The stationary-face argument of Section 1 proves (26).

## 5. New `n=3` Stringer range

For sample size `n=3`, let

\[
 x=p_3(0),\qquad y=p_3(1),\qquad z=p_3(2).
\]

The ascending-knot Stringer weights are

\[
 c=(1-z,\ z-y,\ y-x,\ x).                                  \tag{27}
\]

At every step vertex the beta--binomial identity gives cap probability
exactly `alpha`.  Corollary 2 therefore proves the cap criterion as soon as
(27) is nondecreasing.  We now locate that range exactly.

Put

\[
 f(s)=1-3s^2+2s^3,\qquad f(y)=\alpha,\qquad
 x=1-\alpha^{1/3},\qquad z=(1-\alpha)^{1/3}.                 \tag{28}
\]

The function `f` is strictly decreasing on `(0,1)`.  The last weight
inequality is `y<=2x`.  If `2x>=1` it is automatic; otherwise (28) reduces it
to

\[
 \alpha-f(2x)=x(-3+15x-17x^2)\ge0.                          \tag{29}
\]

The first root in `(0,1)` is

\[
 x_*=\frac{15-\sqrt{21}}{34},
\]

and `alpha=(1-x)^3`.  Hence (29) holds exactly through the endpoint
(\alpha_*=(1-x_*)^3) stated above.

The middle inequality is `y>=(x+z)/2`.  Write `r=alpha^(1/3)=1-x`; for
`alpha<=1/2`, `z>=r` and `z^3+r^3=1`.  Exact simplification gives

\[
 4\left\{f\left(\frac{x+z}{2}\right)-\alpha\right\}
 =3(z-r)(z^2+r^2-1)\ge0,                                   \tag{30}
\]

because `z^2+r^2>=z^3+r^3=1`.  Monotonicity of `f` proves the claim.  Finally,
the first weight inequality is `y<=2z-1`, and

\[
 \alpha-f(2z-1)=(1-z)(17z^2-19z+5)>0.                      \tag{31}
\]

Here \(\alpha\le\alpha_*<1/2\) gives (z>3/4), while the quadratic in (31) is
positive and increasing on `[3/4,1]`.  Thus all three inequalities in (27)
hold for (0<\alpha\le\alpha_*).  The simplex-cap criterion in
[`ORDERED-SIMPLEX-CAP.md`](ORDERED-SIMPLEX-CAP.md) then transfers (26) to
pointwise Stringer--Gaffke domination and hence to distribution-free
Stringer coverage at sample size three.

The endpoint is sharp for this pointwise comparison.  If
`alpha>alpha_*`, then `x<x_*`, the right side of (29) is strictly negative,
and monotonicity of `f` gives `y>2x`.  Thus the last two entries of (27)
satisfy `c_2>c_3`.  For the ordered knot family `k(h)=(0,0,h,1)`, write

\[
 \Phi(h)=\Pr\{k(h)\mathbin{\cdot}D>c\mathbin{\cdot}k(h)\}.
\]

At `h=1`, the beta--binomial identity gives `Phi(1)=alpha`.  Conditional on
`D_2+D_3=c_2+c_3`, symmetry gives conditional mean
`E(D_2 | D_2+D_3=c_2+c_3)=(c_2+c_3)/2`.  The cap derivative (4) therefore
gives

\[
 \Phi'(1)=\frac{g_\alpha}{2}(c_3-c_2)<0,                 \tag{32}
\]

where `g_alpha>0` is the beta density at `c_2+c_3`.  Hence
`Phi(h)>alpha` for `h<1` sufficiently close to one.  The cap criterion then
gives a sample with binomial Stringer strictly below the Gaffke limit for
every `alpha>alpha_*`.  This does not assert Stringer undercoverage beyond
the endpoint; it proves only that this pointwise comparison route cannot be
extended there.

## 6. Mathematical yield and remaining boundary

A literature search updated on 2026-09-02 located the standard Dirichlet
average/B-spline framework of
[Carlson](https://doi.org/10.1016/0021-9045(91)90006-V), the related
integrated-tail and unique-crossing work of
[Diaconis--Perlman](https://doi.org/10.1214/lnms/1215457557) and
[Yu](https://doi.org/10.1214/17-AAP1304), and the Dirichlet convex-order
results of [Letac--Piccioni](https://doi.org/10.3150/15-BEJ765). Those works
provide nearby distributional or stochastic-order tools, but no searched
source states (3), (9), or (26). This is a qualified account of the searched
literature, not an unqualified priority claim; the detailed search boundary
is recorded in [`NOVELTY.md`](../audit/NOVELTY.md).

The result replaces a dimension-specific subdivision proof at `n=3` by a
weight-agnostic geometric theorem: for four uniform Dirichlet coordinates,
monotone barycentric weights force every ordered-knot cap below a step cap.
The key reusable idea is the section centroid.  At a stationary cap its
coordinates reproduce the threshold weights, so a local value comparison at
the section becomes a global face-induction principle.

The two-versus-two inequality is the first case in which the section is a
quadrilateral rather than a simplex.  Its proof also illustrates a more
targeted use of exact computation than exhaustive cap triangulation: one
centroid-order polynomial is used as a nonnegative constraint multiplier,
leaving a globally Bernstein-positive residual.  This separates the
conceptual reduction from 1,200 finite rational signs and makes the
certificate independently regenerable.

The five-coordinate local identities are recorded in
[`FIVE-COORDINATE-VERTEX-BARRIER.md`](FIVE-COORDINATE-VERTEX-BARRIER.md).
They do not currently prove the corresponding global maximum theorem.
Stationary-face induction needs domination by the largest *active* step cap;
an inactive first-coordinate prefix on a repeated-lowest-knot face does not
suffice. Closing this implication in five and six coordinates is a separate
structural problem; see
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).
