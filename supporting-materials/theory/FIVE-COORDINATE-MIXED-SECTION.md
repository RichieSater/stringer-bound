# The five-coordinate distinct-knot section barrier

## Status and mathematical yield

This note proves the complete **distinct-knot** local barrier for five uniform
Dirichlet coordinates.  Every regular section with five distinct ordered
knots and an ordered section centroid is no larger than an active step cap.
For each of the two mixed crossing types, the stronger conclusion holds that
the section cap is no larger than the first step cap.

The proof is ordinary geometry plus two exact rational polynomial
certificates.  In the two-below/three-above case, lexicographic division by an
endpoint-order polynomial directly produces a Bernstein-positive quotient and
remainder.  In the complementary three-below/two-above case, a sparse
363-by-363 exact rational system produces a Bernstein-positive multiplier,
and all 14,700 residual coefficients are checked.  Neither certificate uses
floating point.

This note isolates the open distinct-knot stratum.  The companion
[`FIVE-COORDINATE-VERTEX-BARRIER.md`](FIVE-COORDINATE-VERTEX-BARRIER.md)
classifies repeated-top families and certifies first-prefix inequalities.
Their composition into global vertex maximization remains conjectural on
repeated-lowest faces, as explained in [`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md). Keeping the two layers separate makes clear which identities
are generic and which arise only when the terminal centroid comparison
degenerates on a repeated-top face.

The reusable idea is a certificate architecture: use the final centroid-order
polynomial both to compactify the scale and as the nonnegative constraint in
an exact Bernstein identity.  The first mixed case needs no fitted multiplier
at all; polynomial division chooses it canonically.  The second shows how a
small index witness can determine a much larger exact rational certificate
without storing unchecked numerical coefficients.

## 1. Section centroids and the local target

Let

\[
 D=(D_0,\ldots,D_4)\sim\operatorname{Dirichlet}(1,1,1,1,1),
 \qquad x_0<\cdots<x_4,
\]

and let `q` be a regular cutting level.  Write

\[
 p=\Pr\{x\mathbin{\cdot}D>q\},\qquad
 \gamma_i=\mathbb E(D_i\mid x\mathbin{\cdot}D=q),\qquad
 M_r=\sum_{i<r}\gamma_i.                                    \tag{1}
\]

The step cap at the strict gap `x_(r-1)<x_r` has value

\[
 V_r(\gamma)=I_{M_r}(r,5-r).                                \tag{2}
\]

The local maximum barrier asks for

\[
 p\le \max_{r=1,\ldots,4}V_r(\gamma)                        \tag{3}
\]

whenever `gamma_0 <= ... <= gamma_4`.  The stationary-face argument in
[`TETRAHEDRAL-VERTEX-BARRIER.md`](TETRAHEDRAL-VERTEX-BARRIER.md) explains why
(3), including its block-face analogues, would imply a global vertex theorem.
This note establishes (3) whenever all five knots are distinct.

For centered knots `z_i=x_i-q`, the divided-difference cap formula is

\[
 p=\sum_{i:z_i>0}\frac{z_i^4}{\prod_{j\ne i}(z_i-z_j)}.       \tag{4}
\]

Translation of all knots differentiates oppositely to the cutting level.
Consequently the section density and centroid satisfy

\[
 f(0)=\sum_{i=0}^4\frac{\partial p}{\partial z_i},
 \qquad
 \gamma_i=\frac{1}{f(0)}\frac{\partial p}{\partial z_i}.      \tag{5}
\]

Equations (4)--(5) are the common input to both mixed certificates.

## 2. Two knots below and three above

Translation and positive scaling put this crossing in the form

\[
 z=(-1,-u,tv,tw,t),qquad
 0<u<1,quad 0<v<w<1,quad t>0.                              \tag{6}
\]

### Theorem 1 (two-below/three-above barrier)

For (6), if `gamma_4 >= gamma_3`, then

\[
 \boxed{p\le I_{\gamma_0}(1,4)=1-(1-\gamma_0)^4.}            \tag{7}
\]

In particular, (7) holds when the whole section centroid is ordered.

### Compactification by the endpoint order

Exact simplification of (5) gives

\[
 \gamma_4-\gamma_3=\frac{(w-1)R(u,v,w,t)}{\Delta_0},
 \qquad \Delta_0>0,                                         \tag{8}
\]

where

\[
\begin{aligned}
R={}&t^5vw^2(u+1)(u^2+1)\\
 &+t^4u(u^2+u+1)w(2vw+2v+w)\\
 &+t^3u^2(u+1)(vw^2+4vw+v+2w^2+2w)\\
 &+t^2u^3(2vw+2v+w^2+4w+1)-u^4.                            \tag{9}
\end{aligned}
\]

The hypothesis and `w<1` imply `R<=0`.  Moreover,

\[
 R=(t^2u^3-u^4)+Q,                                         \tag{10}
\]

where every power-basis coefficient of `Q` is positive.  Hence
`t^2<=u<1`, so all four parameters in (6) lie in the unit cube.

### Canonical division certificate

Put

\[
\begin{aligned}
F={}&t^2u^2vw+t^2uvw+t^2vw+t u^2vw+t u^2v+t u^2w\\
   &+tuvw+tuv+tuw+u^2v+u^2w+u^2,\\
B={}&4(t+1)(tv+1)(tw+1)F,\\
U={}&(t+u)(tv+u)(tw+u).
\end{aligned}                                                \tag{11}
\]

Direct exact algebra using (4)--(5) gives

\[
 I_{\gamma_0}(1,4)-p=\frac{tT(u,v,w,t)}{B^4U}.              \tag{12}
\]

The integer polynomial `T` has tensor degree `(11,9,9,22)` and 6,895
nonzero power terms.  Set `G=-R`.  Lexicographic multivariate division over
the rationals in variable order `(u,v,w,t)` gives

\[
 T=S G+H.                                                    \tag{13}
\]

This is not a fitted ansatz: `S` and `H` are the exact quotient and remainder
of the deterministic division.  Their complete sign audit is

| polynomial | tensor degree | power terms | Bernstein positive | zero | negative |
|---|---:|---:|---:|---:|---:|
| `S` | `(7,15,22,54)` | 28,236 | 151,230 | 10,690 | 0 |
| `H` | `(3,16,24,59)` | 35,697 | 80,067 | 21,933 | 0 |

The smallest positive Bernstein coefficients are respectively
`64/36840513835125` and `16/4284919345485303525`.  Tensor Bernstein basis
functions are nonnegative on `[0,1]^4`; hence `S,H>=0`.  Orderedness gave
`G>=0`, so (13) gives `T>=0`, and (12) proves Theorem 1.

## 3. Three knots below and two above

The complementary mixed crossing normalizes to

\[
 z=(-1,-u,-v,tw,t),\qquad
 0<v<u<1,quad 0<w<1,quad t>0.                              \tag{14}
\]

### Theorem 2 (three-below/two-above barrier)

For (14), if `gamma_4 >= gamma_3`, then the same first-step inequality (7)
holds.

### Exact scale bound

Equations (4)--(5) again yield

\[
 \gamma_4-\gamma_3=\frac{(w-1)R_*(u,v,w,t)}{\Delta_*},
 \qquad \Delta_*>0.                                         \tag{15}
\]

Here `R_*` has tensor degree `(3,3,4,7)` and 87 power terms.  Its full
expansion is not conceptually useful; the generator derives it directly from
(4).  Since `w<1`, the hypothesis gives `R_*<=0`.

The scale is compact for an exact rather than numerical reason.  Substitute
`t=2+rho`.  In degree `(3,3,4)` tensor Bernstein form with respect to
`(u,v,w)` and the ordinary power basis in `rho`, the resulting polynomial has
366 positive coefficients, 274 zeros, and no negative coefficient.  The
smallest positive coefficient is `1/18`, and the `rho^0` slice has a positive
coefficient.  Every Bernstein basis function is strictly positive in the
open `(u,v,w)` cube.  Thus

\[
 R_*(u,v,w,t)>0\quad\text{for }t\ge2,                       \tag{16}
\]

contradicting `R_*<=0`.  Therefore `t<2`.

### Sparse exact multiplier certificate

The target margin has the exact form

\[
 I_{\gamma_0}(1,4)-p=\frac{tT_*(u,v,w,t)}{\Delta},
 \qquad \Delta>0,                                           \tag{17}
\]

where `T_*` has tensor degree `(6,6,14,19)` and 4,015 power terms.  Put
`t=2 tau`, so `0<tau<1`, and set

\[
 \widetilde T=T_*(u,v,w,2\tau),qquad
 \widetilde G=-R_*(u,v,w,2\tau).                            \tag{18}
\]

The certificate constructs a tensor-Bernstein polynomial `S` of degree
`(3,3,10,12)`.  Its 363 nonzero positions and 363 target coefficient rows
specify a square rational linear system.  The checker derives every matrix
entry from the exact Bernstein product formula, solves the system over the
rationals, and confirms that it is nonsingular and that all 363 resulting
coefficients are positive.  The remaining 1,925 coefficients of `S` are zero.
It then checks the identity

\[
 \widetilde T=S\widetilde G+H                              \tag{19}
\]

against every coefficient, not only the defining rows.  The degree
`(6,6,14,19)` expansion of `H` has 7,526 positive coefficients, 7,174 zeros,
and no negative coefficient.  Therefore `S,H>=0` on the cube.  Equation
(15) gives `widetilde G>=0`, so (19) proves `widetilde T>=0`; (17) proves
Theorem 2.

The support and row indices are only a witness for the exact linear system.
No decimal multiplier value is stored or trusted.  The checked certificate
records SHA-256 digests of the complete exact multiplier and residual.

## 4. Every distinct-knot crossing

### Corollary 3 (distinct-knot local barrier)

Every regular five-coordinate section with distinct ordered knots and an
ordered section centroid satisfies (3).

**Proof.** Theorems 1 and 2 handle the two mixed crossings.  If only the
largest knot is above the level, put

\[
 a_i=\frac{x_4-q}{x_4-x_i},\qquad i<4.
\]

The upper cap is a simplex and its section is its base, so

\[
 p=\prod_{i<4}a_i,qquad M_4=\frac14\sum_{i<4}a_i.
\]

AM--GM gives `p<=M_4^4=I_(M_4)(4,1)`.  If only the smallest knot is below
the level, put `b_i=(q-x_0)/(x_i-x_0)` for `i>0`.  Then

\[
 p=1-\prod_{i>0}b_i,qquad
 \gamma_0=1-\frac14\sum_{i>0}b_i,qquad
 \gamma_i=b_i/4.
\]

Knot order gives `b_1>=...>=b_4`, while centroid order gives the reverse.
All four would have to be equal, forcing the four upper knots to coincide.
That is impossible in the distinct-knot stratum.  These four cases exhaust
the regular crossings.  ∎

Corollary 3 means that a five-coordinate obstruction cannot occur at an
interior stationary cap with distinct knots.  The only possible remaining
obstructions lie on knot-collision faces, a strictly smaller and structurally
different problem because the collapsed block masses have laws
`Dirichlet(m_1,...,m_k)` with positive integer parameters.  The companion
note gives first-prefix inequalities on those families. An active-prefix
inequality on repeated-lowest faces remains to be proved.

## 5. Verification boundary

The public generator recomputes all of the following from (4):

1. the cap and centroid coordinates used in both mixed crossings;
2. every positive denominator and both endpoint-order factorizations;
3. the exact compactness implications `t<1` and `t<2`;
4. the 6,895- and 4,015-term comparison numerators;
5. the first crossing's exact quotient-remainder division;
6. all 263,920 Bernstein coefficients in that division certificate;
7. the second crossing's 363-by-363 exact rational solve; and
8. all 2,288 multiplier and 14,700 residual coefficients in the second
   certificate.

The compact JSON artifact records degrees, term counts, coefficient extrema,
and SHA-256 digests.  SymPy performs exact polynomial arithmetic, FLINT solves
the rational system, and Python `Fraction` performs every Bernstein conversion
and sign comparison.  Reproduce the artifact byte-for-byte with

```console
make five-coordinate-mixed-section-check
```

The ordinary proof consists of the geometric reduction, the scale bounds, the
endpoint crossings, and the implication from Bernstein positivity.  The
machine component certifies the large exact identities and coefficient signs.
Repeated-knot Dirichlet-block barriers are not consequences of continuity
from the distinct stratum. The separate first-prefix inequalities, with
their own adjacent-order identities and exact certificates, are recorded in
[`FIVE-COORDINATE-VERTEX-BARRIER.md`](FIVE-COORDINATE-VERTEX-BARRIER.md).
They do not yet prove the active-prefix inequalities on every collision face.
