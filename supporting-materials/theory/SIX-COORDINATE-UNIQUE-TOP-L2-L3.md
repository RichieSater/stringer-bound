# The last two six-coordinate unique-top barriers

## Status and theorem

This note proves the two unique-top local barriers with two and three knots
below the section. Together with
[SIX-COORDINATE-UNIQUE-TOP.md](SIX-COORDINATE-UNIQUE-TOP.md) and
[SIX-COORDINATE-REPEATED-TOP.md](SIX-COORDINATE-REPEATED-TOP.md), this
supplies first-prefix calculations for all nine canonical non-elementary
families. The global six-coordinate theorem remains conjectural because
the active-prefix inequalities on repeated-lowest faces are unresolved; see
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).

Let

\[
 D\sim\operatorname{Dirichlet}(1,1,1,1,1,1)
\]

and consider a regular section at zero with one of the ordered knot vectors

\[
 (-1,-u,ta,tb,tc,t),                                      \tag{1}
\]

where `1 >= u > 0`, `0 < a <= b <= c < 1`, or

\[
 (-1,-u,-v,ta,tb,t),                                      \tag{2}
\]

where `1 >= u >= v > 0`, `0 < a <= b < 1`. Write

\[
 p=\Pr\{x\mathbin{\cdot}D>0\},\qquad
 \gamma_i=\mathbb E(D_i\mid x\mathbin{\cdot}D=0).
\]

**First-prefix theorem on the open strata.** If the knots are strictly
ordered and the section centroid is ordered, then in both (1) and (2)

\[
 p\le
 \frac{5\gamma_0}{1+(10974/4651)\gamma_0}
 \le 1-(1-\gamma_0)^5.                                   \tag{3}
\]

The result extends to limits within the stated polynomial domains by
continuity. Vanishing order factors do not automatically allow every
ordered collision point. The first inequality strengthens the first beta
prefix bound, which is an active gap only when the smallest knot is unique.

## 1. A sharp low-degree beta minorant

For `0 < g <= 1/6`, the inequality

\[
 \frac{5g}{1+cg}\le 1-(1-g)^5
\]

is equivalent to

\[
 c\ge h(g):=
 -\frac{g^3-5g^2+10g-10}
 {g^4-5g^3+10g^2-10g+5}.                                \tag{4}
\]

Its derivative factors as

\[
 h'(g)=
 \frac{(g^2-5g+5)(g^4-5g^3+15g^2-20g+10)}
 {(g^4-5g^3+10g^2-10g+5)^2}>0                           \tag{5}
\]

on `[0,1/6]`. Therefore the least admissible coefficient is

\[
 h(1/6)=\frac{10974}{4651}.                              \tag{6}
\]

Equality at zero is understood by continuity. Since ordered barycentric
coordinates sum to one, an ordered section centroid always has
`0 <= gamma_0 <= 1/6`, so this minorant applies throughout the required
domain.

This is the main algebraic simplification. Direct comparison with
`1-(1-gamma_0)^5` produces much larger target polynomials. The
linear-fractional lower bound keeps the exact implication needed for the
beta barrier while reducing the finite sign problem substantially.

## 2. Interpolation reduction

For distinct knots `y_0,...,y_5`, the upper cap at zero is the divided
difference

\[
 p=\sum_{y_i>0}\frac{y_i^5}{\prod_{j\ne i}(y_i-y_j)},      \tag{7}
\]

with the confluent value used at collisions. Its density at the section is
obtained by replacing `y_i^5` with `5y_i^4`, and differentiation gives

\[
 \gamma_i=\frac{\partial p/\partial y_i}{f(0)}.            \tag{8}
\]

Write `p=a/b` and `gamma_0=n/d` after positive denominators are chosen.
The numerator of the first difference in (3) is

\[
 5nb-ad-\frac{10974}{4651}an.                             \tag{9}
\]

Removing positive denominators and monomials from (9) gives a target
polynomial `T_L`. Similarly, every adjacent centroid difference has a
positive denominator and, after its obvious order factor is removed, a gap
polynomial `G_(L,i)`. On the open distinct-knot stratum, orderedness is
therefore exactly the system

\[
 G_{L,i}\ge0,\qquad i=0,\ldots,4.                         \tag{10}
\]

For `L=2`, use

\[
 u=U,\quad a=ABC,\quad b=BC,\quad c=C,\quad t=q.           \tag{11}
\]

The target has tensor degree `(7,4,8,12,14)` and 2,779 nonzero power
coefficients. For `L=3`, use

\[
 u=U,\quad v=UV,\quad a=AB,\quad b=B,\quad t=q.            \tag{12}
\]

The target has tensor degree `(10,5,6,12,15)` and 4,687 nonzero power
coefficients. These polynomials and all five gap polynomials are derived
from (7)--(9), rather than stored as unverified input.

## 3. Projective charts for `L = 2`

The terminal gap first gives the global scale bound `q <= 1`. Split at
`q=U`.

* On `q <= U`, write `q=UZ`.
* On `U <= q`, write `U=qZ` and compare `Z` with `C`.

If `C <= Z`, put

\[
 C=Zd,\qquad q=ZR.
\]

The terminal gap gives `q <= Z`, so this is a unit chart.

If `Z <= C`, put

\[
 Z=Cd,\qquad q=(C/4)R.
\]

Here the terminal gap gives `q <= C/4`. The position of `d` relative to
`AB` and `B` gives three exhaustive regions:

\[
\begin{array}{c|c|c}
\text{region}&\text{substitution}&\text{scale consequence}\\ \hline
d\le AB&d=ABe&R\le4A^4B^3e^5\\
AB\le d\le B&d=Be,\ A=ef&R\le4B^3e^4\\
B\le d&B=de&R\le4d^3 .
\end{array}                                                \tag{13}
\]

Every implication in (13) is an exact coefficientwise-nonnegative
Bernstein identity for the terminal gap. Scaling `R` by the corresponding
right side compactifies each singular corner into a unit cube.

## 4. Projective charts for `L = 3`

The terminal gap gives `q <= 2`. Again split at `q=U`; in the second chart
write `U=qZ` and compare `Z` with `B`.

If `Z <= B`, put

\[
 Z=Bd,\qquad q=BR.
\]

The terminal gap gives `q <= B`. Comparing `A` with `dV` and `d` produces

\[
\begin{array}{c|c|c}
\text{region}&\text{substitution}&\text{scale consequence}\\ \hline
d\le A&d=Ae&R\le8A^3e^4\\
dV\le A\le d&A=de,\ V=ef&R\le4d^3\\
A\le dV&A=dVe&R\le4d^3 .
\end{array}                                                \tag{14}
\]

If `B <= Z`, put

\[
 B=Zd,\qquad q=4ZR.
\]

For `Z <= 1/2`, the terminal gap proves `q <= 4Z`; for `Z >= 1/2`, this
follows from `q <= 2`. The remaining three order regions are

\[
 V=Ade,\qquad V=de,\ A=ef,\qquad d=Ve.                    \tag{15}
\]

Thus (14)--(15) also replace the simultaneous small-scale limit by finitely
many compact unit charts.

## 5. Finite sign conclusion

On each terminal box, one of the following exact implications is used:

1. every local tensor-Bernstein coefficient of `T_L` is nonnegative;
2. every coefficient of `-G_(L,4)` is nonnegative and the relative
   interior is incompatible with orderedness;
3. for an explicit nonnegative integer `mu`, every coefficient of
   `T_L-mu*G_(L,4)` is nonnegative.

The exact rational covers settling the ordered-domain implication are:

| family and chart | boxes |
|---|---:|
| `L=2`, `q <= U` | 37 |
| `L=2`, `C <= Z` | 78 |
| `L=2`, middle region of (13) | 54 |
| `L=2`, upper region of (13) | 30 |
| `L=3`, all seven projective charts | 56 |

The lower region of (13) has 11 further boxes. Its target degree is
`(10,45,37,19,53)`, making a full exact-rational tensor conversion
impractical. This one cover starts with the exact integer power
coefficients, rounds every rational Bernstein weight to a binary64 center,
and records its exact rational rounding error. Exact rational forward-error
bounds are then propagated through every tensor conversion, degree
elevation, dyadic subdivision, and integer-multiplier subtraction.
Underflow and overflow are rejected. Each stored lower center remains
strictly positive after the corresponding absolute error bound is
subtracted; the minimum certified margin is greater than `0.339`.
Accordingly this is a rigorous dyadic enclosure, not an unqualified
binary64 sign computation.

The exact scale identities, the 255 exact-rational terminal boxes, and the
11 rigorously enclosed boxes exhaust (10). They prove `T_L >= 0` and hence
the first inequality in (3). Equations (4)--(6) prove the second.

## 6. Consequence and reusable idea

The companion calculations cover the unique-top `L=4` first-prefix family
and the six repeated-top families. The two elementary crossings are
analytic. If the active-prefix collision inequalities are established,
the section-centroid maximum principle will give the following conjecture
for every nondecreasing six-coordinate probability vector `c` and every
ordered knot vector `x`:

\[
 \Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}
 \le\max_{1\le r<6} I_{C_r}(r,6-r),
 \qquad C_r=\sum_{i<r}c_i.                                \tag{16}
\]

The new transferable mechanism has two parts: order the ratios inside the
central knot cluster to obtain projective scale charts, and replace a
high-degree beta target by a sharp low-degree minorant tailored to the
smallest centroid coordinate. Neither step depends on the
Clopper--Pearson weights. Those weights enter only afterward, when (16) is
conditionally applied to the `n=5` Stringer problem.

## Reproduction

Run

```text
make six-coordinate-unique-top-l2-l3-check
```

The derivation and verification are in
`computations/python/six_coordinate_unique_top_l2_l3.py`; the compact
subdivision trees are in
`computations/python/six_coordinate_unique_top_l2_l3_witness.py`; the
generated result is
`computations/certificates/six-coordinate-unique-top-l2-l3-certificate.json`.
The tests reject altered multipliers and incomplete trees, including a
mutation of the rigorously enclosed binary64 chart.
