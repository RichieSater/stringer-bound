# The \(L=4\) unique-top barrier in six coordinates

## Status and theorem

Let

\[
D=(D_0,\ldots,D_5)\sim\operatorname{Dirichlet}(1,\ldots,1)
\]

and consider the regular section with normalized knots

\[
(-1,-u,-v,-w,ta,t),
\qquad 1\geq u\geq v\geq w>0,\quad 0<a<1,\quad t>0.       \tag{1}
\]

Write

\[
p=\Pr\{-D_0-uD_1-vD_2-wD_3+taD_4+tD_5>0\}
\]

and

\[
\gamma_i=\mathbb E(D_i\mid
 -D_0-uD_1-vD_2-wD_3+taD_4+tD_5=0).
\]

**First-prefix theorem on the open stratum.** If the knots in (1) are
strictly ordered and the section centroid is ordered,
`gamma_0<=...<=gamma_5`, then

\[
p\leq I_{\gamma_0}(1,5)=1-(1-\gamma_0)^5.                 \tag{2}
\]

The conclusion extends to limits of that polynomial domain by continuity,
not automatically to every ordered collision section after a divided-out
order factor vanishes. This is the
unique-top family with four coordinates below the section. The companion
[`SIX-COORDINATE-UNIQUE-TOP-L2-L3.md`](SIX-COORDINATE-UNIQUE-TOP-L2-L3.md)
proves the other two unique-top rows. Together with the six repeated-top
families, these results establish the first-prefix inequalities on their
stated domains. They do not establish the active-prefix collision
inequalities needed for the conjectural global six-coordinate theorem; see
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).

## 1. Interpolation reduction

For distinct knots `y_0,...,y_5`, let

\[
\mathcal D_k(y)=
\sum_{i:y_i>0}\frac{y_i^k}{\prod_{j\ne i}(y_i-y_j)}.       \tag{3}
\]

The cap probability in (1) is `p=mathcal D_5`, its section density is
`5 mathcal D_4`, and differentiation with respect to a knot gives the
corresponding section-centroid coordinate.  Substituting (1), clearing the
positive denominators, and removing common positive factors gives a
comparison polynomial `T` and five adjacent-gap polynomials
`G_0,...,G_4` such that

\[
I_{\gamma_0}(1,5)-p\quad\hbox{has the sign of }T,           \tag{4}
\]

while orderedness is exactly `G_i>=0` in the relative interior.

Use ordered coordinates

\[
u=U,\qquad v=UV,\qquad w=UVW,\qquad a=A,\qquad t=q.        \tag{5}
\]

The degrees and power-term counts are

| polynomial | degree in `(U,V,W,A,q)` | power terms |
|---|---:|---:|
| `T` | `(21,14,7,23,30)` | 87,047 |
| `G_0` | `(4,2,1,5,6)` | 74 |
| `G_1` | `(5,3,1,5,6)` | 74 |
| `G_2` | `(5,4,2,5,6)` | 74 |
| `G_3` | `(7,5,3,7,8)` | 218 |
| `G_4` | `(9,6,3,6,10)` | 572 |

The terminal gap also gives the first compactness bound.  Substitute
`q=3+rho` in `-G_4` and use the tensor-Bernstein basis in `(U,V,W,A)` and
the power basis in `rho`.  Its 21,560 coefficients consist of 10,279
positive and 11,281 zero coefficients.  Hence `G_4>=0` forces `q<3` on a
regular ordered section.

## 2. The first scale split

Split according to `q<=U` or `U<=q`.  In the first chart put `q=UZ`.
After a positive monomial is removed, `T` has degree
`(14,14,7,23,30)`.  All of its 1,339,200 tensor-Bernstein coefficients are
nonnegative: 776,265 are positive and 562,935 vanish.  This proves (2) on
the first chart without case subdivision.

In the second chart put

\[
U=qZ.                                                        \tag{6}
\]

The remaining singularity is the relative size of `A` and `Z`.  Two
projective charts, `Z=AB` and `A=ZB`, cover it.

## 3. The chart `Z<=A`

Set `Z=AB`.  The terminal gap successively implies

\[
q\leq2A,\qquad R:=\frac q{2A}\leq B,
\qquad S:=\frac RB\leq B,
\qquad C:=\frac SB\leq6B.                                  \tag{7}
\]

Each implication is coefficientwise in an exact Bernstein expansion of
`-G_4` after the corresponding shifted substitution.  The preceding bound
puts each new ratio in `[0,1]`: first `R<=1`, then `S<=1`, and finally
`C<=1`.  Thus the compact shifted boxes cover every possible violation of
the next inequality.  The positive/zero coefficient counts in (7) are, in
order,

\[
(9380,1820),\quad(12852,252),\quad(14868,252),
\quad(16884,252).                                           \tag{8}
\]

There are no negative coefficients.  Writing `C=6BD` turns (6)--(7) into

\[
q=12AB^3D,qquad U=12A^2B^4D.                              \tag{9}
\]

In the resulting cube `(A,B,V,W,D)`, the target has degree
`(36,61,14,7,14)`.  A 53-node, 27-leaf rational cover closes eight leaves
directly by `T>=0`, eight by `-G_4>0`, and eleven by

\[
T-\mu G_4\geq0                                             \tag{10}
\]

with an explicit nonnegative integer `mu`.  On the last leaves,
`G_4>=0` and (10) imply `T>=0`; a leaf with `-G_4>0` contains no ordered
relative-interior point.

## 4. The chart `A<=Z`

Set `A=ZB`.  The central four knots are then

\[
qZ(-1,-V,-VW,B).                                           \tag{11}
\]

This common scale makes the correct projective decomposition visible.
Compare `B` first with `V`, and when `B<=V` compare `D=B/V` with `W`.
The three resulting regions are

\[
B\geq V,qquad VW\leq B\leq V,qquad B\leq VW.             \tag{12}
\]

They have the following cube coordinates.

* For `B>=V`, write `V=BC`.  The terminal gap gives `q<=2Z`; after
  `q=2ZR`, a 17-node, 9-leaf cover closes five leaves by `T` and four by
  `-G_4`.
* For `VW<=B<=V`, write `B=VD` and `W=DC`.  On `Z<=1/2`, the terminal gap
  gives `q<=4Z`.  With `Z=z/2` and `q=2zR`, a 19-node, 10-leaf cover closes
  six leaves by `T` and four by `-G_4`.  On `Z>=1/2`, splitting at
  `V=1/2` gives two boxes on which every target coefficient is nonnegative.
* For `B<=VW`, write `B=VWC`.  On `Z<=1/3`, the terminal gap gives
  `q<=9Z`.  With `Z=z/3` and `q=3zR`, a 53-node, 27-leaf cover closes
  thirteen leaves by `T` and fourteen by `-G_4`.  The range `Z>=1/3` is
  covered by four direct target boxes and a final 5-node, 3-leaf split; its
  closures are two by `T` and one by `-G_4`.

The three scale bounds in this paragraph again follow from exact mixed
Bernstein/power expansions.  Their positive/zero counts are respectively

\[
(6356,1484),\qquad(5816,2024),\qquad(9460,4260),            \tag{13}
\]

with no negative coefficients.  Equations (11)--(13) cover the entire
`A<=Z` chart, including their interfaces by continuity.

## 5. Consequence and mathematical yield

The two first-scale charts and the two `A`--`Z` projective charts exhaust
the polynomial domain `G_i>=0`. Every leaf either proves `T>=0` or excludes
an ordered relative-interior point, so (4) proves (2).

The useful mechanism is the resolution of a four-knot central cluster.  A
direct compact cube hides the fact that the positive inner knot must be
compared successively with the first and last negative ratios.  The three
regions in (12) expose a common cluster scale, after which the terminal
centroid inequality yields linear scale bounds. The companion \(L=2,3\)
proofs reuse this projective pattern and combine it with a sharp
linear-fractional minorant of the beta tail.

The exact calculation is implemented in
[`six_coordinate_unique_top.py`](../computations/python/six_coordinate_unique_top.py),
with its finite partitions in
[`six_coordinate_unique_top_witness.py`](../computations/python/six_coordinate_unique_top_witness.py)
and its output in
[`six-coordinate-unique-top-certificate.json`](../computations/certificates/six-coordinate-unique-top-certificate.json).
It is replayed by

```sh
make six-coordinate-unique-top-check
```
