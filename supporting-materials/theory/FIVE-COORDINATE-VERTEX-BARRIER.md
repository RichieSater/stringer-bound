# Five-coordinate first-prefix inequalities and the collision conjecture

## Status and thesis

This note records exact first-prefix inequalities in five coordinates. The
active-prefix local statement and the global monotone cap theorem below are
conjectural on repeated-lowest-knot faces; see [`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md). Let

\[
 D=(D_0,\ldots,D_4)\sim\operatorname{Dirichlet}(1,1,1,1,1),
 \qquad x_0\le\cdots\le x_4,
\]

let `q` be a regular level, and write

\[
 p=\Pr\{x\mathbin{\cdot}D>q\},\qquad
 \gamma_i=\mathbb E(D_i\mid x\mathbin{\cdot}D=q),\qquad
 M_r=\sum_{i<r}\gamma_i.                                    \tag{1}
\]

The unresolved active-prefix conjecture is

\[
 \boxed{
 \gamma_0\le\cdots\le\gamma_4
 \quad\Longrightarrow\quad
 p\le\max_{r:\,x_{r-1}<x_r}I_{M_r}(r,5-r).}                \tag{2}
\]

Statement (2) requires every knot-collision face, not merely a first-prefix
bound there. If (2) is established, it implies the following global
conjecture: for every nondecreasing probability vector
`c=(c_0,...,c_4)`, every ordered `x`, and `C_r=sum_(i<r)c_i`,

\[
 \boxed{
 \Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}
 \le \max_{1\le r\le4}I_{C_r}(r,5-r).}                     \tag{3}
\]

The available first-prefix calculations have two parts. The distinct-knot
identities use only the last centroid comparison, so those identities extend
to collisions that leave the largest knot unique. This does not identify
the first prefix with an active gap when the smallest knot repeats.  If the largest knot is repeated, only three
nontrivial configurations survive.  In each, the adjacent nontrivial
centroid comparison both bounds the remaining scale and acts as a multiplier
constraint in an exact Bernstein certificate.

The large sign layers are computer-assisted but exact.  Polynomial identities
are derived over the rationals, two multipliers are reconstructed by exact
FLINT solves of orders 111 and 38, and every Bernstein sign is checked with
`Fraction` arithmetic.  Floating point is not part of the verification.

## 1. Conditional implication from the active-prefix inequality

For fixed `c`, put

\[
 \Phi_c(x)=\Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}.
\]

At a regular section, differentiation under the integral gives

\[
 \frac{\partial\Phi_c}{\partial x_i}
 =f_x(c\mathbin{\cdot}x)(\gamma_i-c_i),                    \tag{4}
\]

where `f_x` is the density of `x dot D`.  On a knot face, collect equal
consecutive knots into blocks.  Stationarity with respect to each block knot
says that the conditional mass of that block equals the corresponding sum of
the `c_i`.  Symmetry within a block makes every individual centroid
coordinate equal to the average of the `c_i` in that block.  Consecutive
averages of a nondecreasing vector are nondecreasing, so the hypothesis of
(2) holds.

Every strict gap is a boundary between whole blocks.  At such a gap,
`M_r=C_r`; if (2) is established, it bounds the stationary value by the
right side of (3). If a
face has no interior maximum, its maximum is on a proper knot face.  Induction
over the finite face lattice, followed by continuity at singular sections,
would prove (3). The first-prefix certificates below do not establish (2)
when the least knot repeats. No arbitrary Dirichlet-block theorem is needed:
only the
positive-integer block types whose masses sum to five must be checked.

## 2. The unique-largest-knot faces

The companion distinct-stratum proof is
[`FIVE-COORDINATE-MIXED-SECTION.md`](FIVE-COORDINATE-MIXED-SECTION.md).
Its two mixed cases normalize to

\[
 (-1,-u,tv,tw,t),\qquad (-1,-u,-v,tw,t),                    \tag{5}
\]

and use only `gamma_4>=gamma_3`.  Their exact identities and Bernstein
certificates remain valid on the closed parameter boundaries as long as
`w<1`, that is, as long as the largest knot stays unique.  This includes
collisions among the lower knots and among nonlargest upper knots.

If only the largest knot is above the level, the cap is a simplex and AM--GM
gives `p<=I_(M_4)(4,1)`.  If only the smallest knot is below the level, knot
order reverses the four positive interpolation ratios while centroid order
requires them to increase.  Hence all four upper knots coincide; that case
belongs to the repeated-largest endpoint treated below.  Therefore only a
repeated largest knot remains.

## 3. Exhausting the repeated-largest faces

A regular section with a repeated largest knot has one of the following
nontrivial forms after translation and positive scaling:

\[
\begin{array}{ccl}
 \text{three below, top pair:}
   &(-1,-u,-us,t,t),&0<u,s<1,\ t>0,\\
 \text{two below, top pair:}
   &(-1,-u,tv,t,t),&0<u,v<1,\ t>0,\\
 \text{two below, top triple:}
   &(-1,-u,t,t,t),&0<u<1,\ t>0.
\end{array}                                                  \tag{6}
\]

Lower-knot collisions are included by closure.  If one knot is below and four
are above, orderedness forces the four upper knots to coincide.  With two
blocks of masses one and four, the cap is exactly
`I_(gamma_0)(1,4)`.  These cases exhaust all multiplicities and crossing
positions.

For each line of (6), direct differentiation of the divided-difference cap
formula gives

\[
 I_{\gamma_0}(1,4)-p=\frac{tT}{\Delta},\qquad \Delta>0.      \tag{7}
\]

The relevant adjacent centroid order gives a polynomial `G>=0`.  The next
three sections certify `T>=0` from that single constraint.

## 4. Three below and a repeated top pair

Here `gamma_3=gamma_4`, and the nontrivial order condition is
`gamma_3>=gamma_2`.  Its numerator is a degree `(4,2,5)` polynomial
`G(u,s,t)`.  On substituting `t=3/2+rho`, the mixed expansion of `-G` in the
tensor Bernstein basis for `(u,s)` and the power basis for `rho` has 89
positive coefficients, one zero, and no negative coefficient.  Its
`rho^0` slice is strictly positive in the open square.  Thus orderedness
forces

\[
 t<\frac32.                                                  \tag{8}
\]

Put `t=(3/2)tau`.  The target `T` has degree `(12,6,15)` and 369 power terms.
A sparse support-and-row witness defines a 111-by-111 rational system for a
Bernstein multiplier `S` of degree `(8,4,10)`.  The checker derives and solves
the system exactly.  All 111 supported coefficients are positive and the
other 384 are zero.  It then checks every coefficient of

\[
 T=S G+H.                                                     \tag{9}
\]

The degree `(12,6,15)` residual has 727 positive coefficients, 729 zeros, and
no negative coefficient.  Hence `S,H>=0` on the cube; `G>=0` proves (7) is
nonnegative.

## 5. Two below and a repeated top pair

Here the relevant condition is again `gamma_3=gamma_4>=gamma_2`.  Exact
simplification gives

\[
 \gamma_3-\gamma_2=\frac{(v-1)R}{\Delta_0},\qquad
 G=-R,\quad \Delta_0>0.                                     \tag{10}
\]

Because `v<1`, orderedness gives `R<=0`.  Moreover,

\[
 R=t^2u^3-u^4+Q,                                             \tag{11}
\]

where every one of the 19 power coefficients of `Q` is positive.  It follows
that `t^2<u<1`, so the original `(u,v,t)` variables already lie in the unit
cube.

The target has degree `(11,9,18)` and 910 power terms.  Deterministic
lexicographic division by `G` gives

\[
 T=S G+H.                                                     \tag{12}
\]

The quotient `S` has degree `(7,22,50)`; its 9,384 Bernstein coefficients
split into 8,733 positive values and 651 zeros.  The remainder `H` has degree
`(3,24,55)`; its 5,600 coefficients split into 4,424 positive values and
1,176 zeros.  There are no negative coefficients in either expansion, so
(12) proves the desired first-step bound.

## 6. Two below and a repeated top triple

Now `gamma_2=gamma_3=gamma_4`, and orderedness supplies
`gamma_2>=gamma_1`.  Its numerator is a degree `(3,4)` polynomial `G(u,t)`.
The mixed Bernstein/power expansion of `-G(u,2/3+rho)` has 19 positive
coefficients, one zero, and no negative coefficient, with strict positivity
at `rho=0` in the open interval.  Hence

\[
 t<\frac23.                                                  \tag{13}
\]

After `t=(2/3)tau`, the target has degree `(11,14)` and 105 power terms.  A
38-by-38 exact rational system constructs a degree `(8,10)` multiplier with
38 positive supported coefficients and 61 zeros.  The identity (9) is then
checked in all 180 target positions.  The residual has 73 positive
coefficients, 107 zeros, and no negative coefficient. This proves the
first-prefix inequality on this block family, but not (2) on its
repeated-lowest subfaces.

## 7. Verification boundary

The ordinary proof comprises:

1. the section-centroid derivative and stationary-face induction;
2. the face exhaustion separating a unique largest knot from (6);
3. the implication from each adjacent centroid order to `G>=0`;
4. the three exact scale bounds; and
5. the implication from Bernstein positivity in (9) and (12).

The public generator
[`five_coordinate_block_faces.py`](../computations/python/five_coordinate_block_faces.py)
starts from the divided-difference cap formula and recomputes the cap,
centroids, denominators, comparison numerators, compactness expansions,
lexicographic division, both rational systems, and every coefficient sign.
The index-only witnesses in
[`five_coordinate_block_witness.py`](../computations/python/five_coordinate_block_witness.py)
contain no approximate coefficient.  The compact artifact
[`five-coordinate-block-faces-certificate.json`](../computations/certificates/five-coordinate-block-faces-certificate.json)
records degrees, counts, extrema, and SHA-256 digests.  Reproduce it with

```console
make five-coordinate-block-faces-check
```

The distinct/unique-top identities are independently regenerated by

```console
make five-coordinate-mixed-section-check
```

Together these commands check the first-prefix polynomial inequalities.
They do not establish the missing active-prefix implication needed for (3).

## 8. Mathematical yield

The conjecture asks whether every monotone barycentric threshold on the
four-dimensional uniform simplex is maximized at a step vector. The weights
are arbitrary; Clopper--Pearson weights would be one application. The
first-prefix certificates settle the distinct-lowest cases but leave an
active-prefix obligation on repeated-lowest faces.

The proof also exposes the next reusable pattern.  A terminal centroid
comparison controls distinct or unique-top sections.  When that comparison
degenerates because the top knot repeats, move left to the next nontrivial
block boundary.  Its order polynomial again both compactifies the scale and
serves as the multiplier constraint.  Whether this hierarchy continues in
six and higher coordinates is now a precise structural question.

What remains dependent on computation is explicit: two exact sparse
multipliers, one exact quotient--remainder division, and their finite rational
sign tables.  The face classification, global induction, and composition of
those certificates establish the stated first-prefix bounds. Their
composition into (2)--(3) remains unresolved.
