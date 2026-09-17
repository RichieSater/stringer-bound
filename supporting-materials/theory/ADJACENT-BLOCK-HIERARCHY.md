# A finite adjacent-block hierarchy for ordered simplex caps

## Status and thesis

This note turns the section-centroid method into a dimension-independent,
finite proof architecture.  Every knot-collision face is indexed by an
ordered integer composition of the number of Dirichlet coordinates.  On that
face, the relevant section centroid is a Dirichlet-block centroid, and the
global step-maximum theorem follows once one local beta-cap inequality is
proved for each composition and crossing.

The reduction theorem below is proved. It is not by itself a proof of the
local barrier in arbitrary dimension. The five- and six-coordinate global
applications remain conjectural: first-prefix inequalities do not supply
the required active-prefix inequality on repeated-lowest faces. See
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md). The nine six-coordinate families organize the available
polynomial calculations, not a completed global theorem.

## 1. Knot faces are integer compositions

Let

\[
 D=(D_0,\ldots,D_{N-1})\sim\operatorname{Dirichlet}(1,\ldots,1),
 \qquad x_0\le\cdots\le x_{N-1}.
\]

Collect equal consecutive knots into blocks `B_1,...,B_k`, with block sizes

\[
 \lambda=(\lambda_1,\ldots,\lambda_k),\qquad
 \lambda_b=|B_b|,\qquad \sum_b\lambda_b=N.                 \tag{1}
\]

Thus `lambda` is an ordered composition of `N`.  If the distinct block knots
are `y_1<...<y_k`, put

\[
 W_b=\sum_{i\in B_b}D_i.
\]

Dirichlet aggregation gives

\[
 W=(W_1,\ldots,W_k)\sim
 \operatorname{Dirichlet}(\lambda_1,\ldots,\lambda_k).      \tag{2}
\]

For a regular section `y dot W=q`, write

\[
 p=\Pr\{y\mathbin{\cdot}W>q\},\qquad
 \Gamma_b=\mathbb E(W_b\mid y\mathbin{\cdot}W=q).           \tag{3}
\]

Conditional exchangeability inside a block gives

\[
 \gamma_i=\frac{\Gamma_b}{\lambda_b}\quad(i\in B_b).        \tag{4}
\]

Consequently, the full section centroid is nondecreasing exactly when

\[
 \frac{\Gamma_1}{\lambda_1}\le\cdots\le
 \frac{\Gamma_k}{\lambda_k}.                               \tag{5}
\]

At the boundary after block `b`, define

\[
 s_b=\lambda_1+\cdots+\lambda_b,qquad
 M_b=\Gamma_1+\cdots+\Gamma_b.                              \tag{6}
\]

The beta step cap attached to that boundary is
`I_(M_b)(s_b,N-s_b)`.

## 2. The composition criterion

**Theorem (finite composition criterion).**  Fix `N>=2`.  Suppose that for
every ordered composition `lambda` of `N`, every strictly increasing block
knot vector `y`, and every regular section whose block-centroid averages
satisfy (5),

\[
 p\le \max_{1\le b<k} I_{M_b}(s_b,N-s_b).                  \tag{7}
\]

Suppose the same statement holds at singular levels by continuity.  Then for
every nondecreasing probability vector
`c=(c_0,...,c_(N-1))` and every ordered knot vector `x`,

\[
 \Pr\{x\mathbin{\cdot}D>c\mathbin{\cdot}x\}
 \le \max_{1\le r<N}I_{C_r}(r,N-r),
 \qquad C_r=\sum_{i<r}c_i.                                 \tag{8}
\]

Every term on the right of (8) is attained by its corresponding step vector.

**Proof.**  Normalize the least and greatest distinct knots to zero and one.
On a face of type `lambda`, vary a block knot `y_b` while retaining the
fixed barycentric threshold.  Differentiation of the cap gives

\[
 \frac{\partial\Phi_c}{\partial y_b}
 =f(q)\left(\Gamma_b-\sum_{i\in B_b}c_i\right),             \tag{9}
\]

where `f(q)>0` is the section density.  At a relative-interior stationary
point, (9), together with the total-mass and threshold identities, makes
every block centroid mass equal its corresponding `c`-mass.  Formula (4)
then identifies the individual centroid coordinates with consecutive block
averages of the nondecreasing vector `c`, proving (5).  At each strict block
gap, (6) gives `M_b=C_(s_b)`, so (7) bounds the stationary value by (8).

A nonstationary maximum lies on a proper knot face.  Induction down the
finite face lattice ends at a step vector.  There,
`sum_(i<r)D_i` has the `Beta(r,N-r)` law and the cap equals the corresponding
term in (8).  Continuity covers nonregular levels and nonstrict weights.
This proves both the bound and attainment.  ∎

There are `C(N-1,k-1)` compositions with `k` blocks, each having `k-1`
possible crossing gaps.  Thus a literal face-by-face implementation has

\[
 \sum_{k=2}^N(k-1){N-1\choose k-1}=(N-1)2^{N-2}             \tag{10}
\]

local strata.  The adjacent-block hierarchy below compresses this raw list.

## 3. Moving left when the top block grows

Two crossing types are elementary in every dimension.

* If only the largest knot lies above the section, the cap is a simplex and
  AM--GM bounds it by the last beta step cap.
* If only the smallest knot lies below the section, write the positive
  interpolation ratios as `b_1,...,b_(N-1)`.  Knot order makes these ratios
  nonincreasing, whereas centroid order makes them nondecreasing.  They are
  therefore all equal, the upper knots form one block, and the cap is exactly
  the first beta step cap.

For every other crossing, let `L` be the number of coordinates below the
section and let `m` be the multiplicity of the largest knot.  The open
unique-top stratum has

\[
 m=1,\qquad 2\le L\le N-2,                                 \tag{11}
\]

giving `N-3` mixed families.  If a closed-form or polynomial proof on each
family remains valid when lower and nonlargest upper knots collide, it covers
every face whose top block is a singleton.

When the terminal comparison degenerates because the top knot has
multiplicity `m>=2`, move left to the boundary immediately preceding that
top block.  The remaining canonical pairs are

\[
 2\le m\le N-2,
 \qquad 2\le L\le N-m.                                    \tag{12}
\]

For fixed `m`, (12) contains `N-m-1` families, so the repeated-top list has

\[
 \sum_{m=2}^{N-2}(N-m-1)={N-2\choose2}                     \tag{13}
\]

members.  Together, (11)--(13) leave

\[
 (N-3)+{N-2\choose2}=\frac{N(N-3)}2                        \tag{14}
\]

canonical non-elementary families. A certificate on the closed ordered
parameter domain covers its lower collision faces only if its target is an
active-prefix bound and its hypotheses remain valid when order factors
vanish.  This
is the precise sense in which the proof moves left: the active order
constraint is the comparison between the top block average and the block
immediately to its left.

Equation (14) is a proof-obligation count, not a theorem that those
inequalities are true.  In four coordinates it gives two families; in five
coordinates it gives five. The tetrahedral certificates establish the
active-prefix bounds; the five-coordinate first-prefix certificates leave
the repeated-lowest obligation.

## 4. The six-coordinate first-prefix families

For `N=6`, translation and positive scaling put the nine canonical families
in the following forms.  All unnamed positive or negative knots are strictly
ordered; their collision limits belong to the same closed family.

| top multiplicity `m` | coordinates below `L` | normalized centered knots | first-prefix calculation |
|---:|---:|---|---|
| 1 | 2 | `(-1,-u,ta,tb,tc,t)` | computationally certified on the stated polynomial domain |
| 1 | 3 | `(-1,-u,-v,ta,tb,t)` | computationally certified on the stated polynomial domain |
| 1 | 4 | `(-1,-u,-v,-w,ta,t)` | computationally certified on the stated polynomial domain |
| 2 | 2 | `(-1,-u,ta,tb,t,t)` | computationally certified on the stated polynomial domain |
| 2 | 3 | `(-1,-u,-v,ta,t,t)` | computationally certified on the stated polynomial domain |
| 2 | 4 | `(-1,-u,-v,-w,t,t)` | computationally certified on the stated polynomial domain |
| 3 | 2 | `(-1,-u,ta,t,t,t)` | computationally certified on the stated polynomial domain |
| 3 | 3 | `(-1,-u,-v,t,t,t)` | computationally certified on the stated polynomial domain |
| 4 | 2 | `(-1,-u,t,t,t,t)` | computationally certified on the stated polynomial domain |

The five-coordinate certificates suggest the stronger common target

\[
 p\le I_{\gamma_0}(1,5)=1-(1-\gamma_0)^5.                 \tag{15}
\]

The exact theorem in
[`SIX-COORDINATE-REPEATED-TOP.md`](SIX-COORDINATE-REPEATED-TOP.md)
proves (15) in all six repeated-top rows.  It uses a common low-degree
intercept formula for `42`, `33`, and `24`, and exact confluent-interpolation
division identities for `23` and `22`; the `32` row is closed by a quadratic
scale constraint and two projective charts.  The exact theorem in
[`SIX-COORDINATE-UNIQUE-TOP.md`](SIX-COORDINATE-UNIQUE-TOP.md) proves (15)
for the unique-top row with `L=4`, using a central-cluster projective
decomposition. The companion theorem in
[`SIX-COORDINATE-UNIQUE-TOP-L2-L3.md`](SIX-COORDINATE-UNIQUE-TOP-L2-L3.md)
proves a stronger sharp-minorant bound for the `L=2` and `L=3` rows. These
are first-prefix bounds. When the smallest knot repeats, that prefix is
not a block boundary, so these results do not establish (7). The global
six-coordinate monotone cap statement remains conjectural.

## 5. Interface with the `n = 5` Stringer problem

The separate Clopper--Pearson calculation in
[`N5-MONOTONE-RANGE.md`](N5-MONOTONE-RANGE.md) proves that the six Stringer
threshold weights are nondecreasing exactly for

\[
 0<\alpha\le\alpha_5=0.3141689898050253288\ldots .          \tag{16}
\]

At each step vector, beta--binomial calibration makes the cap exactly
`alpha`. If the missing active-prefix inequalities are established, the
composition criterion will yield pointwise Stringer--Gaffke domination and
distribution-free Stringer coverage at `n=5` for every nominal confidence
of at least

\[
 1-\alpha_5=68.5831010194974671\ldots\%.                    \tag{17}
\]

Above `alpha_5`, the terminal two-knot edge already violates pointwise
domination. Thus (17) is a necessary endpoint, and would be sharp if the
conjectural positive comparison is proved; it is
not an undercoverage statement beyond the endpoint.

## 6. Mathematical yield

The composition formulation separates reusable geometry from dimension-
specific algebra.  Stationarity and face induction are proved once.  A new
dimension then asks for a finite list of local Dirichlet-block inequalities,
with an explicit `N(N-3)/2` canonical count rather than an informal promise
to inspect every collision.

It also explains the pattern seen in the five-coordinate proof.  A terminal
centroid comparison controls unique-top sections.  When the top knot repeats,
the comparison at the next block boundary supplies both the compactness
constraint and the multiplier inequality.  Increasing top multiplicity moves
that boundary left one block at a time.  In the final `32` row the same
comparison yields the sharper law `t<=3(u/t)^2`; projective charts then turn
the singular small-scale corner into compact cubes.

Finally, the hierarchy cleanly identifies the verification boundary. The
composition theorem, elementary crossings, canonical count, `n=5` weight
range, and the stated first-prefix polynomial inequalities have their
separate proofs. The active-prefix inequalities on repeated-lowest faces
remain unresolved. The `L=2,3` calculation combines exact rational
Bernstein covers with one rigorously enclosed binary64 chart whose
forward-error bounds are exact rational numbers.
