# Six-coordinate repeated-top first-prefix inequalities

## Status and theorem

Let

\[
D=(D_0,\ldots,D_5)\sim\operatorname{Dirichlet}(1,\ldots,1),
\qquad x_0\leq\cdots\leq x_5,
\]

and consider a regular section `x dot D=0`.  Write

\[
p=\Pr\{x\mathbin{\cdot}D>0\},\qquad
\gamma_i=\mathbb E(D_i\mid x\mathbin{\cdot}D=0).
\]

**First-prefix statement on the certified domains.** Suppose that the section centroid is ordered,
`gamma_0<=...<=gamma_5`, and the largest knot is repeated.  Then

\[
p\leq I_{\gamma_0}(1,5)=1-(1-\gamma_0)^5                 \tag{1}
\]

on each of the following six canonical families:

\[
\begin{array}{c|c}
\text{type}&\text{normalized knots}\\ \hline
42&(-1,-u,-v,-w,t,t)\\
33&(-1,-u,-v,t,t,t)\\
24&(-1,-u,t,t,t,t)\\
23&(-1,-u,ta,t,t,t)\\
22&(-1,-u,ta,tb,t,t)\\
32&(-1,-u,-v,ta,t,t).
\end{array}                                                   \tag{2}
\]

The statement includes limits within the polynomial domains used below.
At a collision where a divided-out order factor vanishes, centroid
orderedness alone need not imply the remaining polynomial constraint; such
faces require a separate argument. Moreover, a first-prefix bound does not
give an active-prefix bound when the least knot repeats. See
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md). Its ordinary reductions and exact rational polynomial
certificates are replayed by

```sh
make six-coordinate-repeated-top-check
```

These calculations address the six repeated-top first-prefix families.
The companion unique-top calculations address the other three families.
The global six-coordinate monotone cap theorem does not yet follow: its
active-prefix hypotheses on repeated-lowest faces remain unresolved.

## 1. One positive knot block

The `42`, `33`, and `24` families admit a common formula.  Scale the top knot
to one, write the `L=6-m` lower knots as

\[
-\frac{1-r_i}{r_i},\qquad
0<r_0\leq\cdots\leq r_{L-1}<1,
\]

and let `m` be the top multiplicity.  If `h_j` denotes the complete
homogeneous symmetric polynomial of degree `j`, the cap probability is

\[
P(r)=\left(\prod_{i=0}^{L-1}r_i\right)
h_{m-1}(1,1-r_0,\ldots,1-r_{L-1}).                         \tag{3}
\]

Put `P_i=partial P/partial r_i` and

\[
Q=\sum_i r_iP_i.
\]

Differentiating the cap with respect to its knots gives

\[
\gamma_i=\frac{r_i^2P_i}{Q}\quad(0\leq i<L),
\qquad
\gamma_{\rm top}
=\frac{\sum_i r_i(1-r_i)P_i}{mQ}.                          \tag{4}
\]

Consequently, the last required centroid comparison is the polynomial
inequality

\[
G_{\rm raw}:=
\sum_i r_i(1-r_i)P_i-mr_{L-1}^2P_{L-1}\geq0.              \tag{5}
\]

The desired margin has the same sign as

\[
T_{\rm raw}
=Q^5-(Q-r_0^2P_0)^5-PQ^5.                                 \tag{6}
\]

The order-preserving coordinates

\[
r_i=\prod_{j=i}^{L-1}z_j,
\qquad 0\leq z_j\leq1,                                    \tag{7}
\]

turn both expressions into polynomials on a cube.  After their greatest
positive monomial factors are removed, their exact sizes are small:

| type | `m` | degree of `G` | terms of `G` | degree of `T` | terms of `T` |
|---|---:|---:|---:|---:|---:|
| `42` | 2 | `(2,2,2,2)` | 15 | `(9,9,9,9)` | 505 |
| `33` | 3 | `(3,3,3)` | 20 | `(14,14,14)` | 660 |
| `24` | 4 | `(4,4)` | 15 | `(19,19)` | 210 |

For each row the exact calculation constructs

\[
T=S G+H,                                                    \tag{8}
\]

where every tensor-Bernstein coefficient of `S` and `H` is nonnegative.
The sparse support of `S` is fixed in the repository, but every rational
coefficient is reconstructed by solving a square system exactly.  The full
sign counts are:

| type | positive coefficients of `S` | positive/zero coefficients of `H` | smallest positive coefficient of `S` | smallest positive coefficient of `H` |
|---|---:|---:|---:|---:|
| `42` | 415 | 9505 / 495 | `54/49` | `1/326592` |
| `33` | 549 | 2588 / 787 | `1576/38115` | `234361/6358393737696` |
| `24` | 126 | 154 / 246 | `57/455455` | `3025/315490896` |

There are no negative coefficients.  Since Bernstein basis functions are
nonnegative on the cube, (5) and (8) imply (6), and hence (1).

The useful point is not merely that three cases can be computed.  Formula
(3) replaces high-degree confluent limits by a single expression valid for
every one-positive-block face and every dimension.  It is therefore a
reusable starting point for the corresponding rows of the hierarchy beyond
six coordinates.

## 2. The mixed `23` and `22` faces

The remaining two proved rows have an additional positive knot.  They are
derived uniformly from confluent divided differences.  For distinct block
knots `y_b`, multiplicities `m_b`, and positive blocks `b`, let

\[
\mathcal D_k(y;m)=
\sum_{b:\,y_b>0}\frac{1}{(m_b-1)!}
\left.
\frac{d^{m_b-1}}{dz^{m_b-1}}
\frac{z^k}{\prod_{c\ne b}(z-y_c)^{m_c}}
\right|_{z=y_b}.                                           \tag{9}
\]

Then `p=mathcal D_5`, the section density is `f=5 mathcal D_4`, and

\[
\Gamma_b=\frac1f\frac{\partial p}{\partial y_b}            \tag{10}
\]

is the section-centroid mass of block `b`.  Its individual centroid value is
`Gamma_b/m_b`.  Substituting the two rows in (2), clearing denominators, and
removing positive factors gives an order polynomial `G` and a target
polynomial `T` such that orderedness implies `G>=0` and (1) is equivalent to
`T>=0`.

Order-preserving substitutions are

\[
\begin{array}{c|c|c}
\text{type}&\text{substitution}&\text{scale confinement}\\ \hline
23&u=U,\ a=A,\ t=\tfrac23\tau&t<\tfrac23,\\
22&u=U,\ a=AB,\ b=B,\ t=\tau&t<1.
\end{array}                                                \tag{11}
\]

The confinement follows directly from `G>=0`: after writing
`t=c+rho`, the mixed Bernstein/power coefficients of `-G` are all
nonnegative and are strict in the relative interior.  There are respectively
92 positive and 34 zero coefficients for `23`, and 196 positive and 140 zero
coefficients for `22`.

Exact lexicographic division on the cube again gives `T=S G+H`.  Rather than
materializing millions of rational Bernstein fractions, the calculation
uses the equivalent multihomogenization

\[
\widetilde R(z)=
\prod_i(1+z_i)^{d_i}
R\!\left(\frac{z_1}{1+z_1},\ldots,
         \frac{z_k}{1+z_k}\right).                         \tag{12}
\]

The coefficient of `z^nu` in (12) is the tensor-Bernstein coefficient of
`R` multiplied by `product_i binomial(d_i,nu_i)`, a positive integer.  Thus
the two arrays have identical signs.  The exact division data are:

| type | degree/terms of `G` | degree/terms of `T` | degree/terms of `S` | positive/zero coefficients of `S` | degree/terms of `H` | positive/zero coefficients of `H` |
|---|---|---|---|---:|---|---:|
| `23` | `(5,2,6)` / 32 | `(19,11,28)` / 2969 | `(14,38,108)` / 18539 | 59659 / 4106 | `(4,40,114)` / 14220 | 19615 / 3960 |
| `22` | `(5,1,3,6)` / 50 | `(19,11,22,33)` / 29801 | `(14,24,62,113)` / 305193 | 2524950 / 168300 | `(4,25,65,119)` / 291265 | 850461 / 179139 |

Again, every omitted count is zero: neither quotient nor remainder has a
negative coefficient.  This proves (1) in both mixed families.

## 3. The final `32` face

It remains to treat

\[
 (-1,-u,-v,ta,t,t),\qquad 1\ge u\ge v>0,\quad 0<a<1.
                                                               \tag{13}
\]

Write `u=U`, `v=UV`, `a=A`, and use `q` for the physical scale `t`.
Confluent divided differences give four adjacent centroid-gap polynomials
`G_0,...,G_3` and a comparison polynomial `T`.  All cleared denominators
are positive; on the relative interior, orderedness is `G_i>=0`, and
(1) is equivalent to `T>=0`.
Their degrees are

| polynomial | degree in `(U,V,A,q)` | power terms |
|---|---:|---:|
| `T` | `(26,13,17,35)` | 40,732 |
| `G_0` | `(5,2,4,7)` | 102 |
| `G_1` | `(6,3,4,7)` | 102 |
| `G_2` | `(7,4,5,8)` | 169 |
| `G_3` | `(8,4,4,9)` | 213 |

The last gap first confines the scale.  Substituting `q=1+rho` in `-G_3`
gives 1,347 positive and 903 zero mixed Bernstein/power coefficients, with
none negative.  Hence an ordered regular section has `q<1`.

Split next according to `q<=U` or `U<=q`.  In the first chart put `q=UZ`.
After positive monomial factors are removed, `T` has degree
`(19,13,17,35)`.  A dyadic cover has 17 nodes and nine leaves: seven leaves
have `T>=0` coefficientwise, while the other two have `-G_3>=0` strictly in
their relative interiors and therefore contain no ordered section.  This
settles the first chart.

The second chart contains the main structural reduction.  Put `U=qZ`.  The
exact identity obtained from `G_3` satisfies

\[
 -G_3(3Z^2+\rho,V,A,Z)\ge0
 \quad(0\le Z,V,A,\rho\le1),                              \tag{14}
\]

because all 2,125 tensor-Bernstein coefficients are nonnegative: 1,630 are
positive and 495 vanish.  It follows that `G_3>=0` forces

\[
 q\le 3Z^2.                                                \tag{15}
\]

Write `q=3Z^2R`.  Two projective charts cover the remaining cube: `A=ZB`
when `A<=Z`, and `Z=AB` when `Z<=A`.  Their exact dyadic decompositions are:

| chart | degree of `T` | nodes / leaves | leaf closures |
|---|---:|---:|---|
| `A=ZB` | `(65,17,13,19)` | 127 / 64 | 35 by `T`, 7 by `T-mu G_3`, 18 by `-G_3`, 4 by `-G_2` |
| `Z=AB` | `(65,62,13,19)` | 149 / 75 | 22 by `T`, 14 by `T-mu G_3`, 7 by `T-mu G_2`, 19 by `-G_3`, 13 by `-G_2` |

Every `mu` in the table is an explicit nonnegative integer, and every
coefficient check is exact.  On a multiplier
leaf, `T-mu G_i>=0` together with `G_i>=0` gives `T>=0`.  A leaf certified by
`-G_i` has no ordered point in its relative interior.  The leaf paths and
multipliers form a complete prefix-free cover, and continuity supplies the
interfaces and collision limits.  This proves (1) for (13).

## 4. What the theorem changes

Exact identities establish the first-prefix comparisons on the six stated
polynomial domains. Separate projective arguments establish the three
unique-top comparisons. These do not eliminate the active-prefix
obstruction on repeated-lowest faces.

The certificates reveal a structural split.  One positive
block has the low-degree formula (3), whereas an intervening positive knot
is naturally handled by confluent interpolation and division modulo the
active centroid constraint.  The `32` proof adds a third reusable mechanism:
an active centroid gap can impose a quadratic scale law, after which
projective charts resolve the singular corner that defeats a direct cube
certificate.

The `n=5` confidence continuum down to `68.58310101949746...%` remains
conditional on the active-prefix collision inequalities. The exact
Clopper--Pearson weight range and terminal-edge obstruction are independent
of that unresolved implication.
