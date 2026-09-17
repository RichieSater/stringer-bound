# The sharp monotone-weight range at `n = 4`

## Status

Let `p_j=p_4(j)` be the one-sided Clopper--Pearson upper limits for four
trials.  This note proves that the ascending-knot Stringer threshold weights

\[
 c=(1-p_3,\ p_3-p_2,\ p_2-p_1,\ p_1-p_0,\ p_0)             \tag{1}
\]

are nondecreasing exactly through the terminal critical level

\[
 \alpha_4=r_4^4=0.3214292869970077\ldots,                   \tag{2}
\]

where `r_4` is the unique root in `(1/2,7/8)` of

\[
 49r^3-79r^2+41r-7=0.                                      \tag{3}
\]

Combined with the conjectural five-coordinate global cap bound, this would
give pointwise Stringer--Gaffke domination and distribution-free Stringer
coverage at `n=4` for every nominal confidence of at least

\[
 1-\alpha_4=0.6785707130029922\ldots .                       \tag{4}
\]

The terminal ordered-knot edge proves failure of pointwise comparison above
this endpoint. Sharpness from below remains conditional on the
active-prefix collision inequalities in [`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md). The weight range itself is proved; the obstruction
does not assert Stringer undercoverage.

## 1. The four spacing inequalities

Write

\[
\begin{aligned}
 d_0&=2p_3-p_2-1, &d_1&=2p_2-p_1-p_3,\\
 d_2&=2p_1-p_0-p_2, &d_3&=2p_0-p_1.
\end{aligned}                                               \tag{5}
\]

The vector (1) is nondecreasing exactly when all four `d_i` are nonnegative.
The terminal-edge theorem gives

\[
 d_3\ge0\quad\Longleftrightarrow\quad\alpha\le\alpha_4.    \tag{6}
\]

Indeed, putting `r=alpha^(1/4)` and `p_0=1-r`, equality `p_1=2p_0` gives

\[
 (2r-1)^3(7-6r)=r^4
 =-(r-1)(49r^3-79r^2+41r-7).                               \tag{7}
\]

It remains to show that none of the other three spacings changes sign first.
Binomial-tail symmetry gives

\[
 p_j(1-\alpha)=1-p_{3-j}(\alpha),                            \tag{8}
\]

and therefore

\[
 d_0(\alpha)=-d_3(1-\alpha),\qquad
 d_1(\alpha)=-d_2(1-\alpha).                                \tag{9}
\]

Thus only `d_2` requires a separate analysis.

## 2. Exact elimination for the middle spacing

Put `x=p_0`, `y=p_1`, and `q=2y-x`.  The relevant binomial lower tails are

\[
 F_0(x)=(1-x)^4,\quad
 F_1(y)=(1-y)^3(1+3y),\quad
 F_2(q)=(1-q)^2(1+2q+3q^2).                                \tag{10}
\]

We have `F_0(x)=F_1(y)=alpha`, and `d_2=q-p_2`.  If `q>=1`, then `d_2>=0`
automatically.  If `q<1` and `d_2=0`, then

\[
 F_0(x)=F_1(y)=F_2(2y-x).                                  \tag{11}
\]

Eliminating `y` from the two polynomial differences in (11) gives

\[
 144x^2(x-1)^6P(x)=0,                                      \tag{12}
\]

where

\[
\begin{aligned}
P(x)={}&43209x^8-432090x^7+1912545x^6-4833544x^5\\
 &+7534696x^4-7294848x^3+4152976x^2-1188352x+110592.
\end{aligned}                                               \tag{13}
\]

Exact Sturm counts show that `P` has exactly two roots in `(0,1)`, one in

\[
 (411/2500,329/2000)=(0.1644,0.1645)                        \tag{14}
\]

and one in

\[
 (139/200,6951/10000)=(0.695,0.6951).                       \tag{15}
\]

The second interval is an algebraic branch on which `q>1`, not a possible
zero of `d_2` in the nontrivial region.  Precisely, comparison with `q=1`
gives the identity

\[
 F_0(x)-F_1\left(\frac{1+x}{2}\right)
 =\frac{(x-1)^3(19x-11)}{16}.                              \tag{16}
\]

Since `F_1` is strictly decreasing, along `F_0(x)=F_1(y)` we have
`2y-x<1` exactly when `x<11/19`.

Now suppose `alpha<=alpha_4`.  The exact bracket
`7529/10000<r_4<753/1000` gives

\[
 x=1-\alpha^{1/4}\ge1-r_4>247/1000.                        \tag{17}
\]

In the only nontrivial region `q<1`, equations (16)--(17) place `x` in
`(247/1000,11/19)`.  Neither root interval (14)--(15) meets that interval, so
`d_2` has no zero there.  At the boundary `q=1`, it equals `1-p_2>0`.
Consequently

\[
 d_2>0\qquad(0<\alpha\le\alpha_4).                          \tag{18}
\]

For later use, `d_2<0` throughout `1/2<=alpha<1`.  Indeed,
`alpha>=1/2` gives

\[
 x\le1-2^{-1/4}<4/25<411/2500,                              \tag{19}
\]

so the same resultant excludes a zero.  The sign near one follows directly
by setting `epsilon=1-alpha`: the defining tails give

\[
 p_0\sim\epsilon/4,\qquad
 p_1\sim(\epsilon/6)^{1/2},\qquad
 p_2\sim(\epsilon/4)^{1/3},                                \tag{20}
\]

and hence `d_2=2p_1-p_0-p_2<0` for sufficiently small positive `epsilon`.
Continuity and the absence of zeros prove the claim.

## 3. Completion and Stringer consequence

Because `alpha_4<1/2`, equations (9), (18), and (6) give

\[
 d_0>0,\qquad d_1>0,\qquad d_2>0,\qquad d_3\ge0
 \quad\text{for }0<\alpha\le\alpha_4.                       \tag{21}
\]

Thus (1) is nondecreasing.  Its prefix sums are
`C_r=1-p_{4-r}`.  At every step vector the beta--binomial identity gives

\[
 I_{C_r}(r,5-r)=\alpha.                                     \tag{22}
\]

The conjectural five-coordinate monotone cap bound in
[`FIVE-COORDINATE-VERTEX-BARRIER.md`](FIVE-COORDINATE-VERTEX-BARRIER.md)
would bound every ordered-knot cap by `alpha`. Its active-prefix collision
obligation remains unresolved, so the continuum of pointwise comparisons
and the resulting coverage continuum are conditional.

When `alpha>alpha_4`, equation (6) reverses the last two weights.  The explicit
terminal edge has positive outward derivative and produces an ordered sample
with Stringer below the Gaffke limit. Hence (2) is a necessary endpoint for
pointwise comparison and the proved exact endpoint for weight monotonicity.

## 4. Verification boundary

The binomial-tail identities, symmetry, sign propagation, asymptotics, and
conditional composition with a global cap theorem are ordinary mathematics.  The public
script
[`n4_weight_monotonicity.py`](../computations/python/n4_weight_monotonicity.py)
uses exact integer polynomial arithmetic to recompute the resultant (12),
uses exact Sturm counts on the rational intervals (14)--(15), checks (16),
and brackets the root of (3).  Its compact artifact is
[`n4-weight-monotonicity-certificate.json`](../computations/certificates/n4-weight-monotonicity-certificate.json).
Reproduce it with

```console
make n4-weight-monotonicity-check
```

No floating-point sign enters the proof. The displayed decimal in (2) is
certified by the exact rational interval recorded in the artifact.
