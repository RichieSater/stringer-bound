# The sharp monotone-weight range at `n = 5`

## Status

Let `p_j=p_5(j)` be the one-sided Clopper--Pearson upper limits for five
trials.  This note proves that the ascending-knot Stringer threshold weights

\[
 c=(1-p_4,\ p_4-p_3,\ p_3-p_2,\ p_2-p_1,\ p_1-p_0,\ p_0)              \tag{1}
\]

are nondecreasing exactly through the terminal critical level

\[
 \alpha_5=r_5^5=0.3141689898050253288\ldots,                          \tag{2}
\]

where `r_5` is the unique root in `(1/2,9/10)` of

\[
 129r^4-271r^3+209r^2-71r+9=0.                                      \tag{3}
\]

Equivalently, the candidate nominal-confidence endpoint is

\[
 1-\alpha_5=0.6858310101949746711\ldots .                             \tag{4}
\]

The weight theorem is proved. The six-coordinate global cap bound remains
conjectural because its active-prefix collision inequalities are unresolved;
see [`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md). Thus domination and coverage throughout (2) are conditional. The
existing fixed-level certificate independently proves the conclusion at
90%, 95%, and 99%.

## 1. The five spacing inequalities

Write

\[
\begin{aligned}
 d_0&=2p_4-p_3-1,       &d_1&=2p_3-p_2-p_4,\\
 d_2&=2p_2-p_1-p_3,     &d_3&=2p_1-p_0-p_2,\\
 d_4&=2p_0-p_1.
\end{aligned}                                                        \tag{5}
\]

The vector (1) is nondecreasing exactly when every `d_i` is nonnegative.
The terminal-edge theorem gives

\[
 d_4\ge0\quad\Longleftrightarrow\quad\alpha\le\alpha_5.             \tag{6}
\]

Indeed, put `r=alpha^(1/5)`, so `p_0=1-r`.  Substitution of
`p_1=2p_0` into its defining binomial tail gives

\[
 (2r-1)^4(9-8r)-r^5
 =-(r-1)(129r^4-271r^3+209r^2-71r+9).                                \tag{7}
\]

Exact Sturm counting gives one root of the quartic in `(1/2,9/10)`, and

\[
 7932/10000<r_5<7933/10000.                                          \tag{8}
\]

Binomial-tail symmetry says

\[
 p_j(1-\alpha)=1-p_{4-j}(\alpha),                                    \tag{9}
\]

and hence

\[
 d_0(\alpha)=-d_4(1-\alpha),\qquad
 d_1(\alpha)=-d_3(1-\alpha),\qquad
 d_2(\alpha)=-d_2(1-\alpha).                                       \tag{10}
\]

It remains to control `d_3` on the low and reflected high ranges and `d_2`
below one half.

## 2. The next-to-terminal spacing

Put `x=p_0` and `y=p_1`.  If `d_3=0`, then `p_2=2y-x` and the three common
binomial lower tails satisfy

\[
 F_0(x)=F_1(y)=F_2(2y-x)=\alpha,                                     \tag{11}
\]

where `F_j(z)=P(Bin(5,z)<=j)`.  Eliminating `y` exactly gives

\[
 400000x^2(x-1)^{12}P_{11}(x)=0,                                    \tag{12}
\]

with

\[
\begin{aligned}
P_{11}(x)={}&682456x^{11}-8871928x^{10}+53222168x^9
 -193600856x^8\\
&+472259515x^7-807107237x^6+980133686x^5
 -838362386x^4\\
&+487990767x^3-179494185x^2+35660800x-2508800.
\end{aligned}                                                        \tag{13}
\]

An exact Sturm count shows that `P_11` has only one root in `(0,1)`, and it
lies in

\[
 1356/10000<x<1357/10000.                                            \tag{14}
\]

For `alpha<=alpha_5`, equations (8) and `x=1-alpha^(1/5)` give

\[
 x>2067/10000,                                                        \tag{15}
\]

so (12)--(14) exclude a zero.  The sign is fixed by the exact rational
witness `alpha=(3/4)^5`, where `p_0=1/4`.  Direct tail evaluations give

\[
 p_1>46/100,\qquad p_2<65/100,
\]

and therefore `d_3>1/50`.  It follows that

\[
 d_3>0\qquad(0<\alpha\le\alpha_5).                                   \tag{16}
\]

The same resultant controls the reflected range.  The bracket (8) gives
`alpha_5<1/3`; hence `beta>=1-alpha_5` implies `beta>2/3`.  Since
`F_0(1/10)=(9/10)^5<2/3`, the corresponding `p_0(beta)` is below `1/10`,
again outside (14).  At the exact witness `beta=(19/20)^5`, where
`p_0=1/20`, rational tail comparisons give

\[
 p_1<182/1000,\qquad p_2>344/1000,
\]

and thus `d_3<-3/100`.  Therefore

\[
 d_3(\beta)<0\qquad(1-\alpha_5\le\beta<1).                            \tag{17}
\]

By (10), equation (17) proves `d_1(alpha)>0` throughout the desired range.

## 3. The central spacing

Now put `x=p_1` and `y=p_2`.  If `d_2=0`, then `p_3=2y-x` and exact
elimination from

\[
 F_1(x)=F_2(y)=F_3(2y-x)                                             \tag{18}
\]

gives

\[
 1600000x^6(x-1)^6A_5(x)B_8(x)=0,                                   \tag{19}
\]

where

\[
 A_5(x)=8x^5-30x^4+40x^3-20x^2+1=2F_1(x)-1                          \tag{20}
\]

and

\[
\begin{aligned}
B_8(x)={}&1597696x^8-14379264x^7+57099936x^6-128866064x^5\\
 &+177551071x^4-148524550x^3+69152775x^2-13296000x-304000.
\end{aligned}                                                        \tag{21}
\]

Thus `A_5(x)=0` is exactly `alpha=1/2`.  The other factor has only one root
in `(0,1)`, lying in

\[
 8209/10000<x<821/1000.                                              \tag{22}
\]

That root is not on the admissible branch.  At a genuine equality in (18),
`q=2y-x=p_3<1`.  Because `F_2` is strictly decreasing, `q<1` forces

\[
 F_1(x)-F_2((1+x)/2)>0.                                               \tag{23}
\]

Exact simplification gives

\[
 F_1(x)-F_2((1+x)/2)
 =\frac{(x-1)^3(67x^2-39x-8)}{16}.                                   \tag{24}
\]

The positive root of the quadratic in (24) is below `3/4`, so (23) forces
`x<3/4`, contradicting (22).  Hence `d_2` has no zero in `(0,1/2)`.

At `alpha=(3/4)^5`, exact tail signs yield

\[
 p_1<463/1000,\qquad p_2>648/1000,\qquad p_3<813/1000,
\]

so `d_2>1/50`.  Consequently

\[
 d_2>0\qquad(0<\alpha<1/2),\qquad d_2(1/2)=0.                        \tag{25}
\]

## 4. Completion and Stringer consequence

For `0<alpha<=alpha_5`, equations (6), (10), (16), (17), and (25) give

\[
 d_0>0,\quad d_1>0,\quad d_2>0,\quad d_3>0,\quad d_4\ge0.             \tag{26}
\]

Thus (1) is nondecreasing exactly through (2).  Its prefix sums are
`C_r=1-p_(5-r)`, so every step cap has value

\[
 I_{C_r}(r,6-r)=\alpha.                                               \tag{27}
\]

The conjectural six-coordinate monotone cap bound would bound every ordered
knot cap by `alpha` in (27), implying pointwise Stringer--Gaffke domination
and coverage at every nominal confidence of at least (4). This continuum
remains conditional on the collision inequalities. Above `alpha_5`, the
last two weights reverse and the terminal edge yields a pointwise cap
violation. Thus (2) is a necessary comparison endpoint and the proved exact
endpoint for weight monotonicity.

## 5. Verification boundary

The binomial-tail identities, symmetry, continuity argument, domain
exclusions, and conditional composition with a global cap theorem are
ordinary mathematics.
The public script
[`n5_weight_monotonicity.py`](../computations/python/n5_weight_monotonicity.py)
recomputes both resultants over the integers, performs exact Sturm counts on
rational intervals, checks every boundary identity and sign witness, and
brackets the root of (3).  Its compact artifact is
[`n5-weight-monotonicity-certificate.json`](../computations/certificates/n5-weight-monotonicity-certificate.json).
Reproduce it with

```console
make n5-weight-monotonicity-check
```

No floating-point sign enters the proof.  The displayed decimals in
(2)--(4) are enclosed by exact rational intervals recorded in the artifact.

## 6. Mathematical yield

This calculation identifies the exact `n=5` monotone-weight range. A proof
of the active-prefix collision inequalities would turn it into the
confidence continuum starting at `68.58310101949746...%`. The continuum
is not established by the weight calculation alone.

The elimination also reveals a stable hierarchy among the Clopper--Pearson
spacings.  The terminal spacing is the first obstruction; the next spacing
has a single algebraic zero well outside the relevant domain, and the central
spacing can vanish only at symmetry level `alpha=1/2` on its admissible
branch.  This pattern is a concrete target for an all-`n` weight-monotonicity
theorem independent of the harder simplex geometry.
