# The ordered-simplex cap problem behind all-`n` Stringer

## Status

This note isolates a dimension-free geometric statement that would prove
ordinary binomial Stringer conservative at conventional confidence levels
for every sample size.  The reduction and all vertex equalities below are
exact.  The complete cap inequality remains **open**.  One entire region is
now proved for every sample size at nominal confidence of at least
`75%`; see
[`ONE-CAP-COMPARISON.md`](ONE-CAP-COMPARISON.md). A certified two-cap collar
strictly extends that region through `n=200` at 90%, 95%, and 99%; see
[`TWO-CAP-COMPARISON.md`](TWO-CAP-COMPARISON.md). Separate fixed-dimension
certificates prove the complete comparison through `n=7` at 90%, 95%, and
99%; see
[`N7-CONVENTIONAL.md`](N7-CONVENTIONAL.md).
The four-coordinate weight-independent theorem proves that every cap with
nondecreasing barycentric weights is bounded by its largest step cap. It
extends the `n=3` Stringer result to every
`0<alpha<=((19+sqrt(21))/34)^3`; see
[`TETRAHEDRAL-VERTEX-BARRIER.md`](TETRAHEDRAL-VERTEX-BARRIER.md).

The analogous global theorems in five and six coordinates remain
conjectural. The local identities and the nine six-coordinate families give
their stated cap bounds, but on repeated-lowest-knot faces these are not yet
the active-prefix bounds required by stationary-face induction. See
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).
The exact Clopper--Pearson weight-order intervals through
`alpha_4=0.3214292869970077...` and
`alpha_5=0.3141689898050253288...` are proved separately, but their
conversion to full Stringer coverage intervals is conjectural:
[`N4-MONOTONE-RANGE.md`](N4-MONOTONE-RANGE.md) and
[`N5-MONOTONE-RANGE.md`](N5-MONOTONE-RANGE.md).
Independent direct certificates still prove the three conventional levels.

Nothing in this note should be cited as a general-`n` coverage theorem.

## 1. Pointwise comparison as a cap at its barycentric threshold

Let `p_j=p_n(j)` be the binomial Clopper--Pearson factors and define the
ascending-knot Stringer weights

\[
 c_0=1-p_{n-1},\qquad
 c_j=p_{n-j}-p_{n-j-1}\ (1\le j<n),\qquad
 c_n=p_0.                                               \tag{1}
\]

They are positive and sum to one.  For ordered sample knots

\[
 x_0\le\cdots\le x_{n-1}\le x_n=1,
\]

the binomial Stringer value is their `c`-barycenter

\[
 s_c(x)=\sum_{i=0}^n c_i x_i.                           \tag{2}
\]

If `D~Dirichlet(1,...,1)`, define the cap functional

\[
 \Phi_c(x)=\Pr_D\{x\mathbin{\cdot}D>s_c(x)\}.           \tag{3}
\]

The valid Gaffke endpoint is no larger than Stringer exactly when
`Phi_c(x)<=alpha`.  The functional is invariant under a common positive
affine transformation of the knots.  Except for the all-equal case, one may
therefore normalize to

\[
 0=x_0\le x_1\le\cdots\le x_{n-1}\le x_n=1.            \tag{4}
\]

This is an ordered `(n-1)`-simplex.  Proving

\[
 \Phi_c(x)\le\alpha\quad\text{throughout (4)}           \tag{5}
\]

would give pointwise binomial Stringer--Gaffke domination, hence ordinary
binomial Stringer coverage, for that `n` and `alpha`.  At
`alpha<exp(-1)`, the all-`n` Poisson-over-binomial factor theorem would
transfer the result to the audit-table Poisson factors.

## 2. Every nontrivial vertex has cap probability exactly `alpha`

The vertices of (4) are

\[
 v^{(r)}=(\underbrace{0,\ldots,0}_{r},
           \underbrace{1,\ldots,1}_{n+1-r}),
 \qquad r=1,\ldots,n.                                  \tag{6}
\]

Put `C_r=sum_(i<r)c_i`.  Equation (1) telescopes to

\[
 C_r=1-p_{n-r}.                                         \tag{7}
\]

At `v^(r)`, the random convex combination in (3) is
`1-Y_r`, where

\[
 Y_r=\sum_{i=0}^{r-1}D_i\sim\operatorname{Beta}(r,n+1-r),
\]

and the threshold is `1-C_r`.  Consequently,

\[
\begin{aligned}
 \Phi_c(v^{(r)})
 &=\Pr\{Y_r<C_r\}\\
 &=\Pr\{\operatorname{Bin}(n,C_r)\ge r\}\\
 &=\Pr\{\operatorname{Bin}(n,p_{n-r})\le n-r\}
 =\alpha.                                               \tag{8}
\end{aligned}
\]

The second equality is the beta--binomial identity; the third swaps
successes and failures; the last is the defining Clopper--Pearson equation.
Thus Stringer's factors place the barycentric threshold so that **every
nontrivial ordered-simplex vertex is exactly tight**.

The central geometric conjecture can now be stated without audit notation:

> For the barycentric vector (1) at conventional `alpha`, the cap functional
> (3) attains its maximum over the ordered knot simplex at a vertex.

By (8), that maximum would be `alpha`.  This statement is stronger than the
coverage conjecture because it asserts pointwise domination of a separately
valid confidence limit.

The conjecture is proved for four Dirichlet coordinates and nondecreasing
weights. The proof compares every stationary section with a beta cap at an
active boundary rank, then inducts over knot faces. Its only nonsimplex
section is a quadrilateral. Five-coordinate mixed-section identities and
six-coordinate projective charts prove local inequalities, but the
active-prefix implication on repeated-lowest-knot faces remains unresolved;
see [`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).
The local algebra alone does not establish either higher global theorem.

## 3. The region already removed from the open problem

When `s_c(x)>=x_(n-1)`, only the knot at one is above the threshold.  The
one-cap lemma proves

\[
 \Phi_c(x)\le
 \max_{1\le r\le n}\left(\frac{C_r}{C_r+c_n}\right)^r. \tag{9}
\]

For the weights (1), the `r=n` term in (9) is exactly `alpha`.  The analytic
binomial-tail comparison in `ONE-CAP-COMPARISON.md` proves every `r<n` term
at most `alpha` for all `n` whenever `0<alpha<=1/4`. Hence (5) is already
established on this complete region of (4) without a sample-size cutoff.
An older exact calculation through `n=200` remains as an independent finite
regression.

The adjacent two-upper-knot region also admits a dimension-free bound. If
`x_(n-2)<=s_c(x)<=x_(n-1)<1` and

\[
q=\frac{1-x_{n-1}}{1-s_c(x)},
\]

then `TWO-CAP-COMPARISON.md` proves

\[
\Phi_c(x)\le \frac{1-(1-q)^n}{q}
\max_{1\le r<n}
\left(\frac{C_r}{C_r+c_n+c_{n-1}(1-q)}\right)^r.     \tag{9a}
\]

Each term has no interior maximum as `q` varies over a collar ending at one.
Exact factor certificates therefore prove complete zero-uplift collars
through `n=200`: `q>=5/12`, `1/3`, and `1/4` at 90%, 95%, and 99%,
respectively. This does not settle the residual two-cap region.

## 4. An exact adjacent-transfer formula for the remaining regions

The following identity turns a possible vertex-maximization proof into a
one-dimensional single-crossing problem.

Fix two adjacent knots `x_i<=x_(i+1)`, hold all other knots fixed, and
preserve their contribution to the threshold.  Write

\[
 q=c_i+c_{i+1},\qquad
 a=\frac{c_i}{q},\qquad b=\frac{c_{i+1}}q=1-a,
\]

and parameterize

\[
 m=ax_i+bx_{i+1},\qquad d=x_{i+1}-x_i,
\]

so that

\[
 x_i=m-bd,\qquad x_{i+1}=m+ad.                          \tag{10}
\]

The Stringer threshold is constant as `d` varies.  Put

\[
 R=D_i+D_{i+1},\qquad W=\frac{D_{i+1}}R.
\]

Dirichlet neutrality gives `W~Unif(0,1)`, independent of `R` and the
remaining coordinates.  If `O=sum_(k notin {i,i+1})x_kD_k`, and `t` is the
fixed threshold, define

\[
 Z=\frac{O+Rm-t}{R}.                                    \tag{11}
\]

The law of `Z` does not depend on `d`.  Conditional on everything except
`W`, the cap probability is therefore

\[
 g(d)=\mathbb E\!\left[\operatorname{clip}
          \left(a+\frac Zd,0,1\right)\right],\qquad d>0.\tag{12}
\]

At differentiability points, or in Stieltjes form when atoms are present,

\[
 g'(d)=-\frac1{d^2}
 \int_{-ad}^{bd} z\,dF_Z(z).                            \tag{13}
\]

Equations (12)--(13) are exact.  They expose the remaining proof target:
show that the truncated first moment in (13) changes sign at most from
positive to negative over the admissible gap interval.  Then `g'` changes
at most from negative to positive, so `g` has no interior maximum and its
maximum is at an endpoint.  Repeating such adjacent endpoint moves, with a
careful termination argument on the ordered simplex, would reduce (5) to
the vertex equalities (8).

That single-crossing statement has **not** been proved.  In particular, no
generic claim is made for arbitrary barycentric weights or arbitrary laws
of `Z`; the distribution in (11) retains special Dirichlet and ordering
structure that a successful proof must use.

## 5. Guardrails and exact next steps

The conventional-confidence restriction cannot simply be deleted.  For
example, at `n=2` and `alpha=0.70`, the nonterminal expression in (9) is
approximately `0.73467`, above `alpha`, even though the paper's separate
`n=2` theorem still proves ordinary Stringer coverage.  Thus pointwise
Gaffke domination and coverage are genuinely different claims.

The structural questions are:

1. prove an active-prefix cap bound on the repeated-lowest-knot faces in
   five and six coordinates, or obtain the proposed global inequalities
   without that induction;
2. retain the proved weight-order intervals separately from their
   conjectural coverage consequences; and
3. use the terminal-edge and one-/two-cap results as boundary information
   for the remaining higher-dimensional cap regions.

More unconstrained grid search over taint distributions does not address
these geometric targets.
