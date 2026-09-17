# Systematic PPS: exact failure and a sharp design safeguard

## Mathematical status and placement

**Status.** The results in this note are proved. The finite enumerations and
displayed factor comparisons are exactly certified by
computations/python/systematic_pps.py and
computations/certificates/systematic-pps-certificate.json.

**One-sentence thesis.** A single random start makes the average taint in a
systematic probability-proportional-to-size (PPS) sample design-unbiased, but
does not make its $n$ taints independent: ordinary Stringer can consequently
have severely subnominal coverage, while exact finite-frame randomization
inversion gives a design-valid, identity-aware Stringer overlay.

The principal new ideas are:

1. a two-item-per-cell periodic population giving an exact two-point law for
   the complete systematic sample and disproving binomial Stringer for every
   $n\ge 2$ and every confidence level;
2. a diagonal impossibility theorem showing that every design-valid rule
   seeing only the taint vector must be at least
   $1-\alpha+\alpha y$ on the constant sample $(y,\ldots,y)$; and
3. a design-valid overlay combining that sharp taint-only floor with a
   deterministic completion bound based on sampled item identities and book
   weights; and
4. an exact least-favourable-completion inversion that uses the entire known
   phase frame, is never larger than either closed-form fallback, and reduces
   to finitely many rational linear programs.

Exact phase enumeration and factor bracketing are routine verification, not
the conceptual contribution. The reusable conclusion is that the effective
distribution-free sample size of a one-start systematic grid can collapse to
one; independent starts, not the number of grid points, restore independent
replication.

These results form the companion manuscript
`../paper/systematic-pps.tex`. The separation is mathematically warranted:
the i.i.d. Stringer paper is organized around uniform-simplex cap geometry,
whereas the companion fixes an ordered finite population and studies its
random-start phase distribution. The finite-frame inversion, disjoint-phase
coalition formula, exact dynamic program, complexity boundary, and
independent-start law now provide a common structural theory beyond the
initial design counterexample. The two manuscripts remain linked by the
ordinary Stringer calculation, but neither proof architecture depends on the
other.

## 1. Finite-population estimand and exact design

Let an ordered finite population consist of items $i=1,\ldots,N$, with
positive recorded amounts $w_i$ and taints $t_i\in[0,1]$. Write

$$
W=\sum_{i=1}^N w_i,
\qquad
\theta=\frac{1}{W}\sum_{i=1}^N w_i t_i.
$$

Thus $W\theta$ is total overstatement and $\theta$ is the overstatement
fraction of recorded amount. Put the items consecutively on the book-value
line $[0,W)$, and let $\tau(u)=t_i$ on item $i$'s interval.

For intended sample size $n$, let $h=W/n$, draw one random start
$R\sim\operatorname{Unif}[0,h)$, and observe

$$
T_k=\tau(R+kh),\qquad k=0,\ldots,n-1.
$$

This is normalized one-start systematic PPS, or systematic MUS. Endpoints of
item intervals have probability zero and may be assigned by any fixed rule.
The exact integer-book-unit analogue chooses the random start uniformly from
the integer positions in one interval. Items larger than $h$ may be hit more
than once; none of the counterexamples below has such an item.

### Proposition 1 (exact design unbiasedness)

For every ordered finite population,

$$
\mathbb E_R\!\left[\bar T(R)\right]=\theta,
\qquad
\bar T(R)=\frac1n\sum_{k=0}^{n-1}T_k.
$$

**Proof.** Translation in the $k$-th systematic cell gives

$$
\mathbb E_R[\bar T(R)]
=\frac1n\sum_{k=0}^{n-1}\frac1h
  \int_0^h\tau(r+kh)\,dr
=\frac1W\int_0^W\tau(u)\,du
=\theta.\qquad\square
$$

This identity proves only unbiasedness of the average. It supplies neither
independence nor a binomial count law.

## 2. A periodic counterexample for every binomial case

Fix $n\ge2$, $0\le y<1$, and $0<q<1$. In each systematic cell put a
leading item of weight $qh$ and taint one, followed by an item of weight
$(1-q)h$ and taint $y$. Repeat this same two-item pattern in all $n$
cells. Then

$$
(T_0,\ldots,T_{n-1})=
\begin{cases}
(1,\ldots,1),&\text{with probability }q,\\
(y,\ldots,y),&\text{with probability }1-q,
\end{cases}
\qquad
\theta=y+q(1-y).
$$

Let $p_n(j)$ be the binomial upper confidence factor at tail level
$\alpha$. In particular,

$$
p_n(0)=1-\alpha^{1/n},\qquad p_n(n)=1.
$$

On a constant sample $(y,\ldots,y)$, telescoping the Stringer increments
gives

$$
\operatorname{SB}_{\rm B}(y,\ldots,y)
=p_n(0)+\{1-p_n(0)\}y
=1-\alpha^{1/n}(1-y).
$$

### Theorem 2 (universal binomial failure family)

For every $n\ge2$ and every $\alpha\in(0,1)$, there is an ordered finite
population with all item weights below one systematic interval for which
ordinary binomial Stringer has coverage strictly below $1-\alpha$.

**Proof.** Since $\alpha^{1/n}>\alpha$, the interval

$$
1-\alpha^{1/n}<q<1-\alpha
$$

is nonempty. Choose $q$ there and use the periodic population. On the low
phase,

$$
\theta-\operatorname{SB}_{\rm B}(y,\ldots,y)
=\{q-(1-\alpha^{1/n})\}(1-y)>0.
$$

That phase has probability $1-q>\alpha$. Hence noncoverage is greater than
$\alpha$. Every item weight is either $qh<h$ or $(1-q)h<h$, so the
failure does not depend on duplicate hits or certainty-item conventions.
$q$ may be chosen rational, say $a/m$; taking $h=m$ gives the same
construction with integer book amounts $a$ and $m-a$ in every cell.
$\square$

For this family the high phase has binomial bound one, so its exact coverage
is $q$, not merely at most $q$.

### Poisson factors

For the usual Poisson confidence factors, the all-zero bound is

$$
\operatorname{SB}_{\rm P}(0,\ldots,0)
=\frac{\lambda_0}{n}=\frac{-\log\alpha}{n}.
$$

Consequently, whenever

$$
\frac{-\log\alpha}{n}<1-\alpha,
$$

choose $q$ strictly between these quantities and take $y=0$. The all-zero
phase has bound below $\theta=q$ and probability $1-q>\alpha$, so ordinary
Poisson Stringer also undercovers. The condition holds at 90% confidence for
$n\ge3$, at 95% for $n\ge4$, and at 99% for $n\ge5$. Final capping at one
does not change the failing all-zero phase. Again $q$ may be chosen rational,
so the counterexample has an exact integer-book-unit version.

## 3. Exact small-population census

At 95% confidence, the script exhausts every ordered binary-taint population
of unit-book-value items through total size nine, for every factorization
$N=nm$ defining an integer systematic interval $m$. It checks 2,188
population-design pairs and 4,444 random-start phases with exact 80-bit
dyadic brackets whose signs are established by integer arithmetic.

There is no failure through $N=8$. The first occurs at $N=9$, $n=3$,
$m=3$, with ordered taints

$$
(1,1,0\mid1,1,0\mid1,1,0).
$$

The three random starts return respectively three ones, three ones, and three
zeros. The target is $2/3$, while

$$
p_3(0)=1-0.05^{1/3}\approx0.6315968501<\frac23.
$$

Thus exact coverage is $2/3$, against nominal coverage $0.95$. This finite
census is corroboration and a minimal witness within its stated search class;
Theorem 2 supplies the general proof.

## 4. A practical-size exact example and uplift

Take total recorded amount 1,000, sample size 100, and interval 10. In each of
100 consecutive cells put:

* one item of recorded amount 1 and taint 1; then
* one item of recorded amount 9 and taint 0.

The target overstatement fraction is $100/1000=0.1$. One of the ten integer
random starts selects all 100 taint-one items; the other nine starts select
all 100 taint-zero items. The exact coverage of both ordinary binomial and
ordinary Poisson Stringer is therefore $0.1$ at each displayed level. For the
Poisson result, the certificate proves the high phase covers by the exact
inequality
$\Pr\{\operatorname{Pois}(10)\le10\}>0.10$, which implies
$\lambda_{100}>10$.

| confidence | all-zero binomial | all-zero Poisson | exact coverage | safe output | binomial uplift | Poisson uplift |
|---:|---:|---:|---:|---:|---:|---:|
| 90% | 0.0227628 | 0.0230259 | 0.1000000 | 0.1000000 | 0.0772372 | 0.0769741 |
| 95% | 0.0295130 | 0.0299573 | 0.1000000 | 0.1000000 | 0.0704870 | 0.0700427 |
| 99% | 0.0450074 | 0.0460517 | 0.1000000 | 0.1000000 | 0.0549926 | 0.0539483 |

The displays are rounded; the JSON certificate contains rigorous lower and
upper rational endpoints for every factor and uplift.

## 5. The sharp taint-only design floor

The counterexample also identifies an unavoidable cost of design robustness.
Call a rule *taint-only* if its output depends on the observed taint vector
but not on sampled item identities, book weights, or population ordering.

### Theorem 3 (necessary diagonal value)

If a taint-only rule $V_\alpha:[0,1]^n\to[0,1]$ has coverage at least
$1-\alpha$ under one-start systematic PPS for every ordered finite
population, then for every $0\le y<1$,

$$
V_\alpha(y,\ldots,y)\ge 1-\alpha+\alpha y.
$$

**Proof.** Suppose the displayed inequality fails. Choose

$$
\frac{V_\alpha(y,\ldots,y)-y}{1-y}<q<1-\alpha
$$

and use the periodic population. Its target $y+q(1-y)$ is larger than
$V_\alpha(y,\ldots,y)$, so the low phase is a noncoverage event. Its
probability $1-q$ is greater than $\alpha$, a contradiction. $\square$

This theorem is deliberately limited to taint-only rules. Item identities and
book weights contain additional information and permit a smaller bound.

## 6. A design-valid Stringer overlay

Define

$$
M_\alpha=1-\alpha+\alpha\bar T.
$$

### Theorem 4 (randomization bound)

For every ordered finite population under the design of Section 1,

$$
\Pr_R\{M_\alpha<\theta\}\le\alpha.
$$

**Proof.** Put $X=1-\bar T$. Proposition 1 gives
$\mathbb E_R X=1-\theta$. The failure event is

$$
X>\frac{1-\theta}{\alpha}=\frac{\mathbb E_R X}{\alpha}.
$$

Markov's inequality bounds its probability by $\alpha$. $\square$

Theorem 3 shows that $M_\alpha$ is exactly sharp on every constant sample
among all taint-only design-valid rules. In this distribution-free sense, one
random start may carry only one observation's worth of information even when
the grid contains many points.

Now let $S$ be the set of distinct sampled item identities and define the
completion bound

$$
C(S)=\frac1W\left(
  \sum_{i\in S}w_it_i+\sum_{i\notin S}w_i
\right).
$$

Because every unobserved taint is at most one, $C(S)\ge\theta$
deterministically. Therefore

$$
H_\alpha=\min\{M_\alpha,C(S)\}
$$

also has coverage at least $1-\alpha$: whenever $M_\alpha\ge\theta$, both
arguments of the minimum are at least $\theta$.

Let $S_\alpha^{\rm cap}$ denote the ordinary binomial- or Poisson-factor
Stringer output capped at one. The pre-specified reporting rule

$$
\boxed{
S_\alpha^{\rm sys}
=\max\{S_\alpha^{\rm cap},H_\alpha\}
}
$$

has design coverage at least $1-\alpha$, is never below the familiar capped
Stringer output, and uses book-weight information only to reduce the sharp
taint-only fallback. No post-sample method selection is involved: the maximum
is the defined procedure.

In the 100-draw all-zero phase, the sampled taint-zero items comprise 900 of
the 1,000 recorded units. Hence $M_{0.05}=0.95$, $C(S)=0.10$,
$H_{0.05}=0.10$, and the safe output is exactly 0.10. Its uplift is about
0.070487 over binomial Stringer and 0.070043 over Poisson Stringer. On the
all-one phase, both design bounds equal one.

## 7. Exact finite-frame minimax inversion

For a fixed known frame, combine random-start intervals that induce the same
item-hit multiplicities into atoms $a=1,\ldots,K$, with probabilities $p_a$.
For a complete taint vector $u$, write $Y_a(u)$ for its grid mean in atom
$a$. If the observed phase audits the distinct set $S$, returns $x=u_S$, and
has mean $y$, define

$$
B^{\rm ff}_\alpha(S,x)
=\max\left(
\left\{\theta(u):u_S=x,\ 0\le u\le1,\
\sum_a p_a\mathbf1\{Y_a(u)\le y\}>\alpha\right\}\cup\{0\}
\right).
$$

This is a minimax completion bound for the specified grid-mean ordering: it
maximizes the target over every compatible population for which the observed
mean is not in a lower randomization tail of mass at most $\alpha$. It does
not assert optimality over other test orderings or all confidence procedures.

For the true population $t$, noncoverage implies
$F_t(Y_R(t))\le\alpha$, where $F_t$ is the exact CDF of the grid mean over
random starts. The elementary discrete-rank inequality
$\Pr\{F_t(Y_R(t))\le\alpha\}\le\alpha$ therefore proves design validity.
Moreover, every compatible completion has target at most $C(S)$. For any
completion counted in the maximum, Markov's inequality applied to $1-Y$
gives

$$
\alpha<F_u(y)\le\frac{1-\theta(u)}{1-y}
$$

when $y<1$. Hence

$$
B^{\rm ff}_\alpha(S,x)\le\min\{M_\alpha,C(S)\}=H_\alpha.
$$

Replacing $H_\alpha$ by $B^{\rm ff}_\alpha$ in the safe overlay preserves
design coverage and capped ordinary Stringer while never increasing the
reported result.

The computation is exact. For every inclusion-minimal phase coalition $A$
with $\sum_{a\in A}p_a>\alpha$, maximize $\theta(u)$ subject to $u_S=x$,
$0\le u\le1$, and $Y_a(u)\le y$ for $a\in A$. These are rational linear
programs, and their maximum is exactly $B^{\rm ff}_\alpha$. The solver
enumerates all minimal coalitions, uses rational simplex rather than a
floating-point MILP certificate, checks matching feasible primal and dual
objectives, expands aggregated item variables to a full least-favourable
taint vector, and checks the witness. Exponential search limits fail closed.
Full proofs and the algorithm boundary are in
`theory/SYSTEMATIC-MINIMAX.md`.

For a practical exact improvement, take 2,000 unit items, $n=100$, and
interval 20. The 20 equally likely phases select disjoint sets of 100 items.
After one all-zero phase at 95% confidence, both $M_{0.05}$ and $C(S)$ equal
0.95. A lower-tail mass strictly above 0.05 requires at least two zero-mean
phases, forcing 200 taints to zero; setting the remaining 1,800 to one shows

$$
B^{\rm ff}_{0.05}=0.90.
$$

The certificate enumerates all $\binom{20}{2}=190$ minimal coalitions and
returns this exact witness, lowering the previous valid fallback by five
percentage points.

More generally, when the phase atoms partition the frame, the finite-frame
bound is

$$
B^{\rm ff}_\alpha(j,y)=1-(1-y)
\min\{p(A):j\in A,\ p(A)>\alpha\}.
$$

For equal phases this is
$1-(\lfloor m\alpha\rfloor+1)(1-y)/m$. Rational unequal phases reduce to a
minimum-excess subset sum: an exact denominator-cleared dynamic program uses
$O(KD)$ time and $O(D)$ memory, while a reduction from SUBSET SUM proves
weak NP-hardness under binary encoding. Full proofs and the specialized
implementation are in `theory/SYSTEMATIC-DISJOINT-PHASES.md`.

## 8. Independent starts

If the design uses $r$ independent random starts and forms the $r$ grid
averages $Y_1,\ldots,Y_r$, then the $Y_\ell$ are i.i.d. variables in
$[0,1]$ with common mean $\theta$. Any proved distribution-free bounded-mean
rule may therefore be applied to these $r$ averages. In particular, the
valid Gaffke limit gives a design-valid multiple-start bound; at $r=1$ its
sharp endpoint is $1-\alpha+\alpha Y_1=M_\alpha$. This distinguishes the
number of independent starts from the number of correlated points within
each grid.

There is also a sharp finite-frame formula at an all-zero observation. For
$m$ equal disjoint phases, if the $r$ starts visit $d$ distinct phases, put

$$
q=\max\{d,\lfloor m\alpha^{1/r}\rfloor+1\}.
$$

Exact inversion of the mean of the $r$ grid averages gives
$B^{\rm ff}_{\alpha,r}=1-q/m$. The lower-tail probability of a candidate
with $q$ zero phases is exactly $(q/m)^r$. The unequal rational analogue is
computed by the same dynamic program using an exact product-tail threshold.

## 9. Literature position and verification boundary

Hoogduin, Hall, Tsay, and Pierce, “Does Systematic Selection Lead to
Unreliable Risk Assessments in Monetary-Unit Sampling Applications?”,
*Auditing: A Journal of Practice & Theory* 34(4) (2015), 85–107,
[doi:10.2308/ajpt-51081](https://doi.org/10.2308/ajpt-51081), already reports
simulation evidence that evaluating systematic MUS with a simple-random-
sample distribution can create material risk-assessment error. Berger,
Chiodini, and Zenga, “Bounds for monetary-unit sampling in auditing: an
adjusted empirical likelihood approach,” *Statistical Papers* 62 (2021),
2739–2761,
[doi:10.1007/s00362-020-01209-w](https://doi.org/10.1007/s00362-020-01209-w),
likewise treats systematic PPS without replacement as an important practical
MUS design. The present result is therefore not positioned as the first
warning about systematic MUS.

The qualified mathematical distinction, based on sources located through
September 1, 2026, is the exact all-$n$, all-confidence binomial failure
family, the sharp taint-only diagonal obstruction, and the closed-form
design-valid overlay. No predecessor for that combination was located. The
full text of the 2015 article was not openly accessible during the search, so
priority language must remain qualified.

The proof boundary is explicit:

* ordinary mathematics proves design unbiasedness, the counterexample family,
  the diagonal lower bound, and validity of $M_\alpha$, $C(S)$,
  $H_\alpha$, and $S_\alpha^{\rm sys}$;
* Python enumerates rational random-start phases, exhausts the stated small
  binary class, regenerates the practical-size failure example, and enumerates
  exact minimal coalitions for the finite-frame inversion;
* exact rational simplex solves the finite-frame least-favourable-completion
  LPs and the implementation checks a full witness population;
* a separate denominator-cleared subset-sum implementation checks the
  disjoint-phase formulas, reconstructs an exact minimum-mass coalition, and
  fails closed at an explicit scaled-state limit;
* the existing exact factor routines certify the binomial and Poisson
  comparisons without floating-point sign decisions; and
* the multiple-start Gaffke statement uses the separately cited external
  validity theorem for that bounded-mean limit.

No result here covers negative taints, ratio estimands with uncertain book
amounts, stratified combinations, successive PPS, nonrandom starts,
post-selection among procedures, or global optimality over test orderings.
The multiple-start finite-frame theorem is limited to phase-partitioned
all-zero observations. Other cases require their own estimands and design
arguments.
