# Disjoint-phase finite-frame inference: formulas, complexity, and replication

## Mathematical status

**Status.** The closed forms, dynamic-program recurrence, subset-sum
reduction, and independent-start theorem below are proved by ordinary
mathematics. The Python implementation performs the rational dynamic program
exactly and regenerates finite corroborating examples.

**One-sentence thesis.** When the atoms of a systematic-PPS random-start
design partition the frame, exact finite-frame inversion reduces to a
minimum-excess subset sum: equal phases have a closed form, rational unequal
phases admit a pseudo-polynomial exact dynamic program, exact evaluation is
weakly NP-hard, and independent starts have an explicit all-zero bound.

The results specialize the grid-mean inversion in
[`SYSTEMATIC-MINIMAX.md`](SYSTEMATIC-MINIMAX.md). They do not assert
optimality among other test orderings or confidence procedures.

## 1. Phase-partition frames

Let the random start have atoms $a=1,\ldots,K$, with probabilities

$$
\pi_a>0,\qquad \sum_{a=1}^K\pi_a=1,
$$

and phase means $Y_a(u)\in[0,1]$. Call the frame a *phase-partition frame*
when every item belongs to the support of exactly one phase atom. Exact
systematic-PPS design unbiasedness then gives

$$
\theta(u)=\sum_{a=1}^K\pi_aY_a(u). \tag{1}
$$

The partition property lets the unobserved phase means be varied
independently: assigning every taint in phase $a$ the common value $z_a$
makes $Y_a=z_a$ without changing another phase. In the observed phase
$j$, the audited taints remain fixed and determine $Y_j=y$.

For a candidate completion $u$, write

$$
F_u(y)=\sum_{a=1}^K\pi_a\mathbf 1\{Y_a(u)\le y\}.
$$

The finite-frame bound retains exactly the candidates satisfying
$F_u(y)>\alpha$.

## 2. Unequal- and equal-phase formulas

### Theorem 1 (unequal disjoint phases)

For an observed phase $j$, define

$$
\rho_j(\alpha)
=\min\left\{
 \sum_{a\in A}\pi_a:
 j\in A,\quad \sum_{a\in A}\pi_a>\alpha
\right\}. \tag{2}
$$

Then the exact finite-frame bound is

$$
B^{\rm ff}_\alpha(j,y)
=1-(1-y)\rho_j(\alpha). \tag{3}
$$

**Proof.** For any feasible completion let
$A=\{a:Y_a(u)\le y\}$. The observed phase belongs to $A$, and feasibility
gives $\pi(A)>\alpha$. By (1),

$$
\theta(u)
=\sum_{a\in A}\pi_aY_a(u)
 +\sum_{a\notin A}\pi_aY_a(u)
\le y\pi(A)+1-\pi(A)
=1-(1-y)\pi(A)
\le1-(1-y)\rho_j(\alpha).
$$

Choose a minimizing coalition $A_\star$ in (2). Retain the observed taints
in phase $j$, assign taint $y$ to every item in each other phase of
$A_\star$, and assign taint one to every phase outside $A_\star$. Its
phase means are $y$ on $A_\star$ and one elsewhere, so its lower-tail
mass is $\pi(A_\star)>\alpha$, and (1) gives equality in (3). $\square$

### Corollary 2 (equal disjoint phases)

If the $m$ phases are equiprobable, put

$$
q_{m,\alpha}=\lfloor m\alpha\rfloor+1.
$$

Then

$$
B^{\rm ff}_\alpha(j,y)
=1-\frac{q_{m,\alpha}}{m}(1-y). \tag{4}
$$

**Proof.** A coalition has mass greater than $\alpha$ exactly when its
cardinality is at least $q_{m,\alpha}$. Such a coalition can always be
chosen to contain the observed phase. Hence
$\rho_j(\alpha)=q_{m,\alpha}/m$, and Theorem 1 applies. $\square$

For $m=20$, $\alpha=1/20$, and $y=0$, formula (4) gives $q=2$ and
$B^{\rm ff}=0.90$. Thus the 2,000-item example in the manuscript is one
instance of a complete formula rather than an isolated LP optimum.

## 3. Exact rational dynamic program

Suppose all phase probabilities are rational. Clear denominators and write

$$
\pi_a=\frac{c_a}{D},\qquad
c_a\in\mathbb Z_{>0},\qquad \sum_a c_a=D.
$$

The smallest integer mass strictly above the one-start threshold is

$$
h=\lfloor\alpha D\rfloor+1. \tag{5}
$$

The observed phase contributes $c_j$. For the other phases, define the
reachable-mass sets recursively by

$$
R_0=\{0\},\qquad
R_{k+1}=R_k\cup(R_k+c_{a_k}), \tag{6}
$$

where $a_1,\ldots,a_{K-1}$ list the phases other than $j$. Let

$$
s_\star=\min\{s\in R_{K-1}:c_j+s\ge h\}. \tag{7}
$$

Then

$$
\rho_j(\alpha)=\frac{c_j+s_\star}{D}. \tag{8}
$$

Every subset of optional phases contributes a state in (6), and every state
comes from such a subset. Condition (5) is exactly the strict inequality
$(c_j+s)/D>\alpha$. Therefore the first crossing in (7) is the coalition
of minimum probability, proving (8). A Boolean array and one back-pointer per
reachable sum compute the coalition in $O(KD)$ time and $O(D)$ memory.
The dependence on $D$, rather than on its bit length, is pseudo-polynomial.

For $r$ independent starts and an all-zero observation, replace (5) by

$$
h_r=\min\left\{h\in\{0,\ldots,D\}:
 \left(\frac hD\right)^r>\alpha\right\}, \tag{9}
$$

and require every distinct observed phase in the initial coalition. The same
recurrence then computes the unequal-phase all-zero bound exactly. The code
tests (9) by integer exponentiation, so no algebraic root or floating-point
comparison enters the result.

## 4. Weak NP-hardness

### Theorem 3 (exact evaluation is weakly NP-hard)

Under binary encoding of rational weights and $\alpha$, exact computation
of $B^{\rm ff}_\alpha$ is weakly NP-hard even when $n=1$, the phases are
disjoint, and the observed taint is zero.

**Proof.** Reduce from **SUBSET SUM**. Given positive integers
$a_1,\ldots,a_k$ and a target $T$ satisfying
$1\le T\le A:=\sum_i a_i$, construct a one-draw systematic-PPS frame with
book weights

$$
w_0=1,\qquad w_i=2a_i\quad(1\le i\le k),
\qquad W=1+2A, \tag{10}
$$

and set

$$
\alpha=\frac{2T}{W}. \tag{11}
$$

With $n=1$, each item is a disjoint phase having probability $w_i/W$.
Observe item $0$ with taint zero. Any candidate's lower tail at zero is the
total probability of its zero-taint phases, which must include phase $0$.
For a subset $I\subseteq\{1,\ldots,k\}$, the corresponding coalition has
integer weight

$$
1+2\sum_{i\in I}a_i.
$$

By (11), its probability is strictly above $\alpha$ exactly when

$$
1+2\sum_{i\in I}a_i>2T,
$$

which, by integrality, is equivalent to
$\sum_{i\in I}a_i\ge T$. Consequently, if

$$
s_\star=\min\left\{
 \sum_{i\in I}a_i:\sum_{i\in I}a_i\ge T
\right\},
$$

then Theorem 1 at $y=0$ gives

$$
B^{\rm ff}_\alpha=1-\frac{1+2s_\star}{W}. \tag{12}
$$

There is a subset summing exactly to $T$ if and only if the exact output in
(12) equals $1-(1+2T)/W$. Thus an exact evaluator decides SUBSET SUM. The
construction has polynomial encoding length. Together with the
pseudo-polynomial algorithm in Section 3, this establishes the stated weak
NP-hardness. $\square$

The reduction explains why the general coalition enumeration cannot be
replaced by a polynomial-time exact algorithm solely from the disjointness
assumption. It also identifies the correct tractable parameter: the scaled
total phase mass $D$.

## 5. Independent-start all-zero theorem

### Theorem 4 (equal phases and independent starts)

Take $r\ge1$ independent random starts on $m$ equiprobable disjoint
phases. Let $d$ distinct phases occur among the $r$ observed starts, and
suppose every observed taint is zero. Invert the lower-tail rank of the mean
of the $r$ grid averages. Put

$$
q_{m,r,\alpha,d}
=\max\left\{
d,\ \left\lfloor m\alpha^{1/r}\right\rfloor+1
\right\}. \tag{13}
$$

Then the exact finite-frame bound is

$$
B^{\rm ff}_{\alpha,r}
=1-\frac{q_{m,r,\alpha,d}}{m}. \tag{14}
$$

**Proof.** For a candidate completion, let $Z$ be the set of phases whose
means are zero. Nonnegative taints imply that the mean of $r$ independently
selected phase means is at most the observed value zero exactly when every
selected phase lies in $Z$. Hence its exact lower-tail probability is

$$
\Pr\{\overline Y_r=0\}
=\left(\frac{|Z|}{m}\right)^r. \tag{15}
$$

Compatibility forces the $d$ observed phases into $Z$. Feasibility of
the inverted confidence set requires the probability in (15) to exceed
$\alpha$, so

$$
|Z|\ge
\max\left\{d,\lfloor m\alpha^{1/r}\rfloor+1\right\}.
$$

Every zero-mean phase contributes zero to the target and every other phase
contributes at most $1/m$, giving the upper bound in (14). It is attained
by setting exactly $q_{m,r,\alpha,d}$ phases, including all observed ones,
to zero and every remaining phase to one. $\square$

For $m=20$ and $\alpha=0.05$, distinct all-zero starts give exact bounds
$0.90$, $0.75$, and $0.60$ for $r=1,2,3$, respectively. The result
quantifies the distinction between more points within one dependent grid and
more independently randomized grids.

## 6. Implementation and verification boundary

The implementation is
[`systematic_disjoint_dp.py`](../computations/python/systematic_disjoint_dp.py).
It:

* clears rational denominators and reduces to positive integer phase masses;
* finds the strict one- or multiple-start threshold exactly;
* performs the $0$--$1$ subset-sum recurrence with exact back-pointers;
* reconstructs and checks a minimum-mass coalition;
* accepts a complete rational systematic-PPS frame, enumerates its atoms, and
  refuses the specialized formula unless their supports partition every
  item;
* fails closed before allocation if the scaled total exceeds
  `max_scaled_total`; and
* exposes JSON interfaces for one-start arbitrary-mean and independent-start
  all-zero calculations.

For example, the frame-aware interface accepts the same observed-phase data
as the general LP solver:

```json
{
  "weights": ["1", "1", "1", "1"],
  "sample_size": 2,
  "observed_phase_index": 0,
  "observed_taints": {"0": "0", "2": "0"},
  "alpha": "1/4"
}
```

```sh
uv run --frozen python \
  supporting-materials/computations/python/systematic_disjoint_dp.py \
  --input frame.json --out result.json
```

The tests exhaust every coalition on several small unequal probability
vectors, compare the equal and unequal formulas against the independent
rational-LP solver, exercise exact-boundary strictness and the complexity
limit, and compare the multiple-start formula with the generic DP. The
committed
[`systematic-disjoint-dp-certificate.json`](../computations/certificates/systematic-disjoint-dp-certificate.json)
records the 20-phase formulas, an unequal-phase fixture, and yes/no finite
instances of the hardness reduction.

Ordinary mathematics, not the certificate, proves Theorems 1, 3, and 4. The
certificate establishes only its listed finite rational calculations. No
claim here covers overlapping phase supports, a nonzero multiple-start
observation, dependent starts, negative taints, or a data-dependent choice of
test ordering.
