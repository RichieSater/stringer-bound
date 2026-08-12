# A clean reformulation, and the n = 2 kernel

> **Status update: the n = 2 wedge inequality stated below is now a
> theorem — see [`N2-PROOF.md`](N2-PROOF.md) for the complete proof
> (all F, all confidence levels). This note is kept as the derivation
> of the reformulation and the record of the partial results that led
> there. Separate Gaffke/simplex-cap arguments now prove the cases
> n = 3 through n = 5 at 90%, 95%, and 99%; see
> [`N3-CONVENTIONAL.md`](N3-CONVENTIONAL.md),
> [`N4-CONVENTIONAL.md`](N4-CONVENTIONAL.md), and
> [`N5-CONVENTIONAL.md`](N5-CONVENTIONAL.md). A five-dimensional
> certificate also proves n = 6 at 95%; see
> [`N6-CONVENTIONAL.md`](N6-CONVENTIONAL.md).**

*Everything in this note is elementary algebra plus case analysis; each
claim states its verification status. Nothing here is a proof of the
conjecture beyond the stated special cases.*

## The complement substitution

Let \(U_i = 1 - T_i \in [0,1]\) ("clean fraction" of a sampled dollar
unit), so the \(U_i\) are i.i.d. with mean \(w = 1 - \mu\), and let
\(u_{(1)} \le \dots \le u_{(n)}\) be their increasing order statistics.
Substituting \(t_{(j)} = 1 - u_{(n+1-j)}\) into the Stringer bound and
telescoping (using \(p_n = 1\)):

\[
\mathrm{SB} \;=\; 1 \;-\; \sum_{j=1}^{n} e_j\, u_{(j)},
\qquad e_j := p_{j} - p_{j-1} \; > 0,
\qquad \sum_j e_j = 1 - p_0 ,
\]

(both orderings reverse together: the \(j\)-th largest taint pairs with
the \(j\)-th smallest clean fraction). The positive factor increments
\(e_j\) encode the confidence-factor structure, and

\[
P(\mathrm{SB} < \mu) \;=\; P\Big( \sum_i e_i\, U_{(i)} > w \Big).
\]

**The conjecture, reformulated.** For \(U_1,\dots,U_n\) i.i.d. on
\([0,1]\) with mean \(w\), and the specific weights \(e_i\) above:
\(P\big(\sum_i e_i U_{(i)} > w\big) \le \alpha\) for every distribution
and every \(w \in [0,1]\), when \(\alpha \le 1/2\).

This is a pure weighted-order-statistic exceedance inequality; the
auditing content is entirely absorbed into the weights. (Identity
machine-checked against the direct definition for n = 2, 5, 9;
`stringer.py` conventions.)

## n = 2: the wedge inequality

For \(n = 2\): \(p_0 = 1-\beta\), \(p_1 = \gamma\), \(p_2 = 1\) with
\(\beta = \sqrt{\alpha}\), \(\gamma = \sqrt{1-\alpha}\), so
\(\beta^2 + \gamma^2 = 1\), and

\[
\mathrm{SB} = 1 - \big( A\,\min(U_1,U_2) + B\,\max(U_1,U_2) \big),
\qquad A = \beta + \gamma - 1,\; B = 1 - \gamma,\; A + B = \beta .
\]

**Wedge inequality (now proved).** For \(U_1, U_2\)
i.i.d. on \([0,1]\) with mean \(w\):

\[
P\big( A \min(U_1,U_2) + B \max(U_1,U_2) > w \big) \;\le\; \beta^2 = \alpha .
\]

The following were the partial checks that preceded the full proof:

- **\(w \ge \beta\)** (i.e. \(\mu \le p_0\); also covers the trivial
  regime): failure requires \(\min U > (w-B)/A\) (from
  \(\max U \le 1\)), so by independence and Markov's inequality
  \(P \le \big(\tfrac{wA}{w-B}\big)^2\), and \(wA/(w-B) \le \beta
  \Leftrightarrow w \ge \beta\). Proved.
- **Two-point distributions** \(U \in \{0, u\}\) with
  \(P(U=0)=q\) (equivalently taints \(\{1, v\}\)): exact case analysis.
  Failure probability is \((1-q)^2\) for \(q \in (1-\beta, \gamma)\)
  and \(1-q^2\) for \(q > \gamma\); both are \(< \alpha\) on their
  ranges and both tend to \(\alpha\) at the endpoints
  (\((1-q)^2 \to \beta^2\) as \(q \downarrow 1-\beta\);
  \(1-q^2 \to 1-\gamma^2\) as \(q \downarrow \gamma\)). So the supremum
  over this family is exactly \(\alpha\), approached, never attained.
  Proved (and machine-checked against `coverage_exact`).
- **General \(0<w<\beta\)**: now proved in
  [`N2-PROOF.md`](N2-PROOF.md). Rectangle bounds on the covering event,
  a survival-integral budget, and a concave potential rule out every
  strict failure of the wedge inequality. Naive Markov/union bounds alone
  remain insufficient; they give \(O(\beta)\), rather than \(\beta^2\),
  on part of the range.

**A cautionary note.** A tempting shortcut — "failure forces the small
order statistic below a threshold whose Markov bound closes the
argument" — fails: the implication "failure \(\Rightarrow T_{\min} <
(\mu - p_1)/(1-p_1)\)" is simply false (counterexample: both taints
equal, moderately sized). We record this because it is the first thing
one tries and it silently gives a wrong theorem for all \(n\).

## Why this case matters

The \(n = 2\) wedge inequality was the smallest instance not settled by
the prior literature. Every ingredient is explicit—two i.i.d. variables,
a linear boundary, and exact constants with
\(\beta^2 + \gamma^2 = 1\)—yet the proof still requires information about
the entire survival function. The resulting rectangle-budget-potential
argument is now a concrete starting point for work at larger \(n\).
