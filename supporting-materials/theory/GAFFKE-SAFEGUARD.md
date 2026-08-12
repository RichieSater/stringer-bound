# An all-sample-size Stringer--Gaffke safeguard

## Purpose

The ordinary Stringer coverage conjecture remains unresolved for arbitrary
sample size at conventional audit confidence levels.  This note records a
simple way to make the reported upper bound distribution-free valid **now**,
without assuming that conjecture and without discarding the familiar Stringer
calculation.

For observed taints `t=(t_1,...,t_n)`, report

\[
  \operatorname{SB}^{\rm safe}_{M,\alpha}(t)
  =\max\{\operatorname{SB}_{M,\alpha}(t),G_\alpha(t)\},
  \tag{1}
\]

where `M` is either the binomial or Poisson factor convention and

\[
G_\alpha(t)
=Q_{1-\alpha}\!\left(D_0+\sum_{i=1}^n t_iD_i\right),
\qquad
(D_0,\ldots,D_n)\sim\operatorname{Dirichlet}(1,\ldots,1).
\tag{2}
\]

Vlassis and Thomas proved Gaffke's finite-sample validity conjecture for
independent nonnegative variables.  Its inversion for observations in
`[0,1]` gives (2), the one-sided Gaffke--Learned-Miller--Thomas upper limit.
See:

- N. Vlassis and P. S. Thomas, *An exact distribution-free test for means of
  nonnegative random variables*, arXiv:2607.08415 (2026).
- J. Ming, A. Ramdas, Y. Shen, R. Wang, and I. Waudby-Smith, *Gaffke's
  confidence interval for the mean of bounded data is inadmissible but
  asymptotically efficient*, arXiv:2607.18661 (2026), Section 5.

## Coverage theorem

**Theorem.** Let `T_1,...,T_n` be independent `[0,1]`-valued observations
with common mean `mu`.  For every `n>=1`, every `alpha in (0,1)`, and either
Stringer factor convention,

\[
 \Pr\{\operatorname{SB}^{\rm safe}_{M,\alpha}(T_1,\ldots,T_n)
       \ge \mu\}\ge 1-\alpha.
 \tag{3}
\]

**Proof.** The validity theorem for (2) gives
`Pr{G_alpha(T)>=mu}>=1-alpha`.  Pointwise,
`SB^safe_{M,alpha}>=G_alpha`, so the event in the preceding display is
contained in the event in (3).  No coverage property of the ordinary
Stringer component is used.  `square`

This result is deliberately not described as a solution of the Stringer
conjecture.  It validates the procedure in (1), not the ordinary Stringer
rule in cases where Stringer is smaller than `G_alpha`.

## When the safeguard leaves Stringer unchanged

If the observed Stringer value is at least `G_alpha(t)`, (1) reports exactly
the ordinary Stringer value.  The procedure must nevertheless be specified
in advance as the maximum in (1); the coverage theorem is for that complete
rule.

The exact comparisons already proved in this repository give a stronger
uniform statement:

- for `n=3,4,5` and `alpha in {0.10,0.05,0.01}`, binomial Stringer is at
  least `G_alpha` for **every** sample, so the safeguard is mathematically
  identical to ordinary binomial Stringer;
- at the same levels and sample sizes, Poisson Stringer is no smaller than
  binomial Stringer, so the Poisson safeguard is also identical to ordinary
  Poisson Stringer.

There is also a larger, sample-checkable region. At every
`0<alpha<=1/4` and every sample size, if the **binomial** Stringer value
is at least the largest observed taint, then it is at least `G_alpha(t)`. The
safeguard therefore has zero uplift on that sample for both factor
conventions. The proof combines a dimension-free one-upper-knot cap lemma
with analytic binomial-tail comparisons; see
[`ONE-CAP-COMPARISON.md`](ONE-CAP-COMPARISON.md). An older set of 59,700
exact Clopper--Pearson vertex checks through `n=200` remains as an
independent regression. This is a pointwise identity for the pre-specified
safeguarded rule, not a conditional-coverage claim for ordinary Stringer.

Outside these proved identity regions, (1) acts as a validated floor. It can
be adopted without taking a position on the unresolved general-`n`
comparison.

## Exact-sign computation of the floor

The implementation in
[`gaffke.py`](../computations/python/gaffke.py) does not treat a floating-point
B-spline value as a certificate.  Suppose the knots
`a_0,...,a_n=(t_1,...,t_n,1)` and the evaluation point `x` are rational.
The uniform-Dirichlet upper tail has the divided-difference representation

\[
 \Pr_D\!\left\{\sum_{i=0}^n a_iD_i>x\right\}
 =[a_0,\ldots,a_n](u-x)_+^n.                 \tag{4}
\]

For sorted knots, write `c_i^(0)=(a_i-x)_+^n`.  The ordinary divided-
difference update is

\[
c_i^{(r)}=
\frac{c_{i+1}^{(r-1)}-c_i^{(r-1)}}{a_{i+r}-a_i}.
\tag{5}
\]

When the denominator vanishes because a block of knots is repeated, the
confluent value is

\[
c_i^{(r)}={n\choose r}(a_i-x)_+^{n-r}.       \tag{6}
\]

Equations (4)--(6) use only rational arithmetic.  The final value
`c_0^(n)` is therefore the exact upper-tail probability.  The implementation
computes the same confluent divided difference through an equivalent
compressed residue formula.  If `m_a` is the multiplicity of a distinct
knot `a`, its contribution when `a>x` is

\[
[y^{m_a-1}]
(a+y-x)^n\prod_{b\ne a}(a-b+y)^{-m_b}.       \tag{7}
\]

Knots at or below `x` contribute zero.  Truncating every power series at
degree `m_a-1` preserves exact rational arithmetic and avoids building a
full quadratic table for the many repeated zeros typical of audit samples.

The code uses a normalized B-spline only to propose a narrow location for
the quantile.  It snaps two surrounding points to a `2^-48` grid and checks
with (4) that

\[
 \Pr_D\{T_D>\ell\}\ge\alpha,
 \qquad
 \Pr_D\{T_D>u\}\le\alpha.                   \tag{8}
\]

If either check fails, the dyadic bracket expands and is checked again.
Thus `u` is a rational upper enclosure of `G_alpha(t)`.  The command reports
that upper endpoint, its dyadic representation, and the certified bracket
width.  Decimal command-line inputs are interpreted as exact rationals.

The Stringer component is evaluated at high precision.  Formal validity of
the returned safeguard does not depend on rounding that component: the
reported maximum is at least the exactly sign-checked upper enclosure `u`.

## Reproduction

From the repository root:

```sh
uv run --frozen python \
  supporting-materials/computations/python/gaffke.py \
  --n 100 --alpha 0.05 --method poisson --taints 1,0.4,0.1
```

Zero taints may be omitted; `--n` remains the full sample size.  The output
separates the Stringer component, the Gaffke bracket, the maximum, and any
uplift.  The regression tests include beta-law special cases, repeated
knots, exact tail-sign checks, Bernoulli/Clopper--Pearson agreement, and the
certified no-uplift cases at `n=3,4,5`, together with representative
all-sample-size one-cap cases.

## Audit-use boundary

Equation (3) is a statistical guarantee under the independent-observations,
common-mean model used in the manuscript.  It does not by itself address
finite-population sampling without replacement, negative taints, audit-unit
selection errors, or the conversion from a taint-rate upper bound to a
book-value misstatement conclusion.  Those design and implementation issues
must be handled by the applicable professional standard and the engagement's
sampling plan.
