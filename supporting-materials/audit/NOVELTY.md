# Novelty and priority verification

An adversarial literature pass (~20 targeted searches, OpenAlex citation
tracing, primary-source reading including the full Bimpeh thesis, the
2021 Statistical Papers MUS survey, the Dutch statistical-auditing essay
collection, and the 2024 Balakrishnan–Rychlik monograph) tested each of
the paper's novelty claims. Initial date: 2026-08-07; updated 2026-08-11
for the new results through \(n=5\), the direct Poisson coverage ranges,
and the analytic all-sample-size one-cap comparison.
Priority statements remain
qualified by the scope and date of this search.

## Claim 1 — first unconditional finite-sample proof at n = 2: HOLDS, with a caveat now addressed

No valid prior proof exists; the literature through 2021 uniformly calls
the problem open. **Caveat found and addressed**: Bimpeh's thesis §5.3
explicitly *claims* the n ≤ 2 result via P̄₂ = 2√(1−α) − (1−α) ≥ 1−α.
Since our Proposition refutes his containment inequality at exactly
n = 1, 2 (continuous counterexamples with CP < P̄ₙ), his claimed proof
is invalid, and the paper now frames the theorem as the first *correct*
finite-sample proof, stating his claim explicitly. The refutation is
load-bearing for both claims and is machine-checked.

## Claim 2 — the wedge inequality is new: HOLDS

- **Rychlik school** (sharp order-statistic/L-statistic bounds;
  Balakrishnan–Rychlik 2024 monograph): covers expectations and
  distribution functions of single order statistics — no sharp bound on
  P(L-statistic ≥ t) under a mean constraint exists there.
- **Closest relative**: Hoeffding & Shrikhande (1955) solved the
  equal-weights case (sharp bound for P(X₁+X₂ ≥ t), nonnegative iid,
  given mean); Łuczak–Mieczkowska–Šileikis (arXiv:1602.03547) extended
  to general k. The wedge inequality has unequal min/max weights except
  at the interior value α = 16/25. Both now cited.
- **Valid-mean-bound genre** (Anderson 1969; Gaffke 2005;
  Learned-Miller–Thomas; Phan et al. 2021; Bentkus–van Zuijlen 2003;
  Vlassis–Thomas 2026 proof of Gaffke's conjecture): none proves or
  implies the wedge inequality; Phan et al. cite the Stringer
  literature and list its coverage as unknown.
- **Trapdoor check (performed, negative)**: could n = 2 conservatism
  follow from a known-valid bound pointwise below SB?
  - Anderson's bound: no — SB fails to dominate it (gap −0.0253 at
    sample (0,1), α = 0.05).
  - Gaffke's bound: at α = 0.05, SB *does* pointwise dominate it
    (minimum gap 0, ties exactly on the two extremal families:
    equal-taint samples and samples with one taint = 1) — so at
    standard confidence an alternative route would exist *given* the
    2026 validity proof plus a proof of the (itself new) domination;
    but at α = 0.7 the domination fails by −0.032 while our theorem
    still holds. The theorem is strictly stronger and self-contained.
    Check script: `computations/python/gaffke_domination_check.py`.

## Claim 3 — exact conventional-level theorems from n = 3 through n = 5: NO PREDECESSOR FOUND

The 2026 proof of Gaffke's test by Vlassis–Thomas and its bounded-mean
inversion described by Learned-Miller–Thomas and Ming et al. supply a valid
comparison limit. No source found in the review compares that limit
pointwise with Stringer, proves the resulting uniform-simplex cap
inequality, or establishes Stringer coverage at \(n=3\) and 90%, 95%, and
99% confidence, much less the corresponding \(n=4\) and \(n=5\) results. The repository's
contribution is the Stringer-specific comparison and exact rational
Bernstein proofs; validity of the Gaffke limit is credited to those sources.
The manuscript does not make an unqualified priority claim for these
results.

## Claim 4 — all-n Poisson-over-binomial factor transfer: DERIVED FROM A CLASSICAL INEQUALITY

The probability comparison is due to Anderson–Samuels (1967). The paper's
contribution is to specialize it to the factor roots and use summation by
parts to compare the complete Stringer bounds for every observed sample.
This should be presented as a new audit-specific consequence of a classical
inequality, not as a new probability inequality.

## Claim 5 — direct Poisson simultaneous-band ranges: NO PREDECESSOR FOUND

The proof combines a standard randomized probability integral transform,
the correctly shifted survival-band constraints, and Bolshev's classical
uniform-order-statistic recursion. No reviewed Stringer source applies that
corrected event to the untruncated Poisson audit factors or obtains the
level-dependent ranges \(n\le8\), \(n\le11\), and \(n\le20\). The exact
rational treatment of the Poisson limits and event probabilities is also new
in the reviewed literature. Because the ingredients are classical and the
priority search was not exhaustive, the manuscript states the theorem
without an unqualified first-priority claim.

## Claim 6 — all-n one-cap comparison at confidence at least \(1-e^{-2}\): NO PREDECESSOR FOUND

The dimension-free cap lemma and the reduction of its
Clopper--Pearson vertices to binomial upper tails are Stringer-specific.
The proof then uses two credited probability results: Anderson--Samuels
(1967) fixed-mean binomial monotonicity and Pinelis's (2021) sharp lower
bound for a binomial exceeding its expectation. No reviewed source combines
these ingredients, proves the resulting factor inequality for all sample
sizes, or derives the zero-uplift region for the Stringer--Gaffke safeguard.
The manuscript presents the result without an unqualified priority claim.

## Claim 7 — first identification of Bimpeh's error: HOLDS, strongly

The thesis has exactly one recorded citation (a 2019 paper with no
engagement with ch. 5); no published criticism, correction, or erratum
exists; the 2002 Bimpeh–Horgan conference abstract contains none of the
chapter-5 material. The error went unnoticed for 18 years; this
repository is the thesis's first substantive engagement in the
literature.

## Claim 8 — first exact explicit counterexample parameters: HOLDS

The only published counterexamples remain Pap–van Zuijlen 1995
(uniform-type taints, existence via recursions, α > 1/2, n > 16); every
later source refers back to them. No prior three-point examples with two
nonzero taint values, no certified rational-arithmetic coverage, and no
violations published in the 30–37% confidence band. Since Pap–van Zuijlen
1996 proved
*asymptotic* conservatism for α ≤ 1/2, certified finite-n violations
approaching the 50% threshold from below sharpen the known picture.

## Bibliography actions taken

Added to the manuscript: Anderson 1969; Balakrishnan–Rychlik 2024;
Gaffke 2005; Hoeffding–Shrikhande 1955; Łuczak–Mieczkowska–Šileikis
2016; Learned-Miller–Thomas 2020; Ming et al. 2026; Pinelis 2021;
Vlassis–Thomas 2026.
Already present: Bickel 1992; Bimpeh 2008;
de Jager–Pap–van Zuijlen 1997; NRC Panel 1989; Pap–van Zuijlen 1995,
1996; Phan et al. 2021; Stringer 1963.
