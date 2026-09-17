# Openness verification

**Verdict (literature search and repository results updated 2026-09-17): the
general i.i.d. finite-sample conservatism
conjecture for the Stringer bound at confidence levels ≥ 50% — including
95% — remains open.** The present repository proves the i.i.d. case \(n=2\) at
every confidence level, proves the exact range
\(0<\alpha\le((19+\sqrt{21})/34)^3\) at \(n=3\), and gives independent
direct certificates at \(n=4,5,6,7\) for 90%, 95%, and 99%. The proposed
\(n=4\) and \(n=5\) continua are conjectural: the five- and six-coordinate
global cap arguments do not yet establish the active-prefix inequality on
repeated-lowest-knot faces. See
[`ORDERED-COLLISION-OBLIGATION.md`](../theory/ORDERED-COLLISION-OBLIGATION.md). Thus \(n=8\)
is the first unresolved binomial-factor sample size at all three levels.
For the Poisson factors used in audit practice, a separate exact
simultaneous-band argument proves coverage through \(n=8\), \(n=11\), and
\(n=20\), respectively. This does not resolve alpha outside the stated \(n=3\) interval,
other confidence levels at \(n=4,5,6,7\), or arbitrary \(n\). No prior general i.i.d. proof or i.i.d. counterexample at
α ≤ 1/2 was found in the reviewed literature. The
strongest partial results, and where a resolution could have hidden, are
recorded below with sources.

## What is proved

| Result | Source | Scope |
|---|---|---|
| Asymptotic conservatism for α ∈ (0, 1/2]; asymptotic **anti**-conservatism for α ∈ (1/2, 1). Threshold exactly α = 1/2 | Pap & van Zuijlen, *Statistica Neerlandica* 50 (1996) 367–389 ([1994 report full text](http://web.archive.org/web/20251019123811/https://repository.ubn.ru.nl/bitstream/handle/2066/100177/1/100177.pdf)) | asymptotic, a.s. expansions, binomial factors |
| Bickel's expansion: SB = T̄ + c(F) z₁₋α/√n + o(n^{-1/2}), with c²(F) ≥ Var(T), equality iff F is supported on ≤ 2 points | Bickel, *Int. Statist. Rev.* 60 (1992) 197–209; made rigorous in Pap–vZ 1996 | asymptotic |
| P(SB ≥ μ) ≥ (1−α)^{n+1} under conditions on F | Bickel 1992 | finite n, weak bound |
| Conjecture holds **exactly** for taints supported on {0,1} (SB reduces to Clopper–Pearson); Stringer coefficients are minimal with this property | de Jager, Pap & van Zuijlen, *Comput. Math. Appl.* 33 (1997) 37–54 | finite n, {0,1} supports |
| Binomial Stringer is distribution-free conservative at \(n=2\) for every α, sharply; it pointwise dominates Gaffke for every \(0<\alpha\le((19+\sqrt{21})/34)^3\) at \(n=3\), and at \(n=4,5,6,7\) for α = 0.10, 0.05, and 0.01 | [`N2-PROOF.md`](../theory/N2-PROOF.md), [`TETRAHEDRAL-VERTEX-BARRIER.md`](../theory/TETRAHEDRAL-VERTEX-BARRIER.md), and the independent [`n=4`](../theory/N4-CONVENTIONAL.md), [`n=5`](../theory/N5-CONVENTIONAL.md), [`n=6`](../theory/N6-CONVENTIONAL.md), [`n=7`](../theory/N7-CONVENTIONAL.md) direct certificates | finite n; all levels at \(n=2\), nominal confidence at least \(66.6314788\ldots\%\) at \(n=3\), and 90%, 95%, 99% at \(n=4,5,6,7\); the proposed \(n=4,5\) continua remain conjectural |
| Poisson audit factors dominate binomial Clopper–Pearson factors coordinatewise whenever nominal confidence exceeds \(1-e^{-1}\approx63.2\%\); hence the Poisson Stringer bound pointwise dominates the binomial version | Anderson–Samuels (1967), specialized to the Stringer formula in [`POISSON-DOMINATION.md`](../theory/POISSON-DOMINATION.md) | every n; factor and bound comparison, not a general coverage proof |
| A corrected simultaneous survival band proves Poisson-factor Stringer coverage for every \(n\le8\) at 90%, every \(n\le11\) at 95%, and every \(n\le20\) at 99% confidence | [`POISSON-SIMULTANEOUS-BAND.md`](../theory/POISSON-SIMULTANEOUS-BAND.md) and exact rational certificate | finite ranges; arbitrary continuous or atomic taint distributions |
| Two scalar paths satisfying the corrected simultaneous-band probability and terminal-factor condition give distribution-free valid Poisson rules for every sample size and confidence level: full scaling and an anchored path that leaves the no-error factor unchanged; capping the calibrated factors at one gives pointwise no-larger valid variants; for nominal confidence at least \(2/e\), the full-scale multiplier two is valid uniformly in \(n\), while a 200-term finite-prefix refinement gives uniform choices 1.525906, 1.436135, and 1.320081 at 90%, 95%, and 99% | [`POISSON-BAND-CALIBRATION.md`](../theory/POISSON-BAND-CALIBRATION.md); exact representative and uniform-multiplier certificates | every `n` and level, but modified calibrated procedures rather than proof of ordinary Stringer; the path must be pre-specified; the uniform multipliers are sufficient, not claimed optimal |
| The pre-specified maximum of either Stringer calculation and the valid Gaffke limit has coverage at least \(1-\alpha\) | [`GAFFKE-SAFEGUARD.md`](../theory/GAFFKE-SAFEGUARD.md) and the Vlassis–Thomas validity theorem | every n and level, but a safeguarded procedure rather than proof of ordinary Stringer |
| At every nominal confidence level of at least \(75\%\), the safeguard has zero uplift for every sample size whenever binomial Stringer is at least the largest observed taint | analytic proof in [`ONE-CAP-COMPARISON.md`](../theory/ONE-CAP-COMPARISON.md); 59,700 exact nonterminal vertex checks through \(n=200\) retained as an independent regression | pointwise one-cap region only; not conditional coverage or a full ordinary-Stringer theorem |
| Through \(n=200\), the safeguard has zero uplift in a complete adjacent two-cap collar: \(q=(1-t_{\max})/(1-\mathrm{SB})\) at least \(5/12\), \(1/3\), or \(1/4\) at 90%, 95%, or 99% | dimension-free formula and endpoint-maximization proof in [`TWO-CAP-COMPARISON.md`](../theory/TWO-CAP-COMPARISON.md); 59,700 exact endpoint checks | pointwise collar only; the residual two-cap region and general ordinary-Stringer coverage remain open |
| Any i.i.d.-valid upper rule at tail \(\eta\) has SRSWOR noncoverage at most \(\eta/\{(N)_n/N^n\}\); choosing \(\eta\le\alpha (N)_n/N^n\) gives finite-population coverage at least \(1-\alpha\) | [`SRSWOR-CONDITIONING.md`](../theory/SRSWOR-CONDITIONING.md) | exact conditioning bridge for an adjusted rule under uniform sampling of distinct indices; not ordinary unadjusted Stringer, systematic PPS, or successive PPS |
| Under one-start systematic PPS, ordinary binomial Stringer fails on some ordered population for every \(n\ge2\) and every confidence level; ordinary Poisson Stringer fails whenever \((-\log\alpha)/n<1-\alpha\). The sharp taint-only floor and item-aware completion overlay give a separate design-valid rule. On a fixed known frame, exact grid-mean randomization inversion gives a design-valid bound no larger than either fallback. Phase-partitioned frames have exact unequal/equal formulas, a pseudo-polynomial DP with a weak-NP-hardness boundary, and an exact independent-start all-zero formula | [`SYSTEMATIC-PPS.md`](../theory/SYSTEMATIC-PPS.md), [`SYSTEMATIC-MINIMAX.md`](../theory/SYSTEMATIC-MINIMAX.md), [`SYSTEMATIC-DISJOINT-PHASES.md`](../theory/SYSTEMATIC-DISJOINT-PHASES.md), and exact finite certificates | analytic failure, safeguard, inversion, dominance, phase-formula, complexity, and independent-start theorems for their stated fixed-population designs; the finite census, 100-draw failure outputs, 2,000-item rational-LP optimum, and DP fixtures have exact certificates |
| Claimed certification at n ≤ 11, α = 0.05 via the lower bound P̄ₙ = P(U_{i:n} ≤ p_n(i) ∀i) (Bolshev recursion). **Reassessed in this repository**: the containment CP ≥ P̄ₙ (his eq. 5.16) rests on an off-by-one in the band constraints, dropping F(t_{n:n}) ≥ α^{1/n}; hand counterexamples with continuous F at n = 1, 2 and an exact atomic one at n = 5 show that P̄ₙ is not the stated coverage bound. The corrected containment probability is ≤ 1−α for every continuous F. See [`BIMPEH-GAP.md`](BIMPEH-GAP.md). | Bimpeh, PhD thesis, DCU 2008, ch. 5 ([full text](https://doras.dcu.ie/600/1/YawThesis.PDF)); reassessment: this repository | the cited argument does not establish the claimed finite-sample range |
| Exact counterexamples at **sub-50% confidence**: uniform-type taints, α > 1/2, n ≳ 17 | Pap & van Zuijlen, *Comput. Math. Appl.* 29 (1995) 51–59 (exact recursions, not Monte Carlo) | finite n, α > 1/2 only |

This repository adds (see `theory/N2-PROOF.md`,
`theory/TETRAHEDRAL-VERTEX-BARRIER.md`,
	`theory/FIVE-COORDINATE-VERTEX-BARRIER.md`,
	`theory/N4-MONOTONE-RANGE.md`,
	`theory/SIX-COORDINATE-UNIQUE-TOP-L2-L3.md`,
	`theory/N5-MONOTONE-RANGE.md`,
`theory/N3-CONVENTIONAL.md`, `theory/N4-CONVENTIONAL.md`,
`theory/N5-CONVENTIONAL.md`,
`theory/N6-CONVENTIONAL.md`,
`theory/N7-CONVENTIONAL.md`,
`theory/POISSON-SIMULTANEOUS-BAND.md`,
`theory/POISSON-BAND-CALIBRATION.md`,
`theory/SYSTEMATIC-PPS.md`, `theory/SYSTEMATIC-MINIMAX.md`,
`theory/SYSTEMATIC-DISJOINT-PHASES.md`, and
`computations/certificates/`): a
distribution-free proof at \(n=2\); a four-coordinate
section-centroid maximum theorem yielding the stated \(n=3\) interval;
local five- and six-coordinate inequalities whose global step-normal
consequences remain conjectural; independent fixed-level proofs at
\(n=4,5,6,7\) for
α = 0.10, 0.05, and 0.01; direct exact
Poisson-factor guarantees through the level-dependent ranges above; the exact
systematic-PPS failure, design-safe overlay, exact finite-frame inversion,
disjoint-phase formulas and complexity boundary, and independent-start
all-zero theorem; and 33 exact rational i.i.d.
counterexamples supported on three points (two nonzero taint values and
zero) at selected nominal confidence levels from 30% to 37% and sample
sizes 50, 100, 200, and 400. Numerical searches at α = 0.05 found smallest
coverage agreeing with \(1-\alpha\) to search precision, approached on the
{0,1}-boundary (\(v_1\to1\)); this is evidence, not an infimum proof.
The all-sample-size scalar calibration, the Gaffke safeguard, and its
certified one-cap identity region do not change the openness verdict. The
first two validate different pre-specified reporting rules; the last controls
only part of the ordinary-Stringer sample space.
The systematic-PPS failure does not resolve the i.i.d. conjecture because it
uses a different joint sampling law; it instead proves that the i.i.d. claim
cannot be transferred to that design.

## Where a resolution could hide (checked, none found)

- **Modified-bound literature**: Lucassen–Moors–van Batenburg 1996
  (simulation only; their orderings shown asymptotically non-competitive in
  Pap & van Zuijlen, *Publ. Math. Debrecen* 57 (2000) 163–183,
  [PDF](https://publi.math.unideb.hu/paper/631/download/10_5486_PMD_2000_2282.pdf));
  Bimpeh–Horgan 2002 (simulation); Clayton–McMullen 2007 offset bound;
  Dworin–Grimlund; Fienberg–Neter–Leitch 1977 (different bound);
  Higgins–Nandram 2009; Rohrbach 1993; Meeden–Sargent stepwise-Bayes
  (explicitly not a frequentist proof). None contains a finite-sample proof.
- **Pap–vZ 2000 "my dollar right or wrong" randomized bound** is exactly
  valid for all n but is a different (randomized) estimator.
- **Bentkus–van Zuijlen program** (2003–2007) built provably valid
  Hoeffding-type bounds *around* the problem; generally weaker than SB.
- **Modern replacements** (Waudby-Smith–Ramdas 2024 betting martingales;
  Shekhar et al., arXiv:2305.06884) cite the conjecture as open and route
  around it. Phan–Thomas–Learned-Miller (ICML 2021,
  [arXiv:2106.03163](https://arxiv.org/pdf/2106.03163)) states it
  explicitly: coverage for α < 0.5 "is unknown".
- **Gaffke's bounded-mean interval** is now known to be finite-sample valid
  (Vlassis–Thomas, arXiv:2607.08415). It agrees with Clopper–Pearson on
  Bernoulli samples. This repository uses pointwise comparison with it to
  prove the three stated cases from \(n=3\) through \(n=7\). No corresponding all-\(n\)
  domination theorem was found.
- **Simplex slicing and barycentric-cut literature** supplies the classical
  cap-volume formula and the general characterization of critical
  halfspace-volume directions by section centroids. Webb's central-simplex
  extremum, its 2026 stability refinement, and recent barycentric-cut work
  solve different geometric problems. The September 8 search found no
  theorem maximizing the cap through an arbitrary ordered barycentric point
  over its matching order cone. The repository proves that statement in
  simplex dimension three; the proposed extensions to dimensions four and
  five have the active-prefix obligation recorded above.
- **arXiv metadata search for "Stringer bound": zero papers** (checked via
  API, 2026-08).
- Dutch survey (Hendriks & Kloosterman, [Essaybundel Statistical Auditing
  2021](https://steekproeven.eu/wp-content/uploads/2021/01/Essaybundel-Statistical-Auditing-pdf.pdf)):
  "De publicaties leiden niet tot een analytisch bewijs voor de Stringer
  Bound."

## Remaining mathematical target

The first unresolved binomial-factor sample size at all three conventional
levels is \(n=8\). A dimension-free geometric comparison would have direct
consequences for ordinary audit sample sizes.
[`ORDERED-SIMPLEX-CAP.md`](../theory/ORDERED-SIMPLEX-CAP.md) records the exact
tight-vertex reduction, the proved one- and two-cap regions, and an open
adjacent-transfer identity.
[`ALL-N-POISSON-PROGRAM.md`](../theory/ALL-N-POISSON-PROGRAM.md) records a
separate reduction to two explicitly unproved high-quantile inequalities.
