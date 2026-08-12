# Openness verification

**Verdict (updated 2026-08-12): the general finite-sample conservatism
conjecture for the Stringer bound at confidence levels ≥ 50% — including
95% — remains open.** The present repository proves the case \(n=2\) at
every confidence level and gives exact computer-assisted proofs from
\(n=3\) through \(n=5\) for 90%, 95%, and 99% confidence, plus \(n=6\)
at 95%. Thus \(n=6\) is the first unresolved binomial-factor sample size
at 90% and 99%, while \(n=7\) is the first unresolved size at 95%.
For the Poisson factors used in audit practice, a separate exact
simultaneous-band argument proves coverage through \(n=8\), \(n=11\), and
\(n=20\), respectively. This does not resolve every confidence level at
\(n=3\) through \(n=6\), nor arbitrary \(n\). No prior general proof or
counterexample at α ≤ 1/2 was found in the reviewed literature. The
strongest partial results, and where a resolution could have hidden, are
recorded below with sources.

## What is proved

| Result | Source | Scope |
|---|---|---|
| Asymptotic conservatism for α ∈ (0, 1/2]; asymptotic **anti**-conservatism for α ∈ (1/2, 1). Threshold exactly α = 1/2 | Pap & van Zuijlen, *Statistica Neerlandica* 50 (1996) 367–389 ([1994 report full text](http://web.archive.org/web/20251019123811/https://repository.ubn.ru.nl/bitstream/handle/2066/100177/1/100177.pdf)) | asymptotic, a.s. expansions, binomial factors |
| Bickel's expansion: SB = T̄ + c(F) z₁₋α/√n + o(n^{-1/2}), with c²(F) ≥ Var(T), equality iff F is supported on ≤ 2 points | Bickel, *Int. Statist. Rev.* 60 (1992) 197–209; made rigorous in Pap–vZ 1996 | asymptotic |
| P(SB ≥ μ) ≥ (1−α)^{n+1} under conditions on F (hypotheses unverified — JSTOR original not yet read) | Bickel 1992 | finite n, weak bound |
| Conjecture holds **exactly** for taints supported on {0,1} (SB reduces to Clopper–Pearson); Stringer coefficients are minimal with this property | de Jager, Pap & van Zuijlen, *Comput. Math. Appl.* 33 (1997) 37–54 | finite n, {0,1} supports |
| Binomial Stringer is distribution-free conservative at \(n=2\) for every α, sharply; from \(n=3\) through \(n=5\), it pointwise dominates the valid Gaffke bounded-mean upper limit at α = 0.10, 0.05, and 0.01; the same domination holds at \(n=6\), α = 0.05 | this repository: [`N2-PROOF.md`](../theory/N2-PROOF.md), [`N3-CONVENTIONAL.md`](../theory/N3-CONVENTIONAL.md), [`N4-CONVENTIONAL.md`](../theory/N4-CONVENTIONAL.md), [`N5-CONVENTIONAL.md`](../theory/N5-CONVENTIONAL.md), [`N6-CONVENTIONAL.md`](../theory/N6-CONVENTIONAL.md), and exact certificates | finite n; all levels at \(n=2\), 90%, 95%, and 99% at \(n=3,4,5\), and 95% at \(n=6\) |
| Poisson audit factors dominate binomial Clopper–Pearson factors coordinatewise whenever nominal confidence exceeds \(1-e^{-1}\approx63.2\%\); hence the Poisson Stringer bound pointwise dominates the binomial version | Anderson–Samuels (1967), specialized to the Stringer formula in [`POISSON-DOMINATION.md`](../theory/POISSON-DOMINATION.md) | every n; factor and bound comparison, not a general coverage proof |
| A corrected simultaneous survival band proves Poisson-factor Stringer coverage for every \(n\le8\) at 90%, every \(n\le11\) at 95%, and every \(n\le20\) at 99% confidence | [`POISSON-SIMULTANEOUS-BAND.md`](../theory/POISSON-SIMULTANEOUS-BAND.md) and exact rational certificate | finite ranges; arbitrary continuous or atomic taint distributions |
| The pre-specified maximum of either Stringer calculation and the valid Gaffke limit has coverage at least \(1-\alpha\) | [`GAFFKE-SAFEGUARD.md`](../theory/GAFFKE-SAFEGUARD.md) and the Vlassis–Thomas validity theorem | every n and level, but a safeguarded procedure rather than proof of ordinary Stringer |
| At every nominal confidence level of at least \(75\%\), the safeguard has zero uplift for every sample size whenever binomial Stringer is at least the largest observed taint | analytic proof in [`ONE-CAP-COMPARISON.md`](../theory/ONE-CAP-COMPARISON.md); 59,700 exact nonterminal vertex checks through \(n=200\) retained as an independent regression | pointwise one-cap region only; not conditional coverage or a full ordinary-Stringer theorem |
| Claimed certification at n ≤ 11, α = 0.05 via the lower bound P̄ₙ = P(U_{i:n} ≤ p_n(i) ∀i) (Bolshev recursion). **Reassessed in this repository**: the containment CP ≥ P̄ₙ (his eq. 5.16) rests on an off-by-one in the band constraints, dropping F(t_{n:n}) ≥ α^{1/n}; hand counterexamples with continuous F at n = 1, 2 and an exact atomic one at n = 5 show that P̄ₙ is not the stated coverage bound. The corrected containment probability is ≤ 1−α for every continuous F. See [`BIMPEH-GAP.md`](BIMPEH-GAP.md). | Bimpeh, PhD thesis, DCU 2008, ch. 5 ([full text](https://doras.dcu.ie/600/1/YawThesis.PDF)); apparently never journal-published; reassessment: this repository | the cited argument does not establish the claimed finite-sample range |
| Exact counterexamples at **sub-50% confidence**: uniform-type taints, α > 1/2, n ≳ 17 | Pap & van Zuijlen, *Comput. Math. Appl.* 29 (1995) 51–59 (exact recursions, not Monte Carlo) | finite n, α > 1/2 only |

This repository adds (see `theory/N2-PROOF.md`,
`theory/N3-CONVENTIONAL.md`, `theory/N4-CONVENTIONAL.md`,
`theory/N5-CONVENTIONAL.md`,
`theory/N6-CONVENTIONAL.md`,
`theory/POISSON-SIMULTANEOUS-BAND.md`, and
`computations/certificates/`): a
distribution-free proof at \(n=2\); exact computer-assisted proofs from
\(n=3\) through \(n=5\) for α = 0.10, 0.05, and 0.01; an exact proof at
\(n=6\), α = 0.05; direct exact
Poisson-factor guarantees through the level-dependent ranges above; and 33 exact rational
counterexamples supported on three points (two nonzero taint values and
zero) at selected nominal confidence levels from 30% to 37% and sample
sizes 50, 100, 200, and 400. Numerical searches at α = 0.05 found smallest
coverage agreeing with \(1-\alpha\) to search precision, approached on the
{0,1}-boundary (\(v_1\to1\)); this is evidence, not an infimum proof.
The all-sample-size safeguard and its certified one-cap identity region do
not change the openness verdict because the former validates a different,
pre-specified maximum rule and the latter controls only part of the sample
space.

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
  prove the three stated cases from \(n=3\) through \(n=5\), and the 95%
  case at \(n=6\). No corresponding all-\(n\)
  domination theorem was found.
- **arXiv metadata search for "Stringer bound": zero papers** (checked via
  API, 2026-08).
- Dutch survey (Hendriks & Kloosterman, [Essaybundel Statistical Auditing
  2021](https://steekproeven.eu/wp-content/uploads/2021/01/Essaybundel-Statistical-Auditing-pdf.pdf)):
  "De publicaties leiden niet tot een analytisch bewijs voor de Stringer
  Bound."

## Open flags needing a deeper read

1. Bickel 1992 original (JSTOR): exact hypotheses of the (1−α)^{n+1}
   theorem.
2. Pap–vZ 1995 full text (Elsevier): exact (n, α) table of their
   counterexample.
3. Bimpeh ch. 5: the F̂-band containment step for atomic F (tie handling);
   this repository re-derives the recursion and cross-checks P̄ₙ against
   exact two-atom coverage (`bolshev.py`).
4. Pyke 1994, one-sided minimax KS bands (*J. Appl. Prob.* 31A, 291–308):
   possibly the right lens for larger \(n\).
5. Extension of the binomial Gaffke/simplex-cap certificate to 90% and 99%
   at \(n=6\), and to 95% at \(n=7\); a
   dimension-free geometric comparison would have direct consequences for
   ordinary audit sample sizes. [`ORDERED-SIMPLEX-CAP.md`](../theory/ORDERED-SIMPLEX-CAP.md)
   records the exact tight-vertex reduction and an open adjacent-transfer
   identity. [`ALL-N-POISSON-PROGRAM.md`](../theory/ALL-N-POISSON-PROGRAM.md)
   records a separate reduction to two explicitly unproved high-quantile
   inequalities.
