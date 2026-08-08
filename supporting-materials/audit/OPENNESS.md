# Openness verification

**Verdict (verified 2026-08-07): the finite-sample conservatism conjecture
for the Stringer bound at confidence levels ≥ 50% — including 95% — is
open.** No proof and no counterexample exists in the literature at α ≤ 1/2.
The strongest partial results, and where a resolution could have hidden,
are recorded below with sources.

## What is proved

| Result | Source | Scope |
|---|---|---|
| Asymptotic conservatism for α ∈ (0, 1/2]; asymptotic **anti**-conservatism for α ∈ (1/2, 1). Threshold exactly α = 1/2 | Pap & van Zuijlen, *Statistica Neerlandica* 50 (1996) 367–389 ([1994 report full text](http://web.archive.org/web/20251019123811/https://repository.ubn.ru.nl/bitstream/handle/2066/100177/1/100177.pdf)) | asymptotic, a.s. expansions, binomial factors |
| Bickel's expansion: SB = T̄ + c(F) z₁₋α/√n + o(n^{-1/2}), with c²(F) ≥ Var(T), equality iff F is supported on ≤ 2 points | Bickel, *Int. Statist. Rev.* 60 (1992) 197–209; made rigorous in Pap–vZ 1996 | asymptotic |
| P(SB ≥ μ) ≥ (1−α)^{n+1} under conditions on F (hypotheses unverified — JSTOR original not yet read) | Bickel 1992 | finite n, weak bound |
| Conjecture holds **exactly** for taints supported on {0,1} (SB reduces to Clopper–Pearson); Stringer coefficients are minimal with this property | de Jager, Pap & van Zuijlen, *Comput. Math. Appl.* 33 (1997) 37–54 | finite n, {0,1} supports |
| Claimed certification at n ≤ 11, α = 0.05 via the lower bound P̄ₙ = P(U_{i:n} ≤ p_n(i) ∀i) (Bolshev recursion). **Refuted in this repository**: the containment CP ≥ P̄ₙ (his eq. 5.16) rests on an off-by-one in the band constraints, dropping F(t_{n:n}) ≥ α^{1/n}; hand counterexamples with continuous F at n = 1, 2 and an exact atomic one at n = 5, and the corrected containment probability is ≤ 1−α for every continuous F, so the method can certify nothing. See [`BIMPEH-GAP.md`](BIMPEH-GAP.md). His Table 5.1 values reproduce exactly; they are just not coverage bounds | Bimpeh, PhD thesis, DCU 2008, ch. 5 ([full text](https://doras.dcu.ie/600/1/YawThesis.PDF)); apparently never journal-published; refutation: this repository | **invalid — the finite-sample record reduces to n = 1 (Bickel), {0,1} supports, and single-atom supports** |
| Exact counterexamples at **sub-50% confidence**: uniform-type taints, α > 1/2, n ≳ 17 | Pap & van Zuijlen, *Comput. Math. Appl.* 29 (1995) 51–59 (exact recursions, not Monte Carlo) | finite n, α > 1/2 only |

This repository adds (see `computations/certificates/`): exact rational
two-atom counterexamples at α = 0.7 (n = 50) and α = 0.68 (n = 100),
margin-certified; and screening evidence that over two- and three-atom
supports at α = 0.05, n ≤ 100, the infimum of coverage is exactly 1−α,
approached on the {0,1}-boundary (v₁ → 1) — consistent with the
de Jager–Pap–van Zuijlen minimality theorem.

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
   possibly the right lens for n ≥ 12.
