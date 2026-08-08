# Independent review record: the n = 2 theorem

Two independent AI referees were tasked to break
[`theory/N2-PROOF.md`](../theory/N2-PROOF.md) (adversarial brief:
attack every step, construct counterexamples, label findings
FATAL / GAP / COSMETIC). Author review and agreement among AI systems
are not substitutes for independent human mathematical review; this
record documents what was checked and what changed.

## Referee 1 (step-by-step attack) — verdict: theorem stands

- **FATAL: none.** Every proof step verified independently, including:
  the factor values and telescoping reduction; the w = β Markov step
  (requires A > 0, proved); the rectangle lemma's inclusion–exclusion
  and closed-event conventions; the left-limit mechanics of bound
  family (iii) under atoms of F; h's monotonicity structure (critical
  point exactly at g = B, so "h ≥ 0 on [B, β]" is rigorous); the
  concavity/endpoint argument for the potential (checked symbolically
  and at α as extreme as 10⁻¹² in 60-digit arithmetic); every endpoint
  and orientation of the tail change of variables (exact bijection in
  both w ≤ B and w > B subcases, no double counting); exhaustiveness
  and each strict link of the three-case endgame.
- **GAP (fixed):** the sharpness clause "approached but not attained"
  was asserted without proof. Fixed by a forced-equality analysis: an
  attaining F would satisfy the weak rectangle hypothesis, the budget
  chain then pins F = (1−β)δ₀ + βδ_{u*}, whose failure event has
  probability 0 — contradiction. Now a proved paragraph
  ("Non-attainment") in the proof file.
- **COSMETIC (all fixed):** "equivalent to" → "implied by" for the
  strict-inequality form (with the ε-scaling remark showing actual
  equivalence); h₊(1) := 0 now defined where previously used
  implicitly; the wing-family sharpness computation restated precisely;
  the reduction identity check now runs at three values of α.
- Independent adversarial search (structured boundary scans within
  10⁻⁸ of the critical masses q = β, q = B, plus 6000 random
  multistarts with Nelder–Mead polishing, α from 0.001 to 0.999):
  no violation of the theorem or any intermediate claim;
  min[P(< w) − γ²] ≥ +6×10⁻¹⁰, converging to 0 exactly along the
  predicted extremal families.

## Referee 2 (independent reconstruction) — verdict: CORRECT, with an independent second proof

Brief: verify or refute the wedge inequality from scratch — own
numerics, own proof attempt — before reading the repo's proof, then
compare.

- **Independent numerics** (own code, no repo code:
  [`computations/independent/wedge_num.py`](../computations/independent/wedge_num.py),
  [`wedge_num2.py`](../computations/independent/wedge_num2.py)): four
  methods (exhaustive 2-atom with analytic masses ~10⁶ configs per
  (α,w); structured 3-atom grids over all critical points; Frank–Wolfe
  on 241-point moment-polytope grids with 40 restarts; random-support
  hill-climbing). Sup of the closed failure event = exactly α, attained
  at the two predicted families, never exceeded; strict event
  approaches α, never attains it.
- **Independent second proof, different architecture**: direct
  decomposition P(fail) = g² + 2∫G(h(u))dF(u) with g = G(w/β), a
  mass-counting tail bound, and a linear majorant making the bound
  linear in the high-part mean m — closed by the inequality
  J(g) = A(β+g) − 2B(1−g) ≥ 0 with J(B) = β²+γ²−1 = 0. sympy confirms
  J(g) ≥ 0 is **algebraically equivalent** to the repo proof's
  potential bound k̃(g) ≥ β: two structurally different proofs hinge on
  the same tight inequality with the same equality points g ∈ {B, β},
  matching the two extremal families. Chain verified against 16,000
  random distributions, zero violations.
- **Repo proof re-derivation**: every step reconstructed and confirmed
  (reduction; degenerate regimes; rectangle lemma; families (i)–(iii)
  including the left-limit and h₊(1) traps; potential identities;
  change of variables in both subcases; all three strictness cases;
  non-attainment).
- **Cosmetic findings (both fixed)**: (&#42;) needed the w > 0
  restriction stated at the display; the non-attainment paragraph
  claimed "a₀ = 0 is forced," which was unjustified and unnecessary —
  the equality links force G = β a.e. regardless of a₀.

## Combined conclusion

Two independent adversarial reviews, one containing a full independent
proof by a different method, concur: the theorem and its sharpness
statement are correct. This does not replace human peer review; it
raises the bar an error would have had to slip past.

## Post-review verifier status

`n2_proof_check.py` (four layers: reduction identity at three α,
symbolic identities, potential bound, ~200k-distribution stress test)
passes in full after the fixes; log in
`computations/certificates/n2-proof-check.log`.
