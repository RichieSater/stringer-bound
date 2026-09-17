# Supporting mathematics and computations

## Principal article

[Ordered simplex caps and the three-observation Stringer bound](paper/stringer.pdf)
([source](paper/stringer.tex)) proves the four-coordinate cap theorem, its
sharp three-observation Stringer comparison, and an all-dimensional
terminal-edge obstruction. Its proof uses the following material:

| Mathematical object | Derivation and exact calculation |
|---|---|
| Tetrahedral section inequality and repeated-upper boundary | [theory](theory/TETRAHEDRAL-VERTEX-BARRIER.md), [Python](computations/python/tetrahedral_vertex_barrier.py), [rational coefficients](computations/certificates/tetrahedral-vertex-barrier-certificate.json) |
| Terminal-edge derivative and critical levels | [theory](theory/TERMINAL-EDGE-BIFURCATION.md), [Python](computations/python/terminal_edge_bifurcation.py), [enclosures](computations/certificates/terminal-edge-bifurcation-certificate.json) |

Ordinary mathematics proves the stationary-face reduction, the
Clopper--Pearson weight inequalities, and the limiting obstruction. Python
with SymPy and rational arithmetic derives and checks the polynomial
identities and Bernstein signs. The statistical coverage implication uses
the cited validity theorem for Gaffke's mean test. No Lean, GAP, Rocq, or
other formal proof is part of this article.

```sh
make sync
make tetrahedral-vertex-barrier-check terminal-edge-bifurcation-check
make paper
```

Run these commands from the repository root. `make paper` compiles both
manuscripts and rejects layout and cross-reference warnings. Python 3.12
and the frozen `uv.lock` environment are used for all Python calculations.

## Higher-dimensional question

The global five- and six-coordinate monotone cap statements are
**conjectural**. The missing implication concerns repeated-lowest knot
blocks: their first-coordinate beta bounds need not be bounds at strict knot
gaps. [The remaining obligation](theory/ORDERED-COLLISION-OBLIGATION.md)
explains the distinction precisely. The local identities and
Clopper--Pearson weight-order calculations do not by themselves prove the
global statements or the complete confidence intervals at sample sizes four
and five.

- [Five-coordinate local calculations](theory/FIVE-COORDINATE-VERTEX-BARRIER.md).
- [Six-coordinate repeated-top calculations](theory/SIX-COORDINATE-REPEATED-TOP.md).
- [Six-coordinate unique-top calculations](theory/SIX-COORDINATE-UNIQUE-TOP.md)
  and [the other unique-top cases](theory/SIX-COORDINATE-UNIQUE-TOP-L2-L3.md).
- Proved weight-order calculations for [four trials](theory/N4-MONOTONE-RANGE.md)
  and [five trials](theory/N5-MONOTONE-RANGE.md).

Most local calculations use exact rational arithmetic. One six-coordinate
chart uses binary64 centers with exact rational forward-error bounds.
These statements describe local inequalities, not a completed
higher-dimensional maximum principle.

## Separate results retained in the research corpus

These results are outside the shortened principal article; they have not
been discarded or folded into its theorem.

- [The two-observation proof](theory/N2-PROOF.md) covers every confidence
  level and gives sharpness and nonattainment.
- Independent direct fixed-level cap calculations are in
  [N3](theory/N3-CONVENTIONAL.md), [N4](theory/N4-CONVENTIONAL.md),
  [N5](theory/N5-CONVENTIONAL.md), [N6](theory/N6-CONVENTIONAL.md), and
  [N7](theory/N7-CONVENTIONAL.md). These do not rely on the conjectural
  higher-dimensional monotone cap theorem. The largest calculations also
  use Singular and Arb; their individual notes state the exact propositions.
- [Poisson factor comparison](theory/POISSON-DOMINATION.md),
  [simultaneous bands](theory/POISSON-SIMULTANEOUS-BAND.md), and
  [calibrated bounds](theory/POISSON-BAND-CALIBRATION.md).
- [Gaffke safeguard](theory/GAFFKE-SAFEGUARD.md),
  [one-cap comparison](theory/ONE-CAP-COMPARISON.md),
  [two-cap collars](theory/TWO-CAP-COMPARISON.md), and
  [decision-scale examples](audit/DECISION-BENCHMARK.md).
- [Without-replacement conditioning](theory/SRSWOR-CONDITIONING.md).
- [Exact low-confidence counterexamples](computations/certificates/certificate-summary.json)
  and [the earlier band-containment issue](audit/BIMPEH-GAP.md).

The [systematic-PPS manuscript](paper/systematic-pps.pdf) is mathematically
separate: it studies a fixed ordered population and random-start inference.
Its proofs and computations are described in
[SYSTEMATIC-PPS](theory/SYSTEMATIC-PPS.md),
[SYSTEMATIC-MINIMAX](theory/SYSTEMATIC-MINIMAX.md), and
[SYSTEMATIC-DISJOINT-PHASES](theory/SYSTEMATIC-DISJOINT-PHASES.md).

## Checks and versions

`make test` runs the regression suite. `make public-corpus-check` checks
tracked public text; its mutation tests ensure excluded process language
and duplicate disclosures are rejected. `make claim-manifest-check` checks
source and evidence links; it is not a mathematical proof checker. The
[manifest](claim-evidence.json) records the separate scopes and limitations.

`make reproduce` runs the wider computational program, including expensive
higher-dimensional local calculations. Passing those calculations does not
resolve the active-prefix obligation above.

The v1.1.0 source package accompanies the
[GitHub release](https://github.com/RichieSater/stringer-bound/releases/tag/v1.1.0)
and the [Zenodo series](https://doi.org/10.5281/zenodo.21850819).
[Archived v1.0.0](https://doi.org/10.5281/zenodo.21850820) predates this article
and is preserved unchanged.
