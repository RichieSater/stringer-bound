# The Stringer Bound Conjecture

[![DOI](https://zenodo.org/badge/1327724172.svg)](https://doi.org/10.5281/zenodo.21850819)
[![Verify paper and certificates](https://github.com/RichieSater/stringer-bound/actions/workflows/verify.yml/badge.svg)](https://github.com/RichieSater/stringer-bound/actions/workflows/verify.yml)

**Repository:** [github.com/RichieSater/stringer-bound](https://github.com/RichieSater/stringer-bound) · **Archived v1.0.0:** [doi:10.5281/zenodo.21850820](https://doi.org/10.5281/zenodo.21850820)

> **Longstanding finite-sample validity question for a bound introduced in
> 1963.** For every
> taint distribution \(F\) on \([0,1]\), every sample size \(n\), and
> standard confidence levels such as \(95\%\), the Stringer bound covers the
> population mean taint with probability at least the nominal level:
> \(P(\mathrm{SB} \ge \theta) \ge 1-\alpha\).

The Stringer bound is widely used to evaluate monetary-unit samples
(Stringer 1963). Its conjectured *conservatism* in finite samples has not
been proved in full generality. Bickel (1992) proved asymptotic
conservatism; Pap and van Zuijlen showed via higher-order expansions that
asymptotic conservatism holds only above a confidence threshold, leaving the
finite-sample question at standard levels like 95% open. This repository is
a systematic, certificate-backed attempt to resolve it: either an exact
counterexample (a taint distribution and sample size where coverage provably
dips below nominal) or structural evidence for the conjecture.

## The bound

Monetary-unit sampling draws \(n\) dollar units i.i.d. from the population;
the taint of a drawn unit is \(T \in [0,1]\) with distribution \(F\), and
the parameter is \(\theta = E[T]\). With observed taints arranged in
nonincreasing order, \(t_{(1)} \ge \dots \ge t_{(n)}\), the Stringer bound
at nominal confidence
\(1-\alpha\) is

```
SB = p_0 + sum_{j>=1} (p_j - p_{j-1}) * t_(j)
```

where \(p_j\) is the upper \(1-\alpha\) confidence limit for a binomial
proportion after \(j\) errors in \(n\) trials (Clopper–Pearson), or its
Poisson analogue \(m_j/n\) (the AICPA "MUS factors"). Every claim in this
repository states which variant it uses.

## Status

**Exact finite-sample guarantees now reach \(n=5\). The \(n=2\) bound is
proved conservative for every taint distribution and every confidence level
([`N2-PROOF.md`](supporting-materials/theory/N2-PROOF.md)). From \(n=3\)
through \(n=5\), exact computer-assisted proofs establish conservatism at 90%,
95%, and 99% confidence
([`N3-CONVENTIONAL.md`](supporting-materials/theory/N3-CONVENTIONAL.md),
[`N4-CONVENTIONAL.md`](supporting-materials/theory/N4-CONVENTIONAL.md),
and [`N5-CONVENTIONAL.md`](supporting-materials/theory/N5-CONVENTIONAL.md)).
At those three levels, \(n=6\) is the first unresolved sample size. The
general-\(n\) conjecture at 95% remains open.**

**A validated all-sample-size reporting safeguard is also available.**
Pre-specify the reported upper bound as the maximum of the familiar
Stringer calculation and the finite-sample-valid Gaffke bounded-mean limit.
That maximum has distribution-free coverage for every \(n\) without assuming
the unresolved Stringer conjecture. It returns the ordinary Stringer value
whenever Stringer is larger; in the certified \(n=3,4,5\) cases, equality
with ordinary Stringer holds for every sample. The implementation uses exact
rational tail-sign checks rather than trusting a floating-point quantile;
see [`GAFFKE-SAFEGUARD.md`](supporting-materials/theory/GAFFKE-SAFEGUARD.md).

Openness was checked against the literature (details and sources in
[`OPENNESS.md`](supporting-materials/audit/OPENNESS.md)): no general proof
or counterexample at confidence levels ≥ 50% was found. The asymptotic
threshold is exactly \(\alpha = 1/2\) (Pap–van Zuijlen 1996); the only
published finite-sample counterexamples found in the review are at sub-50%
confidence; and the strongest previously claimed finite-sample range
(Bimpeh 2008: all \(F\) at \(n \le 11\), \(\alpha = 0.05\)) rests on an
index shift in its containment argument — see below.

Established so far in this repository, each backed by an exact-arithmetic
certificate or a written proof:

- **The \(n = 2\) case is proved in full**
  ([`N2-PROOF.md`](supporting-materials/theory/N2-PROOF.md)): for every
  \(F\) and every \(\alpha\in(0,1)\). The supremum noncoverage probability
  is exactly \(\alpha\), approached but never attained.

- **Conventional confidence levels are certified from \(n=3\) through \(n=5\)**
  ([`N3-CONVENTIONAL.md`](supporting-materials/theory/N3-CONVENTIONAL.md),
  [`N4-CONVENTIONAL.md`](supporting-materials/theory/N4-CONVENTIONAL.md),
  [`N5-CONVENTIONAL.md`](supporting-materials/theory/N5-CONVENTIONAL.md)):
  at 90%, 95%, and 99% confidence, the binomial Stringer bound pointwise
  dominates the recently validated Gaffke bounded-mean upper limit for all
  three sample sizes. The proofs reduce the comparisons to uniform-simplex cap
  inequalities and certify them with exact rational Bernstein coefficients.
  The symbolic derivations and every sign certificate are regenerated from
  source by the top-level reproduction command. These are theorems at the
  three listed levels, not search evidence.

- **All-sample-size Poisson-factor domination at practical confidence
  levels**
  ([`POISSON-DOMINATION.md`](supporting-materials/theory/POISSON-DOMINATION.md)):
  whenever nominal confidence exceeds \(1-e^{-1}\approx63.2\%\), the
  Poisson factors used in audit practice dominate the binomial
  Clopper--Pearson factors coordinatewise for every \(n\). A
  summation-by-parts argument then proves that the complete Poisson
  Stringer bound is pointwise at least the binomial version on every
  observed sample. This replaces the former finite-range numerical
  comparison at 90% and 95% with an analytic result; it does not by itself
  resolve general-\(n\) coverage. It does transfer the proved binomial
  guarantees at \(n=2\), and from \(n=3\) through \(n=5\) for 90%, 95%, and
  99%, to the Poisson-factor bound.

- **All-sample-size safeguarded reporting rule**
  ([`GAFFKE-SAFEGUARD.md`](supporting-materials/theory/GAFFKE-SAFEGUARD.md)):
  for either factor convention, report
  \(\max\{\mathrm{Stringer},\mathrm{Gaffke}\}\). Finite-sample validity of
  the Gaffke component proves coverage of this complete rule for every
  sample size and every confidence level. This does not prove the ordinary
  Stringer conjecture. It is a drop-in statistical floor that preserves the
  familiar calculation whenever Stringer is already larger. The command-line
  implementation certifies a rational dyadic bracket for the Gaffke endpoint
  by exact confluent divided differences, including repeated taints and
  zero-heavy audit samples. A methodology-facing use and scope note is in
  [`PRACTICE-SAFEGUARD.md`](supporting-materials/audit/PRACTICE-SAFEGUARD.md).

- **Two-point lemma** (proof in
  [`two_point_lemma.py`](supporting-materials/computations/python/two_point_lemma.py)):
  a distribution with a single nonzero taint value can never under-cover —
  the coverage event contains the Clopper–Pearson coverage event. So
  \(\{v_1 > v_2 > 0\}\) plus an atom at 0 is the smallest support that can
  carry a counterexample.
- **The bound genuinely fails at low confidence in finite samples**:
  thirty-three exact rational counterexamples with rational factor-interval
  certificates — at selected nominal levels from 30% through 37% and
  sample sizes \(n=50,100,200,400\), using binomial factors.
- **A previous finite-sample certification is reassessed**
  ([`BIMPEH-GAP.md`](supporting-materials/audit/BIMPEH-GAP.md)): Bimpeh's
  (2008) coverage lower bound \(\bar P_n\) — the basis of the belief that
  the 95% conjecture was settled for \(n \le 11\) — rests on an
  off-by-one in the confidence-band constraints. Hand-verifiable
  counterexamples with *continuous* \(F\) at \(n = 1, 2\), an exact
  atomic one at \(n = 5\) (coverage \(31/32 < \bar P_5\)), and a proof
  that the corrected containment probability never exceeds \(1-\alpha\),
  so that containment argument cannot establish the desired guarantee.
  His Table 5.1 itself reproduces
  exactly ([`bolshev.py`](supporting-materials/computations/python/bolshev.py)).
  The cited argument therefore does not add a proven finite-sample range.
- **At 95% no violation was found in the reported searches**: over
  distributions with two nonzero taint values (\(n \le 100\)) and three
  nonzero taint values (\(n \le 25\)), each allowing an atom at zero, the
  smallest coverage found agrees with
  \(1-\alpha\) to numerical precision at every searched \(n\), approached
  as \(v_1 \to 1\) and never crossed — consistent with
  the de Jager–Pap–van Zuijlen minimality theorem on \(\{0,1\}\)
  supports. Single-value populations already reach exact coverage
  \(0.9500302\) (\(n=10\), \(v=1\), \(q=0.85\)).

**Next steps**: determine whether the Gaffke/simplex-cap comparison extends
to \(n=6\), where direct Bernstein certification becomes substantially
larger; pursue the two dimension-free weighted-exponential and
Dirichlet-Poissonization inequalities isolated in
[`ALL-N-POISSON-PROGRAM.md`](supporting-materials/theory/ALL-N-POISSON-PROGRAM.md);
and develop a
certified
branch-and-bound or atoms-reduction argument for ordinary audit sample sizes.
The immediate mathematical target is \(n=6\) at 90%, 95%, and 99%, not more
unstructured grid search. Before journal submission or practice-facing
reliance, the bounded review protocol in
[`HUMAN-REVIEW-PACKET.md`](supporting-materials/audit/HUMAN-REVIEW-PACKET.md)
should be completed by an independent human mathematical statistician.

## Method

Float64 screening over parametric families of taint distributions
(grid + Nelder–Mead refinement), followed by exact recertification of every
candidate: multinomial weights in rational arithmetic, confidence factors
enclosed between dyadic rationals with binomial-CDF endpoint signs evaluated
exactly by integer arithmetic, and rational interval propagation through
every coverage comparison. Nothing is claimed from screening output; only
[`certify.py`](supporting-materials/computations/python/certify.py) verdicts
count.

## Repository layout

```text
supporting-materials/
├── claim-evidence.json       machine-readable claim-to-evidence map
├── theory/
│   ├── N2-PROOF.md           complete n=2 proof
│   ├── N3-CONVENTIONAL.md    exact n=3 proof at 90%, 95%, and 99%
│   ├── N4-CONVENTIONAL.md    exact n=4 proof at 90%, 95%, and 99%
│   ├── N5-CONVENTIONAL.md    exact n=5 proof at 90%, 95%, and 99%
│   ├── POISSON-DOMINATION.md all-n practical-level factor comparison
│   ├── GAFFKE-SAFEGUARD.md   all-n valid reporting floor and exact computation
│   └── ALL-N-POISSON-PROGRAM.md exact reductions for the open all-n target
└── computations/python/
    ├── stringer.py           numerical factors for searches + exact-sign
    │                         dyadic binomial factor intervals
    ├── coverage.py           rational interval certification + float screener
    ├── two_point_lemma.py    proof + machine check: single-value supports
    │                         cannot under-cover
    ├── derive_n3_bernstein_formulas.py
    │                         symbolic derivation of the n=3 certificate
    ├── n3_gaffke_certificate.py
    │                         exact rational n=3 sign certificate
    ├── derive_n4_bernstein_structure.py
    │                         exact n=4 residual and tetrahedral structure
    ├── n4_gaffke_certificate.py
    │                         exact rational n=4 sign certificate
    ├── derive_n5_bernstein_structure.py
    │                         exact n=5 residual, face-ideal, and four-simplex structure
    ├── n5_gaffke_certificate.py
    │                         directed-dyadic n=5 sign certificate
    ├── gaffke.py             exact-sign Gaffke endpoint + safeguarded report
    ├── all_n_poisson_reductions.py
    │                         exact algebra and rejected-shortcut check for all-n route
    ├── search_two_value.py   screening search over {v1 > v2 > 0} supports
    └── certify.py            exact recertification of screening candidates
```

## Reproduction

The canonical environment is pinned in `pyproject.toml` and `uv.lock`
(Python 3.12, `uv` 0.11.14). Tectonic 0.17.0 builds the manuscript.
The single top-level verification command is:

```sh
make reproduce
```

This runs the unit tests, the \(n=2\) symbolic proof checker, the exact
algebra checks for the explicitly open all-\(n\) Poisson route, source
regeneration of the \(n=3\) through \(n=5\) Bernstein structures and exact
sign certificates, the generated counterexample-table check,
claim-to-evidence link validation, and the manuscript build. Individual
search and recertification commands
remain documented in
[`supporting-materials/README.md`](supporting-materials/README.md).

For a directly usable safeguarded audit calculation (zero taints may be
omitted while `--n` remains the full sample size):

```sh
uv run --frozen python \
  supporting-materials/computations/python/gaffke.py \
  --n 100 --alpha 0.05 --method poisson --taints 1,0.4,0.1
```

## Openness

Whether the finite-sample conjecture at 95% is genuinely open — and whether
a resolution hides in the modified-Stringer-bound literature — is tracked in
[`supporting-materials/audit/OPENNESS.md`](supporting-materials/audit/OPENNESS.md).

## AI disclosure

Anthropic Claude and OpenAI Codex assisted with mathematical exploration,
code, literature searches, and verification workflows. The author directed
and reviewed the work and accepts responsibility for the contents.

## License

Code is MIT-licensed under [`LICENSE`](LICENSE).
