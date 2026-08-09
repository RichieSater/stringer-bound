# The Stringer Bound Conjecture

[![DOI](https://zenodo.org/badge/1327724172.svg)](https://doi.org/10.5281/zenodo.21850819)

**Repository:** [github.com/RichieSater/stringer-bound](https://github.com/RichieSater/stringer-bound) · **Archived v1.0.0:** [doi:10.5281/zenodo.21850820](https://doi.org/10.5281/zenodo.21850820)

> **Conjecture (finite-sample conservatism, open since 1963).** For every
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

**Partial resolution: the \(n = 2\) case is proved
([`N2-PROOF.md`](supporting-materials/theory/N2-PROOF.md)) — to our
knowledge, the first finite-sample proof of Stringer conservatism at any
\(n \ge 2\), for all taint distributions and all confidence levels. The
general-\(n\) conjecture at 95% remains open.**

Openness was checked against the literature (details and sources in
[`OPENNESS.md`](supporting-materials/audit/OPENNESS.md)): no general proof
or counterexample at confidence levels ≥ 50% was found. The asymptotic
threshold is exactly \(\alpha = 1/2\) (Pap–van Zuijlen 1996); the only
published finite-sample counterexamples live at sub-50% confidence; and
the strongest finite-sample partial result (Bimpeh 2008: certification for
all \(F\) at \(n \le 11\), \(\alpha = 0.05\)) turns out to have an
essential gap — see below.

Established so far in this repository, each backed by an exact-arithmetic
certificate or a written proof:

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

- **The \(n = 2\) case is proved in full**
  ([`N2-PROOF.md`](supporting-materials/theory/N2-PROOF.md)): via the
  substitution \(U = 1-T\) the bound becomes \(\mathrm{SB} = 1 - (A\min U
  + B\max U)\) with \(A+B = \sqrt{\alpha}\), and conservatism reduces to
  a wedge-exceedance inequality proved by combining three families of
  rectangle bounds with the mean budget \(\int G = w\) and a concavity
  argument for a potential function with equal endpoint values. Valid
  for every \(F\) and every \(\alpha \in (0,1)\), so the certified
  low-confidence violations do not occur at \(n=2\).
  Sharp: the supremum of the failure probability is exactly \(\alpha\),
  approached, never attained. The algebra is symbolically checked and the
  statement is stress-tested in four computational layers
  ([`n2_proof_check.py`](supporting-materials/computations/python/n2_proof_check.py)).

**Next steps**: extend the rectangle-budget-concavity method to
\(n \ge 3\) (the reformulation \(\mathrm{SB} = 1 - \sum_j e_j u_{(j)}\)
is in place), a certified branch-and-bound over the parameter space with two
nonzero taint values for larger \(n\), and an atoms-reduction argument.

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
└── computations/python/
    ├── stringer.py           numerical factors for searches + exact-sign
    │                         dyadic binomial factor intervals
    ├── coverage.py           rational interval certification + float screener
    ├── two_point_lemma.py    proof + machine check: single-value supports
    │                         cannot under-cover
    ├── search_two_value.py   screening search over {v1 > v2 > 0} supports
    └── certify.py            exact recertification of screening candidates
```

## Reproduction

Requires Python 3.9+ with `mpmath`, `numpy`, `scipy`.

```sh
cd supporting-materials/computations/python
python3 two_point_lemma.py --alpha 0.05 --n-max 40      # lemma checks
python3 search_two_value.py --alpha 0.7 --n 50 --out c.json   # known violation
python3 certify.py c.json                                # exact confirmation
python3 search_two_value.py --alpha 0.05 --n 2 30 --range --out c95.json
python3 certify.py c95.json
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
