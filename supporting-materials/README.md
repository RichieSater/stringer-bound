# Supporting materials: verification guide

**Canonical repository:** [github.com/RichieSater/stringer-bound](https://github.com/RichieSater/stringer-bound) · **Archived release:** [doi:10.5281/zenodo.21850820](https://doi.org/10.5281/zenodo.21850820)

All computations are Python 3.9+ with `mpmath`, `numpy`, `scipy`. Every
certified binomial counterexample is decided by rational arithmetic;
float64 appears in screening searches and in separately labeled numerical
Poisson checks.

## Layers

| Layer | Role | Trust status |
|---|---|---|
| `stringer.py` | Numerical factors for searches; for certification, dyadic binomial-factor intervals whose endpoint CDF signs are evaluated exactly with integers | proof-essential |
| `coverage.py` / `coverage_exact` | Exact multinomial weights and rational propagation of factor intervals through every Stringer comparison | proof-essential |
| `two_point_lemma.py` | Written proof that single-value supports cannot under-cover, plus machine checks | proof-essential |
| `certify.py` | The only source of claims: exact coverage, exact nominal comparison, margin certificate | proof-essential |
| `search_two_value.py`, `search_multi_value.py` | float64 screening only | heuristic |
| `bolshev.py` | Reproduces Bimpeh's Table 5.1 and demonstrates that his coverage bound (5.16) is not a coverage bound (`audit/BIMPEH-GAP.md`) | proof-essential |
| `bimpeh_continuous_check.py` | MC corroboration of the hand counterexamples to (5.16) with continuous F | corroboration only |

## Commands

```sh
cd computations/python

# Lemma: single-value supports; high-precision Poisson comparison
python3 two_point_lemma.py --alpha 0.05 --n-max 40

# Known finite-sample violation at low confidence (machinery true-positive)
python3 search_two_value.py --alpha 0.7 --n 50 --out /tmp/c.json
python3 certify.py /tmp/c.json    # expect CONFIRMED lines

# The conjecture at 95% over two-value supports
python3 search_two_value.py --alpha 0.05 --n 2 30 --range --out /tmp/c95.json
python3 certify.py /tmp/c95.json  # expect "nothing to certify" if no dips

# Richer supports
python3 search_multi_value.py --alpha 0.05 --m 3 --n 10 20 --out /tmp/c3.json

# Exact-sign, interval-propagation, and table-generation regression tests
python3 -m unittest discover -s ../tests -v
```

Certified run logs are committed under `computations/certificates/`.
`summarize_certificates.py` regenerates the rows used in the manuscript:

```sh
python3 summarize_certificates.py \
  --out ../certificates/certificate-summary.json
```

## Certificate semantics

A `CONFIRMED` line from `certify.py` states: for the printed rational taint
distribution, the exact rational coverage is below the exact nominal level.
For each binomial factor, the code locates adjacent dyadic endpoints and
evaluates the sign of the binomial CDF minus \(\alpha\) at both endpoints
with integer arithmetic. It then propagates those rational intervals through
each Stringer-bound comparison. If an interval overlaps the exact rational
mean, certification stops and requests a finer dyadic grid. The numerical
root locator affects speed only; its proposed bracket is checked exactly and
an exact grid bisection is the fallback.

Poisson-factor comparisons are high-precision numerical checks, not formal
interval certificates; see
`computations/certificates/poisson-domination-standard-levels.log`.

## What screening output does NOT establish

A screening minimum equal to nominal (slack \(\sim 10^{-13}\)) does not
prove the infimum equals the nominal level; it is evidence subject to grid
and optimizer limitations. Negative results (no candidate found) bound only
the families and sample sizes actually searched.
