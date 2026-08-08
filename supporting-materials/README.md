# Supporting materials: verification guide

All computations are Python 3.9+ with `mpmath`, `numpy`, `scipy`. Every
number that supports a claim comes from exact arithmetic; float64 appears
only in screening searches whose output is never cited directly.

## Layers

| Layer | Role | Trust status |
|---|---|---|
| `stringer.py` | Confidence factors by bisection on provably monotone tails; each factor carries a rigorous bracket width (\(10^{-40}\)) | proof-essential |
| `coverage.py` / `coverage_exact` | Coverage as an exact `Fraction` by full multinomial enumeration, with a minimum-margin certificate | proof-essential |
| `two_point_lemma.py` | Written proof that single-value supports cannot under-cover, plus machine checks | proof-essential |
| `certify.py` | The only source of claims: exact coverage, exact nominal comparison, margin certificate | proof-essential |
| `search_two_value.py`, `search_multi_value.py` | float64 screening only | heuristic |
| `bolshev.py` | Reproduces Bimpeh's Table 5.1 and demonstrates that his coverage bound (5.16) is not a coverage bound (`audit/BIMPEH-GAP.md`) | proof-essential |
| `bimpeh_continuous_check.py` | MC corroboration of the hand counterexamples to (5.16) with continuous F | corroboration only |

## Commands

```sh
cd computations/python

# Lemma: single-value supports cannot under-cover; Poisson factors dominate
python3 two_point_lemma.py --alpha 0.05 --n-max 40

# Known finite-sample violation at low confidence (machinery true-positive)
python3 search_two_value.py --alpha 0.7 --n 50 --out /tmp/c.json
python3 certify.py /tmp/c.json    # expect CONFIRMED lines

# The conjecture at 95% over two-value supports
python3 search_two_value.py --alpha 0.05 --n 2 30 --range --out /tmp/c95.json
python3 certify.py /tmp/c95.json  # expect "nothing to certify" if no dips

# Richer supports
python3 search_multi_value.py --alpha 0.05 --m 3 --n 10 20 --out /tmp/c3.json
```

Certified run logs are committed under `computations/certificates/`.

## Certificate semantics

A `CONFIRMED` line from `certify.py` states: for the exact rational taint
distribution printed, the exact rational coverage is below the exact
nominal level, and the smallest \(|SB - \theta|\) over all count vectors
exceeds \((n+1)\) times the largest factor-bracket width, so no comparison
was decided by numerical error. Trusted inputs: Python, mpmath's arithmetic
on `mpf`, and the hardware.

## What screening output does NOT establish

A screening minimum equal to nominal (slack \(\sim 10^{-13}\)) does not
prove the infimum equals the nominal level; it is evidence subject to grid
and optimizer limitations. Negative results (no candidate found) bound only
the families and sample sizes actually searched.
