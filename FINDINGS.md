# Findings for practitioners: stress-testing the Stringer bound

**Project links:** [GitHub repository](https://github.com/RichieSater/stringer-bound) · [archived v1.0.0 on Zenodo (doi:10.5281/zenodo.21850820)](https://doi.org/10.5281/zenodo.21850820)

*Plain-language summary of what this repository establishes about the
monetary-unit-sampling evaluation method used across the audit profession.
Exact claims are backed by rational certificates or written proofs; search
results are labeled as numerical evidence. See the README for the technical
statements.*

## 1. The general finite-sample guarantee remains open

The Stringer bound is commonly expected to cover the population mean at
least as often as its nominal confidence level, but the claim has not been
proved for general sample size and general taint distributions. A previous
finite-sample certification (Bimpeh 2008: reliability up to n = 11 at 95%)
contains an index shift in its confidence-band argument; the cited argument
therefore does not establish the claimed result
([details](supporting-materials/audit/BIMPEH-GAP.md)). Proven cases include
n = 1, n = 2, populations supported on {0,1}, and populations with one
nonzero taint value
([proof](supporting-materials/computations/python/two_point_lemma.py)).

## 2. No violation was found in the reported 90% and 95% searches

Grid-plus-optimizer searches over populations with two or three nonzero
taint values (n up to 100 in the reported standard-confidence searches)
found no coverage below nominal. The smallest values found were close to
nominal in the examined parameter ranges, with near-minimizers
occurring as the largest taint approaches 1. These are numerical search
results, not an exhaustive proof. At n = 2, the written theorem does prove
that the binomial bound has no uniform safety margin above nominal.

## 3. Poisson-factor searches show a positive numerical cushion

The factors in professional MUS tables are Poisson-based. At the 90% and
95% confidence levels and sample sizes reported here, high-precision checks
find that they dominate the binomial factors
([log](supporting-materials/computations/certificates/poisson-domination-standard-levels.log)).
Smallest coverage found in
searches over populations with two nonzero taint values:

| n | smallest found at nominal 95% | margin | smallest found at nominal 90% | margin |
|---|---|---|---|---|
| 10 | 0.9716 | +2.16% | 0.9270 | +2.70% |
| 25 | 0.9589 | +0.89% | 0.9182 | +1.82% |
| 50 | 0.9573 | +0.73% | 0.9090 | +0.90% |
| 75 | 0.9553 | +0.53% | 0.9145 | +1.45% |
| 100 | 0.9549 | +0.49% | 0.9093 | +0.93% |

(The smallest coverage found for the binomial-factor version of the same
searches agrees numerically with nominal at every listed n.) The Poisson
values show a positive numerical cushion in these searches; the table does
not establish a global worst case.

## 4. At selected low assurance levels the bound provably fails

For several low-confidence parameter choices, including levels relevant to
low-risk strata and limited-assurance procedures, the Stringer bound's
coverage falls **below** its stated level. These examples are certified in
exact rational arithmetic (no simulation error):

| stated confidence | sample size | smallest certified coverage | certificate |
|---|---|---|---|
| 30% | 50 | 29.815% | [certified-alpha0.7-n50.json](supporting-materials/computations/certificates/certified-alpha0.7-n50.json) |
| 32% | 100 | 30.956% | [certified-alpha0.68-n100.json](supporting-materials/computations/certificates/certified-alpha0.68-n100.json) |
| 35% | 200 | 33.729% | [certified-alpha0.65-n200.json](supporting-materials/computations/certificates/certified-alpha0.65-n200.json) |
| 35% | 400 | 32.630% | [certified-alpha0.65-n400.json](supporting-materials/computations/certificates/certified-alpha0.65-n400.json) |
| 37% | 200 | 36.079% | [certified-alpha0.63-n200.json](supporting-materials/computations/certificates/certified-alpha0.63-n200.json) |
| 37% | 400 | 36.777% | [certified-alpha0.63-n400.json](supporting-materials/computations/certificates/certified-alpha0.63-n400.json) |

The rows are generated from the version-2 certificate records by
[`summarize_certificates.py`](supporting-materials/computations/python/summarize_certificates.py);
the machine-readable output is
[`certificate-summary.json`](supporting-materials/computations/certificates/certificate-summary.json).

These examples are consistent with the proven asymptotic threshold at 50%
(Pap–van Zuijlen 1996), but they do not locate a finite-sample threshold.
At exactly 50% confidence, no failure was found through n = 300 in the
families searched. The certified examples show that low-confidence use
requires care; they do not show failure at every confidence level below
50% or for every sample size.

## 5. Conservatism is proved at n = 2

This repository now contains a complete proof that the Stringer bound
is conservative at **n = 2, for every possible error population and
every confidence level**
([N2-PROOF.md](supporting-materials/theory/N2-PROOF.md)). The theorem shows
that the certified low-confidence failures cannot occur at n = 2. Its
machinery—a change of variables that turns the bound into a weighted
order-statistic inequality, followed by a rectangle-and-budget
argument—suggests one possible route to larger n, but whether that approach
extends remains open. The repository's exact machinery is designed so that
future finite-support counterexamples can be checked without relying on
simulation.
