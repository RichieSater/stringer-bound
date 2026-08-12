# Findings for practitioners: stress-testing the Stringer bound

**Project links:** [GitHub repository](https://github.com/RichieSater/stringer-bound) · [archived v1.0.0 on Zenodo (doi:10.5281/zenodo.21850820)](https://doi.org/10.5281/zenodo.21850820)

*Plain-language summary of what this repository establishes about the
monetary-unit-sampling evaluation method used across the audit profession.
Exact claims are backed by rational certificates or written proofs; search
results are labeled as numerical evidence. See the README for the technical
statements.*

## 1. Exact conventional-level guarantees now reach n = 5

For samples of size **n = 3, n = 4, and n = 5**, the binomial-factor
Stringer bound is proved conservative for every taint distribution on [0,1]
at **90%, 95%, and 99% nominal confidence**. The proofs show, sample by
sample, that Stringer
is no smaller than a valid distribution-free upper confidence limit for a
bounded mean. Uniform-simplex cap inequalities are verified with exact
rational Bernstein certificates; these are computer-assisted proofs, not
numerical searches
([n = 3 proof](supporting-materials/theory/N3-CONVENTIONAL.md),
[n = 4 proof](supporting-materials/theory/N4-CONVENTIONAL.md),
and [n = 5 proof](supporting-materials/theory/N5-CONVENTIONAL.md)).

The same guarantees hold for the Poisson-factor version used in professional
MUS tables because those factors pointwise dominate the binomial factors at
all three confidence levels. Together with the all-level theorem at n = 2,
this makes **n = 6 the first unresolved sample size at 90%, 95%, and 99%**.

This is a meaningful finite-sample advance, but n = 5 is far smaller than an
ordinary audit sample. It does not yet justify a change in audit guidance.

A separate corrected simultaneous-band theorem proves ordinary
**Poisson-factor** Stringer coverage farther: every n <= 8 at 90%, every
n <= 11 at 95%, and every n <= 20 at 99%. These are exact finite ranges,
not search results. The sufficient band event falls short at the next n in
each row, which limits that proof route rather than demonstrating Stringer
undercoverage
([proof and certificate](supporting-materials/theory/POISSON-SIMULTANEOUS-BAND.md)).

## 2. A finite-sample-valid reporting safeguard is available for every n

An audit methodology group can evaluate a reporting rule that does not wait
for the open general-n conjecture: before observing the sample, define the
reported upper bound as

```text
max(ordinary Stringer, validated Gaffke bounded-mean limit).
```

Under the paper's independent [0,1]-taint model, this complete rule has at
least nominal coverage for every sample size and either factor convention.
It can only increase, never decrease, a Stringer result. The Gaffke endpoint
is computed with exact rational tail-sign checks rather than an uncertified
floating-point spline calculation
([technical proof](supporting-materials/theory/GAFFKE-SAFEGUARD.md),
[practice note](supporting-materials/audit/PRACTICE-SAFEGUARD.md)).

The safeguard is provably identical to ordinary Stringer on every sample at
n = 3, 4, and 5 at 90%, 95%, and 99%. A new one-cap theorem gives a much
larger sample-checkable identity region: at every nominal confidence level
of at least 75%, and **every sample size**,
there is **zero uplift whenever binomial Stringer is at least the largest
observed taint**. The proof combines
a dimension-free cap lemma with analytic binomial-tail comparisons. An older
set of 59,700 exact Clopper--Pearson vertex checks through n = 200 remains as
an independent regression
([proof and checks](supporting-materials/theory/ONE-CAP-COMPARISON.md)).

This is a candidate statistical control for independent human proof review,
professional-standards analysis, model-risk validation, and software
testing. It is not represented as adopted guidance. The zero-uplift result
describes the output of the pre-specified safeguard; it is not a
conditional-coverage claim for ordinary Stringer.

## 3. The general finite-sample guarantee remains open

The Stringer bound is commonly expected to cover the population mean at
least as often as its nominal confidence level, but the claim has not been
proved for arbitrary sample size and arbitrary taint distributions. At 90%,
95%, and 99%, the first open case is now n = 6; coverage for the larger sample
sizes used in practice remains unresolved.

A previous finite-sample certification (Bimpeh 2008: reliability through
n = 11 at 95%) contains an index shift in its confidence-band argument. The
cited argument therefore does not establish that range
([details](supporting-materials/audit/BIMPEH-GAP.md)). Proven cases now
include n = 1; n = 2 at every confidence level; n = 3 through n = 5 at the
three levels above; populations supported on {0,1}; and populations with one
nonzero taint value
([proof](supporting-materials/computations/python/two_point_lemma.py)).

## 4. No violation was found in the reported 90% and 95% searches

Grid-plus-optimizer searches over populations with two or three nonzero
taint values (n up to 100 in the reported standard-confidence searches)
found no coverage below nominal. The smallest values found were close to
nominal in the examined parameter ranges, with near-minimizers occurring as
the largest taint approaches 1. These are numerical search results, not an
exhaustive proof. The exact n = 2 through n = 5 results above do not
depend on these searches.

## 5. Practical-level Poisson factor domination is proved for every n

The factors in professional MUS tables are Poisson-based. Whenever nominal
confidence exceeds \(1-e^{-1}\approx63.2\%\), they dominate the binomial
factors coordinatewise for **every sample size**. A summation-by-parts
argument then proves that the complete Poisson-factor Stringer bound is at
least the binomial-factor bound on every observed sample
([proof](supporting-materials/theory/POISSON-DOMINATION.md)). This includes
90%, 95%, and 99% confidence. It transfers every established binomial
coverage guarantee in that range to the Poisson version, but it does not
alone settle general-n coverage because the binomial conjecture remains open
there.

The smallest coverage found in searches over populations with two nonzero
taint values was:

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

## 6. At selected low assurance levels the bound provably fails

For several low-confidence parameter choices, the Stringer bound's coverage
falls **below** its stated level. These are not conventional
financial-statement-audit confidence levels. The examples are certified in
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
(Pap–van Zuijlen 1996), but they do not locate a finite-sample threshold. At
exactly 50% confidence, no failure was found through n = 300 in the families
searched. The certified examples show that low-confidence use requires care;
they do not show failure at every confidence level below 50% or for every
sample size.

## 7. The n = 2 theorem is sharp

The Stringer bound is conservative at **n = 2, for every possible taint
population and every confidence level**
([N2-PROOF.md](supporting-materials/theory/N2-PROOF.md)). Its worst-case
coverage infimum equals the nominal level, although no population attains the
infimum. Thus even where coverage is proved, the binomial bound has no
uniform positive safety margin above nominal.

The proof uses a change of variables that turns the bound into a weighted
order-statistic inequality, followed by a rectangle-and-budget argument.
The n = 3 through n = 5 proofs use a different, geometric route. The one-cap
theorem now controls one complete geometric region at every sample size for
nominal confidence of at least 75%.
Determining whether the remaining simplex-cap regions extend to n = 6, and whether they admit a
dimension-free argument, are the repository's next mathematical targets.
