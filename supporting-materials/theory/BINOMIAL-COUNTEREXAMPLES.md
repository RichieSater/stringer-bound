# Exact finite-support Stringer counterexamples

These results concern independent observations with a common distribution
on `[0,1]` and binomial Clopper--Pearson Stringer factors. They do not
establish failure at conventional confidence levels or under a different
sampling design.

## One nonzero taint value

If the distribution is supported on `{0,v}`, where `0<v<=1`, write
`q=P(T=v)` and let `K` be the number of nonzero observations. Then
`K` is binomial with parameters `n,q`, and

\[
 \operatorname{SB}=(1-v)p_n(0)+v p_n(K)\ge v p_n(K),
 \qquad \mu=qv.
\]

The event `p_n(K)>=q` has probability at least `1-alpha` by
Clopper--Pearson coverage. Hence such a distribution cannot undercover.
The independent implementation is
[`two_point_lemma.py`](../computations/python/two_point_lemma.py).

## Certified counterexamples

The distributional parameters are exact dyadic rationals. Full multinomial
enumeration computes each coverage as an exact rational number. Factor
intervals have endpoints on the grid of denominator `2^80`, and integer
binomial-CDF signs certify their endpoints. The displayed coverage values
are rounded summaries, not the exact rational values used to decide failure.

| nominal confidence | sample size | smallest certified coverage | certified examples |
|---:|---:|---:|---:|
| 30% | 50 | 29.815% | 3 |
| 32% | 100 | 30.956% | 9 |
| 35% | 200 | 33.729% | 1 |
| 35% | 400 | 32.630% | 1 |
| 37% | 200 | 36.079% | 6 |
| 37% | 400 | 36.777% | 13 |

The 33 certificates concern three-point distributions with two distinct
nonzero taint values. The table is generated from
[`certificate-summary.json`](../computations/certificates/certificate-summary.json)
by [`summarize_certificates.py`](../computations/python/summarize_certificates.py).
The source certificate files are named `certified-alpha*-n*.json` in the
same certificate directory. These examples establish existence at the six
listed parameter pairs, not at every confidence level below 50%.

For example, the third record of
[`certified-alpha0.7-n50.json`](../computations/certificates/certified-alpha0.7-n50.json)
uses `alpha=0.7`, `n=50`, and

\[
 \begin{aligned}
 v_1&=1,&q_1&=8984394653678481/2^{55},\\
 v_2&=8384891858230937/2^{54},&q_2&=5642023711000999/2^{53},\\
 &&q_0&=4476307521281491/2^{55}.
 \end{aligned}
\]

Its coverage is approximately `0.2988790828`, strictly below `0.30`;
the certificate contains the exact rational value and the strictly positive
comparison margins. It is a representative example, not the minimum
coverage in the first table row.

## Numerical searches

Reported grid-plus-optimizer searches at 50%, 90%, and 95% confidence found
no violation. The two-nonzero-value searches covered `n<=100` at 90% and
95% and `n<=300` at 50%; a three-nonzero-value search covered `n<=25` at
95%. The smallest coverage found agreed with nominal to within `1e-10`.
These are numerical search results, not exhaustive optimization or a
coverage theorem. Relevant records are
[`search-alpha0.05-n35-100.log`](../computations/certificates/search-alpha0.05-n35-100.log),
[`search-practical-poisson-90-95.log`](../computations/certificates/search-practical-poisson-90-95.log),
and [`search-alpha0.5-large-n.log`](../computations/certificates/search-alpha0.5-large-n.log).

The corresponding Poisson-factor searches found coverage above nominal.
At 95%, the reported minima were approximately `0.9716` at `n=10`,
`0.9573` at `n=50`, and `0.9549` at `n=100`. This finite search does not
establish a global minimum or a general Poisson coverage guarantee.

## Reproduction

Run `make certificate-summary-check`. The table-display regression is in
`test_exact_certification.py`; the exact enumeration and factor-enclosure
checks are in `coverage.py`, `certify.py`, and `stringer.py`.
