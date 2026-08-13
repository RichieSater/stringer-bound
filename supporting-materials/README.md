# Supporting materials: verification guide

**Canonical repository:** [github.com/RichieSater/stringer-bound](https://github.com/RichieSater/stringer-bound) · **Archived v1.0.0:** [doi:10.5281/zenodo.21850820](https://doi.org/10.5281/zenodo.21850820)

The current development revision postdates archived v1.0.0; GitHub is the
canonical source until the next archival release.

The canonical reproducibility environment requires Python 3.12 and pins
every Python dependency in the root
`pyproject.toml` and `uv.lock`. Every
regeneration of the `n=6` structure additionally requires Singular for exact
generic rational-function face-ideal reductions.
The certificate was developed with Singular 4.4.1; CI fixes Ubuntu 24.04,
installs its packaged Singular build, and byte-compares the regenerated
artifact. Every certified binomial counterexample is decided by rational
arithmetic;
float64 appears in screening searches and in separately labeled numerical
Poisson checks. The all-sample-size Poisson-versus-binomial factor theorem
at nominal confidence above \(1-e^{-1}\) is a written analytic result, not a
numerical certificate; see `theory/POISSON-DOMINATION.md`. The guarantees
from \(n=3\) through \(n=6\) at 90%, 95%, and 99% are exact computer-assisted
theorems: their symbolic formulas and rational sign certificates are
regenerated from source scripts; see `theory/N3-CONVENTIONAL.md`,
`theory/N4-CONVENTIONAL.md`, `theory/N5-CONVENTIONAL.md`, and
`theory/N6-CONVENTIONAL.md`.
For the Poisson factors used in practice, a separate exact corrected-band
certificate proves coverage for every \(n\le8\) at 90%, every \(n\le11\) at
95%, and every \(n\le20\) at 99%; see
`theory/POISSON-SIMULTANEOUS-BAND.md`. The same band lemma yields an
all-sample-size scalar-calibrated Poisson rule, with exact representative
multiplier certificates in `theory/POISSON-BAND-CALIBRATION.md`. Its anchored
variant fixes the ordinary zero-taint factor and scales only error increments.
For arbitrary sample size, `theory/GAFFKE-SAFEGUARD.md` defines a
distribution-free valid reporting rule obtained by taking the maximum of
Stringer and the valid Gaffke limit. Its implementation brackets the Gaffke
quantile with exact rational tail-sign checks.
Within that rule, `theory/ONE-CAP-COMPARISON.md` proves a directly checkable
zero-uplift region for every sample size whenever nominal confidence is at
least \(75\%\): the safeguard is ordinary Stringer whenever
the binomial Stringer value is at least the largest observed taint. This is a sample-wise
comparison, not a general coverage theorem for ordinary Stringer.

## Layers

| Layer | Role | Trust status |
|---|---|---|
| `stringer.py` | Numerical factors for searches; exact dyadic binomial-factor intervals with integer CDF signs; exact dyadic Poisson-limit intervals with rational exponential enclosures | proof-essential |
| `coverage.py` / `coverage_exact` | Exact multinomial weights and rational propagation of factor intervals through every Stringer comparison | proof-essential |
| `two_point_lemma.py` | Written proof that single-value supports cannot under-cover, plus machine checks | proof-essential |
| `theory/N3-CONVENTIONAL.md` | Complete reduction of the \(n=3\) theorem to simplex-cap inequalities | proof-essential |
| `derive_n3_bernstein_formulas.py` | Symbolically derives and regenerates the 42 middle-region Bernstein formulas from the cap volume | proof-essential |
| `n3_gaffke_certificate.py` | Exact rational interval proof of every coefficient sign at 90%, 95%, and 99% | proof-essential |
| `theory/N4-CONVENTIONAL.md` | Complete reduction of the \(n=4\) theorem to four simplex-cap regions | proof-essential |
| `derive_n4_bernstein_structure.py` | Exactly derives the residual polynomials, boundary factorizations, tetrahedral substitutions, and structural zeros for \(n=4\) | proof-essential |
| `n4_gaffke_certificate.py` | Exact rational interval proof of every nonzero \(n=4\) Bernstein coefficient sign at 90%, 95%, and 99% | proof-essential |
| `theory/N5-CONVENTIONAL.md` | Complete reduction of the \(n=5\) theorem to five four-dimensional simplex-cap regions | proof-essential |
| `derive_n5_bernstein_structure.py` | Exactly derives the residual factorizations, four-simplex substitutions, and face-ideal structural-zero proofs for \(n=5\) | proof-essential |
| `n5_gaffke_certificate.py` | Integer-directed dyadic interval proof of every nonzero \(n=5\) Bernstein coefficient sign at 90%, 95%, and 99% | proof-essential |
| `theory/N6-CONVENTIONAL.md` | Complete reduction of the \(n=6\) theorems at 90%, 95%, and 99% to six five-dimensional simplex-cap regions | proof-essential |
| `derive_n6_bernstein_structure.py` | Exactly derives the residual factorizations and five-simplex structure and proves generic face-ideal vanishing orders with Singular | proof-essential |
| `n6_gaffke_certificate.py` | Integer-directed dyadic proof of face specialization, the triangulation chain, and every nonstructural \(n=6\) Bernstein sign at all three levels | proof-essential |
| `theory/POISSON-DOMINATION.md` | Written proof that practical-level Poisson factors dominate binomial factors for every sample size | proof-essential |
| `theory/POISSON-SIMULTANEOUS-BAND.md` | Written randomized-survival-band proof of the direct Poisson coverage ranges | proof-essential |
| `poisson_band_certificate.py` | Exact rational Poisson-limit signs and Bolshev boundary-crossing probabilities at 90%, 95%, and 99% | proof-essential |
| `theory/POISSON-BAND-CALIBRATION.md` | Written all-sample-size theorems for full-scale and zero-taint-preserving Poisson calibrations | proof-essential |
| `poisson_band_calibration.py` | Exact adjacent-dyadic brackets for representative multipliers on both calibration paths | proof-essential |
| `theory/GAFFKE-SAFEGUARD.md` | Written all-sample-size coverage argument and exact divided-difference computation of the Gaffke floor | proof-essential |
| `theory/ONE-CAP-COMPARISON.md` | Dimension-free cap lemma and analytic all-sample-size Stringer specialization at confidence at least \(75\%\) | proof-essential |
| `one_cap_all_n_check.py` | Symbolic identities, rational constant checks, and finite numerical regression for the analytic proof | proof support |
| `one_cap_certificate.py` | Independent finite regression with integer-checked CP brackets and 59,700 exact nonterminal vertices through \(n=200\) | corroboration |
| `theory/ORDERED-SIMPLEX-CAP.md` | Exact tight-vertex reduction and adjacent-transfer identity for the still-open all-`n` comparison | research roadmap |
| `theory/ALL-N-POISSON-PROGRAM.md` | Exact weighted-exponential reductions, sharp convexity obstruction, and explicitly conjectural all-`n` route | research roadmap |
| `theory/DIRICHLET-POISSONIZATION.md` | Constrained divided-difference reduction, zero-knot boundary reduction, proof for every two-level profile, and exact obstruction to generic `s`-concave localization | research roadmap |
| `gaffke.py` | Exact-rational Dirichlet tail signs, certified dyadic Gaffke bracket, and safeguarded Stringer report | proof-essential for the implemented floor |
| `all_n_poisson_reductions.py` | Symbolic checks of the ordered-weight, two-exponential, and kernel identities; labels both remaining inequalities open | research support |
| `dirichlet_poissonization.py` | Symbolic checks for the zero-knot and all-two-level reductions plus the exact `s`-affine obstruction; does not claim the general inequality | research support |
| `audit/PRACTICE-SAFEGUARD.md` | Methodology-facing workflow, worked example, record-retention fields, and explicit scope boundaries | implementation guidance |
| `audit/HUMAN-REVIEW-PACKET.md` | Scoped independent-review questions, reproduction record, and sign-off template | review protocol |
| `certify.py` | The only source of claims: exact coverage, exact nominal comparison, margin certificate | proof-essential |
| `search_two_value.py`, `search_multi_value.py` | float64 screening only | heuristic |
| `bolshev.py` | Reproduces Bimpeh's Table 5.1 and demonstrates that his coverage bound (5.16) is not a coverage bound (`audit/BIMPEH-GAP.md`) | proof-essential |
| `bimpeh_continuous_check.py` | MC corroboration of the hand counterexamples to (5.16) with continuous F | corroboration only |

## Commands

From the repository root, the complete core check is:

```sh
make reproduce
```

For individual computations:

```sh
# Exact n=3 theorem: derive formulas, then regenerate all sign certificates
make n3-formula-check
make n3-certificate-check

# Exact n=4 theorem: derive tetrahedral structure, then certify all signs
make n4-structure-check
make n4-certificate-check

# Exact n=5 theorem: derive four-simplex/face-ideal structure, then certify
make n5-structure-check
make n5-certificate-check

# Exact n=6 theorems: generic face ideals, then directed intervals
make n6-structure-check
make n6-certificate-check

# Equivalent CI-schedulable split of the long generic structure check
make n6-structure-data-check \
  n6-face-standard-0-check n6-face-standard-1-check \
  n6-face-standard-2-check n6-face-standard-3-check \
  n6-face-standard-4-check n6-face-standard-5-check \
  n6-face-c6-0-check n6-face-c6-1-check n6-face-c6-2-check \
  n6-face-d6-0-check n6-face-d6-1-check n6-face-d6-2-check
make n6-certificate-001-check n6-certificate-005-check \
  n6-certificate-010-check

# Direct Poisson theorem: exact limits and boundary-crossing probabilities
make poisson-band-certificate-check

# All-n scalar Poisson rule: exact representative multiplier brackets
make poisson-band-calibration-check

# Exact calibrated report for a specified zero-heavy audit sample
uv run --frozen python \
  supporting-materials/computations/python/poisson_band_calibration.py \
  --n 25 --alpha 0.05 --taints 1,0.4,0.1 \
  --out /tmp/calibrated-poisson.json
# The JSON returns two paths; pre-specify one rather than taking their minimum.

# Analytic all-n one-cap theorem at confidence >=75%
make one-cap-all-n-check

# Independent exact finite regression through n=200
make one-cap-certificate-check

# Exact algebra behind the explicitly conjectural all-n Poisson program
make all-n-reduction-check

# Zero-knot, all-two-level, and obstruction checks for Dirichlet poissonization
make dirichlet-poissonization-check

cd supporting-materials/computations/python

# Lemma: single-value supports; high-precision Poisson comparison
uv run --frozen python two_point_lemma.py --alpha 0.05 --n-max 40

# Known finite-sample violation at low confidence (machinery true-positive)
uv run --frozen python search_two_value.py --alpha 0.7 --n 50 --out /tmp/c.json
uv run --frozen python certify.py /tmp/c.json    # expect CONFIRMED lines

# The conjecture at 95% over two-value supports
uv run --frozen python search_two_value.py --alpha 0.05 --n 2 30 --range --out /tmp/c95.json
uv run --frozen python certify.py /tmp/c95.json  # expect "nothing to certify" if no dips

# Richer supports
uv run --frozen python search_multi_value.py --alpha 0.05 --m 3 --n 10 20 --out /tmp/c3.json

# Exact-sign, interval-propagation, and table-generation regression tests
uv run --frozen python -m unittest discover -s ../tests -v

# All-n valid reporting safeguard; zeros may be omitted, but n is the full sample
uv run --frozen python gaffke.py \
  --n 100 --alpha 0.05 --method poisson --taints 1,0.4,0.1
```

Certified run logs are committed under `computations/certificates/`.
`summarize_certificates.py` regenerates the rows used in the manuscript:

```sh
uv run --frozen python summarize_certificates.py \
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

The committed Poisson-factor comparison log is a high-precision regression
check, not the proof of domination. The all-\(n\) theorem for
\(\alpha<e^{-1}\) is proved analytically in
`theory/POISSON-DOMINATION.md`. Outside that confidence range, no general
domination claim is made; the log includes a case where domination fails.

The direct Poisson coverage certificate is different: its limit brackets
and boundary-crossing probabilities are formal exact-arithmetic evidence.
At dyadic endpoints, `stringer.py` encloses the exponential between rational
alternating-series partial sums after power-of-two range reduction. Then
`poisson_band_certificate.py` evaluates Bolshev's recursion with exact
`Fraction` arithmetic. Numerical root finding proposes endpoints only.

The scalar-calibration certificate uses the same exact Poisson-limit
enclosures. Floating point proposes one \(2^{-28}\) multiplier cell on each
path. For full scaling, opposite common endpoints enclose the band
probability. For the anchored path, later boundaries increase with their own
Poisson limit but decrease with the zero-count limit, so the certificate uses
the corresponding mixed endpoint substitutions. Exact signs bracket both
path-minimal multipliers. Each upper endpoint is a rigorous valid choice for
its modified rule, with either untruncated or factorwise-capped Poisson
factors. The written theorem also permits capping the calibrated factors
themselves at one, which preserves the band event and gives a pointwise
no-larger report. A separate analytic corollary proves that the full-scale
multiplier two is valid uniformly in sample size whenever nominal confidence
is at least \(2/e\approx73.6\%\). Retaining 100 exact Poisson crossing terms
reduces the certified uniform choices to 1.53, 1.44, and 1.33 at 90%, 95%,
and 99%, respectively. These results do not prove undercoverage of ordinary
Stringer.

## The all-sample-size reporting safeguard

The mathematical procedure is
`max(Stringer, Gaffke)`. Its coverage follows from the independently valid
Gaffke component and therefore does not assume general-`n` Stringer
conservatism. At `n=3,4,5` and 90%, 95%, and 99%, the exact pointwise
comparison theorems show that this maximum is ordinary Stringer on every
sample, for both binomial and Poisson factors. The same statement holds at
`n=6` at 90%, 95%, and 99%.

For decimal taints, `gaffke.py` interprets the inputs as exact rationals.
It uses a B-spline quantile only as a root-location proposal, evaluates the
Dirichlet-average tail at the surrounding dyadic endpoints by exact
confluent divided differences, and expands the bracket unless both tail
signs are proved. The reported Gaffke value is the upper endpoint. Thus the
validity of the implemented floor does not depend on floating-point signs;
see `theory/GAFFKE-SAFEGUARD.md` for the formulas and scope conditions.

At every nominal confidence level of at least
\(75\%\) and every sample size, the additional one-cap
theorem proves that no Gaffke uplift can occur whenever
binomial Stringer is at least the largest sample taint. The proof reduces the
Dirichlet tail to `max_r (C_r/(C_r+p_0))^r` and proves every
Clopper--Pearson vertex inequality analytically. The older exact calculation
through `n=200` remains an independent finite regression. See
`theory/ONE-CAP-COMPARISON.md`. The condition identifies a region where the
pre-specified safeguard returns Stringer; it must not be recast as a
post-selection or conditional-coverage claim.

## The n = 3 theorem certificate

The \(n=3\) proof is separate from both screening and finite-support
counterexample certification. It first reduces pointwise domination of the
valid Gaffke upper limit to three uniform-simplex cap inequalities. Two are
polynomial nonnegativity problems on triangles and one reduces to a
one-dimensional boundary polynomial.

`derive_n3_bernstein_formulas.py` differentiates the middle cap formula,
performs both triangle substitutions, and derives all 42 degree-five
Bernstein coefficients with SymPy. `n3_gaffke_certificate.py` encloses the
Clopper--Pearson factors on a \(2^{-120}\) dyadic grid, verifies every
binomial-CDF endpoint sign with integer arithmetic, and propagates the
enclosures through every coefficient using exact `Fraction` interval
arithmetic. Floating point is used only for readable decimal summaries.
The two committed JSON files are regenerated byte-for-byte by
`make n3-formula-check n3-certificate-check`.

## The n = 4 theorem certificate

The \(n=4\) proof applies the same pointwise Gaffke comparison to the
four-dimensional uniform simplex. After an affine normalization of the
sample, its ordered-knot domain has four cap regions. The lowest region is
closed analytically by weighted AM--GM. In the other three regions,
nonnegative boundary factors are removed exactly and the remaining
degree-four or degree-six residual is checked on a fixed decomposition into
seven tetrahedra.

`derive_n4_bernstein_structure.py` starts from the cap formulas, verifies
the polynomial factorizations by exact division, performs every tetrahedral
substitution, and identifies structural-zero Bernstein coefficients over
the rational-function field `QQ(b,c,d,e)`. `n4_gaffke_certificate.py`
encloses the four Clopper--Pearson factors on a \(2^{-120}\) dyadic grid,
checks their endpoint signs with integers, and evaluates every remaining
coefficient by exact `Fraction` interval arithmetic. The two committed JSON
files are regenerated byte-for-byte by
`make n4-structure-check n4-certificate-check`.

## The n = 5 theorem certificate

The \(n=5\) proof extends the same pointwise Gaffke comparison to a
four-dimensional ordered-knot domain with five cap regions. Exact polynomial
division removes nonnegative boundary factors. The remaining degree-five,
degree-eight, and degree-nine residuals are checked on a fixed decomposition
into thirteen four-simplices.

`derive_n5_bernstein_structure.py` derives the cap residuals and vertex
formulas symbolically. Rather than classifying near-zero coefficients by a
tolerance, it proves every structural-zero pattern through exact membership
in a power of the relevant affine face ideal over `QQ(b,c,d,e,f)`.
`n5_gaffke_certificate.py` brackets the five Clopper--Pearson factors on a
\(2^{-240}\) grid with integer-checked endpoint signs, then propagates closed
dyadic intervals with integer-directed outward rounding at 256 fractional
bits through every vertex, polynomial, and Bernstein calculation. Every
remaining coefficient has a positive lower endpoint. The sign certificate
records the SHA-256 digest of the structure artifact. Both JSON files are
regenerated byte-for-byte by
`make n5-structure-check n5-certificate-check`.

## What screening output does NOT establish

A screening minimum equal to nominal (slack \(\sim 10^{-13}\)) does not
prove the infimum equals the nominal level; it is evidence subject to grid
and optimizer limitations. Negative results (no candidate found) bound only
the families and sample sizes actually searched.
