# Supporting materials: verification guide

**Canonical repository:** [github.com/RichieSater/stringer-bound](https://github.com/RichieSater/stringer-bound) · **Archived release:** [doi:10.5281/zenodo.21850820](https://doi.org/10.5281/zenodo.21850820)

The canonical reproducibility environment requires Python 3.12 and pins
every dependency in the root
`pyproject.toml` and `uv.lock`. Every
certified binomial counterexample is decided by rational arithmetic;
float64 appears in screening searches and in separately labeled numerical
Poisson checks. The all-sample-size Poisson-versus-binomial factor theorem
at nominal confidence above \(1-e^{-1}\) is a written analytic result, not a
numerical certificate; see `theory/POISSON-DOMINATION.md`. The \(n=3\) and
\(n=4\) guarantees at 90%, 95%, and 99% are exact computer-assisted
theorems: their symbolic formulas and rational sign certificates are
regenerated from source scripts; see `theory/N3-CONVENTIONAL.md` and
`theory/N4-CONVENTIONAL.md`.

## Layers

| Layer | Role | Trust status |
|---|---|---|
| `stringer.py` | Numerical factors for searches; for certification, dyadic binomial-factor intervals whose endpoint CDF signs are evaluated exactly with integers | proof-essential |
| `coverage.py` / `coverage_exact` | Exact multinomial weights and rational propagation of factor intervals through every Stringer comparison | proof-essential |
| `two_point_lemma.py` | Written proof that single-value supports cannot under-cover, plus machine checks | proof-essential |
| `theory/N3-CONVENTIONAL.md` | Complete reduction of the \(n=3\) theorem to simplex-cap inequalities | proof-essential |
| `derive_n3_bernstein_formulas.py` | Symbolically derives and regenerates the 42 middle-region Bernstein formulas from the cap volume | proof-essential |
| `n3_gaffke_certificate.py` | Exact rational interval proof of every coefficient sign at 90%, 95%, and 99% | proof-essential |
| `theory/N4-CONVENTIONAL.md` | Complete reduction of the \(n=4\) theorem to four simplex-cap regions | proof-essential |
| `derive_n4_bernstein_structure.py` | Exactly derives the residual polynomials, boundary factorizations, tetrahedral substitutions, and structural zeros for \(n=4\) | proof-essential |
| `n4_gaffke_certificate.py` | Exact rational interval proof of every nonzero \(n=4\) Bernstein coefficient sign at 90%, 95%, and 99% | proof-essential |
| `theory/POISSON-DOMINATION.md` | Written proof that practical-level Poisson factors dominate binomial factors for every sample size | proof-essential |
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

## What screening output does NOT establish

A screening minimum equal to nominal (slack \(\sim 10^{-13}\)) does not
prove the infimum equals the nominal level; it is evidence subject to grid
and optimizer limitations. Negative results (no candidate found) bound only
the families and sample sizes actually searched.
