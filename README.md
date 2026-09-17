# Ordered simplex caps and Stringer bounds

[![Verify paper and certificates](https://github.com/RichieSater/stringer-bound/actions/workflows/verify.yml/badge.svg)](https://github.com/RichieSater/stringer-bound/actions/workflows/verify.yml)

## Principal manuscript

**[Ordered simplex caps and the three-observation Stringer bound](supporting-materials/paper/stringer.pdf)**
([LaTeX source](supporting-materials/paper/stringer.tex)).

For a uniform tetrahedron and any nondecreasing barycentric probability
vector, the largest cap over ordered normals occurs at a step normal.
A section-centroid argument reduces the result to one exact rational
polynomial inequality and a separate repeated-upper-knot identity.

Consequently, for three independent observations from any distribution on
\([0,1]\), the binomial Stringer bound has coverage at least the nominal
level whenever

\[
 0<\alpha\le\left(\frac{19+\sqrt{21}}{34}\right)^3
 =0.3336852118672717\ldots .
\]

This endpoint is sharp for pointwise domination of Gaffke's valid upper mean
bound, not asserted to be a coverage threshold. A two-knot perturbation gives
an obstruction to pointwise domination in every sample size; its critical
tail levels converge to \(0.28466813704083846\ldots\).

The article contains the geometric proof, its statistical application, and
the obstruction. It does not combine the separate Poisson, finite-frame,
or larger fixed-dimension calculations into the article.

## Mathematical scope

- **Proved with exact rational computation:** the four-coordinate cap
  theorem and the three-observation comparison above.
- **Conjectural:** the global five- and six-coordinate monotone cap
  extensions and the dependent full confidence ranges at sample sizes four
  and five. Existing local calculations bound a first-coordinate beta cap;
  they do not establish the required active-prefix bound on repeated-lowest
  knot faces. See [the precise remaining inequality](supporting-materials/theory/ORDERED-COLLISION-OBLIGATION.md).
- **Separate results:** the two-observation proof, the independent fixed-level
  certificates through sample size seven, and the valid modified procedures
  remain in the [supporting materials](supporting-materials/README.md).
  They are not proofs of the conjectural higher-dimensional extensions.
- **Open:** general finite-sample coverage of ordinary Stringer at
  conventional confidence levels under independent sampling.

The distinct fixed-population sampling problem is treated in
[Exact finite-frame inference under one-start systematic PPS](supporting-materials/paper/systematic-pps.pdf)
([source](supporting-materials/paper/systematic-pps.tex)). Its random-start
inversion argument is separate from the independent-observation cap theorem.

## Reproduce the principal result

Python 3.12, `uv`, and Tectonic are required. Dependencies are pinned in
`uv.lock`.

```sh
make sync
make tetrahedral-vertex-barrier-check terminal-edge-bifurcation-check
make claim-manifest-check public-corpus-check
make paper
```

The tetrahedral check reconstructs the cap and centroid formulas, verifies
the rational multiplier identity, and checks all Bernstein coefficient
signs, including the repeated-upper boundary. The terminal-edge check
verifies the derivative identities and numerical enclosures; the
all-dimensional limiting argument is proved in the text. No proof-assistant
formalization is asserted. `make test` runs the wider repository regression
suite, including mutation tests for the public-text policy.

See the [supporting guide](supporting-materials/README.md) for exact file
locations and the other computations. The [claim manifest](supporting-materials/claim-evidence.json)
distinguishes results in the article from results and conjectures in the
research notes.

## Versions

This is **v1.1.0**. The [GitHub release](https://github.com/RichieSater/stringer-bound/releases/tag/v1.1.0)
contains the manuscript PDFs and the exact source package. The
[Zenodo series DOI](https://doi.org/10.5281/zenodo.21850819) resolves to the
latest archived version; `CITATION.cff` cites this version using that stable
series identifier.

[Archived v1.0.0](https://doi.org/10.5281/zenodo.21850820) predates the
current article and is preserved unchanged.
