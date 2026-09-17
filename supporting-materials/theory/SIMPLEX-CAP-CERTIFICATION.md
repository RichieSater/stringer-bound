# A reusable certification criterion for ordered-simplex cap comparisons

## 1. The mathematical interface

Let `n >= 2`, let

```text
0 = z_0 <= z_1 <= ... <= z_n = 1,
(D_0,...,D_n) ~ Dirichlet(1,...,1),
```

and let `s(z)` be an affine threshold. The common proof obligation behind
the repository's `n=3,...,7` results is

```text
P(sum_i z_i D_i > s(z)) <= alpha                         (1)
```

on the complete ordered-knot simplex. For Stringer's binomial factors,

```text
s(z) = p_n(0)
     + sum_{k=1}^{n-1} [p_n(n-k)-p_n(n-k-1)] z_k.       (2)
```

Proposition `prop:capcriterion` in the manuscript proves that (1) is exactly
the pointwise comparison with the valid Gaffke bounded-mean endpoint. Thus
the reusable object is not a particular triangulation or software package;
it is a certificate that closes (1) for a declared `(n,alpha)` domain.

For pairwise distinct knots, the uniform-Dirichlet cap is the divided
difference

```text
H_n(s;z) = sum_{i:z_i>s} (z_i-s)^n
             / product_{j != i}(z_i-z_j).               (3)
```

The confluent divided difference gives the continuous extension when knots
coincide. On a region where the position of `s` among the knots is fixed,
the active terms in (3) are fixed. After multiplying `alpha-H_n` by
denominators of proved sign, the target becomes polynomial
nonnegativity.

## 2. Soundness theorem for a finite certificate

**Proposition (certificate soundness).** Suppose a declared alpha domain and
factor locus satisfy all of the following.

1. Exact factor enclosures cover the full alpha domain and establish every
   ordering and factor identity used.
2. A finite collection of threshold regions covers the complete ordered-knot
   simplex, including its faces.
3. Exact algebra reduces the cap claim either directly to polynomial
   residuals proportional to `alpha-H_n`, or to a section derivative and a
   finite family of local section barriers whose active-prefix step-cap
   values are at most `alpha`. On a repeated-knot face, only strict knot
   gaps are active; a bound at an inactive coordinate prefix is insufficient.
4. The domain is covered either by finite oriented simplex chains or by an
   exact stationary-face induction using those local barriers.
5. Every polynomial residual used by the selected route is nonnegative on its
   declared parameter domain. Bernstein proofs must check every coefficient;
   any omitted coefficient is proved symbolically to vanish on the factor
   locus or a specified face.
6. Repeated-knot and threshold-boundary cases follow from direct formulas or
   continuity.

Then (1) holds on the complete ordered-knot simplex throughout the declared
alpha domain. With the Stringer threshold (2), the Stringer bound pointwise
dominates Gaffke and inherits its finite-sample coverage under the independent
common-distribution model.

**Proof.** In the direct route, condition 3 reduces (1) to `R >= 0`. Under
any simplex parameterization, the Bernstein basis functions are nonnegative
and sum to one. Condition 5 therefore implies `R >= 0` on each simplex, and
condition 4 transfers this inequality to the full region. In the structural
route, the exact section derivative identifies every stationary cap with its
section centroid. The local barriers in conditions 3 and 5 bound that value
by a calibrated step cap, while condition 4 inducts over the finite knot-face
lattice. Condition 2 covers every nondegenerate ordered knot vector in either
route, and condition 6 extends the result to all remaining faces. The final
assertion is the simplex-cap criterion plus the cited finite-sample validity
theorem for Gaffke. ∎

This proposition separates the reusable ordinary mathematics from the
finite sign arithmetic. It also gives a fail-closed audit of a certificate:
a missing region, an unproved denominator sign, an unexplained zero, or an
uncovered degenerate face is a missing proof obligation rather than a small
numerical detail.

## 3. Challenge and solution layers

The repository exposes the interface in two machine-readable layers:

- [`simplex-cap-challenge.json`](../computations/challenges/simplex-cap-challenge.json)
  contains only the theorem template, declared parameter domains, and eight
  verification obligations. It contains no solution paths, digests, or
  commands.
- [`simplex-cap-solution-index.json`](../computations/certificates/simplex-cap-solution-index.json)
  binds each declared instance to its written reduction, deterministic
  structure data, sign artifact, checker source, hashes, and regeneration
  commands.

Run

```sh
make simplex-cap-interface-check
```

to validate the separation, instance scopes, source existence, and hashes and
to regenerate the solution index byte-for-byte. This integrity check does
not replace the proof-essential instance commands listed in the index.

## 4. Proved direct instances and the structural extension

| instance | alpha domain | conceptual reduction | finite calculation |
|---|---|---|---|
| `n=3` direct | every `0.01 <= alpha <= 0.20` | three cap regions; AM--GM; monotonicity in one knot | 414 contiguous rational alpha cells and exact Bernstein signs |
| `n=4` direct | `alpha=0.10,0.05,0.01` | cap-region decomposition and fixed tetrahedra | exact rational Bernstein signs |
| `n=5` direct | same three levels | five cap regions and face ideals | integer-directed dyadic signs on four-simplices |
| `n=6` | same three levels | reflection, six regions, generic face ideals | directed dyadic signs on a fixed chain of 32 five-simplices |
| `n=7` | same three levels | reflection, seven regions, factor-order and outer-face identities | outward-rounded Arb signs on a fixed chain of 64 six-simplices |

The four-coordinate section-centroid theorem separately proves the full
`n=3` interval `0<alpha<=((19+sqrt(21))/34)^3`. The proposed structural
`n=4` and `n=5` intervals through `alpha_4` and `alpha_5` are
conjectural. Their exact weight-order calculations and local Bernstein
identities do not yet supply obligation 4 on repeated-lowest-knot faces.
The missing active-prefix implication is stated in
[`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).
The direct certificates in the table do not use that implication.

The structural-zero layer is essential. A tiny computed coefficient is not
treated as zero: it is either enclosed strictly away from zero or matched to
an exact ideal-membership, factor-order, or boundary identity. Likewise, a
list of positive coefficients proves nothing until the region chain and its
boundary cancellation have been established.

## 5. What can be reused

The criterion applies to any affine ordered-sample upper bound whose
threshold can be compared with a uniform-Dirichlet mean quantile. A new
instance needs only:

1. its factor locus and parameter domain;
2. the normalized affine threshold;
3. an exhaustive threshold-region decomposition;
4. exact residual identities and denominator signs;
5. a fixed simplex cover; and
6. symbolic-zero and rigorous coefficient-sign layers.

The transferable insight is the separation of geometry from arithmetic.
The cap formula and soundness theorem are dimension-independent; only the
region decomposition and finite coefficient set grow with `n`. This makes
the current `n=8` obstacle precise: it is not an undefined request for more
computation, but the need for a tractable complete region chain and a
structural description of its zero coefficients.

## 6. Verification boundary

Ordinary mathematics proves the normalization, cap formula, transfer to
Gaffke, certificate-soundness proposition, and the exact meaning of every
finite obligation. Symbolic programs derive residuals, substitutions, face
relations, and structural zeros. Python `Fraction`, integer-directed dyadic
arithmetic, or Arb proves the remaining finite signs. Gaffke's finite-sample
validity is an external theorem. No proof assistant is used, and the
challenge/solution split is a machine-readable comparator interface rather
than a Lean formalization.
