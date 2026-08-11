# Independent human review packet

## Status

An independent human proof review remains required before journal submission
or practice-facing reliance.  This packet makes the requested review bounded,
traceable, and reproducible; it is not itself a review or a sign-off.

The reviewer should record the exact Git commit reviewed.  A useful review
must be conducted by a person who did not develop the proofs and who is able
to assess measure-theoretic probability inequalities and computer-assisted
certificates.

## Requested dispositions

For each item below, record one of:

- **verified** -- the stated conclusion follows from the supplied argument;
- **verified subject to cited theorem** -- the internal reduction is sound,
  with an explicitly named external result accepted as an assumption;
- **minor correction** -- exposition or a locally repairable omission;
- **major gap** -- a missing inference that affects the claim; or
- **not reviewed**.

Please distinguish mathematical validity from novelty, journal fit, and
professional-standard adoption.

## Review A: the sharp `n=2` theorem

Primary files:

- [`../theory/N2-PROOF.md`](../theory/N2-PROOF.md)
- manuscript Theorem 1 and its proof
- [`REVIEW-N2.md`](REVIEW-N2.md), which discloses the earlier non-attainment
  gap and its repair

Questions:

1. Are the Clopper--Pearson factors and the reduction
   `SB=1-(A min+B max)` correct for every `alpha in (0,1)`?
2. Do the rectangle bounds cover all endpoints and atomic distributions?
3. Does the boundary-pair substitution justify every survival lower bound?
4. Is the concavity/potential inequality used on its valid domain only?
5. Do the three strictness cases exhaust the possibilities?
6. Do both sharpness families approach noncoverage `alpha`?
7. In the repaired equality analysis, does equality force the stated
   two-level survival function, and does every case place positive mass on
   the wedge boundary?

Requested output: a line-by-line proof report, including any silent
regularity assumption and an explicit disposition of sharpness and
non-attainment separately from conservatism.

## Review B: conventional confidence levels at `n=3,4,5`

Primary files:

- [`../theory/N3-CONVENTIONAL.md`](../theory/N3-CONVENTIONAL.md)
- [`../theory/N4-CONVENTIONAL.md`](../theory/N4-CONVENTIONAL.md)
- [`../theory/N5-CONVENTIONAL.md`](../theory/N5-CONVENTIONAL.md)
- the six derivation/certificate programs and their committed JSON artifacts
- manuscript Theorems~`thm:n3`, `thm:n4`, and `thm:n5`

Questions:

1. Accepting the cited Vlassis--Thomas validity theorem, is pointwise
   domination of the Gaffke endpoint sufficient for the stated coverage?
2. Are the affine normalization and simplex-cap formulas complete on all
   boundary faces?
3. Do the regional decompositions cover the complete ordered-knot domains?
4. Are all removed factors proved nonnegative on their claimed regions?
5. Are structural zeros established symbolically rather than by numerical
   tolerance?
6. Do integer-checked factor brackets and directed rational interval
   operations imply every reported positive Bernstein sign?
7. Does `make reproduce` regenerate the proof-essential artifacts without
   an unrecorded computational dependency?

Requested output: separate dispositions for each sample size.  A reviewer
may verify `n=3` fully and audit the certificate architecture for `n=4,5`,
but should say exactly which level of review was completed.

## Review C: the all-`n` Poisson factor comparison

Primary file:

- [`../theory/POISSON-DOMINATION.md`](../theory/POISSON-DOMINATION.md)
- manuscript Proposition~`prop:poissondomination` and
  Corollary~`cor:poissoncoverage`

Questions:

1. Is the Anderson--Samuels inequality quoted with the correct strictness,
   index, and parameter regime?
2. Does `alpha<e^-1` imply `lambda_j>j+1` for every `j` used?
3. Are the cases `lambda_j<n` and `lambda_j>=n` exhaustive and correct?
4. Does summation by parts prove pointwise domination when taints are in
   `[0,1]` and ordered nonincreasingly?
5. Are the paper's coverage transfers limited to binomial cases already
   proved?

## Review D: the all-sample-size reporting safeguard

Primary files:

- [`../theory/GAFFKE-SAFEGUARD.md`](../theory/GAFFKE-SAFEGUARD.md)
- [`../computations/python/gaffke.py`](../computations/python/gaffke.py)
- [`PRACTICE-SAFEGUARD.md`](PRACTICE-SAFEGUARD.md)
- manuscript Proposition~`prop:safeguard`

Questions:

1. Accepting the cited Gaffke validity theorem and bounded-mean inversion,
   does the pre-specified maximum have the claimed all-`n` coverage?
2. Is the distinction from ordinary Stringer stated unambiguously?
3. Is the divided-difference identity for a uniform-Dirichlet average
   normalized correctly?
4. Does the confluent residue formula handle repeated knots and knots equal
   to the evaluation point?
5. Do the exact endpoint signs prove that the reported dyadic upper endpoint
   is conservative for the target quantile?
6. Are the audit-design limitations prominent enough to prevent an iid
   theorem from being misapplied to a different sampling design?

## Minimum reproduction record

Run from a clean checkout:

```sh
git rev-parse HEAD
uv --version
tectonic --version
make reproduce
uv run --frozen python \
  supporting-materials/computations/python/gaffke.py \
  --n 100 --alpha 0.05 --method poisson --taints 1,0.4,0.1
```

Record the operating system, architecture, command exit status, and any
artifact mismatch.  A passing run supports reproducibility; it does not
replace the mathematical review.

## Suggested sign-off block

```text
Reviewer name:
Affiliation (optional):
Relevant expertise:
Conflict-of-interest statement:
Repository commit:
Review dates:

Review A disposition:
Review B disposition (n=3 / n=4 / n=5):
Review C disposition:
Review D disposition:

Issues found and exact locations:
Corrections rechecked:
Permission to publish this review record: yes / no
Signature or authenticated email reference:
```
