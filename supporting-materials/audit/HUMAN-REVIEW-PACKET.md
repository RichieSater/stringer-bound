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

## Review B: conventional confidence levels at `n=3,4,5,6`

Primary files:

- [`../theory/N3-CONVENTIONAL.md`](../theory/N3-CONVENTIONAL.md)
- [`../theory/N4-CONVENTIONAL.md`](../theory/N4-CONVENTIONAL.md)
- [`../theory/N5-CONVENTIONAL.md`](../theory/N5-CONVENTIONAL.md)
- [`../theory/N6-CONVENTIONAL.md`](../theory/N6-CONVENTIONAL.md)
- the eight derivation/certificate programs and their committed JSON artifacts
- manuscript Theorems~`thm:n3`, `thm:n4`, `thm:n5`, and `thm:n6`

Questions:

1. Accepting the cited Vlassis--Thomas validity theorem, is pointwise
   domination of the Gaffke endpoint sufficient for the stated coverage?
2. Are the affine normalization and simplex-cap formulas complete on all
   boundary faces?
3. For `n=6`, does the oriented relative-chain certificate (facet
   cancellation, outer boundary, and degree one) establish that the 32
   listed simplices cover the complete ordered-knot domain without interior
   multiplicity?
4. Are all removed factors proved nonnegative on their claimed regions?
5. Are structural zeros established symbolically rather than by numerical
   tolerance?
6. For `n=6`, do the certified nonzero normal minors justify specialization
   of every generic rational-function face-ideal identity to the 90%, 95%,
   and 99% factor values?
7. Do integer-checked factor brackets and directed rational interval
   operations imply every reported positive Bernstein sign?
8. Does `make reproduce` regenerate the proof-essential artifacts without
   an unrecorded computational dependency?

Requested output: separate dispositions for each sample size.  A reviewer
may verify `n=3` fully and audit the certificate architecture for `n=4,5,6`,
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

## Review D: the direct Poisson simultaneous-band theorem

Primary files:

- [`../theory/POISSON-SIMULTANEOUS-BAND.md`](../theory/POISSON-SIMULTANEOUS-BAND.md)
- [`../computations/python/poisson_band_certificate.py`](../computations/python/poisson_band_certificate.py)
- [`../computations/python/stringer.py`](../computations/python/stringer.py),
  especially the exact rational Poisson-limit enclosures
- manuscript Lemma~`lem:poissonband` and Theorem~`thm:poissonband`

Questions:

1. Does the randomized survival transform remain uniform for arbitrary
   atomic distributions, and does it give the claimed order-statistic
   inequality for every threshold?
2. Is the integral identity for a generic nondecreasing factor sequence
   correct, including zero taints and sample-point conventions?
3. Is the index in the sufficient event `V_(i:n) <= p(i-1)` correct?
4. Is Bolshev's recursion stated and implemented with the correct boundary
   sequence and sample-size indices?
5. Do the odd/even rational exponential partial sums, range reduction, and
   dyadic endpoint signs rigorously enclose every required Poisson limit?
6. Does use of lower limit endpoints give a valid lower bound on the event
   probability, and are all sample sizes in each stated range checked?
7. Is the below-nominal event probability at the next size clearly
   distinguished from Stringer undercoverage?

## Review E: the all-sample-size scalar Poisson calibration

Primary files:

- [`../theory/POISSON-BAND-CALIBRATION.md`](../theory/POISSON-BAND-CALIBRATION.md)
- [`../computations/python/poisson_band_calibration.py`](../computations/python/poisson_band_calibration.py)
- [`../computations/certificates/poisson-band-calibration-certificate.json`](../computations/certificates/poisson-band-calibration-certificate.json)
- manuscript Theorem~`thm:poissoncalibration`

Questions:

1. Does multiplying every Poisson factor by `kappa` multiply the complete
   Stringer expression by exactly the same scalar?
2. Does Lemma~`lem:poissonband` apply to the scaled factors, including its
   terminal-factor condition?
3. For the factor-capped variant, is
   `min(1,kappa*min(1,lambda_j/n))=min(1,kappa*lambda_j/n)` when
   `kappa>=1`, and does the same terminal condition still hold?
4. Is `Q_(n,alpha)(kappa)` continuous and nondecreasing, and is its defining
   feasible set nonempty after imposing `kappa>=max(1,n/lambda_n)`?
5. At conventional levels, does the Anderson--Samuels comparison
   give each marginal crossing probability at most `alpha/n`, including
   capped boundaries?
6. Do upper Poisson-limit endpoints at the lower `kappa` endpoint and lower
   Poisson-limit endpoints at the upper `kappa` endpoint give the stated
   opposite bounds on the joint event probability?
7. Do those exact signs bracket the band-minimal scalar, and are all table
   entries rounded upward before being described as valid choices?
8. Is the theorem consistently described as validating the modified scalar
   rule rather than ordinary Stringer, and is minimality limited to this
   sufficient-event family?

## Review F: the all-sample-size reporting safeguard

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

## Review G: the one-cap zero-uplift theorem

Primary files:

- [`../theory/ONE-CAP-COMPARISON.md`](../theory/ONE-CAP-COMPARISON.md)
- [`../computations/python/one_cap_all_n_check.py`](../computations/python/one_cap_all_n_check.py)
- [`../computations/python/one_cap_certificate.py`](../computations/python/one_cap_certificate.py)
- [`../computations/certificates/one-cap-certificate.json`](../computations/certificates/one-cap-certificate.json)
- manuscript Proposition~`prop:onecap`

Questions:

1. Is the one-upper-knot Dirichlet cap formula normalized correctly?
2. Does `s>=x_(n-1)` imply the ordered nonnegative constraints on the
   normalized complement vector?
3. Does the prefix-step decomposition describe the complete weighted-budget
   section, including its boundary?
4. Does concavity of `sum log(1+u_i)` yield the stated maximum over prefix
   vertices with the correct inequality direction?
5. For the ascending Stringer weights, are
   `C_r=1-p_n(n-r)` and the analytic `r=n` equality indexed correctly?
6. Is the conversion of each factor inequality to the binomial upper-tail
   comparison at `q_(n,r)` correct?
7. Are the Anderson--Samuels and Pinelis theorems applied with all of their
   hypotheses satisfied, including the two cases split at `n q_(n,r)=r-1`?
8. Does the `psi`-derivative argument make the boundary multiplier
   increasing in `x`, and do the separate `r=1` argument, the exact
   `r=2,3,4` margins, and the uniform `r>=5` bound prove the comparison for
   every `0<alpha<=1/4`?
9. Is the older 59,700-comparison certificate described only as an
   independent finite regression, rather than as the all-`n` proof?
10. Is the result consistently limited to zero uplift of the pre-specified
   safeguard, rather than conditional coverage of ordinary Stringer?

## Minimum reproduction record

Run from a clean checkout:

```sh
git rev-parse HEAD
uv --version
tectonic --version
make reproduce
make poisson-band-calibration-check
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
Review B disposition (n=3 / n=4 / n=5 / n=6):
Review C disposition:
Review D disposition:
Review E disposition:
Review F disposition:
Review G disposition:

Issues found and exact locations:
Corrections rechecked:
Permission to publish this review record: yes / no
Signature or authenticated email reference:
```
