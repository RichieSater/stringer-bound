# Practice note: a finite-sample floor for a Stringer MUS evaluation

## Executive summary

The repository does **not** yet prove that ordinary Stringer is conservative
for arbitrary audit sample size.  It does provide a procedure that an audit
methodology group can evaluate without waiting for that open problem:

> Before observing the sample, define the reported taint-rate upper bound as
> the larger of (i) the ordinary Stringer calculation and (ii) the validated
> one-sided Gaffke bounded-mean limit.

Under the manuscript's independent `[0,1]`-taint model, this complete rule has
finite-sample coverage at least the selected confidence level for every
sample size.  It never lowers a Stringer result.  When Stringer is already
larger, the reported number is unchanged.

This is a candidate statistical control for technical-methodology review,
not a statement that any audit standard has adopted the procedure.

## The rule

For sample taints `t_1,...,t_n`, tail probability `alpha`, and the firm's
chosen binomial or Poisson Stringer factors, compute

\[
U_{\rm safe}=max\{U_{\rm Stringer},U_{\rm Gaffke}\},
\]

where

\[
U_{\rm Gaffke}=Q_{1-\alpha}
\left(D_0+\sum_{i=1}^n t_iD_i\right),
\qquad
(D_0,\ldots,D_n)\sim\operatorname{Dirichlet}(1,\ldots,1).
\]

The proof is short but the procedural distinction matters.  Gaffke's limit
has coverage at least `1-alpha`; the maximum is never below that limit.
Coverage therefore belongs to the **pre-specified maximum rule**, not to an
after-the-fact assertion that ordinary Stringer was valid.

At `n=3,4,5` and 90%, 95%, and 99% confidence, exact pointwise theorems in
this repository prove that the rule is ordinary Stringer on every possible
sample. The same statement holds at `n=6` and 95% confidence. At other
sample-size and confidence-level combinations, the Gaffke component is a valid floor while
the general pointwise comparison and the Stringer conjecture remain open.

At every nominal confidence level of at least
`75%`, an additional analytic theorem supplies a simple
zero-uplift check for **every sample size**: if the
**binomial-factor** Stringer value is at least the largest observed taint,
then the valid Gaffke limit is no larger. The pre-specified safeguard
therefore returns ordinary Stringer on that sample, for either factor
convention. The check is sufficient, not necessary; failure of the condition
says only that the Gaffke endpoint must be computed. It is a pointwise
statement about the safeguard's output, not a conditional-coverage guarantee
for unmodified Stringer. The proof and supporting checks are in
[`ONE-CAP-COMPARISON.md`](../theory/ONE-CAP-COMPARISON.md).

A separate theorem validates the ordinary **Poisson-factor** Stringer bound
for every `n<=8` at 90%, every `n<=11` at 95%, and every `n<=20` at 99%
confidence. It uses a corrected simultaneous survival band rather than the
Gaffke comparison. Consequently, it proves coverage in those ranges but does
not assert that the safeguard's Gaffke component is pointwise inactive there.
Beyond those ranges, the safeguard remains the proved all-sample-size option
under the stated model.

## Reproducible command

From the repository root, with zero taints omitted from the input list but
the full sample size retained in `--n`:

```sh
uv run --frozen python \
  supporting-materials/computations/python/gaffke.py \
  --n 100 --alpha 0.05 --method poisson --taints 1,0.4,0.1
```

The current output is:

```json
{
  "alpha": 0.05,
  "gaffke": 0.05215837212760732,
  "gaffke_bracket_width": 1.7763568394002505e-14,
  "gaffke_certificate_bits": 48,
  "gaffke_lower": 0.05215837212758956,
  "gaffke_upper_dyadic": "3670319144971/70368744177664",
  "governing_bound": "stringer",
  "method": "poisson",
  "n": 100,
  "safeguarded": 0.055104224503892164,
  "stringer": 0.055104224503892164,
  "uplift": 0.0
}
```

For this illustrative sample, the safeguard does not change the Poisson
Stringer result.  That observation is sample-specific, not a proof of the
unmodified rule at `n=100`.

## Why the numerical floor is certifiable

The Gaffke endpoint is a quantile of a uniform-Dirichlet average.  A B-spline
calculation proposes its location, but no floating-point sign is trusted.
The program:

1. interprets decimal taints and `alpha` as exact rationals;
2. evaluates each endpoint tail as an exact confluent divided difference;
3. checks that the lower endpoint has tail at least `alpha` and the upper
   endpoint has tail at most `alpha`;
4. expands the dyadic bracket if either exact check fails; and
5. reports the upper endpoint, which is conservative for the mathematical
   Gaffke quantile.

Repeated taints and the many zero observations typical of MUS are handled
directly.  The implementation computes only the Stringer factor prefix
needed by the observed nonzero taints, so a zero-heavy sample does not require
solving all `n+1` Poisson equations.

## Suggested engagement-file record

If a methodology owner authorizes evaluation of the rule, retain at least:

- the full sample size, confidence level, and factor convention;
- the complete taint list, including a count reconciliation for omitted
  zeros;
- the ordinary Stringer component;
- both endpoints and the width of the certified Gaffke bracket;
- the safeguarded maximum and uplift;
- the repository commit and locked dependency file used; and
- the unedited JSON output.

Any conversion from a taint-rate bound to a projected currency misstatement
should follow the approved MUS design.  Under the usual dollar-weighted
interpretation, the population mean taint is the population overstatement
divided by the relevant recorded amount, but that identity depends on the
sampling frame and taint definition actually used.

## Scope boundaries requiring methodology review

The theorem and command do not, by themselves, validate:

- sampling without replacement or systematic PPS designs that are not
  represented by the independent model;
- negative taints, credits, understatements, or taints outside `[0,1]`;
- contaminated or incorrectly constructed sampling frames;
- stratification, certainty items, or separately evaluated high-value items;
- qualitative audit evidence, tolerable-misstatement decisions, or risk
  assessments; or
- compliance with AICPA, PCAOB, IAASB, governmental, or firm-specific
  requirements.

Before production use, an audit organization should obtain an independent
human proof review, validate the implementation against approved factor
tables and test vectors, decide how the rule interacts with its sampling
designs, and subject the workflow to its normal model-risk and software-change
controls.

## Relationship to the research program

The safeguard creates an all-sample-size valid reporting option, but it does
not end the mathematical program. Proving ordinary Stringer conservative at
conventional levels for all `n` would validate the familiar calculation
without the floor and could show that the floor is never active. The direct
Poisson theorem now reaches `n=8`, `n=11`, and `n=20` at 90%, 95%, and 99%,
respectively, but audit sample sizes can be larger. Extending the
exact pointwise comparison beyond the current `n=6`, 95% frontier, or finding a dimension-free
simplex-cap argument, therefore remains the central theoretical objective.
The one-cap theorem already controls, at nominal confidence of at least
`75%`, the complete region in which the binomial Stringer threshold is
at least the largest observed taint for every sample size. The remaining cap
regions are the next target.
