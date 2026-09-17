# Practice note: finite-sample reporting controls for a Stringer MUS evaluation

## Executive summary

The repository does **not** yet prove that ordinary Stringer is conservative
for arbitrary audit sample size.  It does provide two families of procedures
that an audit methodology group can evaluate without waiting for that open
problem.
The principal safeguard in this note is:

> Before observing the sample, define the reported taint-rate upper bound as
> the larger of (i) the ordinary Stringer calculation and (ii) the validated
> one-sided Gaffke bounded-mean limit.

Under the manuscript's independent `[0,1]`-taint model, this complete rule has
finite-sample coverage at least the selected confidence level for every
sample size.  It never lowers a Stringer result.  When Stringer is already
larger, the reported number is unchanged.

This is a candidate statistical control for technical-methodology analysis,
not a statement that any audit standard has adopted the procedure.

The rule above is for independent taints. A separate proved rule now covers
the exact one-uniform-start systematic-PPS design. It combines capped
ordinary Stringer with a sharp randomization floor and a deterministic
item-weight completion bound. When the complete fixed frame is available, an
exact randomization inversion uses every possible start and is pointwise no
larger than that closed-form fallback; see the design-specific section below.

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

At `n=3` for every nominal confidence of at least `66.6314788...%`,
and at `n=4,5,6,7` for 90%, 95%, and 99%, pointwise theorems in this
repository prove that the rule is ordinary Stringer on every possible
sample. The proposed full confidence intervals at `n=4,5` remain
conjectural; their repeated-lowest-knot face implication is unresolved
([exact obligation](../theory/ORDERED-COLLISION-OBLIGATION.md)).
The direct fixed-level certificates do not use that implication. At other sample-size and
confidence-level combinations, the Gaffke
component is a valid floor while
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

For `2<=n<=200`, a certified two-cap test extends the zero-uplift region when
binomial Stringer lies between the second-largest and largest observed
taints. Put

\[
q=\frac{1-t_{\max}}{1-U_{\rm Stringer,B}}.
\]

The Gaffke endpoint is no larger than binomial Stringer when `q>=5/12` at
90%, `q>=1/3` at 95%, or `q>=1/4` at 99%. This is again a sufficient
sample-wise identity test for the pre-specified safeguard, not a
conditional-coverage statement. The exact formula and 59,700 endpoint
certificates are in
[`TWO-CAP-COMPARISON.md`](../theory/TWO-CAP-COMPARISON.md).

A separate theorem validates the ordinary **Poisson-factor** Stringer bound
for every `n<=8` at 90%, every `n<=11` at 95%, and every `n<=20` at 99%
confidence. It uses a corrected simultaneous survival band rather than the
Gaffke comparison. Consequently, it proves coverage in those ranges but does
not assert that the safeguard's Gaffke component is pointwise inactive there.
Beyond those ranges, the safeguard remains the proved all-sample-size option
under the stated model.

## Separate rule for one-start systematic PPS

Let the ordered population have known positive book weights `w_i`, audited
taints in `[0,1]`, total weight `W`, and one uniform random start in an
interval of length `W/n`. For the `n` equally spaced sample hits, let
`sample_mean` be the mean taint and let `S` be the set of distinct sampled
item identities. Compute

```text
M = 1 - alpha + alpha * sample_mean
C = (sum_{i in S} w_i*t_i + sum_{i not in S} w_i) / W
H = min(M, C)
U_systematic = max(min(1, ordinary Stringer), H)
```

The sample mean is exactly design-unbiased. Markov's inequality proves that
`M` has at least nominal coverage, and no universally design-valid taint-only
rule can be smaller on a constant sample. The completion value `C` is a
deterministic upper bound because every unobserved taint is set to its maximum
one. Therefore `H` and `U_systematic` retain nominal design coverage, while
the last maximum never lowers capped ordinary Stringer.

For a complete known rational frame, the sharper exact calculation is

```text
B_ff = maximum population target over all taint completions that
       (a) match the audited item identities and taints, and
       (b) give the observed grid mean exact lower-tail mass > alpha
U_systematic_ff = max(min(1, ordinary Stringer), B_ff)
```

Exact randomization-rank inversion proves that `B_ff` has nominal design
coverage. It also proves `B_ff <= H` on every observation. The implementation
enumerates every inclusion-minimal phase coalition of probability above
`alpha`, solves its least-favourable-completion LP by exact rational simplex,
checks matching feasible primal and dual objectives, and checks a full witness
population. It does not treat a floating-point MILP gap as a certificate. The
word *minimax* refers to maximizing over compatible
taint completions for this fixed grid-mean ordering, not to global optimality
among confidence rules.

If the phase atoms partition the frame, the LPs collapse to an exact
minimum-mass coalition. For observed phase `j` with mean `y`, compute

```text
rho = min{phase_probability(A): j in A and phase_probability(A) > alpha}
B_ff = 1 - (1-y)*rho.
```

For equal phases, `rho=(floor(m*alpha)+1)/m`. For rational unequal masses,
`systematic_disjoint_dp.py` clears denominators and solves the resulting
subset-sum problem exactly with checked back-pointers. Its scaled-state limit
fails closed independently of the general LP limits.

In the committed 95% benchmark, 2,000 unit items and sample size 100 produce
20 equiprobable disjoint phases. After one all-zero phase, `M=C=H=0.95`, but
the exact 190-coalition calculation gives `B_ff=0.90`. Both ordinary all-zero
Stringer calculations are below 0.03, so the exact-frame overlay reports 0.90.
See [`SYSTEMATIC-MINIMAX.md`](../theory/SYSTEMATIC-MINIMAX.md) and the
regenerated certificate. The phase-partition formula, exact DP, and its
complexity boundary are in
[`SYSTEMATIC-DISJOINT-PHASES.md`](../theory/SYSTEMATIC-DISJOINT-PHASES.md).
If an explicit exact-search limit is reached, the
implementation reports no finite-frame certificate; the closed-form `H`
remains a proved fallback and is pointwise at least the mathematical `B_ff`.

With independent starts, each complete grid average is one replication. In
the equal-disjoint-phase all-zero case, `r` starts visiting `d` distinct
phases give the exact finite-frame value

```text
1 - max(d, floor(m*alpha^(1/r))+1)/m.
```

For nonzero observations, the proved generic route applies a validated
bounded-mean procedure to the independently randomized grid averages; the
displayed finite-frame formula is not asserted outside the all-zero case.

This correction is necessary in the stated distribution-free model:
ordinary binomial Stringer fails on some ordered population for every
`n>=2` and every confidence level. Ordinary Poisson Stringer likewise fails
at 90%, 95%, and 99% on some population whenever `n>=3`, `n>=4`, and `n>=5`.
In the exact 100-draw example, both ordinary conventions have 10% coverage at
all three levels. On its all-zero phase at 95%, the safe output is 10%, about
7.0487 percentage points above binomial Stringer and 7.0043 points above
Poisson Stringer. The proof, exact finite census, and uplift brackets are in
[`SYSTEMATIC-PPS.md`](../theory/SYSTEMATIC-PPS.md).

This design-specific rule assumes one genuinely uniform start, a fixed
ordered frame, known book weights, and taints bounded above by one. It does
not inherit the i.i.d. Gaffke or Poisson-calibration calculations above, and
it does not cover successive PPS, stratified combinations, nonrandom starts,
negative taints, or uncertain book amounts.

## Alternative all-sample-size calibrations within the Poisson factor family

The corrected band also yields a second rule that a methodology group can
evaluate. Let \(\kappa_{n,\alpha}\) be the smallest scalar satisfying
\(\kappa\ge\max\{1,n/\lambda_n\}\) for which

\[
\Pr\left\{V_{i:n}\le
\min\left(1,\frac{\kappa_{n,\alpha}\lambda_{i-1}}n\right),
\ i=1,\ldots,n\right\}\ge1-\alpha,
\]

where \(\lambda_j\) is the ordinary Poisson count limit. Then report

\[
U_{\rm scalar}
=\min\{1,\kappa_{n,\alpha}U_{\rm Stringer,P}\}.
\]

This rule is distribution-free valid for every sample size and every
confidence level under the same i.i.d. model. It keeps the existing Poisson
factor curve and requires only one precomputed multiplier. Exact certificates give these simple
six-decimal valid choices:

| nominal confidence | `n=25` | `n=50` | `n=100` | `n=200` |
|---:|---:|---:|---:|---:|
| 90% | 1.189306 | 1.251208 | 1.286176 | 1.305049 |
| 95% | 1.126246 | 1.195804 | 1.235956 | 1.257979 |
| 99% | 1.027950 | 1.111273 | 1.161122 | 1.189239 |

The multiplier applies to the complete Poisson Stringer result and the
product is capped at one. The theorem also permits each Poisson factor to be
capped at one before forming and scaling the complete expression; a
methodology must document which convention it uses. These values validate
the **modified scalar rule**;
they do not imply that ordinary Stringer fails at those sample sizes. They
are minimal only for this particular corrected-band criterion, to the
certified resolution. See
[`POISSON-BAND-CALIBRATION.md`](../theory/POISSON-BAND-CALIBRATION.md).

If a methodology needs one conservative multiplier that is independent of
sample size, an analytic corollary permits `2` at every `n` whenever nominal
confidence is at least `2/e`, approximately `73.6%`. This includes 90%, 95%,
and 99%. The value `2` is only a simple uniform upper bound; it is not claimed
to be necessary or optimal. The elementary closed-form criterion gives
`1.76`, `1.66`, and `1.52` at 90%, 95%, and 99%. Retaining the first 200
Poisson crossing probabilities and bounding the remainder geometrically
reduces the certified uniform choices to `1.525906`, `1.436135`, and
`1.320081`. The sample-size-specific values in the table are smaller still.

A pointwise tighter implementation uses

\[
\widehat p_j=\min\{1,\kappa_{n,\alpha}\lambda_j/n\}
\]

as the calibrated factor curve and forms the Stringer expression from that
curve. This has the same coverage guarantee and is never larger than
multiplying the complete result and capping only at the end. It also removes
the distinction between first capping or not capping the ordinary Poisson
factors.

The six-decimal multiplier table controls rounding of the multiplier, not
rounding in a third-party Poisson factor table. A production implementation
must compute the underlying factors with a certified upper enclosure or
show that the approved table rounds them upward. The command below performs
the former check with exact rational endpoints.

### A zero-taint-preserving path

The full-scale calibration raises the no-error factor. A second proved path
keeps that factor fixed. With

\[
\bar p_j=\min\{1,\lambda_j/n\},
\]

precompute the smallest admissible \(\eta_{n,\alpha}\ge1\) for which the
corrected band probability based on

\[
d_j=\bar p_0+\eta_{n,\alpha}(\bar p_j-\bar p_0)
\]

is at least \(1-\alpha\), including the terminal-factor condition. Then
either factor convention has a valid anchored report:

\[
\begin{aligned}
U^{\rm A}_{\rm capped}
&=\min\{1,\bar p_0+\eta_{n,\alpha}
(U_{\rm Stringer,\bar P}-\bar p_0)\},\\
U^{\rm A}_{\rm untruncated}
&=\min\{1,p_0^{\rm P}+\eta_{n,\alpha}
(U_{\rm Stringer,P}-p_0^{\rm P})\}.
\end{aligned}
\]

On an all-zero sample, these expressions return the ordinary zero-taint
factor exactly before the final cap. The price is a no-smaller multiplier on
the error increments. Neither anchored rule uniformly dominates the
full-scale rule. A methodology must choose the calibration path before
seeing the sample; selecting the smaller reported result afterward is not
covered by either theorem.

Here too, a methodology can cap every effective anchored factor at one
before forming the Stringer expression. The capped and untruncated base
curves then produce the same adjusted factors. This version has the same
coverage guarantee, preserves the all-zero result, and is pointwise no
larger than capping only the final affine report.

The same exact certificate gives these simple six-decimal valid choices for
the anchored error-increment multiplier:

| nominal confidence | `n=25` | `n=50` | `n=100` | `n=200` |
|---:|---:|---:|---:|---:|
| 90% | 1.820063 | 2.299904 | 2.764607 | 3.219411 |
| 95% | 1.511563 | 1.947538 | 2.367560 | 2.778290 |
| 99% | 1.101565 | 1.480516 | 1.839132 | 2.188420 |

For a fixed factor convention, let `s` be ordinary Stringer and `p_0` its
no-error factor. Before the final cap, the anchored report minus the
full-scale report is

\[
(\eta_{n,\alpha}-\kappa_{n,\alpha})s
-(\eta_{n,\alpha}-1)p_0.
\]

The anchored report is therefore lower on an all-zero sample, but it can be
higher after observed errors raise `s`. The paths are not pointwise ordered;
the methodology choice should reflect the expected sample profile and must
be made before seeing the sample.

The scalar rules and the Stringer--Gaffke maximum solve the same validity
problem differently. The scalar rules are operationally close to existing
factor-table workflows. A nontrivial full-scale path raises every sample
result at a fixed `(n,alpha)`; the anchored path instead preserves
the all-zero result. The Gaffke safeguard is sample-adaptive and often has zero
uplift, but it requires a second bounded-mean calculation. No general
pointwise ordering is claimed. A methodology evaluation should compare both
on representative engagement populations before selecting either rule.

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

For an arbitrary specified sample size and a conventional tail probability
`alpha < exp(-1)`, write the exact calibration certificate to an evidence file
and extract its conservative upper endpoint:

```sh
uv run --frozen python \
  supporting-materials/computations/python/poisson_band_calibration.py \
  --n 25 --alpha 0.05 --taints 1,0.4,0.1 \
  --out /tmp/poisson-kappa.json
jq '{kappa_upper:.case.kappa_upper,
     eta_upper:.zero_anchor_case.eta_upper,report}' /tmp/poisson-kappa.json
```

The command can take longer at large `n` because it brackets every required
Poisson limit and evaluates the joint probability with exact rational
arithmetic. Most decimal fields summarize exact numerators and denominators
stored alongside them; the
`case.kappa_upper.valid_decimal_ceiling_12` and
`zero_anchor_case.eta_upper.valid_decimal_ceiling_12` fields are instead
explicitly rounded upward for direct use. For this example the valid choices
are `1.126245908440` and `1.511562403292`. The ordinary Poisson component is
enclosed above by `0.2204168980155687`; the full-scale and anchored reports
are enclosed above by `0.2482436295408857` and `0.2718737360313792`,
respectively. The JSON also reports rigorous upper enclosures for the
pointwise no-larger calibrated-factor-capped versions; they coincide with
the final-cap values in this three-error example because none of the factors
used is above one. The ordinary component and the final-cap fields use the
untruncated Poisson factors; the JSON records that base convention and labels
both calibrated capping variants explicitly. It also records the exact taint
multiset, including the implied zero count, and every dyadic Poisson-limit
bracket used in the calculation. At 90%, 95%, and 99%, the
`elementary_uniform_full_scale_case` and
`refined_uniform_full_scale_case` fields supply the separately certified
multipliers that are valid without a sample-size lookup.

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
- if a scalar path is used, its pre-specified path name, whether capping is
  applied only at the end or to each calibrated factor, whether the
  multiplier is sample-size-specific or uniform, the exact certified
  multiplier and its applicable certificate fields, the uncapped ordinary
  Poisson result, and the calibrated result;
- if the systematic rule is used, the ordered frame total, random-start
  convention, sampled item identities and book weights, `sample_mean`, `M`,
  `C`, `H`, capped ordinary Stringer, and the final maximum;
- the repository commit and locked dependency file used; and
- the unedited JSON output.

Any conversion from a taint-rate bound to a projected currency misstatement
should follow the approved MUS design.  Under the usual dollar-weighted
interpretation, the population mean taint is the population overstatement
divided by the relevant recorded amount, but that identity depends on the
sampling frame and taint definition actually used.

## Scope boundaries for deployment

The i.i.d. theorem and command at their unadjusted tail do not, by themselves,
validate:

- ordinary unadjusted evaluation under systematic or successive PPS designs.
  The one-start systematic design has the separate rule above; uniform
  sampling of distinct finite-population indices has a different exact
  bridge that requires the
  pre-specified adjustment `eta=alpha*(N)_n/N^n`; see
  [`SRSWOR-CONDITIONING.md`](../theory/SRSWOR-CONDITIONING.md);
- negative taints, credits, understatements, or taints outside `[0,1]`;
- contaminated or incorrectly constructed sampling frames;
- stratification, certainty items, or separately evaluated high-value items;
- qualitative audit evidence, tolerable-misstatement decisions, or risk
  assessments; or
- compliance with AICPA, PCAOB, IAASB, governmental, or firm-specific
  requirements.

The theorem does not by itself authorize production deployment. A deployed
implementation must be tested against the certified factor tables and test
vectors, integrated explicitly with the sampling design actually used, and
covered by the organization's model-risk and software-change controls.

## Relationship to the research program

The safeguard creates an all-sample-size valid reporting option, but it does
not end the mathematical program. Proving ordinary Stringer conservative at
conventional levels for all `n` would validate the familiar calculation
without the floor and could show that the floor is never active. The direct
Poisson theorem now reaches `n=8`, `n=11`, and `n=20` at 90%, 95%, and 99%,
respectively, but audit sample sizes can be larger. Extending the
exact pointwise comparison beyond the current `n=7` frontier, or finding a dimension-free
simplex-cap argument, therefore remains the central theoretical objective.
The one-cap theorem already controls, at nominal confidence of at least
`75%`, the complete region in which the binomial Stringer threshold is
at least the largest observed taint for every sample size. The remaining cap
regions are the next target.
