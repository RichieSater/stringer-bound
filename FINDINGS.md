# Stringer bounds: scope and practical consequences

## Independent observations

The [principal article](supporting-materials/paper/stringer.pdf) concerns
three independent observations from a common distribution on `[0,1]`.
It proves binomial Stringer coverage at every confidence of at least
**66.6314788...%**, using a tetrahedral cap theorem and Gaffke's valid mean
bound. The endpoint is sharp for that pointwise comparison, not asserted to
be the boundary of coverage itself.

The five- and six-coordinate global extensions are **conjectural**. The
existing local computations do not establish the required inequality on
repeated-lowest-knot faces; see the
[exact remaining obligation](supporting-materials/theory/ORDERED-COLLISION-OBLIGATION.md).
In particular, the complete confidence ranges previously proposed at sample
sizes four and five do not follow from those computations. Independent
fixed-level certificates at 90%, 95%, and 99% remain separate results.

Ordinary Stringer coverage at general sample sizes and conventional
confidence levels remains open. Searches without violations are numerical
evidence, not guarantees.

## Valid alternatives are different procedures

The pre-specified maximum of ordinary Stringer and the valid Gaffke bound
has distribution-free coverage at every sample size under independent
sampling. Calibrated Poisson procedures give other valid alternatives.
Neither result proves validity of ordinary Stringer at all sample sizes.
See [the safeguard](supporting-materials/theory/GAFFKE-SAFEGUARD.md) and
[Poisson calibration](supporting-materials/theory/POISSON-BAND-CALIBRATION.md).

## One-start systematic sampling is a different model

The [separate systematic-PPS manuscript](supporting-materials/paper/systematic-pps.pdf)
constructs fixed populations on which ordinary Stringer undercovers.
Dependence among grid points can make a large sample reveal only one random
phase. The same manuscript derives design-valid bounds using sampled item
identities and the known finite frame; these are not consequences of the
i.i.d. cap theorem.

## Supporting results

The two-observation proof, larger fixed-level certificates, exact
low-confidence counterexamples, finite-population conditioning, and
practical-size comparisons remain in the
[supporting guide](supporting-materials/README.md). Each has its own stated
hypotheses and verification method. No guarantee transfers to a different
sampling design without an applicable theorem.
