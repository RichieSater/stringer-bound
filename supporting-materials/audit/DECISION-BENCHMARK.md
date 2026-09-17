# Decision-oriented benchmark for practical-size Stringer calculations

This benchmark makes the statistical results operationally legible without
turning them into audit guidance. It compares fixed, zero-heavy sample
profiles at `n = 25, 100, 200` under three conventional confidence levels.
Every number is regenerated from exact factor brackets, exact Gaffke tail
signs, and certified calibration multipliers.

The cases answer a narrow decision question. If a separately specified
upper-taint-rate threshold is `tau`, two methods can lead to different
numerical comparisons whenever `tau` lies between their reported upper
bounds. The machine-readable artifact gives that complete threshold interval
in exact rational form, its width in basis points, and its illustrative dollar
equivalent per `$10 million` of recorded amount. It does **not** select a
threshold or convert the comparison into an audit conclusion.

## Methods and mathematical status

| Method | Status in these practical-size cases |
|---|---|
| ordinary binomial Stringer | exact reference calculation; general coverage unresolved |
| ordinary Poisson Stringer | exact reference calculation; general coverage unresolved |
| Gaffke bounded-mean endpoint | proved finite-sample valid under the paper's i.i.d. `[0,1]` model, using the cited external validity theorem |
| `max(Poisson Stringer, Gaffke)` | proved finite-sample valid under that model |
| full-scale calibrated Poisson, factorwise capped | proved valid modified rule under that model |
| zero-anchor calibrated Poisson, factorwise capped | proved valid modified rule under that model |

The two calibrated paths are different pre-specified procedures. The
benchmark compares them but does not certify selecting the smaller result
after observing the sample.

## Exact benchmark table

Zero taints are omitted from the profile label but included in `n`. Entries
are conservative upper enclosures, displayed as percentages. The exact
rationals are in
[`audit-decision-benchmark.json`](../computations/certificates/audit-decision-benchmark.json).

| fixed profile | confidence | ordinary Poisson reference | valid Gaffke | valid maximum | valid full-scale | valid zero-anchor |
|---|---:|---:|---:|---:|---:|---:|
| all zero, `n=25` | 95% | 11.9830% | 11.2929% | 11.9830% | 13.4958% | 11.9830% |
| taints `1,.4,.1`, `n=25` | 95% | 22.0417% | 19.6049% | 22.0417% | 24.8244% | 27.1874% |
| all zero, `n=100` | 95% | 2.9958% | 2.9514% | 2.9958% | 3.7026% | 2.9958% |
| taint `.1`, `n=100` | 95% | 3.1706% | 3.0536% | 3.1706% | 3.9187% | 3.4097% |
| taint `1`, `n=100` | 95% | 4.7439% | 4.6560% | 4.7439% | 5.8633% | 7.1346% |
| taints `1,.4,.1`, `n=100` | 95% | 5.5105% | 5.2159% | 5.5105% | 6.8107% | 8.9495% |
| taints `.2,.1,.05,.02,.01,.005`, `n=100` | 95% | 3.6217% | 3.3533% | 3.6217% | 4.4763% | 4.4778% |
| taints `1,.4,.1`, `n=200` | 95% | 2.7553% | 2.6352% | 2.7553% | 3.4660% | 4.9912% |
| taints `1,.4,.1`, `n=100` | 90% | 4.5987% | 4.3913% | 4.5987% | 5.9147% | 8.6502% |
| taints `1,.4,.1`, `n=100` | 99% | 7.5094% | 7.0129% | 7.5094% | 8.7193% | 9.9463% |

## What the benchmark shows

1. **The valid alternatives differ at a decision-relevant scale.** For the
   `n=100`, 95%, three-taint profile, the full-scale calibrated endpoint is
   130.021214 basis points above ordinary Poisson Stringer; the zero-anchor
   endpoint is 343.898891 basis points above it. Those widths correspond to
   conservative ceilings of `$130,021.22` and `$343,898.90` per illustrative
   `$10 million` of recorded amount. A threshold inside either interval would
   separate the corresponding numerical comparisons.

2. **The anchor trades zero-error stability for larger error increments.** On
   every all-zero profile it equals ordinary Poisson Stringer exactly, whereas
   full scaling raises the result. With one or more taints the anchored path
   can be higher; for the `n=100`, 95%, three-taint profile it exceeds the
   full-scale path by 213.877677 basis points. For six small taints, the two
   calibrated results are nearly equal. Neither path uniformly dominates the
   other.

3. **The safeguard preserves the familiar calculation in all ten fixed
   cases.** Exact brackets prove that Poisson Stringer exceeds Gaffke in each
   case, so their valid maximum has zero uplift. This is case-specific
   evidence, not a theorem that ordinary Poisson Stringer dominates Gaffke at
   these sample sizes. The separately valid Gaffke endpoint is lower than
   ordinary Poisson Stringer in every listed case.

4. **Poisson-factor conservatism is visible.** The artifact also reports the
   exact binomial Stringer reference. On the all-zero `n=100`, 95% profile,
   the Poisson result is 4.442732 basis points above both the binomial and
   Gaffke endpoints. Reporting this distinction prevents the factor
   convention from being mistaken for evidence about the unresolved ordinary
   coverage claim.

## Reproduction and exactness boundary

Run serially:

```sh
make audit-benchmark-check
```

The generator
[`audit_decision_benchmark.py`](../computations/python/audit_decision_benchmark.py):

1. uses 96-bit dyadic Clopper--Pearson brackets with integer endpoint signs;
2. reads the exact 64-bit Poisson brackets and adjacent-dyadic calibration
   multipliers from the existing calibration certificate;
3. evaluates Stringer through a nonnegative summation-by-parts formula;
4. encloses each Gaffke quantile on a 64-bit dyadic grid and verifies both
   tail signs exactly;
5. regenerates every rational rate, threshold-switch interval, basis-point
   width, and currency scaling; and
6. compares the resulting JSON byte-for-byte with the committed artifact.

The run is small: it is single-process, takes about two seconds on the
development machine, and peaked below 85 MB of resident memory.

## Scope boundary

The validity labels quantify over independent `[0,1]`-valued observations
with a common mean. They do not establish a bridge to systematic PPS,
sampling without replacement, stratification, certainty selections,
negative taints, item-level caps, or a particular monetary projection. The
separate one-start systematic-PPS failure theorem and design-safe overlay are
in `../theory/SYSTEMATIC-PPS.md`; none of these i.i.d. benchmark columns is
that overlay. The
`$10 million` column is only the identity `rate × recorded amount`; it does
not validate a sampling frame, taint definition, tolerable-misstatement
threshold, or professional-standards conclusion.
