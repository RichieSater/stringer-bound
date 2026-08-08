# Findings for practitioners: stress-testing the Stringer bound

*Plain-language summary of what this repository establishes about the
monetary-unit-sampling evaluation method used across the audit profession.
Every number below is backed by an exact-arithmetic certificate or a
written proof in this repository; see the README for the technical
statements.*

## 1. The guarantee auditors rely on has never been proven — and the one claimed proof is wrong

Since Stringer (1963), MUS evaluation has rested on the belief that the
Stringer bound's true coverage is at least its stated confidence level in
real (finite) samples. Sixty-three years on, there is **no proof** of this
for any realistic sample size. The only claimed finite-sample
certification (Bimpeh 2008: reliability up to n = 11 at 95%) contains an
index error in its confidence-band argument; this repository refutes it
with counterexamples checkable by hand and shows the method could never
have certified any sample size
([details](supporting-materials/audit/BIMPEH-GAP.md)). What is actually
proven today: n = 1, populations whose errors are all 100% overstatements,
and populations with a single distinct error size
([proof](supporting-materials/computations/python/two_point_lemma.py)).

## 2. At standard confidence (95%, 90%) the bound passes every stress test — with zero margin

Exhaustive grid-plus-optimizer searches over two- and three-error-size
populations (n up to 100), evaluated in exact arithmetic: the worst-case
coverage never falls below the stated confidence — but it comes
arbitrarily close. The worst populations are those dominated by 100%
overstatements (taints near 1), where coverage approaches exactly 0.95 at
95% confidence. **Practical reading: the binomial Stringer bound has no
finite-sample safety margin at all in the worst case, but no failure
either.** The conjecture, if true, is perfectly tight.

## 3. The AICPA (Poisson-factor) tables carry a real safety margin

The factors in professional MUS tables are Poisson-based, which are
strictly more conservative than the exact-binomial factors (dominance
checked rigorously per sample size). Worst-case coverage over
two-error-size populations:

| n | worst case at nominal 95% | margin | worst case at nominal 90% | margin |
|---|---|---|---|---|
| 10 | 0.9716 | +2.16% | 0.9270 | +2.70% |
| 25 | 0.9589 | +0.89% | 0.9182 | +1.82% |
| 50 | 0.9573 | +0.73% | 0.9090 | +0.90% |
| 75 | 0.9553 | +0.53% | 0.9145 | +1.45% |
| 100 | 0.9549 | +0.49% | 0.9093 | +0.93% |

(The binomial-factor version of the same searches bottoms out at exactly
the nominal level at every one of these n, at both 95% and 90%.) The
Poisson margin shrinks as n grows: the extra conservatism practitioners
have long observed in MUS is a finite-sample cushion that erodes with
larger samples.

## 4. At low assurance levels the bound provably fails

For evaluations run at low confidence — the regime of low-risk strata and
limited-assurance procedures — the Stringer bound's coverage falls
**below** its stated level, certified in exact rational arithmetic
(no simulation error):

| stated confidence | sample size | true worst-case coverage found | certificate |
|---|---|---|---|
| 30% | 50 | 29.888% | [certified-alpha0.7-n50.json](supporting-materials/computations/certificates/certified-alpha0.7-n50.json) |
| 32% | 100 | 30.956% | [certified-alpha0.68-n100.json](supporting-materials/computations/certificates/certified-alpha0.68-n100.json) |
| 35% | 200 | 33.668% | [certified-alpha0.65-n200.json](supporting-materials/computations/certificates/certified-alpha0.65-n200.json) |
| 35% | 400 | 32.630% | [certified-alpha0.65-n400.json](supporting-materials/computations/certificates/certified-alpha0.65-n400.json) |
| 37% | 200 | 36.079% | [certify-alpha0.63.log](supporting-materials/computations/certificates/certify-alpha0.63.log) (6 confirmed) |
| 37% | 400 | 35.508% | [certify-alpha0.63.log](supporting-materials/computations/certificates/certify-alpha0.63.log) (13 confirmed) |

The shortfall *grows* with sample size at a fixed low confidence level,
and the failure threshold moves toward 50% confidence as samples grow —
consistent with the proven asymptotic threshold at exactly 50%
(Pap–van Zuijlen 1996). At exactly 50% confidence (risk of incorrect
acceptance 50%, a floor some methodologies permit), no failure appears
through n = 300. **Practical reading: never rely on a Stringer evaluation
below 50% confidence; at 50% and above no failure has ever been
exhibited, but no proof protects you either.**

## 5. A first piece of the proof auditors never had

This repository now contains a complete proof that the Stringer bound
is conservative at **n = 2, for every possible error population and
every confidence level**
([N2-PROOF.md](supporting-materials/theory/N2-PROOF.md)) — the first
finite-sample proof at any n ≥ 2 in the method's 63-year history. Two
samples is not an audit, but the proof's machinery (a change of
variables that turns the bound into a weighted order-statistic
inequality, plus a rectangle-and-budget argument) is built to scale,
and it explains *why* the bound survives every stress test at 95% while
failing at low confidence: the failure is strictly a
large-sample-size phenomenon. The general case for audit-sized samples
remains open — and this repository's exact machinery is built so that
any resolution can be verified without trusting simulations.
