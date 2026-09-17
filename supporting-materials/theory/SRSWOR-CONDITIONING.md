# A finite-population conditioning bridge for simple random sampling

## Status and scope

This note proves an exact bridge from any with-replacement, distribution-free
upper confidence bound to **simple random sampling without replacement**
(SRSWOR) from a finite population.  It gives a valid, explicit tail-level
adjustment for Stringer wherever the corresponding i.i.d. Stringer theorem is
already proved, and it applies at every sample size to the valid Gaffke and
calibrated-Poisson procedures in this repository.

The result does **not** cover systematic PPS, successive PPS, stratification,
certainty selections, or a design whose ordered sample is not uniform over
distinct population indices.  It validates an adjusted procedure, not the
ordinary unadjusted Stringer calculation under an arbitrary audit design.

## 1. Exact conditioning theorem

Fix a finite population of `N` taints

\[
 t_1,\ldots,t_N\in[0,1],\qquad
 \mu_N=\frac1N\sum_{i=1}^N t_i,
\]

and let `1<=n<=N`.  Write

\[
 c_{N,n}:=\frac{(N)_n}{N^n}
 =\prod_{k=0}^{n-1}\left(1-\frac{k}{N}\right),       \tag{1}
\]

the probability that `n` independent uniform draws of population indices
are all distinct.

Let `U_eta(x_1,...,x_n)` be any measurable upper bound with the following
with-replacement guarantee at tail probability `eta`:

\[
 \Pr_{X_1,\ldots,X_n\stackrel{\mathrm{iid}}\sim F}
 \{U_\eta(X_1,\ldots,X_n)<\mu(F)\}\le\eta           \tag{2}
\]

for every distribution `F` on `[0,1]`.  No independence claim is made about
the sample in the next theorem.

**Conditioning theorem.**  If `(T_1,...,T_n)` is an ordered SRSWOR sample
from the finite population, then

\[
 \Pr_{\rm WOR}\{U_\eta(T_1,\ldots,T_n)<\mu_N\}
 \le \min\left\{1,\frac{\eta}{c_{N,n}}\right\}.     \tag{3}
\]

Consequently, for a desired tail probability `alpha`, every value
`eta<=alpha*c_(N,n)` for which (2) is proved gives

\[
 \Pr_{\rm WOR}\{U_\eta(T_1,\ldots,T_n)\ge\mu_N\}
 \ge 1-\alpha.                                      \tag{4}
\]

**Proof.**  Let `I_1,...,I_n` be independent uniform indices in
`{1,...,N}` and set `X_j=t_(I_j)`.  Their common distribution is the
empirical population distribution, whose mean is exactly `mu_N`.  Let `A`
be the event that all indices are distinct.  Then `P(A)=c_(N,n)`, and the
conditional law of `(I_1,...,I_n)` given `A` is uniform over all ordered
distinct index tuples; this is precisely ordered SRSWOR.  For the failure
event `E_eta={U_eta(X_1,...,X_n)<mu_N}`, (2) gives `P(E_eta)<=eta`, so

\[
 \Pr_{\rm WOR}(E_\eta)
 =\frac{\Pr(E_\eta\cap A)}{\Pr(A)}
 \le\frac{\Pr(E_\eta)}{c_{N,n}}
 \le\frac{\eta}{c_{N,n}}.
\]

Taking the trivial upper bound one proves (3), and
`eta<=alpha*c_(N,n)` proves (4).  \(\square\)

The argument is deliberately elementary: it uses neither negative
dependence nor a comparison between without- and with-replacement tail
events.  Those comparisons need not hold for a general nonlinear statistic.

## 2. Consequences for the bounds in this repository

Write `SB_eta` for Stringer computed with tail probability `eta`.

* At `n=2`, the i.i.d. theorem holds for every `eta in (0,1)`.  Thus, for
  every `N>=2`,

  \[
  \mathrm{SB}_{\alpha(1-1/N)}                         \tag{5}
  \]

  covers the finite-population mean under SRSWOR with probability at least
  `1-alpha`.

* At `n=3`, the i.i.d. theorem holds throughout
  `0<eta<=alpha_*`, where
  `alpha_*=((19+sqrt(21))/34)^3=0.3336852118672717...`. Therefore

  \[
  \mathrm{SB}_{\alpha(1-1/N)(1-2/N)}                 \tag{6}
  \]

  is SRSWOR-valid whenever the adjusted tail is at most `alpha_*`. In
  particular, (6) works for every `N>=3` whenever the target tail
  `alpha<=alpha_*`, including nominal 90%, 95%, and 99% confidence.

* At `n=4,5,6,7`, the independent direct i.i.d. certificates establish
  the fixed tails `0.10`, `0.05`, and `0.01`. A proved tail
  `eta_0` gives a target SRSWOR tail `alpha` whenever
  `eta_0<=alpha*c_(N,n)`. The proposed `n=4,5` continuous ranges do
  not currently justify setting `eta=alpha*c_(N,n)` at arbitrary adjusted
  tails; see
  [`ORDERED-COLLISION-OBLIGATION.md`](ORDERED-COLLISION-OBLIGATION.md).
  The exact smallest population sizes for selected conventional transfers
  are shown below:

  | `n` | use 95% i.i.d. rule for 90% SRSWOR (`c>=1/2`) | use 99% i.i.d. rule for 90% SRSWOR (`c>=1/10`) | use 99% i.i.d. rule for 95% SRSWOR (`c>=1/5`) |
  |---:|---:|---:|---:|
  | 5 | 17 | 7 | 8 |
  | 6 | 24 | 9 | 12 |
  | 7 | 33 | 12 | 16 |

  These are adjusted procedures: for example, the third column computes the
  99% i.i.d. Stringer rule to guarantee 95% under SRSWOR.  It is not a claim
  that ordinary 95% Stringer is design-valid.

* Any upper bound that satisfies (2) for every `eta` transfers at every
  `N>=n` by setting `eta=alpha*c_(N,n)`.  This includes the finite-sample
  valid Gaffke endpoint and the all-sample-size calibrated Poisson rules
  proved in this repository.  The pre-specified maximum of an arbitrary
  Stringer calculation and the adjusted Gaffke endpoint is therefore also
  SRSWOR-valid.

The ordinary, unadjusted rule at a proved i.i.d. tail `alpha` inherits only
the weaker bound

\[
 \Pr_{\rm WOR}\{U_\alpha\ge\mu_N\}
 \ge 1-\min\{1,\alpha/c_{N,n}\}.                     \tag{7}
\]

Equation (7) is useful as a diagnostic but is not nominal coverage unless
`c_(N,n)=1`.

## 3. Interpretation and reproducibility

The correction is small when collisions would be rare under hypothetical
with-replacement sampling and can be substantial when `n` is not small
relative to `N`.  It is transparent rather than claimed optimal: conditioning
discards all information about how the failure event interacts with repeated
indices.  Direct without-replacement concentration or betting methods may be
less conservative.

`computations/python/srswor_conditioning.py` regenerates exact threshold
sizes and illustrative collision factors in
`computations/certificates/srswor-conditioning-certificate.json`.  The
certificate corroborates the arithmetic; the proof above establishes the
theorem.
