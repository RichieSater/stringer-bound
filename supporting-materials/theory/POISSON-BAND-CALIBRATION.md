# An all-sample-size scalar calibration of Poisson Stringer

> **Status.** This note proves a distribution-free all-sample-size upper
> confidence bound obtained by multiplying every ordinary Poisson Stringer
> factor by one precomputed scalar.  It does not prove that ordinary Stringer
> is conservative outside the ranges already established in the manuscript.

The ordinary Poisson factors are attractive in audit practice because they
are tabulated and because the resulting bound has the familiar Stringer form.
The corrected simultaneous-band argument gives a direct way to retain that
form while obtaining a finite-sample guarantee for every sample size.

## 1. Definition

Fix `n>=1` and a tail probability `0<alpha<1`. Let

\[
 e^{-\lambda_j}\sum_{r=0}^j\frac{\lambda_j^r}{r!}=\alpha,
 \qquad j=0,\ldots,n,                                  \tag{1}
\]

and let `SB_P` be ordinary Poisson Stringer, whose factors are
`lambda_j/n`. Put

\[
 \kappa_0=\max\{1,n/\lambda_n\}.
\]

For `kappa>=kappa_0`, define

\[
 a_i(\kappa)=\min\!\left\{1,
       \frac{\kappa\lambda_{i-1}}n\right\},
 \qquad i=1,\ldots,n,                                  \tag{2}
\]

and

\[
 Q_{n,\alpha}(\kappa)
 =\Pr\{V_{i:n}\le a_i(\kappa),\ i=1,\ldots,n\},       \tag{3}
\]

where `V_(1:n)<=...<=V_(n:n)` are uniform order statistics.  Put

\[
 \kappa_{n,\alpha}
 =\inf\{\kappa\ge\kappa_0:
        Q_{n,\alpha}(\kappa)\ge1-\alpha\}.              \tag{4}
\]

The function in (3) is continuous and nondecreasing. Section 3 below
shows that the set in (4) is nonempty. As the inverse image of
`[1-alpha,1]` under a continuous function, it is closed, so its infimum is
attained.

## 2. Coverage theorem

> **Theorem (band-calibrated Poisson Stringer).**  Let
> `T_1,...,T_n` be independent with an arbitrary common distribution on
> `[0,1]` and mean `mu`.  For every `0<alpha<1`,
>
> \[
>  U_{n,\alpha}
>  =\min\{1,\kappa_{n,\alpha}\operatorname{SB}_{\rm P}\}
>                                                               \tag{5}
> \]
>
> satisfies
>
> \[
>  \Pr\{U_{n,\alpha}\ge\mu\}\ge1-\alpha.             \tag{6}
> \]
>
> More generally, (6) holds with any certified
> `kappa>=kappa_(n,alpha)`.  It also holds if each Poisson factor is
> capped at one before the complete factor-capped Stringer expression is
> multiplied by `kappa`.

**Proof.**  Multiply every Poisson factor in the Stringer formula by
`kappa`.  Linearity gives

\[
 c_0+\sum_{j=1}^n(c_j-c_{j-1})T_{(j)}
 =\kappa\operatorname{SB}_{\rm P},
 \qquad c_j=\frac{\kappa\lambda_j}{n}.                 \tag{7}
\]

The factors are nondecreasing, and `c_n>=1` by the definition of
`kappa_0`. The randomized-survival-band lemma in
[`POISSON-SIMULTANEOUS-BAND.md`](POISSON-SIMULTANEOUS-BAND.md) therefore
gives

\[
 \Pr\{\kappa\operatorname{SB}_{\rm P}\ge\mu\}
 \ge Q_{n,\alpha}(\kappa).                              \tag{8}
\]

Use `kappa=kappa_(n,alpha)` in (8).  Capping the reported upper bound at
one preserves coverage because `mu<=1`.  The randomized probability
integral transform in the band lemma makes the argument valid for both
continuous and atomic taint distributions.

For the factor-capped convention, put
`pbar_j=min(1,lambda_j/n)` and use `c_j=kappa*pbar_j`.  Because `kappa>=1`,

\[
 \min\{1,c_{i-1}\}
 =\min\{1,\kappa\lambda_{i-1}/n\}=a_i(\kappa),
\]

and the definition of `kappa_0` gives `c_n>=1`.  The same band event and
the same proof therefore apply.  `square`

This is an all-sample-size theorem for the **calibrated rule** in (5), not
for ordinary Poisson Stringer.  If `kappa_0=1` and the unadjusted event already satisfies
`Q_(n,alpha)(1)>=1-alpha`, then `kappa_(n,alpha)=1` and the rule is exactly
ordinary Stringer.  The direct theorem in the manuscript establishes this
through `n=8`, `n=11`, and `n=20` at 90%, 95%, and 99% confidence,
respectively.

## 3. Analytic finite upper bounds

The calibration exists without relying on a numerical root. Because the
`lambda_j` are increasing, the elementary choice

\[
 \kappa=\max\{1,n/\lambda_0\}
\]

belongs to the permitted domain, makes every boundary in (2) equal to one,
and hence makes `Q_(n,alpha)(kappa)=1`.

At conventional levels, an explicit bound follows from a Bonferroni
construction. Let

\[
 \theta=\frac\alpha n,
 \qquad
 e^{-\lambda_j(\theta)}
 \sum_{r=0}^j\frac{\lambda_j(\theta)^r}{r!}=\theta,
 \qquad
 \kappa_{\mathrm B}
 =\max\left\{\kappa_0,
       \max_{0\le j<n}\frac{\lambda_j(\theta)}{\lambda_j}\right\}.
                                                               \tag{9}
\]

Then

\[
 \kappa_{n,\alpha}\le\kappa_{\mathrm B}<\infty.         \tag{10}
\]

Indeed, for

\[
 b_i=\min\{1,\lambda_{i-1}(\theta)/n\},
\]

the Anderson--Samuels binomial--Poisson comparison used in
[`POISSON-DOMINATION.md`](POISSON-DOMINATION.md) gives

\[
 \Pr\{V_{i:n}>b_i\}
 =\Pr\{\operatorname{Bin}(n,b_i)\le i-1\}
 \le\theta.                                             \tag{11}
\]

If `lambda_(i-1)(theta)>=n`, the probability on the left is zero; otherwise the
comparison is strict.  The union bound gives

\[
 \Pr\{V_{i:n}\le b_i\text{ for every }i\}
 \ge1-n\theta=1-\alpha.                                 \tag{12}
\]

By (9), `a_i(kappa_B)>=b_i`, so (10) follows.  Thus the elementary choice
`kappa_B` is already a proved all-`n` calibration.  Computing the exact joint
probability in (3), rather than applying the union bound, gives the no-larger
calibration in (4).

## 4. Exact representative calibrations

Bolshev's recursion evaluates (3): with `Q_0=1`,

\[
 Q_m=1-\sum_{i=1}^m {m\choose i}
       (1-a_{m-i+1})^iQ_{m-i}.                          \tag{13}
\]

The committed certificate encloses each `lambda_j` between adjacent
64-bit dyadic rationals and brackets `kappa_(n,alpha)` between adjacent
28-bit dyadic rationals.  At the lower `kappa` endpoint it substitutes the
**upper** lambda endpoints into (13), obtaining an upper bound below
`1-alpha`.  At the upper `kappa` endpoint it substitutes the **lower**
lambda endpoints, obtaining a lower bound at least `1-alpha`.  All
probabilities in these two sign checks are exact rational numbers. The JSON
stores every Poisson-limit endpoint used, together with the endpoint-sign
direction, so the event calculation is self-contained rather than only a
table of final decimals.

The following entries round the certified dyadic upper endpoints upward to
12 decimal places, so every displayed multiplier is itself a valid choice.
The JSON records the same value as `valid_decimal_ceiling_12` alongside the
exact numerator and denominator.

| nominal confidence | `n` | valid `kappa` choice | scalar uplift |
|---:|---:|---:|---:|
| 90% | 25  | 1.189305365086 | 18.9306% |
| 90% | 50  | 1.251207165421 | 25.1207% |
| 90% | 100 | 1.286175295711 | 28.6176% |
| 90% | 200 | 1.305048119277 | 30.5049% |
| 95% | 25  | 1.126245908440 | 12.6246% |
| 95% | 50  | 1.195803422481 | 19.5804% |
| 95% | 100 | 1.235955074430 | 23.5956% |
| 95% | 200 | 1.257978815586 | 25.7979% |
| 99% | 25  | 1.027949459851 | 2.7950% |
| 99% | 50  | 1.111272465438 | 11.1273% |
| 99% | 100 | 1.161121416837 | 16.1122% |
| 99% | 200 | 1.189238607884 | 18.9239% |

The multiplier applies to the complete ordinary Poisson Stringer result,
not separately to the observed taints.  For example, the exact certificate
shows that multiplying the ordinary 95% Poisson result by
`1.235955074430` (and then capping at one) is a valid choice at `n=100`
under the manuscript's i.i.d. model.

These are the smallest multipliers under the sufficient-event criterion
(3), to the certified dyadic resolution.  They are **not** asserted to be
the smallest possible corrections among all valid confidence procedures,
and a value above one is not evidence that ordinary Stringer undercovers.
The displayed multiplier does not by itself control rounding of an external
Poisson factor table; a report must also use exact factors or conservative
upper enclosures, as the command in Section 6 does.

## 5. Relation to prior adjusted-level proposals

Bimpeh (2008, Section 5.4) explored an extended Stringer calculation using
Rom's adjusted significance levels and reported numerical boundary-crossing
probabilities through `n=20`.  The construction here is conceptually
related but logically separate.  It uses the corrected right-side
order-statistic constraint, evaluates its actual joint probability, and
scales one existing Poisson factor curve.  Its coverage proof does not
apply Rom's independent-multiple-testing result to the dependent
order-statistic constraints.

The scalar rule is also separate from the Stringer--Gaffke safeguard.  Both
are valid for every sample size under the model.  The Gaffke rule raises
Stringer only on samples where a separately valid bounded-mean limit is
larger; the scalar rule preserves the traditional factor-table form but can
be more conservative.  No pointwise ordering between the two rules is
claimed here.

## 6. Reproduction and scope

Regenerate the exact table and compare it byte for byte with the committed
artifact:

```sh
make poisson-band-calibration-check
```

For a specified sample, the same program computes a rational upper enclosure
of ordinary Poisson Stringer by using upper factor endpoints in its
nonnegative summation-by-parts form, then applies the certified scalar:

```sh
uv run --frozen python \
  supporting-materials/computations/python/poisson_band_calibration.py \
  --n 25 --alpha 0.05 --taints 1,0.4,0.1 \
  --out /tmp/calibrated-poisson.json
jq '{kappa_upper:.case.kappa_upper,report}' \
  /tmp/calibrated-poisson.json
```

Zero taints may be omitted while `--n` remains the full sample size. Decimal
taints are interpreted as exact rationals. The returned calibrated report is
an upper enclosure of the mathematical rule in (5), so factor rounding
cannot invalidate the guarantee. This command uses untruncated Poisson
factors and records that convention in its JSON; the factor-capped variant
proved in Section 2 is not silently substituted.

The source is
[`poisson_band_calibration.py`](../computations/python/poisson_band_calibration.py),
and the artifact is
[`poisson-band-calibration-certificate.json`](../computations/certificates/poisson-band-calibration-certificate.json).

The theorem assumes independent observations from a common `[0,1]` taint
distribution.  It does not by itself cover systematic PPS sampling,
sampling without replacement, negative taints, or the professional-judgment
steps that convert a taint bound into an audit conclusion.  A methodology
owner would also need to assess the efficiency of the scalar correction
against the Gaffke safeguard and other validated procedures before practice
adoption.

## References

T. W. Anderson and S. M. Samuels, “Some inequalities among binomial and
Poisson probabilities,” *Proceedings of the Fifth Berkeley Symposium on
Mathematical Statistics and Probability*, vol. 1, 1967, pp. 1--12.

Y. Bimpeh, *Statistical Modelling and Inference for Financial Auditing*,
PhD thesis, Dublin City University, 2008,
<https://doras.dcu.ie/600/>.

D. M. Rom, “A sequentially rejective test procedure based on a modified
Bonferroni inequality,” *Biometrika* 77 (1990), 663--665,
<https://doi.org/10.1093/biomet/77.3.663>.
