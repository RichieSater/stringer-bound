# A simultaneous-band theorem for the Poisson Stringer factors

## Result

Let the Poisson count limits `lambda_j` be defined by

\[
  e^{-\lambda_j}\sum_{k=0}^j\frac{\lambda_j^k}{k!}=\alpha,
  \qquad p_n^{\rm P}(j)=\lambda_j/n.                    \tag{1}
\]

For independent taints with any common distribution on `[0,1]`, the
Poisson-factor Stringer bound has coverage at least `1-alpha` throughout the
following ranges:

| nominal confidence | `alpha` | every sample size |
|---:|---:|---:|
| 90% | 0.10 | `1 <= n <= 8` |
| 95% | 0.05 | `1 <= n <= 11` |
| 99% | 0.01 | `1 <= n <= 20` |

This theorem is independent of the binomial-factor comparison. In
particular, it proves new Poisson cases beyond the pointwise-transfer range
`n<=6` at all three levels.

The proof has two parts.  A distribution-free lemma turns a corrected
simultaneous survival-band event into Stringer coverage.  An exact rational
certificate evaluates the probability of that event at the three stated
levels.

## 1. The sufficient event

For a sample `T_1,...,T_n`, put

\[
 G(x)=\Pr\{T>x\},\qquad
 N_x=\#\{r:T_r>x\}.
\]

Ignoring sample points, which form a Lebesgue-null set in the integral, the
Stringer identity is

\[
 \operatorname{SB}_{\rm P}
 =\int_0^1p_n^{\rm P}(N_x)\,dx,
 \qquad
 \mu=\int_0^1G(x)\,dx.                                  \tag{2}
\]

The argument must allow atoms.  Introduce independent auxiliary uniforms
`R_1,...,R_n` and set

\[
 V_r=\Pr\{T>T_r\}+R_r\Pr\{T=T_r\}.                       \tag{3}
\]

The randomized probability integral transform makes the `V_r` independent
`Unif(0,1)` variables for every distribution of `T`. Every observation above
`x` has `V_r<=G(x)`, while every observation at or below `x` has
`V_r>=G(x)`. Hence, if `N_x=j<n`,

\[
 G(x)\le V_{j+1:n}.                                      \tag{4}
\]

Define

\[
 a_i=\min\{1,p_n^{\rm P}(i-1)\},\qquad i=1,\ldots,n,
\]

and the corrected simultaneous event

\[
 \mathcal E_n=\{V_{i:n}\le a_i, i=1,\ldots,n\}.          \tag{5}
\]

On `E_n`, equations (4)--(5) give
`G(x)<=p_n^P(N_x)` whenever `N_x<n`.  When `N_x=n`, the same conclusion
follows from `G(x)<=1` and `p_n^P(n)>1` at the three levels considered here.
The latter inequality also follows from the established fact
`lambda_j>j+1` for `alpha<e^-1`.  Integrating proves

\[
 \mathcal E_n\subseteq\{\operatorname{SB}_{\rm P}\ge\mu\}.
                                                               \tag{6}
\]

Although `E_n` uses the auxiliary randomizers, the event on the right of
(6) depends only on the original taint sample.  Its probability in the
augmented product space is therefore its ordinary sampling probability.

The proof also applies if every factor is capped at one before integration:
on `E_n`, `G(x)<=min(1,p_n^P(N_x))`.  Capping only the final reported bound
at one is valid for the simpler reason that `mu<=1`.

## 2. Exact probability of the event

For nondecreasing boundaries `0<=a_1<=...<=a_n<=1`, write

\[
 Q_m=\Pr\{V_{i:m}\le a_i, i=1,\ldots,m\},\qquad Q_0=1.
\]

Bolshev's recursion gives

\[
 Q_m=1-\sum_{i=1}^m {m\choose i}
       (1-a_{m-i+1})^iQ_{m-i}.                              \tag{7}
\]

The certificate encloses every needed `lambda_j` between adjacent dyadic
rationals `L_j<U_j` with denominator `2^80`.  The endpoint signs

\[
 e^{-L_j}\sum_{k=0}^jL_j^k/k!>\alpha,
 \qquad
 e^{-U_j}\sum_{k=0}^jU_j^k/k!<\alpha                       \tag{8}
\]

are proved with rational arithmetic.  Specifically, the code divides a
rational argument by a power of two until it is in `[0,1]`, brackets the
exponential between consecutive odd and even partial sums of its alternating
series, and raises the positive bounds back to the required power.  It uses
48 term pairs.  Thus `L_j<lambda_j<U_j` without assuming a correctly rounded
library exponential.

Substituting the smaller rational boundaries

\[
 a_i^-=\min\{1,L_{i-1}/n\}                                  \tag{9}
\]

into (7) gives a rational lower bound for `Pr(E_n)`.  Every comparison below
is exact; the decimals are only summaries.

| `alpha` | largest certified `n` | lower bound for `Pr(E_n)` | margin over `1-alpha` |
|---:|---:|---:|---:|
| 0.10 | 8  | 0.9028501343687063 | 0.0028501343687063 |
| 0.05 | 11 | 0.9532937128380234 | 0.0032937128380234 |
| 0.01 | 20 | 0.9906297230587641 | 0.0006297230587641 |

The certificate separately checks every smaller positive sample size, not
only the endpoint displayed in the table.  Equations (6)--(9) therefore
prove the stated coverage ranges for arbitrary continuous or atomic taint
distributions.

## 3. Scope of the result

At the next sample size, even the rational upper boundaries based on `U_j`
give simultaneous-event probabilities below nominal:

| `alpha` | next `n` | upper bound for `Pr(E_n)` |
|---:|---:|---:|
| 0.10 | 9  | 0.8903656999283113 |
| 0.05 | 12 | 0.9482646529430472 |
| 0.01 | 21 | 0.9899786983569257 |

This does **not** show Stringer undercoverage at those sample sizes.  The
event in (5) is sufficient, not necessary; the integral comparison in (2)
can hold when the pointwise band crosses `G`.  The table identifies the exact
limit of this particular proof route at the three stated levels.

The result assumes independent sampling from a common taint distribution,
as does the manuscript's model.  It does not by itself cover systematic PPS
sampling, sampling without replacement, negative taints, or a change in
professional auditing standards.

## 4. Reproduction

The source is
[`poisson_band_certificate.py`](../computations/python/poisson_band_certificate.py),
and the committed output is
[`poisson-simultaneous-band-certificate.json`](../computations/certificates/poisson-simultaneous-band-certificate.json).
Regenerate and compare it byte for byte with

```sh
make poisson-band-certificate-check
```

The relevant unit tests independently check special cases of Bolshev's
recursion, the rational exponential enclosure, every dyadic endpoint sign,
the three certified frontiers, and the next-size limitations.
