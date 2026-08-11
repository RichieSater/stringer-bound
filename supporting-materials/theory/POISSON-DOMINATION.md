# Poisson factors dominate binomial factors at practical confidence levels

## Statement

Fix a sample size (n\ge 1) and tail probability
(0<\alpha<e^{-1}). Let (p^{\mathrm B}_{n}(j)) be the one-sided
Clopper--Pearson upper factor, defined for (0\le j<n) by

\[
\Pr\{\operatorname{Bin}(n,p^{\mathrm B}_{n}(j))\le j\}=\alpha,
\qquad p^{\mathrm B}_{n}(n)=1.
\]

Let the Poisson audit factor be

\[
p^{\mathrm P}_{n}(j)=\frac{\lambda_j}{n},
\qquad
\Pr\{\operatorname{Pois}(\lambda_j)\le j\}=\alpha,
\qquad 0\le j\le n.
\]

Then

\[
p^{\mathrm P}_{n}(j)\ge p^{\mathrm B}_{n}(j)
\quad (0\le j\le n),
\]

and the Poisson-factor Stringer bound is at least the binomial-factor
Stringer bound for every sample. Equivalently, this holds at every
nominal confidence level greater than
(1-e^{-1}\approx 0.6321205588), including conventional audit levels
such as 90%, 95%, and 99%.

The Poisson factors here are the untruncated confidence factors used in
the Stringer formula. A factor can exceed one.

## Proof of coordinatewise domination

Write

\[
B(j;n,p)=\Pr\{\operatorname{Bin}(n,p)\le j\},\qquad
P(j;\lambda)=\Pr\{\operatorname{Pois}(\lambda)\le j\}.
\]

Anderson and Samuels (1967, Corollary 2.1 and Section 4.2) prove

\[
B(j;n,\lambda/n)<P(j;\lambda)
\quad\text{if }j<\lambda-1\text{ and }\lambda<n.
\]

Their Theorem 4.1 also implies that

\[
P(j;\lambda)<P(0;1)=e^{-1}
\quad\Longrightarrow\quad
\lambda>j+1.
\]

Apply these statements at (lambda=\lambda_j). Because
(P(j;\lambda_j)=\alpha<e^{-1}), we have
(lambda_j>j+1).

* If (lambda_j<n), then
  (B(j;n,\lambda_j/n)<P(j;\lambda_j)=\alpha). The binomial CDF is
  strictly decreasing in (p), so its root satisfies
  (p^{\mathrm B}_{n}(j)<\lambda_j/n=p^{\mathrm P}_{n}(j)).
* If (lambda_j\ge n), then
  (p^{\mathrm P}_{n}(j)\ge1\ge p^{\mathrm B}_{n}(j)).

This covers every (j=0,\ldots,n).

## From factors to the Stringer bound

Let

\[
d_j=p^{\mathrm P}_{n}(j)-p^{\mathrm B}_{n}(j)\ge0
\]

and arrange the observed taints as
(t_{(1)}\ge\cdots\ge t_{(n)}\). Subtracting the two Stringer bounds and
summation by parts give

\[
\begin{aligned}
\operatorname{SB}_{\mathrm P}-\operatorname{SB}_{\mathrm B}
={}&d_0(1-t_{(1)})
  +\sum_{j=1}^{n-1}d_j(t_{(j)}-t_{(j+1)})
  +d_n t_{(n)}\\
\ge{}&0.
\end{aligned}
\]

Every coefficient multiplying a (d_j) is nonnegative. Thus
coordinatewise factor domination implies pointwise domination of the
complete Stringer bound even though domination of successive factor
increments need not hold.

## Scope

This result replaces finite-range numerical comparisons at 90% and 95%
with an analytic statement for every sample size. It transfers any
proved binomial-factor coverage result to the Poisson-factor bound; in
particular, it transfers the paper's (n=2) theorem and its
one-nonzero-taint result throughout the stated confidence range.

It does **not** prove general-(n) conservatism, because pointwise
domination alone cannot supply the still-unproved binomial-factor
coverage guarantee for arbitrary distributions when (n\ge3).

A separate corrected simultaneous-band argument proves direct Poisson
coverage for every `n<=8` at 90%, every `n<=11` at 95%, and every `n<=20`
at 99% confidence; see
[`POISSON-SIMULTANEOUS-BAND.md`](POISSON-SIMULTANEOUS-BAND.md). That result
does not rely on transferring a binomial theorem and likewise does not claim
general-`n` coverage.

## Source

T. W. Anderson and S. M. Samuels, “Some inequalities among binomial and
Poisson probabilities,” *Proceedings of the Fifth Berkeley Symposium on
Mathematical Statistics and Probability*, Vol. 1, University of
California Press, 1967, pp. 1--12.

Public scan: <https://digicoll.lib.berkeley.edu/record/112999>
