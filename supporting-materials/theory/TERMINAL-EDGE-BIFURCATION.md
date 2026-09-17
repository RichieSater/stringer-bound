# The terminal Clopper--Pearson edge bifurcation

## Status and thesis

This note proves an all-sample-size limitation of the pointwise
Stringer--Gaffke comparison.  It does **not** prove undercoverage of the
Stringer bound.

For every sample size `n >= 2`, there is a unique nontrivial critical tail
level \(\alpha_n\).  Above it, the cap calibrated at the terminal step vertex
immediately rises along a two-knot edge, so binomial Stringer cannot
pointwise dominate the Gaffke bounded-mean limit on every ordered sample.
The critical levels converge to

\[
 \alpha_\infty=0.28466813704083846\ldots,
 \qquad \alpha_\infty(1-2\log\alpha_\infty)=1.                 \tag{1}
\]

Consequently, for every \(\alpha>\alpha_\infty\), pointwise domination fails
for all sufficiently large sample sizes.  The ordinary Stringer coverage
question is weaker and remains open at conventional confidence levels from
`n=8` onward.

## 1. Exact edge derivative

Let \(p_0=p_n(0)\) and \(p_1=p_n(1)\) be the one-sided
Clopper--Pearson upper limits.  In the ascending-knot cap formulation, the
last two threshold weights are

\[
 c_{n-1}=p_1-p_0,\qquad c_n=p_0.                               \tag{2}
\]

Consider

\[
 x(\varepsilon)=(0^{\,n-1},1-\varepsilon,1),
 \qquad 0\le\varepsilon<1.                                    \tag{3}
\]

For a uniform Dirichlet vector with `n+1` coordinates, put

\[
 S=D_{n-1}+D_n\sim\operatorname{Beta}(2,n-1),\qquad
 W=\frac{D_{n-1}}S\sim\operatorname{Unif}(0,1).                \tag{4}
\]

The two variables in (4) are independent.  The random score and its moving
threshold are

\[
 x(\varepsilon)\mathbin{\cdot}D=S(1-\varepsilon W),\qquad
 c\mathbin{\cdot}x(\varepsilon)=p_1-\varepsilon(p_1-p_0).
\]

Conditioning on \(W\), differentiating the beta survival probability at
zero, and using \(\mathbb EW=1/2\) gives

\[
 \boxed{
 \left.\frac{d}{d\varepsilon}
 \Pr\{x(\varepsilon)\mathbin{\cdot}D>
          c\mathbin{\cdot}x(\varepsilon)\}
 \right|_{\varepsilon=0}
 =f_{2,n-1}(p_1)\left(\frac{p_1}{2}-p_0\right),}              \tag{5}
\]

where \(f_{2,n-1}(p_1)>0\).  At \(\varepsilon=0\), the cap is the calibrated
terminal step cap and has probability \(\alpha\).  It therefore exceeds
\(\alpha\) for all sufficiently small positive \(\varepsilon\) exactly when

\[
 p_1>2p_0.                                                      \tag{6}
\]

By (2), condition (6) is exactly the reversal
\(c_{n-1}>c_n\) of the last two Stringer threshold weights.  The cap
criterion then gives an ordered sample on which \(\mathrm{SB}<G_\alpha\).
When the sign in (6) is reversed, (5) proves only local decay on this edge;
it does not establish global pointwise domination.

## 2. The unique finite-sample transition

Write

\[
 r=\alpha^{1/n},\qquad p_0=1-r.
\]

At the transition, \(p_1=2p_0\).  Substitution in the defining binomial tail
for \(p_1\) gives

\[
 \boxed{
 (2r-1)^{n-1}\{2n-1-2(n-1)r\}=r^n.}                          \tag{7}
\]

Besides the trivial endpoint \(r=1\), equation (7) has exactly one solution
\(r_n\) in

\[
 \frac12<r_n<1-\frac1{2n}.                                    \tag{8}
\]

Indeed, define

\[
 H_n(r)=\frac{(2r-1)^{n-1}\{2n-1-2(n-1)r\}}{r^n}.
]

On \((1/2,1)\), direct differentiation gives

\[
 \frac{d}{dr}\log H_n(r)
 =\frac{n(2nr-2n+1)}
 {r(2r-1)(2nr-2n-2r+1)}.                                     \tag{9}
\]

The final denominator factor is negative.  Hence \(H_n\) increases up to
\(1-1/(2n)\) and decreases strictly thereafter to \(H_n(1)=1\).  Since
\(H_n(1/2)=0\), the strict maximum is above one and (8) follows.

Let

\[
 \alpha_n=r_n^n.                                               \tag{10}
\]

The binomial lower tail is strictly decreasing in its success probability.
Thus \(H_n(r)>1\) is equivalent to \(p_1>2p_0\), proving the sharp
finite-`n` edge transition:

\[
 \begin{cases}
  \alpha<\alpha_n &: \text{the cap decreases initially on (3)},\\
  \alpha=\alpha_n &: \text{the first derivative vanishes},\\
  \alpha>\alpha_n &: \text{the cap rises above \(\alpha\) on (3)}.
 \end{cases}                                                   \tag{11}
\]

The first values are

| \(n\) | \(\alpha_n\) |
|---:|---:|
| 2 | \(0.36\) |
| 3 | \(0.3336852118672717\ldots\) |
| 4 | \(0.3214292869970077\ldots\) |
| 5 | \(0.3141689898050253\ldots\) |
| 8 | \(0.3032656125408426\ldots\) |
| 10 | \(0.2996033710921845\ldots\) |
| 30 | \(0.2897088973919907\ldots\) |

For `n=2`, \(r_2=3/5\).  For `n=3`,

\[
 r_3=\frac{19+\sqrt{21}}{34},
 \qquad
 \alpha_3=\left(\frac{19+\sqrt{21}}{34}\right)^3,             \tag{12}
\]

so the endpoint in the four-coordinate theorem is the first instance of the
same all-dimensional bifurcation.

## 3. Limiting critical level

Fix \(\alpha\in(0,1)\) and put \(r=\alpha^{1/n}\).  Then

\[
 H_n(\alpha^{1/n})
 \longrightarrow \alpha(1-2\log\alpha).                       \tag{13}
\]

The limiting function has, besides the endpoint \(\alpha=1\), one root in
\((0,1)\), namely (1).  The convergence in (13) is uniform on compact
subintervals around that root; together with the uniqueness in (8), this gives

\[
 \alpha_n\longrightarrow\alpha_\infty.                        \tag{14}
\]

In particular, if \(\alpha>\alpha_\infty\), the right side of (13) is
strictly greater than one (until the trivial root at one).  Hence (6) holds
for every sufficiently large `n`, and pointwise Stringer--Gaffke domination
fails in all those dimensions.

## 4. Exact verification boundary

The derivative calculation, uniqueness argument, and convergence proof are
ordinary mathematics.  The generator
[`terminal_edge_bifurcation.py`](../computations/python/terminal_edge_bifurcation.py)
checks the symbolic identities, brackets every \(\alpha_n\) for
`2 <= n <= 30` by exact rational polynomial signs, and brackets
\(\alpha_\infty\) by exact rational atanh-series bounds for the logarithm.
Its machine-readable output is
[`terminal-edge-bifurcation-certificate.json`](../computations/certificates/terminal-edge-bifurcation-certificate.json).
Run

```bash
make terminal-edge-bifurcation-check
```

to regenerate the certificate and compare it byte for byte.  The finite
table is reproducible numerical information; no finite table or floating
search is used to infer the all-sample-size theorem.

## Mathematical yield

The transition identifies the exact structural assumption behind the
four-coordinate range.  Nondecreasing threshold weights are not merely a
convenient hypothesis: the moment the terminal pair reverses, its calibrated
step cap loses local maximality.  The same mechanism persists in every
sample size and converges to the nontrivial constant (1).

This also separates two questions that can otherwise be conflated.  The
Gaffke comparison is a strong pointwise certificate for Stringer validity,
so its failure marks a limit of that proof mechanism.  Coverage averages over
samples from a common distribution and is strictly weaker; (5)--(14) do not
produce a distribution with Stringer undercoverage.

The reusable object is the likelihood-ratio curve \(H_n\), whose unimodality
turns a multivariate cap failure into a one-variable phase transition.  Any
future all-dimensional positive theorem based on monotone threshold weights
must live below the finite \(\alpha_n\) boundary, while a proof of Stringer
coverage beyond it must use information lost by pointwise Gaffke domination.
