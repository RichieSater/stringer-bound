# Ordered knot collisions and the active-prefix inequality

## Status

The four-coordinate monotone cap theorem is proved. The proposed
five- and six-coordinate global monotone cap theorems remain conjectural:
the available first-prefix inequalities do not establish the active-prefix
inequality on every repeated-lowest-knot face. This is a missing implication,
not a counterexample to either conjecture. The polynomial identities and
coefficient inequalities retain their stated algebraic meanings.

## The necessary local statement

Let `D` be uniform on the simplex with `N` coordinates, let
`x_0<=...<=x_(N-1)`, and, at a regular section `x dot D=q`, put

\[
 p=\Pr\{x\mathbin{\cdot}D>q\},\qquad
 \gamma_i=\mathbb E(D_i\mid x\mathbin{\cdot}D=q),\qquad
 M_r=\sum_{i<r}\gamma_i.
\]

The stationary-face argument requires

\[
 p\le\max_{r:\,x_{r-1}<x_r}I_{M_r}(r,N-r).                 \tag{1}
\]

At a stationary point on a collision face, each centroid coordinate is the
average of the threshold weights `c_i` in its equal-knot block. Thus
`M_r=C_r=sum_(i<r)c_i` at block boundaries, but not necessarily inside a
block. Replacing the maximum in (1) by a maximum over all prefixes weakens
the statement and does not suffice for the global theorem.

The mixed five- and six-coordinate calculations compare `p` with
`I_(gamma_0)(1,N-1)`. This is an active prefix when `x_0<x_1`. If
`x_0=x_1`, the same inequality, even if valid by a closed polynomial
identity, need not be an inequality of the form (1).

## An exact illustration of the distinction

Take

\[
 c=(1/100,9/100,3/20,1/4,1/2),\qquad
 \bar c=(1/20,1/20,3/20,1/4,1/2).
\]

The second vector averages the first two entries of the first. Both vectors
are nondecreasing probability vectors. Their prefix sums agree for `r>=2`,
but the four beta caps for `c` are

\[
 \frac{3940399}{10^8},\quad\frac{523}{10000},\quad
 \frac{13}{256},\quad\frac1{16}.
\]

Their maximum is `1/16`, whereas

\[
 I_{\bar c_0}(1,4)=1-(19/20)^4
 =\frac{29679}{160000}>\frac1{16}.
\]

This example disproves the proposed inference from a first-prefix bound
after block averaging to the desired original-prefix bound. It does not
assert that these vectors occur at a particular section, and is not a
counterexample to a cap theorem.

## Collision limits of polynomial constraints

Another boundary distinction must be retained when a centroid difference is
factored as `a(x) G(x)/Delta(x)`, with `Delta>0` and `a>0` on an open
stratum. Orderedness implies `G>=0` there. On a collision face where
`a=0`, the difference vanishes regardless of the sign of `G`. Therefore a
certificate assuming `G>=0` extends to limits satisfying that constraint,
not automatically to every ordered point of the collision face. A density
argument or a separate face calculation is needed for the latter claim.

## What remains to be proved

For `N=5,6`, establish (1) on repeated-lowest-knot faces, with the actual
block multiplicities, or prove a separate maximum principle excluding those
faces as obstructions. Any use of continuity must also preserve the
centroid-order hypotheses after vanishing factors are removed. Until then,
the global cap bounds and the confidence continua inferred from them are
conditional. The `n=4,5` weight-monotonicity calculations and the
terminal-edge obstructions remain independent results; the separate
fixed-confidence certificates do not use this missing implication.

For `N=4`, the mixed-section certificate instead uses `I_(M_2)(2,2)`.
The boundary after the second knot stays strict throughout a two-versus-two
regular crossing, including collisions within either pair. The repeated-top
case has its own polynomial identity, and the endpoint crossings use their
active first or last prefix. Thus this active-prefix issue does not remove
the four-coordinate theorem or its `n=3` application.
