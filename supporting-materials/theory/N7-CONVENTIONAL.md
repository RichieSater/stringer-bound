# Certified finite-sample conservatism at `n = 7` and conventional confidence levels

This note records the complete reduction and computer-assisted proof boundary
for the manuscript's `n = 7` result. It proves that the binomial-factor
Stringer bound is conservative for every distribution on `[0,1]` when

```text
n = 7,       alpha in {0.10, 0.05, 0.01},
             nominal confidence in {90%, 95%, 99%}.
```

The proof compares Stringer pointwise with the valid one-sided Gaffke
bounded-mean upper limit. Vlassis and Thomas proved finite-sample validity of
Gaffke's test; its inversion for observations in `[0,1]` is the
Dirichlet-average limit described by Learned-Miller--Thomas and Ming et al.
A pointwise no-smaller Stringer limit inherits that distribution-free
coverage guarantee.

Primary sources:

- N. Vlassis and P. S. Thomas, *An Exact Distribution-Free Test for Means of
  Nonnegative Random Variables*, arXiv:2607.08415 (2026).
- J. Ming, A. Ramdas, Y. Shen, R. Wang, and I. Waudby-Smith, *Gaffke's
  confidence interval for the mean of bounded data is inadmissible but
  asymptotically efficient*, arXiv:2607.18661 (2026).

## 1. Reduction to a seven-simplex cap

Write `p_j = p_7(j)` and define the ascending-knot weights

```text
a = 1-p_6,
b = p_6-p_5,
c = p_5-p_4,
d = p_4-p_3,
e = p_3-p_2,
f = p_2-p_1,
g = p_1-p_0,
h = p_0.
```

For ascending sample taints `t1 <= ... <= t7`,

```text
SB = a t1 + b t2 + c t3 + d t4 + e t5 + f t6 + g t7 + h.
```

The eight weights are positive and sum to one. If the smallest taint is one,
both endpoints equal one. Otherwise, affine equivariance reduces the
comparison to

```text
(0,x,y,z,w,u,v,1),       0 <= x <= y <= z <= w <= u <= v <= 1,
s = b x + c y + d z + e w + f u + g v + h.                (1)
```

For `(D0,...,D7) ~ Dirichlet(1,...,1)`, it is enough to prove

```text
V(x,y,z,w,u,v)
  := P_D(xD1+yD2+zD3+wD4+uD5+vD6+D7 > s)
  <= alpha.                                                  (2)
```

For distinct knots `t0=0,t1=x,...,t6=v,t7=1`, the uniform-simplex cap
formula is

```text
V = sum_{i: ti>s} (ti-s)^7 / product_{j!=i}(ti-tj),          (3)
```

with coincident-knot values defined by continuity. It has seven polynomial
pieces, according as `s` lies above `v`, between two successive interior
knots, or below `x`.

Put

```text
q0=h, q1=g+h, q2=f+g+h, q3=e+f+g+h,
q4=d+e+f+g+h, q5=c+d+e+f+g+h.
```

The Clopper--Pearson equations give

```text
alpha = sum_{j=0}^k C(7,j) qk^j (1-qk)^(7-j),  k=0,...,5,
1-alpha = (b+c+d+e+f+g+h)^7.                                (4)
```

## 2. Exact residual factorizations and reflection

In the first six regions, multiply `alpha-V` by the positive common product
of the denominators appearing in (3). In the final region use

```text
s^7 - (b+c+d+e+f+g+h)^7 x y z w u v.                        (5)
```

After substituting the appropriate identity in (4), exact polynomial
division extracts all indicated upper-endpoint and coincident-knot factors.
The remaining residual degrees are

| region | threshold range | residual degree | six-simplices |
|---|---|---:|---:|
| A | `s >= v` | 7 | 1 |
| B | `u <= s <= v` | 12 | 6 |
| C | `w <= s <= u` | 15 | 15 |
| D | `z <= s <= w` | 16 | 20 |
| E | `y <= s <= z` | 15 | 15 |
| F | `x <= s <= y` | 12 | 6 |
| G | `s <= x` | 7 | 1 |

Only the four source residuals A--D need to be stored. Reversing the eight
Dirichlet spacings sends

```text
(a,b,c,d,e,f,g,h) -> (h,g,f,e,d,c,b,a),
(x,y,z,w,u,v) -> (1-v,1-u,1-w,1-z,1-y,1-x),
V -> 1-V,       alpha -> 1-alpha.
```

Consequently, the exact residual identities are

```text
P_G = -R(P_A),       P_F = -R(P_B),       P_E = -R(P_C).    (6)
```

The source derivation is performed in Singular over
`QQ(b,c,d,e,f,g,h)`. Every requested division is followed by an exact
remainder check.

## 3. The fixed 64-simplex chain

Let `V_i`, for `i=0,...,6`, have `6-i` leading zeros and `i` trailing ones.
These are the vertices of the ordered-knot simplex. For `r=1,...,6`, write
`x_r=x,y,z,w,u,v` and
`L_r=s-x_r`. Whenever an edge `V_i V_j` crosses `L_r=0`, define

```text
h_{r,ij} = V_i + tau (V_j-V_i),
tau = L_r(V_i) / [L_r(V_i)-L_r(V_j)].                       (7)
```

Successive slicing by `s=v,u,w,z,y,x`, followed by a fixed triangulation of
each slice, gives `1,6,15,20,15,6,1` six-simplices. The complete ordered
vertex lists and exact rational-function coordinates are stored in
[`n7-gaffke-bernstein-structure.json.gz`](../computations/certificates/n7-gaffke-bernstein-structure.json.gz).

The rigorous numerical layer encloses every affine determinant away from
zero at each certified level, checks every exact vertex-boundary identity,
and proves every other domain and cap-region inequality strictly. It also
certifies all 378 pairs of the 28 named vertices as distinct. The integer
relative-chain check pairs 112 internal five-facets with opposite
orientations. The remaining 224 five-facets lie on the seven outer
ordered-simplex hyperplanes, 32 on each. The sum of the 64 positively
oriented normalized determinant balls contains one and lies strictly
between `1/2` and `3/2`. Hence the relative-chain degree is the integer one.
Because every simplex lies in the ordered domain, the chain covers that
domain without interior multiplicity.

## 4. Generic structural zeros

After mapping each six-simplex affinely to the standard six-simplex, a
polynomial is nonnegative if all coefficients in its fixed-degree Bernstein
expansion are nonnegative. The exact structural-zero totals are

| region | degree | six-simplices | structural zeros |
|---|---:|---:|---:|
| A | 7 | 1 | 1 |
| B | 12 | 6 | 10,462 |
| C | 15 | 15 | 99,507 |
| D | 16 | 20 | 193,352 |
| E | 15 | 15 | 99,507 |
| F | 12 | 6 | 10,462 |
| G | 7 | 1 | 1 |

These are exact identities, not tolerance classifications. For each source
region, the derivation records the exact polynomial identity `N=F P`, where
`N` is the compact cleared cap residual, `F` is the extracted product of
linear factors, and `P` is the stored residual. For an affine face ideal `I`,
an affine coordinate change identifies the polynomial ring with
`K[t,r]` and `I=(r1,...,rc)`. Its associated graded ring is again a polynomial
ring and therefore a domain, so the `I`-adic order is additive on nonzero
products.

Twenty-two of the 26 generic face conditions use only threshold equations
`s=t_i`. For each one, exact restriction checks determine which linear
factors in `N` and `F` belong to `I`. Every factored summand of `N` has order
at least `ord_I(F)+q`, where `q` is the required residual order. Hence
`ord_I(P)>=q`, proving `P in I^q`. The artifact records every term order and
the extracted-factor order. The four remaining conditions are outer faces.
For those, an exact affine change of variables turns the independent
generators into normal coordinates, and sparse multivariate Horner evaluation
in the quotient by all normal monomials of total degree `q` gives a zero
remainder in Singular over `QQ(b,c,d,e,f,g,h)`.

Reflection (6) links these 26 generic statements to all 322 simplex-face
conditions. At each certified weight box, rigorous nonzero normal minors
prove that every geometric face basis and every proof-source ideal basis
specializes with the required rank under the original or reflected weights,
as applicable. The certificate also encloses away from zero every
proof-source pivot used for a restriction or inverse chart, so no generic
proof is specialized across a pole.

## 5. Rigorous Bernstein signs

Every coefficient not covered by the exact face identities is enclosed as
follows:

1. each Clopper--Pearson root is bracketed by adjacent dyadic rationals on a
   `2^-512` grid;
2. each binomial-CDF sign at a bracket endpoint is evaluated with integer
   arithmetic;
3. every weight, vertex, determinant, residual, affine substitution, and
   Bernstein transformation is evaluated with outward-rounded Arb real-ball
   arithmetic at 768-bit precision through `python-flint`; and
4. every nonstructural Bernstein coefficient has a strictly positive lower
   endpoint.

The smallest positive lower endpoint in each region is:

| alpha | A | B | C | D | E | F | G |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 8.29e-10 | 6.54e-19 | 2.41e-25 | 1.13e-27 | 4.97e-25 | 9.07e-17 | 1.12e-3 |
| 0.05 | 2.21e-12 | 6.73e-23 | 5.89e-30 | 3.54e-32 | 9.75e-29 | 5.17e-19 | 7.01e-4 |
| 0.01 | 4.80e-18 | 1.11e-31 | 3.65e-40 | 3.21e-42 | 4.35e-37 | 3.41e-24 | 1.76e-4 |

The combined certificate records the exact rational endpoint represented by
each decimal display.

Thus every residual in (5) and its six predecessors is nonnegative, so (2)
holds. Pointwise domination of the Gaffke upper limit, whose finite-sample
validity is established in the cited work, proves
finite-sample conservatism of the binomial Stringer bound at `n=7` and 90%,
95%, and 99% nominal confidence. The all-`n` factor comparison in
[`POISSON-DOMINATION.md`](POISSON-DOMINATION.md) transfers the conclusion to
the uncapped Poisson-factor Stringer bound.

## 6. Reproduction and trust boundary

The proof has three independently checkable layers.

```bash
# Re-derive A--D, all vertices, the 64 simplices, and every face specification.
make n7-structure-data-check

# Prove all 26 generic ideal-power conditions.
make n7-factor-order-check
make n7-face-a-check n7-face-b-basic-check n7-face-b5-check
make n7-face-c-basic-check n7-face-c6-check n7-face-c8-check
make n7-face-d-basic-check n7-face-d6-check n7-face-d9-check

# Regenerate the three rigorous ball certificates.
make n7-certificate-001-check n7-certificate-005-check \
  n7-certificate-010-check
```

The source-data and four outer-face checks require Singular. The data command
derives the four source residuals independently rather than trusting a cached
expansion. The nine Make targets partition the 26 generic ideal-power checks.
Continuous integration verifies the 22 factor-order records together and
runs the four outer-face Singular checks separately; this changes scheduling,
not the proof. The three sign jobs each
regenerate a complete confidence-level record and compare it exactly with
the corresponding record in the combined certificate.

The combined sign certificate records the SHA-256 digest of the compressed
structure artifact, binding the numerical layer to the exact symbolic
layer. Conventional binary floating-point arithmetic is not used for a
factor-root, face, determinant, region, or Bernstein sign decision. Decimal
strings are only readable renderings of rational endpoints of rigorous Arb
balls.
