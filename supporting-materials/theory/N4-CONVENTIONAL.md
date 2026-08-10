# Certified finite-sample conservatism at `n = 4`

This note supplies the complete reduction and exact-computation details for
the manuscript's `n = 4` result. It proves the binomial-factor Stringer bound
conservative at:

| `alpha` | nominal confidence |
|---:|---:|
| `0.10` | 90% |
| `0.05` | 95% |
| `0.01` | 99% |

As in the `n = 3` theorem, the proof is a pointwise comparison with the valid
one-sided Gaffke bounded-mean upper limit. Vlassis and Thomas proved the
finite-sample validity of Gaffke's test; inversion for observations in
`[0,1]` gives the Dirichlet-average upper limit described by
Learned-Miller--Thomas and Ming et al. Therefore a pointwise no-smaller
Stringer limit inherits its distribution-free coverage guarantee.

Primary sources:

- N. Vlassis and P. S. Thomas, *An Exact Distribution-Free Test for Means of
  Nonnegative Random Variables*, arXiv:2607.08415 (2026).
- J. Ming, A. Ramdas, Y. Shen, R. Wang, and I. Waudby-Smith, *Gaffke's
  confidence interval for the mean of bounded data is inadmissible but
  asymptotically efficient*, arXiv:2607.18661 (2026).

## 1. Reduction to a four-simplex cap

Write `p_j = p_4(j)` and define the ascending-knot weights

```text
a = 1-p_3,
b = p_3-p_2,
c = p_2-p_1,
d = p_1-p_0,
e = p_0.
```

For ascending sample taints `t1 <= t2 <= t3 <= t4`,

```text
SB = a t1 + b t2 + c t3 + d t4 + e.
```

The weights are positive and sum to one. If
`(D0,...,D4) ~ Dirichlet(1,...,1)`, the one-sided Gaffke upper bound is the
`(1-alpha)` quantile of the random convex combination of
`(t1,t2,t3,t4,1)`. Affine equivariance reduces the comparison to the knots

```text
(0,x,y,z,1),       0 <= x <= y <= z <= 1,
s = b x + c y + d z + e.
```

It is enough to prove

```text
V(x,y,z) := P_D(xD1+yD2+zD3+D4 > s) <= alpha.            (1)
```

The uniform-simplex cap formula is

```text
        (1-s)^4
V = ----------------------------,                         s >= z,
    (1-x)(1-y)(1-z)

        (1-s)^4                    (z-s)^4
V = ---------------------------- - ---------------------, y <= s <= z,
    (1-x)(1-y)(1-z)          z(z-x)(z-y)(1-z)

        (1-s)^4                    (z-s)^4
V = ---------------------------- - ---------------------
    (1-x)(1-y)(1-z)          z(z-x)(z-y)(1-z)

                  (y-s)^4
    + ----------------------------------------,           x <= s <= y,
       y(y-x)(z-y)(1-y)

              s^4
V = 1 - --------------,                                  s <= x.
              xyz
```

Boundary values are understood by continuity. The Clopper--Pearson equations
give four identities, one for each possible Bernoulli error count:

```text
alpha = (1-e)^4,
alpha = (1-d-e)^3 [1+3(d+e)],
alpha = (1-c-d-e)^2 [1+2(c+d+e)+3(c+d+e)^2],
1-alpha = (b+c+d+e)^4.                                  (2)
```

## 2. The region `s <= x`: analytic AM--GM

Put `g=b+c+d+e=p_3=(1-alpha)^(1/4)`. At each certified level the exact
factor enclosures prove

```text
b <= g/4,       c <= g/4,       d <= g/4.                (3)
```

Weighted AM--GM gives

```text
s/g = (b/g)x + (c/g)y + (d/g)z + e/g
    >= x^(b/g)y^(c/g)z^(d/g)
    >= (xyz)^(1/4),
```

where the last step uses (3) and `0 <= x,y,z <= 1`. Hence
`s^4 >= (1-alpha)xyz`, so the final cap formula is at most `alpha`.

## 3. The three polynomial regions

For the other regions, clear the positive denominators in the cap formulas.
When `s >= z`, it is enough to prove

```text
P_A = alpha(1-x)(1-y)(1-z) - (1-s)^4 >= 0.              (4)
```

For `y <= s <= z`, use

```text
D_B = (1-x)(1-y)(1-z) z(z-x)(z-y),
P_B = alpha D_B
      -(1-s)^4 z(z-x)(z-y)
      +(z-s)^4(1-x)(1-y).
```

The second identity in (2) gives the exact factorization

```text
P_B = (1-z) R_B,                                        (5)
```

where `R_B` has degree six. Since `1-z >= 0`, it remains to prove
`R_B >= 0`.

For `x <= s <= y`, let

```text
D_C = (1-x)(1-y)(1-z) z(z-x)(z-y) y(y-x),
P_C = alpha D_C
      -(1-s)^4 z(z-x)(z-y)y(y-x)
      +(z-s)^4(1-x)(1-y)y(y-x)
      -(y-s)^4(1-x)(1-z)z(z-x).
```

The third identity in (2) gives

```text
P_C = (1-z)(1-y)(z-y) R_C,                              (6)
```

where `R_C` also has degree six. Every extracted factor in (6) is
nonnegative on the ordered-knot domain, so it remains to prove `R_C >= 0`.
The symbolic derivation script performs the polynomial divisions in (5)--(6)
over the exact ring `QQ[b,c,d,e,x,y,z]` and rejects a nonzero remainder.

## 4. Exact polytope decomposition

The three regions in (4)--(6) are convex polytopes. Define

```text
r0 = (0,0,0),
r1 = (0,0,e/(1-d)),
r2 = (0,e/(1-c-d),e/(1-c-d)),
r3 = (t,t,t),                 t=e/(1-b-c-d),
v001 = (0,0,1),
q2 = (0,(d+e)/(1-c),1),
q3 = (u,u,1),                 u=(d+e)/(1-b-c),
v011 = (0,1,1),
q4 = ((c+d+e)/(1-b),1,1).
```

Direct intersection of the order constraints with the planes `s=z`, `s=y`,
and `s=x` gives:

- region A (`s >= z`): the tetrahedron
  `[r0,r1,r2,r3]`;
- region B (`y <= s <= z`): the convex hull of
  `{v001,r1,r2,q2,r3,q3}`, triangulated as

  ```text
  [r3,r2,r1,v001],
  [r3,q2,q3,v001],
  [r3,q2,r2,v001];
  ```

- region C (`x <= s <= y`): the convex hull of
  `{v011,r2,q2,r3,q3,q4}`, triangulated as

  ```text
  [r3,q2,r2,v011],
  [q3,r3,q4,v011],
  [q3,r3,q2,v011].
  ```

All denominators in these vertices are positive because the five Stringer
weights are positive and sum to one. The decompositions follow by solving the
active triples of linear boundary equations; the listed tetrahedra meet only
on shared faces and cover the corresponding convex hulls.

Map each tetrahedron affinely to the standard tetrahedron

```text
U >= 0, V >= 0, W >= 0, U+V+W <= 1.
```

A polynomial is nonnegative there whenever every coefficient in its
fixed-degree tetrahedral Bernstein expansion is nonnegative. The exact
structure calculation finds:

| polynomial | degree | tetrahedra | structural-zero counts |
|---|---:|---:|---|
| `P_A` | 4 | 1 | 1 |
| `R_B` | 6 | 3 | 10, 5, 10 |
| `R_C` | 6 | 3 | 5, 10, 10 |

The zero coefficients are identities in `QQ(b,c,d,e)`, not numerical
near-zeros. The derivation script represents each coefficient as a sparse
rational function and classifies it as zero only when its exact numerator
cancels. Every other coefficient is then evaluated with rational interval
arithmetic at each of the three confidence levels.

The smallest positive exact lower bounds are:

| `alpha` | `P_A` | `R_B` | `R_C` |
|---:|---:|---:|---:|
| 0.01 | 3.10e-10 | 4.78e-14 | 5.84e-12 |
| 0.05 | 3.47e-7 | 2.66e-10 | 4.60e-9 |
| 0.10 | 7.74e-6 | 1.18e-8 | 7.93e-8 |

Thus (4)--(6) are nonnegative throughout their regions. Together with the
AM--GM region, this proves (1). Pointwise domination of the valid Gaffke
upper limit proves distribution-free finite-sample conservatism of the
binomial Stringer bound at `n=4` and 90%, 95%, and 99% nominal confidence.
The all-`n` factor comparison in `POISSON-DOMINATION.md` transfers the same
conclusion to the Poisson-factor Stringer bound.

## 5. What is exact and how to reproduce it

Run:

```bash
make n4-structure-check
make n4-certificate-check
```

The first command:

1. starts from the four-simplex cap formula;
2. derives the exact residual polynomials `P_A`, `R_B`, and `R_C`;
3. verifies the factorizations (5)--(6) by exact polynomial division;
4. applies every tetrahedral affine map; and
5. derives all structural-zero indices over `QQ(b,c,d,e)`.

It regenerates `n4-gaffke-bernstein-structure.json` byte-for-byte.

The second command:

1. encloses every `n=4` Clopper--Pearson factor on a `2^-120` dyadic grid;
2. checks every binomial-CDF bracket sign with integer arithmetic;
3. verifies the AM--GM margins in (3);
4. evaluates every nonzero Bernstein coefficient with exact `Fraction`
   interval arithmetic; and
5. regenerates `n4-gaffke-certificate.json` byte-for-byte.

Floating-point arithmetic is used only to print decimal summaries, never to
decide a factor sign, structural zero, or Bernstein-coefficient sign.
