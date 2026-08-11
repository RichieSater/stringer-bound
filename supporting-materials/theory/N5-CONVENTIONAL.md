# Certified finite-sample conservatism at `n = 5`

This note gives the complete reduction and exact-computation details for the
manuscript's `n = 5` result. It proves the binomial-factor Stringer bound
conservative for every distribution on `[0,1]` at:

| `alpha` | nominal confidence |
|---:|---:|
| `0.10` | 90% |
| `0.05` | 95% |
| `0.01` | 99% |

As in the `n = 3` and `n = 4` results, the proof compares Stringer
pointwise with the valid one-sided Gaffke bounded-mean upper limit. Vlassis
and Thomas proved the finite-sample validity of Gaffke's test; its inversion
for observations in `[0,1]` is the Dirichlet-average upper limit described
by Learned-Miller--Thomas and Ming et al. A pointwise no-smaller Stringer
limit inherits that distribution-free coverage guarantee.

Primary sources:

- N. Vlassis and P. S. Thomas, *An Exact Distribution-Free Test for Means of
  Nonnegative Random Variables*, arXiv:2607.08415 (2026).
- J. Ming, A. Ramdas, Y. Shen, R. Wang, and I. Waudby-Smith, *Gaffke's
  confidence interval for the mean of bounded data is inadmissible but
  asymptotically efficient*, arXiv:2607.18661 (2026).

## 1. Reduction to a five-simplex cap

Write `p_j = p_5(j)` and define the ascending-knot weights

```text
a = 1-p_4,
b = p_4-p_3,
c = p_3-p_2,
d = p_2-p_1,
e = p_1-p_0,
f = p_0.
```

For ascending sample taints `t1 <= ... <= t5`,

```text
SB = a t1 + b t2 + c t3 + d t4 + e t5 + f.
```

The six weights are positive and sum to one. If
`(D0,...,D5) ~ Dirichlet(1,...,1)`, the one-sided Gaffke upper bound is
the `(1-alpha)` quantile of the random convex combination of
`(t1,...,t5,1)`. Affine equivariance reduces the comparison to

```text
(0,x,y,z,w,1),       0 <= x <= y <= z <= w <= 1,
s = b x + c y + d z + e w + f.
```

It is enough to prove

```text
V(x,y,z,w) := P_D(xD1+yD2+zD3+wD4+D5 > s) <= alpha.    (1)
```

Put

```text
H1 = (1-s)^5 / [(1-x)(1-y)(1-z)(1-w)],
Hw = (w-s)^5 / [w(w-x)(w-y)(w-z)(1-w)],
Hz = (z-s)^5 / [z(z-x)(z-y)(w-z)(1-z)],
Hy = (y-s)^5 / [y(y-x)(z-y)(w-y)(1-y)].
```

The uniform-simplex cap formula is

```text
V = H1,                                      s >= w,
V = H1-Hw,                              z <= s <= w,
V = H1-Hw+Hz,                           y <= s <= z,
V = H1-Hw+Hz-Hy,                        x <= s <= y,
V = 1-s^5/(xyzw),                            s <= x.    (2)
```

Boundary values are understood by continuity. With

```text
q1 = e+f,       q2 = d+e+f,       q3 = c+d+e+f,
g  = b+c+d+e+f = p_4,
```

the Clopper--Pearson equations give

```text
alpha = (1-f)^5,
alpha = (1-q1)^4 (1+4q1),
alpha = (1-q2)^3 (1+3q2+6q2^2),
alpha = (1-q3)^2 (1+2q3+3q3^2+4q3^3),
1-alpha = g^5.                                             (3)
```

## 2. Exact residual factorizations

Let

```text
K = (1-x)(1-y)(1-z)(1-w),
W = w(w-x)(w-y)(w-z),
Z = z(z-x)(z-y),
Y = y(y-x).
```

In the first four regions, multiply `alpha-V` respectively by the
nonnegative denominators `K`, `KW`, `KWZ`, and `KWZY`. Call the resulting
cleared residuals `P_A`, `P_B`, `P_C`, and `P_D`. In the last region the
equivalent residual is

```text
P_E = s^5 - g^5 x y z w.                                  (4)
```

Substitution of the identities in (3), followed by exact polynomial
division, gives

```text
P_B = (1-w) R_B,
P_C = (1-w)(1-z)(w-z) R_C,
P_D = (1-w)(1-z)(1-y)(w-z)(w-y)(z-y) R_D.                 (5)
```

The divisions are carried out over the exact polynomial ring
`QQ[b,c,d,e,f,x,y,z,w]`, and the derivation rejects any nonzero remainder.
Every extracted factor in (5) is nonnegative on the ordered-knot domain.
The polynomials to be certified and their total degrees are therefore

| region | polynomial | degree |
|---|---|---:|
| `s >= w` | `P_A` | 5 |
| `z <= s <= w` | `R_B` | 8 |
| `y <= s <= z` | `R_C` | 9 |
| `x <= s <= y` | `R_D` | 8 |
| `s <= x` | `P_E` | 5 |

## 3. Four-simplex decompositions

Let `V_i=(0^(4-i),1^i)`, for `i=0,...,4`, be the vertices of the ordered
knot simplex. For `r=1,2,3,4`, write `x_r=x,y,z,w` and
`L_r(q)=s(q)-q_r`. Whenever the edge `V_i V_j` crosses `L_r=0`, define

```text
h_{r,ij} = V_i + tau (V_j-V_i),
tau = L_r(V_i) / [L_r(V_i)-L_r(V_j)].                    (6)
```

All denominators in (6) are positive sums of Stringer weights. Direct
intersection of the ordered simplex with the four threshold planes gives
the following triangulations; brackets list the five vertices of a
four-simplex.

```text
A: [V0,h4_01,h4_02,h4_03,h4_04]

B: [h4_04,h3_12,h4_01,h4_03,h4_02]
   [h3_13,h4_04,h3_12,h4_01,h3_14]
   [h3_13,h4_04,h3_12,h4_01,h4_03]

C: [h2_23,h3_04,h3_02,h3_12,h3_03]
   [h2_23,h3_04,h3_02,h3_12,h2_24]
   [h2_23,h3_13,h3_04,h3_12,h3_03]
   [h3_14,h2_23,h3_04,h3_12,h2_24]
   [h3_14,h2_23,h3_13,h3_04,h3_12]

D: [h2_13,h2_04,h2_23,h1_34,h2_03]
   [h2_14,h2_13,h2_04,h2_23,h1_34]
   [h2_24,h2_14,h2_04,h2_23,h1_34]

E: [h1_04,h1_14,h1_24,h1_34,V4]
```

Some labels in (6) represent the same point on a shared threshold
boundary, for example `h4_04=h3_04=h2_04=h1_04`. The listed simplices
meet only on common faces and cover their respective convex regions. This
can be checked by slicing successively at `s=w,z,y,x`; the resulting
numbers of simplices are `1,3,5,3,1`. The exact coordinate formulas for
every listed vertex are recorded in the structure certificate.

## 4. Exact Bernstein signs

Map each listed four-simplex affinely to the standard four-simplex. A
polynomial is nonnegative there if every coefficient in its fixed-degree
Bernstein expansion is nonnegative. The exact structure calculation gives:

| region | degree | four-simplices | structural-zero counts |
|---|---:|---:|---|
| A | 5 | 1 | 1 |
| B | 8 | 3 | 66, 15, 39 |
| C | 9 | 5 | 59, 35, 92, 59, 92 |
| D | 8 | 3 | 15, 39, 66 |
| E | 5 | 1 | 1 |

The zero classification is symbolic, not tolerance-based. If `I(F)` is
the affine ideal of a simplex face and a residual belongs to
`I(F)^(q+1)`, then all Bernstein coefficients whose total index away from
that face is at most `q` vanish. The derivation proves every required
ideal-power membership by exact Groebner reduction over
`QQ(b,c,d,e,f)`. Every coefficient not covered by these identities is
then enclosed numerically, but rigorously, as follows:

1. each Clopper--Pearson root is bracketed by adjacent dyadic rationals on
   a `2^-240` grid;
2. each binomial-CDF sign at a bracket endpoint is evaluated with integer
   arithmetic;
3. polynomial, vertex, affine-map, and Bernstein operations use closed
   dyadic intervals on a `2^-256` grid, with every multiplication and
   division rounded outward by integer division; and
4. every four-simplex affine determinant is enclosed away from zero; and
5. every nonstructural Bernstein coefficient has a strictly positive
   lower endpoint.

The smallest positive certified lower bounds in each region are:

| `alpha` | A | B | C | D | E |
|---:|---:|---:|---:|---:|---:|
| 0.01 | 8.16e-13 | 7.45e-20 | 3.23e-21 | 7.81e-16 | 3.47e-4 |
| 0.05 | 6.70e-9 | 1.89e-14 | 6.24e-16 | 3.57e-12 | 1.38e-3 |
| 0.10 | 3.81e-7 | 4.93e-12 | 1.36e-13 | 1.34e-10 | 2.20e-3 |

Thus every residual in (4)--(5) is nonnegative throughout its region, so
(1) holds. Pointwise domination of the valid Gaffke upper limit proves
distribution-free finite-sample conservatism of the binomial Stringer bound
at `n = 5` and 90%, 95%, and 99% nominal confidence. The all-`n` factor
comparison in `POISSON-DOMINATION.md` transfers the same conclusion to the
Poisson-factor Stringer bound.

## 5. Reproduction and trust boundary

Run:

```bash
make n5-structure-check
make n5-certificate-check
```

The first command regenerates
`n5-gaffke-bernstein-structure.json` byte-for-byte. It derives (3)--(5),
the exact vertex formulas, all structural-zero index sets, and every
face-ideal power proof. The second command regenerates
`n5-gaffke-certificate.json` byte-for-byte, verifies the integer-checked
factor brackets, and proves every remaining coefficient positive with
directed dyadic interval arithmetic. The sign certificate contains a
SHA-256 digest of the structure file, binding the numerical layer to the
exact symbolic layer.

Floating-point arithmetic is used only to render decimal summaries. It is
not used for a factor-root, structural-zero, or coefficient-sign decision.
