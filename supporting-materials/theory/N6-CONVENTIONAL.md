# Certified finite-sample conservatism at `n = 6` and 95% confidence

This note gives the complete reduction and exact-computation boundary for the
manuscript's `n = 6` result. It proves that the binomial-factor Stringer bound
is conservative for every distribution on `[0,1]` when

```text
n = 6,       alpha = 0.05,       nominal confidence = 95%.
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

## 1. Reduction to a six-simplex cap

Write `p_j = p_6(j)` and define the ascending-knot weights

```text
a = 1-p_5,
b = p_5-p_4,
c = p_4-p_3,
d = p_3-p_2,
e = p_2-p_1,
f = p_1-p_0,
g = p_0.
```

For ascending sample taints `t1 <= ... <= t6`,

```text
SB = a t1 + b t2 + c t3 + d t4 + e t5 + f t6 + g.
```

The seven weights are positive and sum to one. If the smallest taint is one,
both endpoints equal one. Otherwise, affine equivariance (subtract the
smallest taint and divide by its distance from one) reduces the comparison to

```text
(0,x,y,z,w,u,1),       0 <= x <= y <= z <= w <= u <= 1,
s = b x + c y + d z + e w + f u + g.                    (1)
```

If `(D0,...,D6) ~ Dirichlet(1,...,1)`, it is enough to prove

```text
V(x,y,z,w,u) := P_D(xD1+yD2+zD3+wD4+uD5+D6 > s)
                <= alpha.                                (2)
```

Define

```text
H1 = (1-s)^6 / [(1-x)(1-y)(1-z)(1-w)(1-u)],
Hu = (u-s)^6 / [u(u-x)(u-y)(u-z)(u-w)(1-u)],
Hw = (w-s)^6 / [w(w-x)(w-y)(w-z)(u-w)(1-w)],
Hz = (z-s)^6 / [z(z-x)(z-y)(w-z)(u-z)(1-z)],
Hy = (y-s)^6 / [y(y-x)(z-y)(w-y)(u-y)(1-y)].
```

The uniform-simplex cap formula is

```text
V = H1,                                           s >= u,
V = H1-Hu,                                   w <= s <= u,
V = H1-Hu+Hw,                                z <= s <= w,
V = H1-Hu+Hw-Hz,                             y <= s <= z,
V = H1-Hu+Hw-Hz+Hy,                          x <= s <= y,
V = 1-s^6/(xyzwu),                                s <= x. (3)
```

Boundary values are understood by continuity. Put

```text
q1=f+g, q2=e+f+g, q3=d+e+f+g, q4=c+d+e+f+g,
h=b+c+d+e+f+g=p_5.
```

The Clopper--Pearson equations give

```text
alpha = (1-g)^6,
alpha = (1-q1)^5 (1+5q1),
alpha = (1-q2)^4 (1+4q2+10q2^2),
alpha = (1-q3)^3 (1+3q3+6q3^2+10q3^3),
alpha = (1-q4)^2 (1+2q4+3q4^2+4q4^3+5q4^4),
1-alpha = h^6.                                             (4)
```

## 2. Exact residual factorizations

In the first five regions, multiply `alpha-V` by the positive common
product of the denominators appearing in (3). Call the cleared residuals
`P_A,...,P_E`; in the last region use

```text
P_F = s^6 - h^6 x y z w u.                                (5)
```

Substitution of (4), followed by exact polynomial division in the
coordinate-and-weight polynomial ring, extracts all coincident-knot and
upper-endpoint factors. On the ordered-knot domain those factors are
nonnegative. The residual polynomials that remain to be certified have
these degrees:

| region | threshold range | degree | extracted factors |
|---|---|---:|---|
| A | `s >= u` | 6 | none |
| B | `w <= s <= u` | 10 | `1-u` |
| C | `z <= s <= w` | 12 | `(1-w)(1-u)(u-w)` |
| D | `y <= s <= z` | 12 | `(1-z)(1-w)(1-u)(w-z)(u-z)(u-w)` |
| E | `x <= s <= y` | 10 | all `1-t` and pair differences for `t in {y,z,w,u}` |
| F | `s <= x` | 6 | none |

The derivation rejects any nonzero polynomial-division remainder.

## 3. Five-simplex decompositions

Let `V_i=(0^(5-i),1^i)`, for `i=0,...,5`, be the vertices of the ordered
knot simplex. For `r=1,...,5`, write `x_r=x,y,z,w,u` and
`L_r(q)=s(q)-q_r`. Whenever the edge `V_i V_j` crosses `L_r=0`, define

```text
h_{r,ij} = V_i + tau (V_j-V_i),
tau = L_r(V_i) / [L_r(V_i)-L_r(V_j)].                     (6)
```

Successive slicing by `s=u,w,z,y,x`, followed by a fixed triangulation of
each slice, gives respectively `1,5,10,10,5,1` five-simplices. The complete
ordered vertex lists are stored in
[`n6-gaffke-bernstein-structure.json`](../computations/certificates/n6-gaffke-bernstein-structure.json).
Every coordinate in (6) is retained there as an exact rational function of
`b,...,g`; no rounded vertex is used. The sign certificate encloses every
affine determinant away from zero at `alpha=0.05`. Exact boundary identities
and directed strict-sign checks place every vertex in its claimed convex
region, and directed coordinate differences certify all 210 pairs of the 21
named vertices as distinct. The certificate also verifies the
triangulation as an oriented relative chain: the 48 internal four-facets
occur twice with opposite orientations, while the 96 unpaired four-facets
lie on the six outer ordered-simplex hyperplanes, 16 on each. The sum of the
32 positively oriented determinant intervals lies strictly between `1/2`
and `3/2`. The chain represents an integer multiple of the generator of
`H_5(Delta, boundary Delta; Z)`, and its determinant sum is that integer
times the normalized volume one. Hence its degree is one. Since each
simplex has positive orientation and lies in `Delta`, the simplices cover
the ordered domain without interior multiplicity; this
coverage is certified rather than inferred from a plotting triangulation.

## 4. Exact Bernstein signs

Map each listed five-simplex affinely to the standard five-simplex. A
polynomial is nonnegative there if every coefficient in its fixed-degree
Bernstein expansion is nonnegative. The exact structure calculation gives:

| region | degree | five-simplices | structural-zero counts by simplex |
|---|---:|---:|---|
| A | 6 | 1 | 1 |
| B | 10 | 5 | 456, 456, 302, 162, 57 |
| C | 12 | 10 | 666, 666, 911, 911, 911, 911, 813, 358, 554, 813 |
| D | 12 | 10 | 813, 358, 554, 813, 666, 666, 911, 911, 911, 911 |
| E | 10 | 5 | 456, 456, 57, 302, 162 |
| F | 6 | 1 | 1 |

The zero classification is symbolic, not tolerance-based. Let `I(F)` be
the affine ideal of a simplex face. In characteristic zero, a polynomial
belongs to `I(F)^q` exactly when it and every ordinary partial derivative of
total order below `q` vanish modulo `I(F)`. The structure derivation checks
that criterion over the generic rational-function field
`QQ(b,c,d,e,f,g)` with exact Singular reductions. It therefore proves the
required vanishing orders without specializing the weights. An exact
rational specialization is used only to choose independent generators from
redundant face equations. The numerical layer then encloses a nonzero normal
minor for each selected generator set at the certified weight intervals;
thus every generic face identity specializes to the certified case.

Every coefficient not covered by those identities is enclosed as follows:

1. each Clopper--Pearson root is bracketed by adjacent dyadic rationals on a
   `2^-320` grid;
2. each binomial-CDF sign at a bracket endpoint is evaluated with integer
   arithmetic;
3. polynomial, vertex, affine-map, determinant, and Bernstein operations use
   closed dyadic intervals on a `2^-384` grid, with multiplication and
   division rounded outward by integer division; and
4. every nonstructural Bernstein coefficient has a strictly positive lower
   endpoint.

At `alpha=0.05`, the smallest positive lower bound in each region is:

| region | A | B | C | D | E | F |
|---|---:|---:|---:|---:|---:|---:|
| lower bound | 1.23e-10 | 1.17e-18 | 6.57e-23 | 4.01e-22 | 1.66e-15 | 9.57e-4 |

Thus every residual in (5) and its five predecessors is nonnegative, so
(2) holds. Pointwise domination of the valid Gaffke upper limit proves
finite-sample conservatism of the binomial Stringer bound at `n=6` and 95%
nominal confidence. The all-`n` factor comparison in
[`POISSON-DOMINATION.md`](POISSON-DOMINATION.md) transfers the conclusion to
the uncapped Poisson-factor Stringer bound.

## 5. Reproduction and trust boundary

Run:

```bash
make n6-structure-check
make n6-certificate-check
```

The first command regenerates the polynomial factorizations, exact vertex
formulas, fixed triangulations, structural-zero index sets, and generic
face-ideal derivative proofs. It requires Singular. The second command
regenerates the directed-dyadic sign certificate and verifies every remaining
coefficient. The sign certificate contains a SHA-256 digest of the structure
file, binding the numerical layer to the exact symbolic layer.

Floating-point arithmetic is used only to render decimal summaries. It is
not used for a factor-root, structural-zero, determinant, or
Bernstein-coefficient sign decision.
