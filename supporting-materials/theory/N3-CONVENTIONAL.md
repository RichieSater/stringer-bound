# Certified finite-sample conservatism at `n = 3`

This note supplies the proof and exact-computation details for the manuscript's
`n = 3` result. It proves the binomial-factor Stringer bound conservative at
the three confidence levels most directly relevant to audit practice:

| `alpha` | nominal confidence |
|---:|---:|
| `0.10` | 90% |
| `0.05` | 95% |
| `0.01` | 99% |

The proof is a pointwise comparison with the one-sided Gaffke bound. Vlassis
and Thomas proved the finite-sample validity of Gaffke's test, and inversion
gives the bounded-mean upper limit described explicitly by Ming et al. Thus a
pointwise no-smaller Stringer limit inherits the coverage guarantee.

Primary sources:

- N. Vlassis and P. S. Thomas, *An Exact Distribution-Free Test for Means of
  Nonnegative Random Variables*, arXiv:2607.08415 (2026).
- J. Ming, A. Ramdas, Y. Shen, R. Wang, and I. Waudby-Smith, *Gaffke's
  confidence interval for the mean of bounded data is inadmissible but
  asymptotically efficient*, arXiv:2607.18661 (2026).

## 1. Reduction to a simplex-cap inequality

Let

```text
a = 1 - p_3(2),
b = p_3(2) - p_3(1),
c = p_3(1) - p_3(0),
d = p_3(0).
```

For sample taints `t1 <= t2 <= t3`, the Stringer bound is

```text
SB = a t1 + b t2 + c t3 + d.
```

The weights are positive and sum to one. If
`D = (D0,D1,D2,D3) ~ Dirichlet(1,1,1,1)`, the one-sided Gaffke upper bound is
the `(1-alpha)` quantile of the random convex combination of the four knots
`(t1,t2,t3,1)`. Consequently, Stringer pointwise dominates Gaffke exactly
when

```text
P_D(D0 t1 + D1 t2 + D2 t3 + D3 > SB) <= alpha.          (1)
```

Both sides are affine equivariant in the knots. Unless `t1 = 1` (the trivial
constant case), subtract `t1` and divide by `1-t1`. It is enough to prove (1)
for

```text
(t1,t2,t3,1) = (0,x,y,1),       0 <= x <= y <= 1,
s = b x + c y + d.
```

Let `V(x,y)` be the probability in (1). The standard uniform-simplex cap
formula gives

```text
             (1-s)^3
V = ---------------------------,                         s >= y,
       (1-x)(1-y)

        1       [ (1-s)^3       (y-s)^3  ]
V = ---------  [ ---------  -  --------- ],             x <= s <= y,
      (1-y)     [   1-x          y(y-x)   ]

                 s^3
V = 1 - ----------------,                               s <= x.
                 xy
```

The formulas agree by continuity on region boundaries.

The Clopper--Pearson defining equations imply the useful identities

```text
(a+b+c)^3                 = alpha,
(b+c+d)^3                 = 1-alpha,
(a+b)^2 (3-2(a+b))        = alpha.                       (2)
```

## 2. The region `s <= x`: an analytic AM--GM bound

Put `g = b+c+d = (1-alpha)^(1/3)`. For each certified level, the exact
factor enclosures verify

```text
b <= g/3,       c <= g/3.                                (3)
```

Weighted AM--GM gives

```text
s/g = (b/g)x + (c/g)y + (d/g)
    >= x^(b/g) y^(c/g)
    >= (xy)^(1/3),
```

where the last inequality uses (3) and `0 <= x,y <= 1`. Hence

```text
s^3 >= (1-alpha)xy,
```

and the third cap formula yields `V <= alpha`.

## 3. The region `s >= y`: a triangular Bernstein certificate

Set

```text
X = 1-x,       Y = 1-y,       L = 1-s = a+bX+cY,
A = 1-b-c = a+d.
```

The region `0 <= Y <= X <= 1` and `L <= Y` is the triangle with vertices

```text
(a/A,a/A),        (1,1),        (1,(a+b)/(1-c)).
```

On this triangle, `V <= alpha` is equivalent to nonnegativity of the cubic

```text
P_A(X,Y) = alpha X Y - (a+bX+cY)^3.                       (4)
```

After the affine map from the standard triangle, every degree-three
Bernstein coefficient of (4) is strictly positive except the coefficient at
`(X,Y)=(1,1)`. That coefficient is exactly zero by the first identity in
(2). The exact interval checker verifies all remaining signs separately for
`alpha = 0.01, 0.05, 0.10`.

## 4. The middle region `x <= s <= y`

Define distances from the threshold to the two middle knots,

```text
rho = s-x >= 0,       tau = y-s >= 0.
```

Solving for the knots gives

```text
x = [d-(1-c)rho+c tau]/A,       y = x+rho+tau.
```

The conditions `x >= 0` and `y <= 1` make the `(rho,tau)` domain the
quadrilateral with vertices

```text
(0,0),   (d/(1-c),0),   (c+d,a+b),   (0,a/(1-b)).         (5)
```

Triangulate (5) along the diagonal from `(0,0)` to `(c+d,a+b)`.

For the middle cap formula, direct differentiation and cancellation give

```text
partial V / partial y = Q(x,y) / [y^2(1-x)(y-x)^2],       (6)
```

where `Q` is a degree-five polynomial. The denominator in (6) is positive in
the open middle region. Under each affine triangle map in (5), all
degree-five Bernstein coefficients of `A^3 Q` are nonnegative. Three are
structural zeros on each triangle and every other coefficient has a strictly
positive exact rational lower bound. Therefore `V` is nondecreasing in `y`
while `x` is held fixed.

Increasing `y` preserves `s-x >= 0` and `y-s >= 0`, because their derivatives
are `c` and `1-c`. It follows that the middle-region volume is no larger than
its limit at `y=1`.

At `y=1`, put `q=1-x`. The middle-region constraint is
`q >= q0 := a/(1-b)`, and the limiting cap volume is

```text
             (a+bq)^2 [3q-(a+bq)(1+q)]
V_1(q) = ----------------------------------.              (7)
                           q^2
```

After mapping `[q0,1]` to `[0,1]`, every degree-four Bernstein coefficient of

```text
alpha q^2 - (a+bq)^2 [3q-(a+bq)(1+q)]                    (8)
```

is positive except the endpoint coefficient at `q=1`. That coefficient is
zero by the third identity in (2). Thus (7) is at most `alpha`, completing
the middle-region proof.

Together, Sections 2--4 prove (1) at all three certified levels. Pointwise
domination of the valid Gaffke bound then proves distribution-free
finite-sample conservatism of the binomial Stringer bound at `n=3` and 90%,
95%, and 99% nominal confidence. The all-`n` Poisson-factor comparison in
`POISSON-DOMINATION.md` transfers the same conclusion to the Poisson-factor
Stringer bound.

## 5. What is exact and how to reproduce it

Run:

```bash
make n3-formula-check
make n3-certificate-check
```

The first command starts from the middle cap formula, differentiates it with
SymPy, performs both affine triangle substitutions, regenerates all 42
Bernstein-coefficient formulas, and compares them byte-for-byte with
`n3-gaffke-bernstein-formulas.json`.

The second command:

1. encloses every `n=3` Clopper--Pearson factor on a `2^-120` dyadic grid;
2. checks the binomial-CDF sign at every bracket endpoint using integer
   arithmetic;
3. propagates the brackets through exact rational interval arithmetic;
4. verifies (3) and every nonzero Bernstein coefficient in (4), (6), and
   (8); and
5. regenerates `n3-gaffke-certificate.json` byte-for-byte.

The smallest positive certified Bernstein lower bound is about
`1.37e-9` (at `alpha=0.01`), many orders of magnitude larger than the factor
bracket widths. Floating-point arithmetic is used only to print readable
decimal summaries, never to decide a sign.
