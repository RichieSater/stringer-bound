"""Checks for the exact reductions in ALL-N-POISSON-PROGRAM.md.

This script does *not* prove either open all-n quantile inequality.  It checks
the algebra that isolates them:

* the ordered-weight Poisson Stringer identity;
* the exact two-exponential quantile-convexity threshold ``4*exp(-3)``;
* the exact equal-weight Hessian in every dimension;
* the three-exponential tilted-simplex curvature reduction, its proved
  repeated-maximum boundary, and its proved equal-smaller symmetry line;
* the gamma-kernel antiderivative identity; and
* a decisive numerical counterexample to the tempting SymPol shortcut.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from mpmath import mp
from sympy import (Function, Rational, diff, exp, factorial, integrate, log,
                   simplify, symbols, series, together)

from stringer import exact_exp_neg_bounds, factor_prefix


def check_ordered_weight_identity() -> None:
    """Symbolically compare the two forms of Poisson Stringer."""
    for n in range(1, 9):
        p = symbols(f"p0:{n + 1}")
        t = symbols(f"t1:{n + 1}")
        direct = p[0] + sum(
            (p[j] - p[j - 1]) * t[j - 1]
            for j in range(1, n + 1)
        )
        ordered = ((1 - t[0]) * p[0]
                   + sum((t[j - 1] - t[j]) * p[j]
                         for j in range(1, n))
                   + t[n - 1] * p[n])
        assert simplify(direct - ordered) == 0


def check_two_exponential_obstruction() -> None:
    """Verify q'(1) and q''(1) by an implicit tail expansion."""
    h, x = symbols("h x", positive=True)
    t = 1 + h
    # Upper tail of E_1+tE_2, continued analytically through t=1.
    tail = (exp(-x) - t * exp(-x / t)) / (1 - t)
    expansion = series(tail, h, 0, 3).removeO()
    expected = (exp(-x) * (1 + x)
                + h * x**2 * exp(-x) / 2
                + h**2 * exp(-x) * (x**3 / 6 - x**2 / 2))
    assert simplify(expansion - expected) == 0

    q1 = x / 2
    q2 = x * (x - 3) / 12
    q_path = x + q1 * h + q2 * h**2 / 2
    implicit = series(
        expansion.subs(x, q_path), h, 0, 3).removeO()
    baseline = exp(-x) * (1 + x)
    assert simplify(implicit - baseline) == 0


def check_two_exponential_global_convexity() -> None:
    """Verify the exact identities behind the global two-weight theorem.

    The written proof in ``TWO-EXPONENTIAL-QUANTILE.md`` supplies the sign
    arguments.  This routine independently checks the nontrivial symbolic
    factorization and the coefficient formula used there.
    """

    z = symbols("z", positive=True)
    ez = exp(z)
    b = (1 - exp(-z)) / z
    a = 1 - exp(-z) - z**2 / (ez - 1)
    h = simplify(z * diff(a, z) / a)
    expected_h = (
        z * (z * ez - ez + 1)**2
        / ((ez - 1) * ((ez - 1)**2 - z**2 * ez))
    )
    assert simplify(h - expected_h) == 0

    margin = 3 - h + log(1 + h * b) - log(4)
    k = (2 * z**2 * ez + z * ez**2 - z
         - 4 * ez**2 + 8 * ez - 4)
    expected_derivative = (
        (ez - 1 - z)**2 * (z * ez - ez + 1)**2 * k * ez
        / (
            (ez - 1)**2
            * (ez**2 - 1 - 2 * z * ez)
            * (z**2 * ez - ez**2 + 2 * ez - 1)**2
        )
    )
    assert simplify(diff(margin, z) - expected_derivative) == 0

    # K(z) has zero coefficients through degree five.  For m>=2 its
    # z^m coefficient, multiplied by m!, is the following integer.  It is
    # positive at m=6,7 and manifestly positive for every m>=8.
    def scaled_k_coefficient(m: int) -> int:
        if m == 0 or m == 1:
            return 0
        return 2 * m * (m - 1) + (m - 8) * 2 ** (m - 1) + 8

    k_series = series(k, z, 0, 25).removeO().expand()
    assert k_series.coeff(z, 0) == 0
    assert k_series.coeff(z, 1) == 0
    for m in range(2, 25):
        assert simplify(
            factorial(m) * k_series.coeff(z, m)
            - scaled_k_coefficient(m)
        ) == 0
    assert [scaled_k_coefficient(m) for m in range(6)] == [0] * 6
    assert scaled_k_coefficient(6) == 4
    assert scaled_k_coefficient(7) == 28
    assert all(scaled_k_coefficient(m) > 0 for m in range(8, 65))

    # These series also check the removable z=0 endpoint used in the proof:
    # h(0)=3 and the comparison margin starts positively at order two.
    assert simplify(series(h, z, 0, 3).removeO()
                    - (3 - z / 2 - z**2 / 60)) == 0
    assert simplify(series(margin, z, 0, 3).removeO()
                    - 3 * z**2 / 40) == 0


def check_equal_weight_hessian_thresholds() -> None:
    """Check the all-dimensional Hessian formula at equal weights.

    The written argument uses the Dirichlet conditional law.  Here we verify
    its density-weighted conditional-variance differentiation exactly for a
    range of symbolic dimensions, along with the Poisson-CDF derivative used
    to prove that the local thresholds increase with dimension.
    """

    c, lam, s1, s2, x = symbols("c lam s1 s2 x", positive=True)
    for m in range(2, 13):
        density = (
            x ** (m - 1) * exp(-x / c)
            / (c**m * factorial(m - 1))
        )
        conditional_variance = (
            (x / c)**2 * (m * s2 - s1**2)
            / (m**2 * (m + 1))
        )
        curvature = simplify(
            -diff(density * conditional_variance, x) / density
        )
        expected = (
            x * (x - c * (m + 1)) * (m * s2 - s1**2)
            / (c**3 * m**2 * (m + 1))
        )
        assert simplify(curvature - expected) == 0

        poisson_cdf = exp(-lam) * sum(
            lam**j / factorial(j) for j in range(m)
        )
        poisson_mass = exp(-lam) * lam ** (m - 1) / factorial(m - 1)
        assert simplify(diff(poisson_cdf, lam) + poisson_mass) == 0


def check_three_exponential_reduction() -> None:
    """Check the exact tail and curvature prefactor in the three-weight note."""

    w, x, z = symbols("w x z", positive=True)

    def bfun(value):
        return (1 - exp(-value)) / value

    partition = (bfun(w) - bfun(z)) / (z - w)
    linear = (z * bfun(w) - w * bfun(z)) / (z - w)
    assert simplify(linear - bfun(z) - z * partition) == 0
    assert simplify(linear - bfun(w) - w * partition) == 0

    a = x / (x + z)
    b = x / (x + w)
    weights = (a, b, Rational(1, 1))
    direct_tail = 0
    for i, weight in enumerate(weights):
        denominator = 1
        for j, other in enumerate(weights):
            if i != j:
                denominator *= weight - other
        direct_tail += weight**2 * exp(-x / weight) / denominator
    expected_tail = exp(-x) * (
        1 + linear * x + partition * x**2
    )
    assert simplify(direct_tail - expected_tail) == 0

    # At fixed original weights, z and w vary proportionally with x.  This
    # total derivative checks the sign and the factor x/Z in equation (10).
    k = Function("k")(z, w)
    weighted_covariance = exp(-x) * x**4 * k
    total_derivative = (
        diff(weighted_covariance, x)
        + z / x * diff(weighted_covariance, z)
        + w / x * diff(weighted_covariance, w)
    )
    density_without_constants = exp(-x) * x**2 * partition
    radial_k = z * diff(k, z) + w * diff(k, w)
    expected_curvature = x / partition * ((x - 4) * k - radial_k)
    assert simplify(
        -total_derivative / density_without_constants
        - expected_curvature
    ) == 0


def check_three_exponential_repeated_max_boundary() -> None:
    """Check the exact repeated-max boundary proof for the open 3D target."""

    s, v = symbols("s v", positive=True)
    i0 = integrate((1 - v)**3 * exp(-s * v), (v, 0, 1))
    i1 = integrate(v * (1 - v)**3 * exp(-s * v), (v, 0, 1))
    r_minus = simplify(s * i1 / i0)

    raw = [
        integrate(v**j * (1 - v) * exp(-s * v), (v, 0, 1))
        for j in range(4)
    ]
    mean = raw[1] / raw[0]
    variance = raw[2] / raw[0] - mean**2
    weighted_centered_square = (
        raw[3] / raw[0]
        - 2 * mean * raw[2] / raw[0]
        + mean**2 * raw[1] / raw[0]
    )
    r_plus = simplify(s * weighted_centered_square / variance)

    # Verify directly that (2,1) and (0,1) are the two generalized
    # eigendirections of (s E[V(R-mu)(R-mu)^T], Cov(R)).  Integrating out U
    # leaves the common triangle normalizer raw[0].
    def triangle_expectation(integrand):
        return integrate(integrand * exp(-s * v), (v, 0, 1)) / raw[0]

    mean_v = mean
    mean_u = triangle_expectation((1 - v)**2 / 2)
    ev2 = raw[2] / raw[0]
    euv = triangle_expectation(v * (1 - v)**2 / 2)
    eu2 = triangle_expectation((1 - v)**3 / 3)
    c_uv = simplify(euv - mean_u * mean_v)
    c_vv = variance
    c_uu = simplify(eu2 - mean_u**2)

    evu2 = triangle_expectation(v * (1 - v)**3 / 3)
    ev2u = triangle_expectation(v**2 * (1 - v)**2 / 2)
    h_uu = simplify(
        s * (evu2 - 2 * mean_u * euv + mean_u**2 * mean_v)
    )
    h_uv = simplify(
        s * (
            ev2u - mean_v * euv - mean_u * ev2
            + mean_u * mean_v**2
        )
    )
    h_vv = simplify(s * weighted_centered_square)

    assert simplify(2 * c_uv + c_vv) == 0
    assert simplify(2 * h_uv + h_vv) == 0
    assert simplify(
        (4 * h_uu + 4 * h_uv + h_vv)
        / (4 * c_uu + 4 * c_uv + c_vv)
        - r_minus
    ) == 0
    assert simplify(h_vv / c_vv - r_plus) == 0

    nfun = s**2 * exp(s) - 2 * s * exp(s) + 2 * exp(s) - 2
    dfun = (
        s**3 * exp(s) - 3 * s**2 * exp(s) + 6 * s * exp(s)
        - 6 * exp(s) + 6
    )
    lfun = (
        s**4 * exp(s) + 2 * s**2 * exp(2 * s)
        + 8 * s**2 * exp(s) + 2 * s**2
        - 12 * s * exp(2 * s) + 12 * s
        + 12 * exp(2 * s) - 24 * exp(s) + 12
    )

    difference_numerator, difference_denominator = together(
        r_plus - r_minus
    ).as_numer_denom()
    assert simplify(
        difference_numerator - s**2 * nfun * lfun * exp(s)
    ) == 0
    gfun = simplify(s**2 * exp(s) * raw[0])
    hden = simplify(s**2 * gfun**2 * variance)
    assert simplify(difference_denominator - gfun * dfun * hden) == 0

    h = simplify(4 - r_minus)
    assert simplify(h - 3 * s * nfun / dfun) == 0

    bfun = (1 - exp(-s)) / s
    partition = (1 - bfun) / s
    positive_tail_factor = 1 + h + partition * h**2
    margin = 3 - h + log(positive_tail_factor) - log(4)
    expected_derivative = (
        -27 * s**2 * nfun**2 * lfun
        / (dfun**4 * positive_tail_factor)
    )
    assert simplify(diff(margin, s) - expected_derivative) == 0

    lseries = series(lfun, s, 0, 19).removeO().expand()
    assert all(lseries.coeff(s, j) == 0 for j in range(8))
    for j in range(8, 19):
        scaled = simplify(factorial(j) * lseries.coeff(s, j))
        expected = (
            j * (j - 1) * (j - 2) * (j - 3)
            + 8 * j * (j - 1) - 24
            + 2 ** (j - 1) * (j**2 - 13 * j + 24)
        )
        assert scaled == expected
        assert scaled > 0


def check_three_exponential_equal_smaller_line() -> None:
    """Check the exact two-equal-smaller-weights proof for the 3D target."""

    s, y = symbols("s y", positive=True)
    moments = [
        integrate(y**j * exp(-s * y), (y, 0, 1))
        for j in range(1, 5)
    ]
    mean = moments[1] / moments[0]
    variance = moments[2] / moments[0] - mean**2
    weighted_centered_square = (
        moments[3] / moments[0]
        - 2 * mean * moments[2] / moments[0]
        + mean**2 * moments[1] / moments[0]
    )
    r_antisymmetric = simplify(s * moments[3] / moments[2])
    r_symmetric = simplify(s * weighted_centered_square / variance)

    cfun = s**2 + 2 * s - 2 * exp(s) + 2
    pfun = s * exp(s) + s - 2 * exp(s) + 2
    qfun = (
        s**3 * exp(s) - s**2 * exp(s) - s**2
        + 4 * s * exp(s) - 4 * s
        - 2 * exp(2 * s) + 4 * exp(s) - 2
    )
    lfun = (
        s**4 * exp(s) + 2 * s**2 * exp(2 * s)
        + 8 * s**2 * exp(s) + 2 * s**2
        - 12 * s * exp(2 * s) + 12 * s
        + 12 * exp(2 * s) - 24 * exp(s) + 12
    )

    eigen_numerator, _ = together(
        r_antisymmetric - r_symmetric
    ).as_numer_denom()
    assert simplify(eigen_numerator + s**2 * cfun * lfun) == 0

    h = simplify(4 - r_symmetric)
    expected_h = s**2 * pfun**2 / ((s - exp(s) + 1) * qfun)
    assert simplify(h - expected_h) == 0
    assert simplify(
        variance + qfun / (s**2 * (s - exp(s) + 1)**2)
    ) == 0

    rfun = (
        s**4 * exp(2 * s) - s**4 * exp(s)
        - 2 * s**3 * exp(2 * s) - 2 * s**3 * exp(s)
        + s**2 * exp(3 * s) + 9 * s**2 * exp(2 * s)
        - 9 * s**2 * exp(s) - s**2
        - 6 * s * exp(3 * s) + 6 * s * exp(2 * s)
        + 6 * s * exp(s) - 6 * s
        + 4 * exp(3 * s) - 12 * exp(2 * s)
        + 12 * exp(s) - 4
    )
    expected_h_derivative = (
        s * cfun * pfun * rfun
        / ((s - exp(s) + 1)**2 * qfun**2)
    )
    assert simplify(diff(h, s) - expected_h_derivative) == 0

    rseries = series(rfun, s, 0, 25).removeO().expand()
    assert all(rseries.coeff(s, j) == 0 for j in range(9))
    for j in range(4, 25):
        scaled = simplify(factorial(j) * rseries.coeff(s, j))
        expected = (
            Rational(2**j, 16)
            * (j**4 - 10 * j**3 + 59 * j**2 - 2 * j - 192)
            + Rational(3**j, 9) * (j**2 - 19 * j + 36)
            - j**4 + 4 * j**3 - 14 * j**2 + 17 * j + 12
        )
        assert scaled == expected
        if j >= 9:
            assert scaled > 0

    zdiag = (1 - (s + 1) * exp(-s)) / s**2
    adiag = 2 * (1 - exp(-s)) / s - exp(-s)
    tail_factor = 1 + adiag * h + zdiag * h**2
    margin = 3 - h + log(tail_factor) - log(4)
    ffun = (
        (-2 * s**2 + 12 * s - 4) * exp(4 * s)
        + (-s**4 - 4 * s**3 - 20 * s**2 - 32 * s + 16)
        * exp(3 * s)
        + (
            2 * s**6 - 2 * s**5 + 18 * s**4 + 24 * s**3
            + 48 * s**2 + 24 * s - 24
        ) * exp(2 * s)
        + (
            -s**6 - 6 * s**5 - 21 * s**4 - 20 * s**3
            - 28 * s**2 + 16
        ) * exp(s)
        + 2 * s**2 - 4 * s - 4
    )
    expected_margin_derivative = (
        cfun**2 * pfun**2
        * (s**2 * exp(s) - (exp(s) - 1)**2) * ffun
        / (
            (s - exp(s) + 1)**2 * qfun**4 * exp(s)
            * tail_factor
        )
    )
    assert simplify(diff(margin, s) - expected_margin_derivative) == 0

    fseries = series(ffun, s, 0, 41).removeO().expand()
    assert all(fseries.coeff(s, j) == 0 for j in range(12))
    assert [factorial(j) * fseries.coeff(s, j) for j in range(12, 15)] == [
        18480, 240240, 480480,
    ]
    for j in range(4, 41):
        scaled = simplify(factorial(j) * fseries.coeff(s, j))
        expected = (
            Rational(4**j, 8) * (-j**2 + 25 * j - 32)
            + Rational(2**j, 32)
            * (
                j**6 - 17 * j**5 + 141 * j**4 - 415 * j**3
                + 866 * j**2 - 192 * j - 768
            )
            - Rational(3**j, 81)
            * (j**4 + 6 * j**3 + 155 * j**2 + 702 * j - 1296)
            - j**6 + 9 * j**5 - 46 * j**4 + 121 * j**3
            - 173 * j**2 + 90 * j + 16
        )
        assert scaled == expected
        if 15 <= j:
            assert scaled < 0

    # Exact rational endpoint checks.  The elementary interval operations
    # below are inclusion-monotone and introduce no floating-point step.
    def iadd(left, right):
        return left[0] + right[0], left[1] + right[1]

    def ineg(value):
        return -value[1], -value[0]

    def isub(left, right):
        return iadd(left, ineg(right))

    def imul(left, right):
        products = (
            left[0] * right[0], left[0] * right[1],
            left[1] * right[0], left[1] * right[1],
        )
        return min(products), max(products)

    def ipow(value, exponent):
        result = (Fraction(1), Fraction(1))
        for _ in range(exponent):
            result = imul(result, value)
        return result

    def idiv(left, right):
        assert not (right[0] <= 0 <= right[1])
        reciprocal = (Fraction(1, right[1]), Fraction(1, right[0]))
        return imul(left, reciprocal)

    def iconstant(value):
        value = Fraction(value)
        return value, value

    def endpoint_intervals(value, exp_bracket):
        value = Fraction(value)
        sint = iconstant(value)
        rint = exp_bracket
        h_numerator = imul(
            imul(ipow(sint, 2), rint),
            ipow(iadd(iconstant(value - 2), imul(
                iconstant(value + 2), rint)), 2),
        )
        h_denominator = imul(
            isub(imul(iconstant(value + 1), rint), iconstant(1)),
            isub(
                iadd(
                    imul(
                        iconstant(value**3 - value**2 + 4 * value + 4),
                        rint,
                    ),
                    ineg(imul(
                        iconstant(value**2 + 4 * value + 2),
                        ipow(rint, 2),
                    )),
                ),
                iconstant(2),
            ),
        )
        h_interval = idiv(h_numerator, h_denominator)

        coefficients = (
            (-2 * value**2 + 12 * value - 4),
            (-value**4 - 4 * value**3 - 20 * value**2
             - 32 * value + 16),
            (2 * value**6 - 2 * value**5 + 18 * value**4
             + 24 * value**3 + 48 * value**2 + 24 * value - 24),
            (-value**6 - 6 * value**5 - 21 * value**4
             - 20 * value**3 - 28 * value**2 + 16),
            (2 * value**2 - 4 * value - 4),
        )
        scaled_f = iconstant(0)
        for exponential_power, coefficient in zip(
                range(4, -1, -1), coefficients):
            scaled_f = iadd(
                scaled_f,
                imul(iconstant(coefficient), ipow(
                    rint, 4 - exponential_power)),
            )

        a_interval = isub(
            imul(iconstant(Fraction(2) / value), isub(
                iconstant(1), rint)),
            rint,
        )
        z_interval = imul(
            iconstant(Fraction(1) / value**2),
            isub(iconstant(1), imul(iconstant(value + 1), rint)),
        )
        return h_interval, scaled_f, a_interval, z_interval

    ten_billion = 10**10
    simple_brackets = {
        Fraction(1): (
            Fraction(3678794411, ten_billion),
            Fraction(3678794412, ten_billion),
        ),
        Fraction(6, 5): (
            Fraction(3011942119, ten_billion),
            Fraction(3011942120, ten_billion),
        ),
    }
    for value, bracket in simple_brackets.items():
        exp_lower, exp_upper = exact_exp_neg_bounds(value)
        assert bracket[0] < exp_lower < exp_upper < bracket[1]

    one = endpoint_intervals(Fraction(1), simple_brackets[Fraction(1)])
    six_fifths = endpoint_intervals(
        Fraction(6, 5), simple_brackets[Fraction(6, 5)],
    )
    assert one[0][1] < Fraction(347, 100)
    assert one[1][0] > Fraction(45, 10**8)
    assert one[1][1] < Fraction(48, 10**8)
    assert six_fifths[0][0] > Fraction(84, 25)
    assert Fraction(-27, 10**7) < six_fifths[1][0]
    assert six_fifths[1][1] < Fraction(-26, 10**7)
    assert six_fifths[2][0] > Fraction(43, 50)
    assert six_fifths[3][0] > Fraction(117, 500)

    x = Fraction(47, 100)
    exp_x_upper = (
        1 + x + x**2 / 2 + x**3 / 6
        + x**4 / (24 * (1 - x / 5))
    )
    assert exp_x_upper == Fraction(17395178081, 10872000000)
    assert exp_x_upper < Fraction(8, 5)


def check_kernel_identity() -> None:
    """Verify the nth-derivative identity for several symbolic n."""
    y = symbols("y", positive=True)
    for n in range(1, 11):
        left = diff(y**n * exp(-Rational(n, 1) / y), y, n) / factorial(n)
        right = (exp(-Rational(n, 1) / y)
                 * sum((Rational(n, 1) / y)**k / factorial(k)
                       for k in range(n + 1)))
        assert simplify(left - right) == 0
        # For y>1, the nth derivative of (y-1)^n/n! is one; for
        # 0<y<1, the positive-part term is identically zero.
        assert diff((y - 1)**n, y, n) / factorial(n) == 1


def sympol_shortcut_counterexample() -> tuple[mp.mpf, mp.mpf, mp.mpf, int]:
    """Return the n=15, alpha=.05 one-full-taint comparison.

    The calculation is deliberately labeled numerical.  Its margin is over
    2e-3, so it is more than enough to rule out this proposed proof route.
    """
    n = 15
    errors = 1
    alpha = Fraction(1, 20)
    with mp.workdps(80):
        poisson_p1 = factor_prefix(
            n, str(float(alpha)), errors, "poisson", 80)[errors][0]
        candidates = []
        nonzero_complements = n - errors
        for k in range(1, nonzero_complements + 1):
            elementary_mean = (mp.mpf(comb(nonzero_complements, k))
                               / comb(n, k))
            candidates.append(
                (mp.mpf(alpha.numerator) / alpha.denominator
                 * elementary_mean) ** (mp.mpf(1) / k))
        index = max(range(len(candidates)), key=candidates.__getitem__)
        sympol_upper = 1 - candidates[index]
        gap = poisson_p1 - sympol_upper
        assert gap < mp.mpf("-0.002")
        return +poisson_p1, +sympol_upper, +gap, index + 1


def main() -> int:
    check_ordered_weight_identity()
    check_two_exponential_obstruction()
    check_two_exponential_global_convexity()
    check_equal_weight_hessian_thresholds()
    check_three_exponential_reduction()
    check_three_exponential_repeated_max_boundary()
    check_three_exponential_equal_smaller_line()
    check_kernel_identity()
    poisson, sympol, gap, active_k = sympol_shortcut_counterexample()
    print("all-n reduction algebra: PASS")
    print("exact two-exponential convexity threshold: 4*exp(-3) =",
          mp.nstr(4 * mp.e**-3, 16))
    print("equal-weight local thresholds: strictly increasing from the "
          "same m=2 value")
    print("three-exponential convexity: reduced exactly to the documented "
          "two-variable inequality")
    print("three-exponential repeated-max boundary: PROVED")
    print("three-exponential equal-smaller symmetry line: PROVED")
    print("SymPol shortcut counterexample (numerical):")
    print("  n=15 alpha=0.05 one full taint")
    print("  active elementary-symmetric order:", active_k)
    print("  Poisson Stringer:", mp.nstr(poisson, 16))
    print("  SymPol upper:", mp.nstr(sympol, 16))
    print("  Stringer - SymPol:", mp.nstr(gap, 16))
    print("higher-dimensional inequality A and general inequality B: "
          "NOT CLAIMED AS PROVED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
