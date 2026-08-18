"""Checks for the exact reductions in ALL-N-POISSON-PROGRAM.md.

This script does *not* prove either open all-n quantile inequality.  It checks
the algebra that isolates them:

* the ordered-weight Poisson Stringer identity;
* the exact two-exponential quantile-convexity threshold ``4*exp(-3)``;
* the exact equal-weight Hessian in every dimension;
* the three-exponential tilted-simplex curvature reduction and its proved
  repeated-maximum, equal-smaller, and infinite-gap boundary families,
  including the boundary-trace identity, both finite symmetry-boundary
  derivative signs, the small-gap and two-large-gap regions, and the positive
  sharp-corner expansion;
* the gamma-kernel antiderivative identity; and
* a decisive numerical counterexample to the tempting SymPol shortcut.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from mpmath import mp
from sympy import (Function, Poly, Rational, cancel, diff, exp, factorial,
                   factor, integrate, limit, log, oo, simplify, symbols,
                   series, together)

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


def check_three_exponential_trace_identity() -> None:
    """Check the boundary-trace reformulation of the 3D threshold."""

    u, v, z, w = symbols("u v z w", positive=True)
    c1, c2, mu1, mu2 = symbols("c1 c2 mu1 mu2", real=True)
    ell = c1 * (u - mu1) + c2 * (v - mu2)
    tilt = z * u + w * v
    centered_offset = c1 * mu1 + c2 * mu2

    divergence = (
        diff(u * ell**2 * exp(-tilt), u)
        + diff(v * ell**2 * exp(-tilt), v)
    )
    expected = (
        (4 - tilt) * ell**2 + 2 * centered_offset * ell
    ) * exp(-tilt)
    assert simplify(divergence - expected) == 0

    # The coordinate-edge fluxes vanish.  On u+v=1, the outward normal
    # times arclength is (1,1) du, so the radial flux is exactly the
    # boundary quadratic form in equation (13a).
    assert simplify((u * ell**2 * exp(-tilt)).subs(u, 0)) == 0
    assert simplify((v * ell**2 * exp(-tilt)).subs(v, 0)) == 0
    edge_flux = simplify(
        ((u + v) * ell**2 * exp(-tilt)).subs(v, 1 - u)
    )
    expected_edge = simplify(
        (ell**2 * exp(-tilt)).subs(v, 1 - u)
    )
    assert simplify(edge_flux - expected_edge) == 0

    bulk, weighted = symbols("bulk weighted", positive=True)
    boundary = 4 * bulk - weighted
    assert simplify(
        4 - weighted / bulk - boundary / bulk
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


def check_three_exponential_axis_transversality() -> None:
    """Check the strict fixed-sum inward derivative at the proved axis."""

    s, u, v = symbols("s u v", positive=True)
    contrast = 2 * u + v
    path_score = u - v

    def triangle_integral(integrand):
        return simplify(integrate(
            integrate(integrand, (u, 0, 1 - v)) * exp(-s * v),
            (v, 0, 1),
        ))

    partition = triangle_integral(1)
    mean = simplify(triangle_integral(contrast) / partition)
    assert mean == 1
    centered = contrast - mean
    mean_derivative = simplify(
        -triangle_integral(path_score * centered) / partition
    )

    denominator = triangle_integral(centered**2)
    denominator_derivative = -triangle_integral(
        path_score * centered**2
    )
    numerator = triangle_integral(s * v * centered**2)
    numerator_derivative = triangle_integral(
        path_score * centered**2
        - s * v * path_score * centered**2
        - 2 * s * v * mean_derivative * centered
    )
    rho = simplify(numerator / denominator)
    rho_derivative = simplify(
        (numerator_derivative * denominator
         - numerator * denominator_derivative)
        / denominator**2
    )

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
    h = simplify(4 - rho)
    assert simplify(h - 3 * s * nfun / dfun) == 0
    h_derivative = -rho_derivative

    a_derivative = -Rational(1, 2) + partition
    z_derivative = -triangle_integral(path_score)
    tail_factor = 1 + h + partition * h**2
    margin_derivative = simplify(
        -h_derivative
        + (
            a_derivative * h + h_derivative
            + z_derivative * h**2
            + 2 * partition * h * h_derivative
        ) / tail_factor
    )
    expected = (
        81 * s**2 * nfun**2 * lfun
        / (2 * dfun**4 * tail_factor)
    )
    assert simplify(margin_derivative - expected) == 0

    axis_derivative = (
        -27 * s**2 * nfun**2 * lfun
        / (dfun**4 * tail_factor)
    )
    assert simplify(margin_derivative
                    + Rational(3, 2) * axis_derivative) == 0


def check_three_exponential_diagonal_transversality() -> None:
    """Check strict fixed-sum transverse concavity on the diagonal.

    The written proof uses the boundary-trace representation.  This routine
    derives the second-order generalized-eigenvalue perturbation from the
    centered bulk and boundary Gram matrices, verifies its exact
    factorization, and checks the coefficient proof for the only remaining
    exponential polynomial.
    """

    s, t = symbols("s t", positive=True)
    moments = symbols("j0:6", positive=True)
    edge_weight = symbols("edge_weight", positive=True)

    # Along (z,w)=(s-epsilon,s+epsilon), use Y=U+V and D=U-V.
    # The entries below are the coefficients through the order needed for
    # the simple largest generalized eigenvalue.
    j0, j1, j2, j3, j4, j5 = moments
    mean_y_zero = j2 / j1
    mean_y_second = j4 / (6 * j1) - j2 * j3 / (6 * j1**2)
    mean_d_first = j3 / (3 * j1)

    k_yy_zero = j3 - j2**2 / j1
    k_yy_second = (
        j5 / 6 - j2 * j4 / (3 * j1)
        + j2**2 * j3 / (6 * j1**2)
    )
    k_dd_zero = j3 / 3
    k_yd_first = j4 / 3 - j2 * j3 / (3 * j1)

    boundary_offset_zero = 1 - mean_y_zero
    boundary_offset_second = -mean_y_second
    l_yy_zero = edge_weight * boundary_offset_zero**2
    l_yy_second = edge_weight * (
        boundary_offset_zero**2 / 6
        + 2 * boundary_offset_zero * boundary_offset_second
    )
    l_dd_zero = edge_weight / 3
    l_yd_first = edge_weight * boundary_offset_zero * (
        Rational(1, 3) - mean_d_first
    )

    h_zero = cancel(l_yy_zero / k_yy_zero)
    off_diagonal = l_yd_first - h_zero * k_yd_first
    other_eigenvalue_gap = l_dd_zero - h_zero * k_dd_zero
    h_second = cancel(
        (
            l_yy_second - h_zero * k_yy_second
            - off_diagonal**2 / other_eigenvalue_gap
        ) / k_yy_zero
    )

    a_zero = j0 + s * j1
    a_second = j2 / 2 + s * j3 / 6
    substitutions = {edge_weight: t}
    for order, moment in enumerate(moments):
        substitutions[moment] = (
            factorial(order) / s ** (order + 1)
            * (
                1 - t * sum(
                    s**power / factorial(power)
                    for power in range(order + 1)
                )
            )
        )

    h_zero_exact = factor(cancel(h_zero.subs(substitutions)))
    h_second_exact = factor(cancel(h_second.subs(substitutions)))
    a_zero_exact = cancel(a_zero.subs(substitutions))
    a_second_exact = cancel(a_second.subs(substitutions))
    partition_exact = substitutions[j1]
    partition_second_exact = substitutions[j3] / 6
    tail_factor = cancel(
        1 + a_zero_exact * h_zero_exact
        + partition_exact * h_zero_exact**2
    )
    tail_factor_second = cancel(
        a_second_exact * h_zero_exact
        + a_zero_exact * h_second_exact
        + partition_second_exact * h_zero_exact**2
        + 2 * partition_exact * h_zero_exact * h_second_exact
    )
    margin_second = factor(cancel(
        -h_second_exact + tail_factor_second / tail_factor
    ))

    exponential_polynomials = {
        7: 8 * (s**2 - 5 * s - 3),
        6: -4 * (s**4 - 8 * s**3 - 62 * s**2 - 52 * s - 54),
        5: -2 * (
            4 * s**6 + 19 * s**5 + 107 * s**4 + 400 * s**3
            + 504 * s**2 + 348 * s + 396
        ),
        4: (
            9 * s**8 + 26 * s**7 + 142 * s**6 + 544 * s**5
            + 1370 * s**4 + 2200 * s**3 + 1760 * s**2
            + 1760 * s + 1560
        ),
        3: -2 * (
            s**10 + 5 * s**9 + 16 * s**8 + 93 * s**7
            + 282 * s**6 + 678 * s**5 + 1090 * s**4
            + 1200 * s**3 + 1100 * s**2 + 1420 * s + 900
        ),
        2: (
            s**10 + 12 * s**9 + 71 * s**8 + 238 * s**7
            + 520 * s**6 + 968 * s**5 + 1240 * s**4
            + 1520 * s**3 + 2232 * s**2 + 2640 * s + 1224
        ),
        1: -2 * (
            3 * s**7 + 22 * s**6 + 71 * s**5 + 163 * s**4
            + 416 * s**3 + 704 * s**2 + 644 * s + 228
        ),
        0: 2 * (
            s**6 + 12 * s**5 + 57 * s**4 + 140 * s**3
            + 184 * s**2 + 128 * s + 36
        ),
    }
    exponential_polynomial_in_t = sum(
        polynomial * t**(-power)
        for power, polynomial in exponential_polynomials.items()
    )

    pfun = s / t + s - 2 / t + 2
    cfun = s**2 + 2 * s - 2 / t + 2
    hden = s - 1 / t + 1
    qfun = (
        s**3 / t - s**2 / t - s**2 + 4 * s / t - 4 * s
        - 2 / t**2 + 4 / t - 2
    )
    expected_margin_second = (
        t * pfun**2 * cfun * exponential_polynomial_in_t
        / (6 * hden**2 * qfun**4 * tail_factor)
    )
    assert factor(cancel(
        margin_second / expected_margin_second
    )) == 1

    # If E(s)=sum_k E_k(s)e^(ks), its nth derivative at zero is
    # sum_k k^n p_k(n), where p_k is obtained from the falling factorials
    # of the coefficients of E_k.  The k=0 polynomial has degree six and
    # therefore drops out in the range n>=18 used below.
    n = symbols("n", integer=True, nonnegative=True)
    coefficient_polynomials = {}
    for base in range(1, 8):
        derived = 0
        for (power,), coefficient in Poly(
                exponential_polynomials[base], s).terms():
            falling = 1
            for offset in range(power):
                falling *= n - offset
            derived += coefficient * falling / Rational(base)**power
        coefficient_polynomials[base] = factor(derived)

    expected_coefficient_polynomials = {
        7: Rational(8, 49) * (n**2 - 36 * n - 147),
        6: -(
            n**4 - 54 * n**3 - 2077 * n**2 - 9102 * n - 69984
        ) / Rational(324),
        5: -Rational(2, 15625) * (
            4 * n**6 + 35 * n**5 + 2065 * n**4 + 36375 * n**3
            + 190771 * n**2 + 858250 * n + 6187500
        ),
        4: (
            9 * n**8 - 148 * n**7 + 2986 * n**6 + 1296 * n**5
            + 180161 * n**4 + 906548 * n**3 + 3124364 * n**2
            + 24620624 * n + 102236160
        ) / Rational(65536),
        3: -Rational(2, 59049) * (
            n**10 - 30 * n**9 + 474 * n**8 - 2781 * n**7
            + 11712 * n**6 + 46719 * n**5 - 67789 * n**4
            + 1269696 * n**3 + 2940822 * n**2
            + 23751036 * n + 53144100
        ),
        2: (
            n**10 - 21 * n**9 + 290 * n**8 - 2394 * n**7
            + 14193 * n**6 - 47813 * n**5 + 108620 * n**4
            - 44556 * n**3 + 340768 * n**2 + 982592 * n
            + 1253376
        ) / Rational(1024),
        1: -2 * (
            3 * n**7 - 41 * n**6 + 266 * n**5 - 882 * n**4
            + 1845 * n**3 - 1565 * n**2 + 1018 * n + 228
        ),
    }
    assert all(
        simplify(
            coefficient_polynomials[base]
            - expected_coefficient_polynomials[base]
        ) == 0
        for base in range(1, 8)
    )

    full_exponential_polynomial = sum(
        polynomial * exp(power * s)
        for power, polynomial in exponential_polynomials.items()
    )
    assert series(
        full_exponential_polynomial, s, 0, 18
    ).removeO().expand() == 0

    def scaled_coefficient(order):
        return simplify(sum(
            base**order * coefficient_polynomials[base].subs(n, order)
            for base in range(1, 8)
        ))

    initial_coefficients = [
        scaled_coefficient(order) for order in range(18, 41)
    ]
    assert min(initial_coefficients) == 61751289600

    def has_positive_shift_coefficients(polynomial, start):
        offset = symbols("offset", nonnegative=True)
        shifted = Poly(polynomial.subs(n, offset + start), offset)
        return all(coefficient > 0 for coefficient in shifted.all_coeffs())

    negative_parts = {
        base: -coefficient_polynomials[base]
        for base in (1, 3, 5, 6)
    }
    for base in (1, 3, 5):
        assert has_positive_shift_coefficients(negative_parts[base], 18)
    assert has_positive_shift_coefficients(coefficient_polynomials[2], 18)
    assert has_positive_shift_coefficients(coefficient_polynomials[4], 18)
    assert has_positive_shift_coefficients(coefficient_polynomials[7], 40)
    assert has_positive_shift_coefficients(negative_parts[6], 82)
    assert all(
        coefficient_polynomials[6].subs(n, order) > 0
        for order in range(41, 82)
    )

    # For a negative base-k contribution, the following polynomial is
    # positive exactly when its ratio to the positive base-7 contribution
    # decreases from n to n+1.
    for base in (1, 3, 5):
        ratio_difference = (
            7 * negative_parts[base]
            * coefficient_polynomials[7].subs(n, n + 1)
            - base * negative_parts[base].subs(n, n + 1)
            * coefficient_polynomials[7]
        )
        assert has_positive_shift_coefficients(ratio_difference, 41)
    ratio_difference_six = (
        7 * negative_parts[6]
        * coefficient_polynomials[7].subs(n, n + 1)
        - 6 * negative_parts[6].subs(n, n + 1)
        * coefficient_polynomials[7]
    )
    assert has_positive_shift_coefficients(ratio_difference_six, 88)

    ratio_135_at_41 = sum(
        negative_parts[base].subs(n, 41) * base**41
        for base in (1, 3, 5)
    ) / (coefficient_polynomials[7].subs(n, 41) * 7**41)
    assert ratio_135_at_41 < Rational(1, 2)
    for order in range(82, 89):
        ratio_six = (
            negative_parts[6].subs(n, order) * 6**order
            / (coefficient_polynomials[7].subs(n, order) * 7**order)
        )
        assert 0 < ratio_six < Rational(1, 10000)


def check_three_exponential_two_large_gap_region() -> None:
    """Check the exact constants in the analytic two-large-gap bound."""

    # Check the coefficient proof of
    # integral_0^1 r^2 exp(-d r) dr <= 4/(1+d)^3.
    d = symbols("d", positive=True)
    hfun = (
        d**5 + 5 * d**4 + 11 * d**3 + 13 * d**2 + 8 * d + 2
        + (2 * d**3 - 6 * d**2 - 6 * d - 2) * exp(d)
    )
    hseries = series(hfun, d, 0, 20).removeO().expand()
    assert all(hseries.coeff(d, order) == 0 for order in range(3))
    assert [
        factorial(order) * hseries.coeff(d, order)
        for order in range(3, 6)
    ] == [22, 70, 88]
    for order in range(6, 20):
        scaled = factorial(order) * hseries.coeff(d, order)
        expected = 2 * order**3 - 12 * order**2 + 4 * order - 2
        assert scaled == expected
        assert expected > 0

    # The full-quadrant complement polynomial has the documented value at
    # z=13.  The two boundary/bulk coordinate comparisons contribute the
    # constants 16 and 16+8=24 in the common trace bound.
    complement_polynomial = 13**3 + 3 * 13**2 + 8 * 13 + 8
    assert complement_polynomial == 2816
    assert 4 * 4 == 16
    assert 16 + 8 == 24

    # e > 163/60 by its series through degree five.  These exact integer
    # inequalities respectively prove theta(13)<1/50 and the final trace
    # bound at z=13 is less than 8/5.
    exponential_partial_sum = sum(
        Fraction(1, factorial(order)) for order in range(6)
    )
    assert exponential_partial_sum == Fraction(163, 60)
    assert exponential_partial_sum > Fraction(19, 7)
    assert 163**13 > 100 * complement_polynomial * 60**13
    assert 49 * 163**13 > 750 * 13**4 * 60**13
    assert 19**7 > 4**5 * 7**7


def check_three_exponential_small_gap_region() -> None:
    """Check the fixed-sum identities and constants in the small-gap theorem."""

    from sympy import cosh, sinh

    p, q = symbols("p q", positive=True)

    def bfun(value):
        return (1 - exp(-value)) / value

    def first_moment(value):
        return (1 - (1 + value) * exp(-value)) / value**2

    bminus = bfun(p - q)
    bplus = bfun(p + q)
    fixed_sum_a = bminus + bplus - exp(-p) * sinh(q) / q
    expected_derivative = (
        first_moment(p - q) - first_moment(p + q)
        - exp(-p) * (q * cosh(q) - sinh(q)) / q**2
    )
    assert simplify(diff(fixed_sum_a, q) - expected_derivative) == 0
    assert simplify(
        diff(q * cosh(q) - sinh(q), q) - q * sinh(q)
    ) == 0

    radius = Fraction(4, 9)
    exponential_lower = sum(
        radius**order / Fraction(factorial(order))
        for order in range(6)
    )
    assert exponential_lower == Fraction(1381403, 885735)
    assert exponential_lower > Fraction(2500, 1603)

    # For k>=7, k! >= 7! 8^(k-7), so the omitted exponential tail is at
    # most (1/7!) sum_{j>=0} 8^(-j) = 1/4410.
    exponential_upper = (
        sum(Fraction(1, factorial(order)) for order in range(7))
        + Fraction(1, 4410)
    )
    assert exponential_upper == Fraction(31967, 11760)
    assert exponential_upper < Fraction(271829, 100000)

    tail_polynomial_lower = (
        Fraction(100) - 139 * Fraction(1603, 2500)
    )
    assert tail_polynomial_lower == Fraction(27183, 2500)
    assert 4 * Fraction(271829, 100000) == Fraction(271829, 25000)
    assert tail_polynomial_lower > Fraction(271829, 25000)


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


def check_three_exponential_infinite_gap_boundary() -> None:
    """Check the reduction of the infinite-gap edge to two exponentials."""

    z, u, w = symbols("z u w", positive=True)
    ez = exp(z)
    bfun = (1 - exp(-z)) / z
    moments = [
        integrate(u**j * exp(-z * u), (u, 0, 1)) / bfun
        for j in range(1, 4)
    ]
    mean = moments[0]
    variance = moments[1] - mean**2
    weighted_centered_square = (
        moments[2] - 2 * mean * moments[1] + mean**2 * moments[0]
    )
    rho_u = simplify(1 + z * weighted_centered_square / variance)
    rho_y = simplify(3 + z * mean)
    h_two = (
        z * (z * ez - ez + 1)**2
        / ((ez - 1) * ((ez - 1)**2 - z**2 * ez))
    )
    assert simplify(rho_u - (4 - h_two)) == 0
    assert simplify(mean - (1 / z - 1 / (ez - 1))) == 0
    expected_gap = (
        -z**2 * ez * (z * ez + z - 2 * ez + 2)
        / ((ez - 1) * (z**2 * ez - (ez - 1)**2))
    )
    assert simplify(rho_y - rho_u - expected_gap) == 0

    b_w = (1 - exp(-w)) / w
    partition = (b_w - bfun) / (z - w)
    tail_linear = bfun + z * partition
    assert simplify(limit(w * partition, w, oo) - bfun) == 0
    assert simplify(limit(tail_linear, w, oo) - bfun) == 0


def check_three_exponential_sharp_corner_expansion() -> None:
    """Check the joint zero-gap/infinite-gap margin expansion."""

    z, t, eps, a, b = symbols("z t eps a b", positive=True)
    # Terms omitted from B(1/t)=t*(1-exp(-1/t)) are flat at t=0.  Eight
    # terms of B(z) are more than enough for every total degree checked.
    bpoly = sum((-z)**k / factorial(k + 1) for k in range(9))
    partition = t * (bpoly - t) / (1 - z * t)
    log_partition = log(t) + log(bpoly - t) - log(1 - z * t)

    c11 = diff(log_partition, z, 2)
    c12 = -t**2 * diff(log_partition, z, t)
    c22 = (
        t**4 * diff(log_partition, t, 2)
        + 2 * t**3 * diff(log_partition, t)
    )
    k11 = partition * c11
    k12 = partition * c12
    k22 = partition * c22

    def radial_derivative(value):
        return z * diff(value, z) - t * diff(value, t)

    h11 = -radial_derivative(k11) / partition
    h12 = -radial_derivative(k12) / (t * partition)
    h22 = -radial_derivative(k22) / (t**2 * partition)
    # Congruence by diag(1,1/t) converts (U,V) to (U,wV).
    covariance = (c11, c12 / t, c22 / t**2)
    weighted_covariance = (h11, h12, h22)

    def scaled_series(value, order=5):
        return series(
            value.subs({z: eps * a, t: eps * b}), eps, 0, order,
        ).removeO().expand()

    c_scaled = tuple(scaled_series(value) for value in covariance)
    h_scaled = tuple(
        scaled_series(value) for value in weighted_covariance
    )
    assert [value.subs(eps, 0) for value in c_scaled] == [
        Rational(1, 12), 0, 1,
    ]
    assert [value.subs(eps, 0) for value in h_scaled] == [
        Rational(1, 12), 0, 3,
    ]

    rho = (
        1 + eps * (a / 2 - 3 * b)
        + eps**2 * (a**2 / 60 + 3 * b**2)
        + eps**3 * (
            -a**2 * b / 60 - Rational(3, 2) * a * b**2
            + 9 * b**3
        )
    )
    characteristic = (
        (h_scaled[0] - rho * c_scaled[0])
        * (h_scaled[2] - rho * c_scaled[2])
        - (h_scaled[1] - rho * c_scaled[1])**2
    )
    assert series(characteristic, eps, 0, 4).removeO().expand() == 0

    z_scaled = eps * a
    t_scaled = eps * b
    b_scaled = sum((-z_scaled)**k / factorial(k + 1)
                   for k in range(9))
    partition_scaled = (
        t_scaled * (b_scaled - t_scaled) / (1 - z_scaled * t_scaled)
    )
    tail_linear = b_scaled + z_scaled * partition_scaled
    h_threshold = 4 - rho
    margin = (
        3 - h_threshold
        + log(1 + tail_linear * h_threshold
              + partition_scaled * h_threshold**2)
        - log(4)
    )
    expected_margin = (
        Rational(3, 40) * a**2 * eps**2
        + (-Rational(3, 20) * a**2 * b
           + Rational(9, 2) * b**3) * eps**3
    )
    assert simplify(
        series(margin, eps, 0, 4).removeO().expand() - expected_margin
    ) == 0


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
    check_three_exponential_trace_identity()
    check_three_exponential_repeated_max_boundary()
    check_three_exponential_axis_transversality()
    check_three_exponential_equal_smaller_line()
    check_three_exponential_diagonal_transversality()
    check_three_exponential_small_gap_region()
    check_three_exponential_two_large_gap_region()
    check_three_exponential_infinite_gap_boundary()
    check_three_exponential_sharp_corner_expansion()
    check_kernel_identity()
    poisson, sympol, gap, active_k = sympol_shortcut_counterexample()
    print("all-n reduction algebra: PASS")
    print("exact two-exponential convexity threshold: 4*exp(-3) =",
          mp.nstr(4 * mp.e**-3, 16))
    print("equal-weight local thresholds: strictly increasing from the "
          "same m=2 value")
    print("three-exponential convexity: reduced exactly to the documented "
          "two-variable inequality")
    print("three-exponential boundary-trace identity: PROVED")
    print("three-exponential repeated-max boundary: PROVED")
    print("three-exponential axis transversality: STRICTLY POSITIVE")
    print("three-exponential equal-smaller symmetry line: PROVED")
    print("three-exponential diagonal transverse second derivative: NEGATIVE")
    print("three-exponential total-gap-at-most-8/9 region: PROVED")
    print("three-exponential min-gap-at-least-13 region: PROVED")
    print("three-exponential infinite-gap boundary: REDUCED TO PROVED 2D CASE")
    print("three-exponential sharp corner: POSITIVE PUNCTURED NEIGHBORHOOD")
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
