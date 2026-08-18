"""Exact checks for the open Dirichlet--Poissonization comparison.

The all-sample-size Poisson program compares

    T_y = sum_i y_i D_i,       Z_y = sum_i y_i E_i,

where ``D`` is uniform on the ``n``-simplex and the ``E_i`` are independent
unit exponentials.  At the normalized threshold used in the research note,
the desired inequality is

    P(Z_y > n) >= P(T_y > 1).

This module does *not* claim that inequality in full.  It records a proved
boundary reduction, a proved structured family, and one obstruction:

* a radial deformation moves any strictly positive profile to a zero-knot
  boundary while keeping the simplex tail fixed and decreasing the
  poissonized tail; and
* the comparison is true for every vector having at most two distinct
  coefficient values.  Its radial endpoint follows from the published
  Anderson--Samuels binomial--Poisson inequality.  This module checks the
  exact algebra and rational instances both on and strictly inside the sum
  boundary; and
* the radial reduction and a one-variable divided-difference argument prove
  the comparison in every dimension when at most two coefficients are
  nonzero, and therefore the complete comparison for ``n=2``; and
* a constrained three-knot secant argument proves the complete comparison
  for ``n=3``; and
* exact scalar bounds plus deterministic real-ball branch certificates prove
  the complete comparison for ``n=4``; and
* a second deterministic branch certificate proves every ``n=5`` profile
  having at most four nonzero coefficients; and
* a recursive knot-insertion bound, exact scalar estimates, and a third
  deterministic branch certificate prove the complete comparison for
  ``n=5``; and
* the ``k=2`` Anderson--Samuels comparison proves a dimension-free convex
  core for profiles with three nonzero coefficients; and
* the same derivative-CDF identity at arbitrary order proves a sparse
  convex core on every coordinate face, including a first all-``n`` region
  for profiles with four nonzero coefficients; and
* an exact monotonicity and endpoint argument proves a dimension-free far
  cap for profiles with three nonzero coefficients; and
* a weighted divided-difference argument and exact endpoint bounds prove a
  dimension-free far cap for profiles with four nonzero coefficients; and
* a one-crossing curvature argument, a quantitative curvature inequality,
  and exact endpoint comparisons prove the full region where the middle knot
  is at most ``n/3``; and
* a derivative comparison, exact analytic tail bounds, and a rigorous Arb
  check of the finite prefix complete every profile with three nonzero
  coefficients; and
* a fully explicit, mean-constrained ``1/n``-concave law shows why generic
  one-dimensional s-concave localization is too broad to prove the
  comparison.

The latter check is formal.  Its probability above one is rational, and its
gamma-smoothed probability is a rational multiple of ``exp(-750/97)``.  The
exponential is enclosed by exact rational alternating-series bounds from
``stringer.exact_exp_neg_bounds``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from stringer import EXACT_EXP_PAIRS, exact_exp_neg_bounds


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "dirichlet-poissonization-certificate.json"
)


def verify_saffine_symbolic_identities() -> None:
    """Derive the two closed forms used in the localization obstruction."""
    from sympy import Rational, diff, exp, factorial, integrate, simplify, symbols

    n = 10
    y = symbols("y", positive=True)
    upper = Rational(97, 75)
    raw_density = (Rational(3, 4) + Rational(75, 388) * y) ** (n - 1)
    normalizer = integrate(raw_density, (y, 0, upper))
    assert normalizer == Rational(95984119, 196608000)
    density = raw_density / normalizer

    mean = simplify(integrate(y * density, (y, 0, upper)))
    assert mean == Rational(729166363, 816359775)
    assert mean < Rational(10, 11)

    tail = simplify(integrate(density, (y, 1, upper)))
    assert tail == Rational(
        3108309643939756140704768,
        6633646218308706152889893,
    )

    # H is C^(n-1) at one.  Integrating density*H^(n)/n! by parts n
    # times therefore leaves this endpoint expression; every lower-endpoint
    # term vanishes because y^n exp(-n/y) is flat at zero.
    h_upper_piece = y**n * exp(-Rational(n, 1) / y) - (y - 1) ** n
    boundary = sum(
        (-1) ** order
        * diff(density, y, order)
        * diff(h_upper_piece, y, n - 1 - order)
        for order in range(n)
    )
    smoothed_minus_tail = simplify(
        boundary.subs(y, upper) / factorial(n)
    )
    expected_smoothed = (
        Rational(
            5573507995079350591862317513,
            5265884111440931688376513,
        )
        * exp(-Rational(750, 97))
    )
    assert simplify(tail + smoothed_minus_tail - expected_smoothed) == 0


def _decimal(value: Fraction, digits: int = 18) -> str:
    value = Fraction(value)
    with localcontext() as context:
        context.prec = digits + 30
        result = Decimal(value.numerator) / Decimal(value.denominator)
        return format(result, f".{digits}f")


def _fraction_record(value: Fraction, digits: int = 18) -> dict[str, str]:
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": _decimal(value, digits),
    }


def _large_fraction_record(value: Fraction, digits: int = 18) -> dict[str, object]:
    """Record a huge exact rational compactly but reproducibly.

    Alternating-series bounds followed by repeated squaring can have tens of
    thousands of decimal digits.  Embedding each integer in JSON would make
    the research certificate needlessly large.  The source regenerates the
    exact value; this record preserves its bit lengths, a canonical binary
    SHA-256 digest, and a readable decimal summary.
    """
    value = Fraction(value)
    digest = hashlib.sha256()
    for integer in (value.numerator, value.denominator):
        sign = b"-" if integer < 0 else b"+"
        magnitude = abs(integer)
        encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
        digest.update(sign)
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return {
        "numerator_bits": abs(value.numerator).bit_length(),
        "denominator_bits": value.denominator.bit_length(),
        "sha256_canonical_binary": digest.hexdigest(),
        "decimal": _decimal(value, digits),
    }


def _binomial_lower_tail(n: int, j: int, p: Fraction) -> Fraction:
    """Exact ``P(Bin(n,p) <= j)`` for rational ``p``."""
    if not 0 <= j <= n:
        raise ValueError("j must lie in {0,...,n}")
    if not 0 <= p <= 1:
        raise ValueError("p must lie in [0,1]")
    return sum(
        Fraction(math.comb(n, ell))
        * p**ell
        * (1 - p) ** (n - ell)
        for ell in range(j + 1)
    )


def equal_block_check(n: int, k: int, lam: Fraction) -> dict[str, object]:
    """Certify one rational equal-block binomial--Poisson comparison.

    The coefficient vector has ``k`` entries equal to ``n/lam`` and
    ``n+1-k`` zero entries.  The condition ``lam>=k`` is exactly the
    coefficient-sum condition.  The Dirichlet tail is a binomial lower tail;
    the exponential tail is a Poisson lower tail.

    This function is only a regression check.  The proof for all real
    ``lam>=k`` is the Anderson--Samuels theorem plus continuity at ``lam=k``.
    """
    lam = Fraction(lam)
    if n < 1 or not 1 <= k <= n:
        raise ValueError("require n>=1 and 1<=k<=n")
    if not Fraction(k) <= lam <= n:
        raise ValueError("require k<=lam<=n")

    p = lam / n
    dirichlet_tail = _binomial_lower_tail(n, k - 1, p)
    polynomial = sum(lam**ell / math.factorial(ell) for ell in range(k))
    exp_lower, exp_upper = exact_exp_neg_bounds(lam, EXACT_EXP_PAIRS)
    poisson_lower = polynomial * exp_lower
    poisson_upper = polynomial * exp_upper
    if poisson_lower < dirichlet_tail:
        raise AssertionError("equal-block comparison was not certified")
    return {
        "n": n,
        "positive_coordinates": k,
        "lambda": _fraction_record(lam),
        "common_positive_coordinate": _fraction_record(Fraction(n, 1) / lam),
        "dirichlet_tail": _fraction_record(dirichlet_tail),
        "poisson_tail_lower_bound": _large_fraction_record(poisson_lower),
        "poisson_tail_upper_bound": _large_fraction_record(poisson_upper),
        "certified_lower_margin": _large_fraction_record(
            poisson_lower - dirichlet_tail
        ),
    }


def verify_radial_symbolic_identities() -> None:
    """Check algebra behind the zero-knot radial reduction.

    The written proof also uses log-concavity and the elementary mode bound
    derived from the convolution identity whose Laplace-transform form is
    checked here.
    """
    from sympy import diff, factor, prod, symbols, together

    t, s = symbols("t s", positive=True)
    c = symbols("c0:4", real=True)
    e = symbols("e0:4", positive=True)
    weights = [1 + t * value for value in c]
    x_t = sum(weight * value for weight, value in zip(weights, e))
    total = sum(e)
    direction = sum(value * exponential for value, exponential in zip(c, e))
    assert factor(x_t - total - t * direction) == 0
    assert factor(sum(weights) - (4 + t * sum(c))) == 0

    # If L is the Laplace transform of a weighted exponential sum, then
    # -L' = L*K, where K is the transform of sum_i exp(-u/a_i).  Inversion
    # gives x*f(x)=(f*k)(x), the identity used to put the mode below the mean.
    a = symbols("a0:4", positive=True)
    laplace = prod(1 / (1 + value * s) for value in a)
    kernel_transform = sum(value / (1 + value * s) for value in a)
    assert factor(together(-diff(laplace, s) - laplace * kernel_transform)) == 0


def verify_two_level_symbolic_identities() -> None:
    """Check the two-level endpoint algebra after radial deformation."""
    from sympy import factor, symbols, together

    a, b, k, ell = symbols("a b k ell", positive=True)
    m = k + ell
    t_star = 1 / (1 - a)
    b_star = 1 + t_star * (b - 1)
    q = (1 - a) / (b - a)
    original_sum = k * a + ell * b

    assert factor(together(b_star - (b - a) / (1 - a))) == 0
    assert factor(together(b_star - 1 / q)) == 0
    endpoint_sum_from_path = m + t_star * (original_sum - m)
    assert factor(together(endpoint_sum_from_path - ell * b_star)) == 0


def verify_two_positive_knot_theorem() -> None:
    """Check the algebra in the all-``n`` two-positive-knot proof."""
    from sympy import diff, exp, factor, log, powsimp, simplify, symbols

    n = symbols("n", integer=True, positive=True)
    u, x, r = symbols("u x r", positive=True)
    upper_piece = u * (exp(-n / u) - (1 - 1 / u) ** n)
    expected_derivative = (
        exp(-n / u) * (1 + n / u)
        - (1 - 1 / u) ** (n - 1) * (1 + (n - 1) / u)
    )
    assert simplify(diff(upper_piece, u) - expected_derivative) == 0

    ratio = (
        exp(-n * x) * (1 + n * x)
        / ((1 - x) ** (n - 1) * (1 + (n - 1) * x))
    )
    expected_log_derivative = (
        -n * x * (n * (n - 1) * x**2 + n * x - 1)
        / ((x - 1) * (1 + n * x) * (1 + (n - 1) * x))
    )
    assert powsimp(
        factor(diff(log(ratio), x) - expected_log_derivative),
        force=True,
    ) == 0

    # The lower bounds in the proof of log R_n(1/(n-1))>0 leave exactly
    # 5r^2/24.
    log_ratio_lower = (
        -1 - r + (r / 2 - r**2 / 8) + 1 + r / 2 + r**2 / 3
    )
    assert simplify(log_ratio_lower - 5 * r**2 / 24) == 0

    # The termwise logarithmic-series comparison uses j*2^j>=j+1.  The
    # base case is equality and this recurrence makes the induction strict.
    j = symbols("j", integer=True, positive=True)
    sequence = j * 2**j - j - 1
    assert sequence.subs(j, 1) == 0
    assert simplify(sequence.subs(j, j + 1) - 2 * sequence - (2 ** (j + 1) + j)) == 0

    # The final endpoint comparison only needs e>2.
    e_lower = sum(Fraction(1, math.factorial(order)) for order in range(3))
    assert e_lower > 2


def verify_n2_three_knot_theorem() -> None:
    """Check the complete ``n=2`` corollary through its stronger theorem."""
    verify_two_positive_knot_theorem()


def _exp_taylor_lower(x: Fraction, order: int) -> Fraction:
    """Positive Taylor lower bound for ``exp(x)`` at rational ``x>=0``."""
    x = Fraction(x)
    if x < 0 or order < 0:
        raise ValueError("require x>=0 and order>=0")
    return sum(x**k / math.factorial(k) for k in range(order + 1))


def _exp_taylor_upper(x: Fraction, order: int) -> Fraction:
    """Taylor upper bound with a geometric majorant for the positive tail."""
    x = Fraction(x)
    if x < 0 or order < 0 or x >= order + 2:
        raise ValueError("require x>=0, order>=0, and x<order+2")
    partial = _exp_taylor_lower(x, order)
    first_omitted = x ** (order + 1) / math.factorial(order + 1)
    ratio = x / (order + 2)
    return partial + first_omitted / (1 - ratio)


def verify_n3_four_knot_theorem() -> None:
    """Check every exact identity and rational bound in the ``n=3`` proof."""
    from sympy import Rational, diff, exp, factor, simplify, symbols

    u = symbols("u", positive=True)
    lower_piece = u**2 * exp(-3 / u)
    upper_piece = lower_piece - (u - 1) ** 3 / u
    expected_first = (2 * u + 3) * exp(-3 / u) - (
        2 * u - 3 + u**-2
    )
    expected_second = (
        (2 + 6 / u + 9 / u**2) * exp(-3 / u) - 2 + 2 / u**3
    )
    expected_third = (27 * exp(-3 / u) - 6) / u**4
    h_upper_piece = u**3 * exp(-3 / u) - (u - 1) ** 3
    assert simplify(h_upper_piece - u * upper_piece) == 0
    assert simplify(diff(upper_piece, u) - expected_first) == 0
    assert simplify(diff(upper_piece, u, 2) - expected_second) == 0
    assert simplify(diff(upper_piece, u, 3) - expected_third) == 0
    assert factor(diff(lower_piece, u, 2)) == (
        (2 * u**2 + 6 * u + 9) * exp(-3 / u) / u**2
    )
    assert simplify(
        expected_second.subs(u, Rational(7, 5))
        - (
            Rational(533, 49) * exp(-Rational(15, 7))
            - Rational(436, 343)
        )
    ) == 0
    assert simplify(expected_second.subs(u, 1) - 17 * exp(-3)) == 0
    assert simplify(
        expected_second.subs(u, 3) - (5 * exp(-1) - Rational(52, 27))
    ) == 0
    assert simplify(
        diff(lower_piece, u).subs(u, Rational(9, 10))
        - Rational(24, 5) * exp(-Rational(10, 3))
    ) == 0
    assert simplify(
        expected_first.subs(u, Rational(3, 2))
        - (6 * exp(-2) - Rational(4, 9))
    ) == 0
    assert simplify(
        expected_first.subs(u, 3) - (9 * exp(-1) - Rational(28, 9))
    ) == 0

    # Taylor bounds used to certify strict convexity through 7/5.
    exp_15_7_lower = _exp_taylor_lower(Fraction(15, 7), 2)
    exp_15_7_upper = _exp_taylor_upper(Fraction(15, 7), 6)
    assert exp_15_7_lower > Fraction(9, 2)
    assert exp_15_7_upper == Fraction(4605295343, 540244208)
    assert exp_15_7_upper < Fraction(341, 40)
    f_second_7_5_lower = (
        Fraction(533, 49) * Fraction(40, 341) - Fraction(436, 343)
    )
    assert f_second_7_5_lower == Fraction(564, 116963)
    assert f_second_7_5_lower > 0

    # Exact sandwich for exp(10/3), used at the 9/10 slope split.
    exp_10_3_lower = _exp_taylor_lower(Fraction(10, 3), 7)
    exp_10_3_upper = _exp_taylor_upper(Fraction(10, 3), 7)
    assert exp_10_3_lower == Fraction(3781751, 137781)
    assert exp_10_3_lower > Fraction(192, 7)
    assert exp_10_3_upper == Fraction(65696017, 2342277)
    assert exp_10_3_upper < Fraction(144, 5)

    # The e-bounds control f'(1), f'(3/2), f'(3), and Q(3/2).
    e_lower = _exp_taylor_lower(Fraction(1), 7)
    e_upper = _exp_taylor_upper(Fraction(1), 6)
    assert e_lower == Fraction(685, 252)
    assert e_lower > Fraction(405, 149)
    assert e_lower > Fraction(8, 3)
    assert e_upper == Fraction(31967, 11760)
    assert e_upper < Fraction(68, 25)
    assert e_upper < Fraction(11, 4)

    # Translate the exponential bounds into the displayed derivative chain.
    assert Fraction(24, 5) * Fraction(5, 144) == Fraction(1, 6)
    assert Fraction(24, 5) * Fraction(7, 192) == Fraction(7, 40)
    f_prime_3_lower = 9 * Fraction(25, 68) - Fraction(28, 9)
    assert f_prime_3_lower == Fraction(121, 612)
    assert f_prime_3_lower > Fraction(7, 40)
    f_prime_3_upper = 9 * Fraction(149, 405) - Fraction(28, 9)
    assert f_prime_3_upper == Fraction(1, 5)
    f_prime_1_lower = 5 * Fraction(4, 11) ** 3
    assert f_prime_1_lower > Fraction(1, 5)
    f_prime_3_2_lower = 6 * Fraction(4, 11) ** 2 - Fraction(4, 9)
    assert f_prime_3_2_lower == Fraction(380, 1089)
    assert f_prime_3_2_lower > Fraction(1, 3)

    # Endpoint checks for Q(u)=u/3-4/15-f(u).
    f_9_10_upper = Fraction(81, 100) * Fraction(7, 192)
    assert f_9_10_upper == Fraction(189, 6400)
    assert Fraction(1, 30) - f_9_10_upper > 0
    assert Fraction(19, 60) - Fraction(81, 256) == Fraction(1, 3840)

    # Algebra in the final secant comparison.
    a, b, c = symbols("a b c", real=True)
    assert factor(
        (3 - b - c)
        - (Rational(8, 5) - b)
        + (c - Rational(7, 5))
    ) == 0


def verify_three_positive_convex_core() -> None:
    """Check the all-``n`` three-positive convex-core identities."""
    from sympy import binomial, exp, factor, simplify, symbols

    n = symbols("n", integer=True, positive=True)
    u = symbols("u", positive=True)
    lam = n / u
    lower_piece = u**2 * exp(-n / u)
    upper_piece = lower_piece - u**2 * (1 - 1 / u) ** n
    poisson_lower_two = exp(-lam) * (1 + lam + lam**2 / 2)
    binomial_lower_two = (
        (1 - 1 / u) ** n
        + n / u * (1 - 1 / u) ** (n - 1)
        + binomial(n, 2) / u**2 * (1 - 1 / u) ** (n - 2)
    )
    assert simplify(
        upper_piece.diff(u, 2) / 2
        - (poisson_lower_two - binomial_lower_two)
    ) == 0
    assert factor(lower_piece.diff(u, 2)) == (
        (2 * u**2 + 2 * n * u + n**2) * exp(-n / u) / u**2
    )

    c_n = n**2 / (2 * (n + 1))
    lambda_at_boundary = simplify(n / c_n)
    assert simplify(lambda_at_boundary - 2 * (n + 1) / n) == 0
    assert simplify(
        lambda_at_boundary
        - lambda_at_boundary / (n + 1)
        - 2
    ) == 0


def verify_sparse_convex_core() -> dict[str, object]:
    """Check the arbitrary-order sparse convex-core identities.

    The written proof uses induction in ``k``.  These exact SymPy checks
    cover symbolic ``n`` through the first six derivative orders, including
    the new four-positive case ``k=3``, and every admissible order for the
    integer prefix ``2<=n<=10``.  They corroborate rather than replace the
    all-parameter induction and Anderson--Samuels argument in the note.
    """
    from sympy import (
        binomial,
        diff,
        exp,
        factorial,
        simplify,
        symbols,
        together,
    )

    n = symbols("n", integer=True, positive=True)
    u = symbols("u", positive=True)

    symbolic_orders: list[int] = []
    for order in range(1, 7):
        poisson_piece = diff(
            u**order * exp(-n / u),
            u,
            order,
        ) / factorial(order)
        poisson_cdf = exp(-n / u) * sum(
            (n / u) ** j / factorial(j)
            for j in range(order + 1)
        )
        assert simplify(poisson_piece - poisson_cdf) == 0

        binomial_piece = diff(
            u**order * (1 - 1 / u) ** n,
            u,
            order,
        ) / factorial(order)
        binomial_cdf = sum(
            binomial(n, j)
            * u ** (-j)
            * (1 - 1 / u) ** (n - j)
            for j in range(order + 1)
        )
        assert simplify(together(binomial_piece - binomial_cdf)) == 0
        symbolic_orders.append(order)

    integer_pairs: list[str] = []
    for dimension in range(2, 11):
        for order in range(1, dimension):
            poisson_piece = diff(
                u**order * exp(-dimension / u),
                u,
                order,
            ) / factorial(order)
            poisson_cdf = exp(-dimension / u) * sum(
                (dimension / u) ** j / factorial(j)
                for j in range(order + 1)
            )
            assert simplify(poisson_piece - poisson_cdf) == 0

            binomial_piece = diff(
                u**order * (1 - 1 / u) ** dimension,
                u,
                order,
            ) / factorial(order)
            binomial_cdf = sum(
                math.comb(dimension, j)
                * u ** (-j)
                * (1 - 1 / u) ** (dimension - j)
                for j in range(order + 1)
            )
            assert simplify(together(binomial_piece - binomial_cdf)) == 0
            integer_pairs.append(f"n={dimension},k={order}")

    k = symbols("k", integer=True, positive=True)
    core_endpoint = n**2 / (k * (n + 1))
    lambda_at_boundary = simplify(n / core_endpoint)
    assert simplify(
        lambda_at_boundary
        - lambda_at_boundary / (n + 1)
        - k
    ) == 0

    return {
        "symbolic_n_derivative_orders_checked": symbolic_orders,
        "integer_parameter_pairs_checked": integer_pairs,
        "integer_parameter_pair_count": len(integer_pairs),
        "threshold_identity": (
            "At u=n^2/(k(n+1)), lambda-lambda/(n+1)=k"
        ),
    }


def verify_three_positive_far_cap() -> None:
    """Check the all-``n`` three-positive far-cap identities and bounds."""
    from sympy import exp, factor, log, powsimp, simplify, symbols

    n = symbols("n", integer=True, positive=True)
    u, lam = symbols("u lambda", positive=True)
    upper_piece = u**2 * (
        exp(-n / u) - (1 - 1 / u) ** n
    )
    scaled_derivative = (
        exp(-lam) * (2 + lam)
        - (1 - lam / n) ** (n - 1)
        * (2 + lam - 2 * lam / n)
    )
    assert simplify(
        upper_piece.diff(u) / u
        - scaled_derivative.subs(lam, n / u)
    ) == 0

    ratio = (
        exp(-lam) * (2 + lam)
        / (
            (1 - lam / n) ** (n - 1)
            * (2 + lam - 2 * lam / n)
        )
    )
    expected_log_derivative = (
        -lam**2 * (lam * n - 2 * lam + 3 * n - 2)
        / (
            (lam + 2)
            * (lam - n)
            * (lam * n - 2 * lam + 2 * n)
        )
    )
    assert powsimp(
        factor(log(ratio).diff(lam) - expected_log_derivative),
        force=True,
    ) == 0

    # The positive logarithmic-series terms at lambda=1 leave a strictly
    # positive rational lower bound for every n>=4.
    log_ratio_one_lower = (
        -1
        + (n - 1) * (1 / n + 1 / (2 * n**2))
        + 2 / (3 * n)
    )
    assert simplify(
        log_ratio_one_lower - (n - 3) / (6 * n**2)
    ) == 0

    # Exact algebra in the general endpoint estimate.
    x_plus_y = n / (n - 1) + n / (2 * (n - 1) ** 2)
    assert simplify(
        x_plus_y
        - (1 + 3 / (2 * (n - 1)) + 1 / (2 * (n - 1) ** 2))
    ) == 0

    # At n=4, e<11/4 implies exp(4/3)<4, while the degree-nine
    # positive Taylor sum gives exp(4)>54.
    e_upper = _exp_taylor_upper(Fraction(1), 6)
    exp_four_lower = _exp_taylor_lower(Fraction(4), 9)
    assert e_upper < Fraction(11, 4)
    assert Fraction(11, 4) ** 4 < 64
    assert exp_four_lower == Fraction(153527, 2835)
    assert exp_four_lower > 54
    assert Fraction(17, 36) > Fraction(25, 54)

    # The induction proving the endpoint estimate for n>=5 starts at this
    # exact degree-four Taylor bound.  The difference below proves that the
    # ratio of consecutive right sides is less than two.
    exp_three_lower = _exp_taylor_lower(Fraction(3), 4)
    assert exp_three_lower == Fraction(131, 8)
    assert exp_three_lower > Fraction(72, 5)
    assert factor(
        2 * (n + 1) ** 3 - n * (n + 2) ** 2
    ) == n**3 + 2 * n**2 + 2 * n + 2


def verify_four_positive_far_cap() -> dict[str, object]:
    """Check the all-``n`` four-positive far-cap proof and constants."""
    from sympy import Rational, factor, simplify, symbols

    a, b, c, d = symbols("a b c d", real=True)
    ga, gb, gc, gd = symbols("g_a g_b g_c g_d", real=True)

    def first_dd(x0, x1, y0, y1):
        return (y1 - y0) / (x1 - x0)

    def second_dd(x0, x1, x2, y0, y1, y2):
        return (
            first_dd(x1, x2, y1, y2)
            - first_dd(x0, x1, y0, y1)
        ) / (x2 - x0)

    q_abc = second_dd(a, b, c, ga, gb, gc)
    q_bcd = second_dd(b, c, d, gb, gc, gd)
    f_abc = second_dd(a, b, c, a * ga, b * gb, c * gc)
    f_bcd = second_dd(b, c, d, b * gb, c * gc, d * gd)
    third_f = (f_bcd - f_abc) / (d - a)
    assert simplify(
        third_f - (d * q_bcd - a * q_abc) / (d - a)
    ) == 0

    n = symbols("n", integer=True, positive=True)
    endpoint_coefficient = (
        1 / (n - 1)
        + n
        + 2
        + (1 + n + n**2 / 2) / 3
    )
    expected_coefficient = (
        n**2 / 6 + 4 * n / 3 + Rational(7, 3) + 1 / (n - 1)
    )
    assert simplify(endpoint_coefficient - expected_coefficient) == 0
    assert simplify(expected_coefficient.subs(n, 4) - Rational(32, 3)) == 0
    assert simplify(expected_coefficient.subs(n, 5) - Rational(161, 12)) == 0
    assert factor(
        2 * expected_coefficient
        - expected_coefficient.subs(n, n + 1)
    ) == (
        (n + 1) * (n**3 + 4 * n**2 - 5 * n + 6)
        / (6 * n * (n - 1))
    )

    # The n=4 endpoint margins use positive Taylor lower bounds and the
    # geometric-tail upper bound for exp(4).
    exp_three_lower = _exp_taylor_lower(Fraction(3), 8)
    exp_eight_thirds_lower = _exp_taylor_lower(Fraction(8, 3), 10)
    exp_four_upper = _exp_taylor_upper(Fraction(4), 12)
    n4_left_endpoint_margin = (
        9 * exp_eight_thirds_lower
        - Fraction(16, 9) * exp_four_upper
        - 32
    )
    n4_right_endpoint_margin = (
        4 * exp_three_lower
        - Fraction(81, 64) * exp_four_upper
        - Fraction(32, 3)
    )
    assert exp_three_lower == Fraction(89641, 4480)
    assert exp_eight_thirds_lower == Fraction(12045015679, 837019575)
    assert exp_four_upper == Fraction(553360529, 10135125)
    assert n4_left_endpoint_margin == Fraction(
        29899236229,
        66496555125,
    )
    assert n4_right_endpoint_margin == Fraction(6461863, 24024000)
    assert n4_left_endpoint_margin > 0
    assert n4_right_endpoint_margin > 0

    # For n>=5, the endpoint coefficient times exp(-n) decreases.  These
    # two Taylor checks give T_5<1/10 and
    # exp(-45/32)/2>T_5, respectively.
    exp_five_lower = _exp_taylor_lower(Fraction(5), 8)
    exp_115_over_32_lower = _exp_taylor_lower(Fraction(115, 32), 5)
    assert exp_five_lower == Fraction(1115309, 8064)
    assert exp_five_lower > Fraction(805, 6)
    assert exp_115_over_32_lower == Fraction(
        24748696103,
        805306368,
    )
    assert exp_115_over_32_lower > Fraction(161, 6)

    return {
        "n4_exp_3_degree_8_lower": _fraction_record(exp_three_lower),
        "n4_exp_8_over_3_degree_10_lower": _fraction_record(
            exp_eight_thirds_lower
        ),
        "n4_exp_4_degree_12_upper": _fraction_record(exp_four_upper),
        "n4_left_endpoint_margin": _fraction_record(
            n4_left_endpoint_margin
        ),
        "n4_right_endpoint_margin": _fraction_record(
            n4_right_endpoint_margin
        ),
        "n_ge_5_exp_5_degree_8_lower": _fraction_record(exp_five_lower),
        "n_ge_5_exp_115_over_32_degree_5_lower": _fraction_record(
            exp_115_over_32_lower
        ),
    }


def verify_n4_five_knot_theorem() -> dict[str, object]:
    """Rigorously certify the complete constrained comparison for ``n=4``.

    The proof combines exact one-variable inequalities with two deterministic
    real-ball branch certificates.  The smaller certificate proves an
    auxiliary constrained second-divided-difference inequality.  The larger
    certificate closes the remaining ordered four-knot boxes after the
    analytic convex-core, central, far-cap, and scalar-comparison regions
    have been removed.
    """
    from flint import arb, ctx, fmpq
    from sympy import (
        Poly,
        Rational,
        diff,
        exp,
        expand,
        factor,
        log,
        simplify,
        symbols,
        together,
    )

    lam = symbols("lambda", positive=True)
    u = symbols("u", positive=True)
    x = symbols("x", nonnegative=True)

    def bernstein_coefficients(poly, variable, lo, hi):
        transformed = Poly(
            expand(poly.subs(variable, lo + (hi - lo) * x)),
            x,
        )
        degree = transformed.degree()
        power = [transformed.nth(index) for index in range(degree + 1)]
        return [
            factor(
                sum(
                    power[order]
                    * Rational(
                        math.comb(index, order),
                        math.comb(degree, order),
                    )
                    for order in range(index + 1)
                )
            )
            for index in range(degree + 1)
        ]

    # Differentiate every piece used later, before applying any inequalities.
    exponential_piece = exp(-4 / u)
    f3_lower = u**3 * exponential_piece
    f3_upper = f3_lower - (u - 1) ** 4 / u
    g2_lower = u**2 * exponential_piece
    g2_upper = g2_lower - (u - 1) ** 4 / u**2
    f1_lower = u * exponential_piece
    f1_upper = f1_lower - (u - 1) ** 4 / u**3
    h0_upper = exponential_piece - (u - 1) ** 4 / u**4
    q_lower_u = exponential_piece * (
        1 + 4 / u + 8 / u**2 + Rational(32, 3) / u**3
    )
    q_upper_u = q_lower_u - 1 + u**-4
    assert simplify(diff(f3_lower, u, 3) / 6 - q_lower_u) == 0
    assert simplify(diff(f3_upper, u, 3) / 6 - q_upper_u) == 0
    assert simplify(
        diff(g2_lower, u, 2) / 2
        - exponential_piece * (1 + 4 / u + 8 / u**2)
    ) == 0
    assert simplify(
        diff(f1_lower, u, 2) / 2
        - 8 * exponential_piece / u**3
    ) == 0

    # The exact n=4 convexity endpoint.  If q=F_(4,3)'''/6, then
    # q'=4/u times the difference of the fourth Poisson and binomial masses.
    # Their ratio is (32/3)exp(-4/u), so q decreases through 19/15.
    exp_eight_thirds_lower = _exp_taylor_lower(Fraction(8, 3), 4)
    assert exp_eight_thirds_lower > Fraction(32, 3)
    core_lambda = Fraction(60, 19)
    core_poisson_polynomial = sum(
        core_lambda**j / math.factorial(j) for j in range(4)
    )
    core_binomial_tail = 1 - Fraction(15, 19) ** 4
    core_exp_threshold = core_poisson_polynomial / core_binomial_tail
    core_exp_upper = _exp_taylor_upper(core_lambda, 7)
    core_endpoint_margin = core_exp_threshold - core_exp_upper
    assert core_poisson_polynomial == Fraction(98719, 6859)
    assert core_binomial_tail == Fraction(79696, 130321)
    assert core_exp_threshold == Fraction(110333, 4688)
    assert core_exp_upper == Fraction(778209455563, 33073254343)
    assert core_endpoint_margin == Fraction(
        825443746875,
        155047416359984,
    )
    assert core_endpoint_margin > 0

    # Central affine minorant q(u)>=(1-u)/6 on [9/16,37/16].
    # Below one, q increases and the affine function decreases, so one exact
    # endpoint check suffices.
    central_u_lower = Fraction(9, 16)
    central_lambda_upper = Fraction(64, 9)
    central_poisson_polynomial = sum(
        central_lambda_upper**j / math.factorial(j) for j in range(4)
    )
    central_endpoint_rhs = (1 - central_u_lower) / 6
    central_exp_threshold = (
        central_poisson_polynomial / central_endpoint_rhs
    )
    central_exp_upper = _exp_taylor_upper(central_lambda_upper, 9)
    central_endpoint_margin = central_exp_threshold - central_exp_upper
    assert central_exp_threshold == Fraction(6531424, 5103)
    assert central_exp_upper == Fraction(
        242470438074164581,
        192208990105125,
    )
    assert central_endpoint_margin == Fraction(
        3541400253951419,
        192208990105125,
    )
    assert central_endpoint_margin > 0

    poisson_three = sum(lam**j / math.factorial(j) for j in range(4))
    poisson_four_mass = exp(-lam) * lam**4 / math.factorial(4)
    binomial_four_mass = lam**4 / 4**4
    q_lambda = exp(-lam) * poisson_three - 1 + lam**4 / 4**4
    assert simplify(q_upper_u.subs(u, 4 / lam) - q_lambda) == 0
    assert simplify(
        -lam**2 / 4 * diff(q_lambda, lam)
        - lam * (poisson_four_mass - binomial_four_mass)
    ) == 0
    assert simplify(
        poisson_four_mass / binomial_four_mass
        - Rational(32, 3) * exp(-lam)
    ) == 0
    central_r = (
        Rational(7, 6) - Rational(2, 3) / lam - lam**4 / 256
    )
    central_factor_polynomial = (
        3 * lam**4
        + 12 * lam**3
        + 48 * lam**2
        + 192 * lam
        - 128
    )
    assert simplify(
        factor(together(central_r))
        + (lam - 4)
        * central_factor_polynomial
        / (768 * lam)
    ) == 0
    assert min(
        Poly(diff(central_factor_polynomial, lam), lam).all_coeffs()
    ) > 0
    assert central_factor_polynomial.subs(lam, Rational(64, 37)) > 0
    central_exp_upper_poly = sum(
        lam**j / math.factorial(j) for j in range(7)
    ) + lam**7 / math.factorial(7) / (1 - lam / 8)
    central_gap = together(
        poisson_three - central_exp_upper_poly * central_r
    )
    central_numerator, central_denominator = central_gap.as_numer_denom()
    assert simplify(
        central_denominator - 3870720 * lam * (8 - lam)
    ) == 0
    central_breaks = [
        Rational(64, 37),
        Rational(2),
        Rational(5, 2),
        Rational(3),
        Rational(7, 2),
        Rational(4),
    ]
    central_bernstein = []
    for lo, hi in zip(central_breaks, central_breaks[1:]):
        central_bernstein.extend(
            bernstein_coefficients(central_numerator, lam, lo, hi)
        )
    assert min(central_bernstein) == Rational(11480818817, 25344)
    assert min(central_bernstein) > 0

    # Global scalar upper bound g_4''(u)/2<=u/4.  For u<=1 it follows from
    # monotonicity of lambda*P(Pois(lambda)<=2) on lambda>=4 and exp(4)>52.
    exp_four_degree_eight = _exp_taylor_lower(Fraction(4), 8)
    assert exp_four_degree_eight > 52
    assert simplify(
        diff(lam * exp(-lam) * (1 + lam + lam**2 / 2), lam)
        + (lam**3 - lam**2 - 2 * lam - 2) * exp(-lam) / 2
    ) == 0
    scalar_lower_monotonicity = Poly(
        (lam**3 - lam**2 - 2 * lam - 2).subs(lam, x + 4),
        x,
    )
    assert min(scalar_lower_monotonicity.all_coeffs()) > 0

    # For u>1, multiply the desired inequality by exp(lambda), retain the
    # degree-eight positive Taylor lower bound, and certify the resulting
    # polynomial by exact Bernstein coefficients.
    q2_r = 1 / lam + 1 - lam**3 / 16 + 3 * lam**4 / 256
    assert simplify(
        (
            1 / lam
            - (diff(g2_upper, u, 2) / 2).subs(u, 4 / lam)
        )
        - exp(-lam)
        * (q2_r * exp(lam) - (1 + lam + lam**2 / 2))
    ) == 0
    # The multiplier of exp(lambda) is positive on (0,4]: its polynomial
    # part decreases from one to zero, and the remaining 1/lambda is
    # positive.  This makes substitution of a Taylor lower bound legitimate.
    q2_polynomial_part = 1 - lam**3 / 16 + 3 * lam**4 / 256
    assert factor(diff(q2_polynomial_part, lam)) == (
        3 * lam**2 * (lam - 4) / 64
    )
    assert q2_polynomial_part.subs(lam, 4) == 0
    exp_degree_eight = sum(
        lam**j / math.factorial(j) for j in range(9)
    )
    q2_upper_polynomial = factor(
        lam * (q2_r * exp_degree_eight - (1 + lam + lam**2 / 2))
    )
    q2_breaks = [Rational(index, 2) for index in range(9)]
    q2_bernstein = []
    for lo, hi in zip(q2_breaks, q2_breaks[1:]):
        q2_bernstein.extend(
            bernstein_coefficients(q2_upper_polynomial, lam, lo, hi)
        )
    assert min(q2_bernstein) == Rational(90499, 640640)
    assert min(q2_bernstein) > 0

    # Algebra for the auxiliary lower bound [b,c,d]g_4>=g_4(d)/d^2.
    # It is obtained by two knot insertions once [b,c,d]f_4>=0 and
    # [b,d]h_4>=0 have been established.
    b, c, d = symbols("b c d", positive=True)
    gb, gc, gd = symbols("g_b g_c g_d")

    def symbolic_dd(nodes, values):
        current = list(values)
        for order in range(1, len(nodes)):
            current = [
                (current[index + 1] - current[index])
                / (nodes[index + order] - nodes[index])
                for index in range(len(current) - 1)
            ]
        return current[0]

    q_bcd = symbolic_dd([b, c, d], [gb, gc, gd])
    q_0bd = symbolic_dd([0, b, d], [0, gb, gd])
    q_00d = gd / d**2
    insertion_four = symbolic_dd([0, b, c, d], [0, gb, gc, gd])
    assert simplify(q_bcd - q_0bd - c * insertion_four) == 0
    insertion_b = (gb / b**2 - gd / d**2) / (b - d)
    assert simplify(q_0bd - q_00d - b * insertion_b) == 0

    # On [0,2], h_4 is increasing; thereafter it is one-turn.  These exact
    # ratio derivatives encode the analytic comparison h_4(d)>=h_4(b) under
    # 2b+d<=4.  The final endpoint ordering is checked just below.
    h_ratio_log_derivative = -1 + 3 / (4 - lam)
    assert simplify(h_ratio_log_derivative - (lam - 1) / (4 - lam)) == 0
    exp_two_upper = _exp_taylor_upper(Fraction(2), 4)
    assert exp_two_upper < 8

    # The endpoint comparison needed after the single-turn argument for h_4.
    # Multiplication by exp(4) reduces h_4(4)>h_4(1) to
    # exp(3){1-81 exp(1)/256}>1.  The displayed Taylor bounds leave a
    # strictly positive rational margin.
    exp_one_upper = _exp_taylor_upper(Fraction(1), 1)
    exp_three_degree_two_lower = _exp_taylor_lower(Fraction(3), 2)
    h_endpoint_margin = (
        exp_three_degree_two_lower
        * (1 - Fraction(81, 256) * exp_one_upper)
        - 1
    )
    assert simplify(
        exp(4) * (h0_upper.subs(u, 4) - exp(-4))
        - (exp(3) * (1 - Rational(81, 256) * exp(1)) - 1)
    ) == 0
    assert exp_one_upper == Fraction(11, 4)
    assert exp_three_degree_two_lower == Fraction(17, 2)
    assert h_endpoint_margin == Fraction(213, 2048)
    assert h_endpoint_margin > 0

    # The scalar pruning test below evaluates f_4 at the lower d-endpoint.
    # Section 11 gives the one-maximum property of f_4'.  The following
    # exact endpoint estimates show that f_4 increases through 3 and that
    # f_4(4)>f_4(3).  Consequently f_4(d)>=f_4(d_0) whenever
    # 1<d_0<3 and d_0<=d<=4, which is precisely the branch-test regime.
    exp_four_thirds_upper = _exp_taylor_upper(Fraction(4, 3), 1)
    exp_four_thirds_lower = _exp_taylor_lower(Fraction(4, 3), 5)
    exp_one_degree_two_upper = _exp_taylor_upper(Fraction(1), 2)
    exp_one_degree_two_lower = _exp_taylor_lower(Fraction(1), 2)
    f_curvature_ratio = (
        4 * exp(-lam) / (3 * (1 - lam / 4) ** 2)
    )
    assert simplify(
        diff(log(f_curvature_ratio), lam) - (lam - 2) / (4 - lam)
    ) == 0
    assert simplify(
        diff(f1_upper, u).subs(u, 3)
        - (Rational(7, 3) * exp(-Rational(4, 3)) - Rational(16, 27))
    ) == 0
    assert exp_four_thirds_upper == Fraction(59, 15)
    assert exp_four_thirds_upper < Fraction(63, 16)
    assert exp_one_degree_two_lower == Fraction(5, 2)
    assert exp_one_degree_two_lower > Fraction(64, 27)
    assert exp_four_thirds_lower == Fraction(13793, 3645)
    assert exp_one_degree_two_upper == Fraction(49, 18)
    f_four_minus_f_three_margin = (
        4 / exp_one_degree_two_upper
        - Fraction(81, 64)
        - 3 / exp_four_thirds_lower
        + Fraction(16, 27)
    )
    assert f_four_minus_f_three_margin == Fraction(
        4159877,
        1167880896,
    )
    assert f_four_minus_f_three_margin > 0
    assert simplify(
        f1_upper.subs(u, 4)
        - f1_upper.subs(u, 3)
        - (
            4 * exp(-1)
            - Rational(81, 64)
            - 3 * exp(-Rational(4, 3))
            + Rational(16, 27)
        )
    ) == 0

    # The auxiliary convex core d<=4/3 follows from the exact mass ratio for
    # f_4''/2.  It is increasing in lambda on [3,4), and its lambda=3 value
    # exceeds one because exp(3)<64/3.
    pmf_two_ratio = f_curvature_ratio
    assert simplify(
        diff(log(pmf_two_ratio), lam) - (lam - 2) / (4 - lam)
    ) == 0
    exp_three_upper = _exp_taylor_upper(Fraction(3), 4)
    assert exp_three_upper == Fraction(817, 40)
    assert exp_three_upper < Fraction(64, 3)

    # Uniform derivative bound for the auxiliary Hermite--Genocchi
    # integrand r=f_4''/2.  For u<=1 its derivative is maximized at lambda=6.
    exp_six_degree_six = _exp_taylor_lower(Fraction(6), 6)
    assert exp_six_degree_six > Fraction(405, 2)
    auxiliary_lower_derivative = (
        exp(-lam) * lam**4 * (lam - 3) / 32
    )
    assert factor(diff(auxiliary_lower_derivative, lam)) == (
        -exp(-lam) * lam**3 * (lam - 6) * (lam - 2) / 32
    )
    r_prime = (
        15 * lam**6 / 2048
        - 3 * lam**5 / 64
        + 9 * lam**4 / 128
        + exp(-lam) * (lam**5 / 32 - 3 * lam**4 / 32)
    )
    r_lambda = (
        exp(-lam) * lam**3 / 8
        - Rational(3, 32) * lam**3 * (1 - lam / 4) ** 2
    )
    assert simplify(
        (diff(f1_upper, u, 2) / 2).subs(u, 4 / lam) - r_lambda
    ) == 0
    assert simplify(-lam**2 / 4 * diff(r_lambda, lam) - r_prime) == 0

    expected_h_prime = 4 / u**2 * (
        exp(-4 / u) - (1 - 1 / u) ** 3
    )
    assert simplify(diff(h0_upper, u) - expected_h_prime) == 0

    # The main Hermite--Genocchi integrand has |q'|<4 globally.  For u>1
    # this is immediate from the difference of two probabilities.  For
    # u<=1, lambda*P(Pois(lambda)=4) is maximized at lambda=5.
    exp_five_degree_four = _exp_taylor_lower(Fraction(5), 4)
    assert exp_five_degree_four > Fraction(3125, 96)
    main_lower_derivative = exp(-lam) * lam**5 / 24
    assert factor(diff(main_lower_derivative, lam)) == (
        -exp(-lam) * lam**4 * (lam - 5) / 24
    )

    previous_precision = ctx.prec
    ctx.prec = 160

    def arb_exact(value: Fraction) -> arb:
        value = Fraction(value)
        return arb(fmpq(value.numerator, value.denominator))

    def arb_interval(lo: Fraction, hi: Fraction) -> arb:
        lo_q = fmpq(lo.numerator, lo.denominator)
        hi_q = fmpq(hi.numerator, hi.denominator)
        return arb((lo_q + hi_q) / 2, (hi_q - lo_q) / 2)

    derivative_stack = [(Fraction(0), Fraction(4), 0)]
    derivative_leaves = 0
    derivative_max_depth = 0
    try:
        while derivative_stack:
            lo, hi, depth = derivative_stack.pop()
            lam_ball = arb_interval(lo, hi)
            value = (
                15 * lam_ball**6 / 2048
                - 3 * lam_ball**5 / 64
                + 9 * lam_ball**4 / 128
                + (-lam_ball).exp()
                * (lam_ball**5 / 32 - 3 * lam_ball**4 / 32)
            )
            target = arb_exact(Fraction(3, 5))
            if value < target and value > -target:
                derivative_leaves += 1
                derivative_max_depth = max(derivative_max_depth, depth)
                continue
            if depth >= 30:
                raise AssertionError(
                    "n=4 auxiliary derivative bound did not close"
                )
            midpoint = (lo + hi) / 2
            derivative_stack.append((lo, midpoint, depth + 1))
            derivative_stack.append((midpoint, hi, depth + 1))
        assert derivative_leaves == 438
        assert derivative_max_depth == 12

        def arb_f1(value: Fraction) -> arb:
            if value == 0:
                return arb(0)
            u = arb_exact(value)
            result = u * (-4 / u).exp()
            if value > 1:
                result -= (u - 1) ** 4 / u**3
            return result

        def arb_f3(value: Fraction) -> arb:
            if value == 0:
                return arb(0)
            u = arb_exact(value)
            result = u**3 * (-4 / u).exp()
            if value > 1:
                result -= (u - 1) ** 4 / u
            return result

        def arb_divided_difference(function, parameters):
            knots = []
            knot = Fraction()
            for parameter in parameters:
                knot += parameter
                knots.append(knot)
            values = [function(knot_value) for knot_value in knots]
            for order in range(1, len(knots)):
                values = [
                    (values[index + 1] - values[index])
                    / arb_exact(knots[index + order] - knots[index])
                    for index in range(len(values) - 1)
                ]
            return values[0]

        def update_transcript(digest, reason, lower, upper):
            digest.update(reason.encode("ascii") + b"|")
            for endpoint in (lower, upper):
                for value in endpoint:
                    digest.update(
                        f"{value.numerator}/{value.denominator},".encode(
                            "ascii"
                        )
                    )
            digest.update(b"\n")

        def tightened_upper(lower, upper, weights):
            result = list(upper)
            for index, weight in enumerate(weights):
                available = 4 - sum(
                    weights[other] * lower[other]
                    for other in range(len(weights))
                    if other != index
                )
                result[index] = min(result[index], available / weight)
            return tuple(result)

        # Auxiliary certificate: [b,c,d]f_4>=0 when b+c+d<=4.
        secondary_weights = (3, 2, 1)
        secondary_stack = [
            (
                (Fraction(0),) * 3,
                (Fraction(4, 3), Fraction(2), Fraction(4)),
                0,
            )
        ]
        secondary_lipschitz_leaves = 0
        secondary_core_leaves = 0
        secondary_calls = 0
        secondary_max_depth = 0
        secondary_digest = hashlib.sha256()
        while secondary_stack:
            lower, upper, depth = secondary_stack.pop()
            secondary_calls += 1
            upper = tightened_upper(lower, upper, secondary_weights)
            if any(
                upper[index] < lower[index] for index in range(3)
            ):
                update_transcript(
                    secondary_digest, "infeasible", lower, upper
                )
                continue
            if sum(upper) <= Fraction(4, 3):
                secondary_core_leaves += 1
                update_transcript(secondary_digest, "core", lower, upper)
                continue

            center = tuple(
                (lower[index] + upper[index]) / 2 for index in range(3)
            )
            weighted_center = sum(
                secondary_weights[index] * center[index]
                for index in range(3)
            )
            if weighted_center > 4:
                scale = Fraction(4, 1) / weighted_center
                center = tuple(value * scale for value in center)
            if center[1] == 0 or center[2] == 0:
                raise AssertionError("unexpected confluent secondary center")
            center_value = arb_divided_difference(arb_f1, center)
            radii = [
                max(
                    center[index] - lower[index],
                    upper[index] - center[index],
                )
                for index in range(3)
            ]
            error = Fraction(1, 5) * sum(
                secondary_weights[index] * radii[index]
                for index in range(3)
            )
            if center_value > arb_exact(error):
                secondary_lipschitz_leaves += 1
                secondary_max_depth = max(secondary_max_depth, depth)
                update_transcript(
                    secondary_digest, "L", lower, upper
                )
                continue
            if depth >= 40:
                raise AssertionError("n=4 secondary subdivision did not close")
            split_index = max(
                range(3),
                key=lambda index: secondary_weights[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            secondary_stack.append((tuple(lower_child), upper, depth + 1))
            secondary_stack.append((lower, tuple(upper_child), depth + 1))

        assert secondary_calls == 24479
        assert secondary_lipschitz_leaves == 12187
        assert secondary_core_leaves == 53
        assert secondary_calls == 2 * (
            secondary_lipschitz_leaves + secondary_core_leaves
        ) - 1
        assert secondary_max_depth == 20
        assert secondary_digest.hexdigest() == (
            "baf76e5da205718ac2f7e7037bde03a4d02e5cacfa9351c14254c24f1ca31dfe"
        )

        # Main certificate for the residual ordered four-knot face.
        main_weights = (4, 3, 2, 1)
        main_stack = [
            (
                (Fraction(0),) * 4,
                (
                    Fraction(1),
                    Fraction(4, 3),
                    Fraction(2),
                    Fraction(4),
                ),
                0,
            )
        ]
        main_counts = {
            "core": 0,
            "far": 0,
            "central": 0,
            "scalar": 0,
            "lipschitz": 0,
        }
        main_calls = 0
        main_max_depth = 0
        main_digest = hashlib.sha256()
        while main_stack:
            lower, upper, depth = main_stack.pop()
            main_calls += 1
            upper = tightened_upper(lower, upper, main_weights)
            if any(
                upper[index] < lower[index] for index in range(4)
            ):
                update_transcript(main_digest, "infeasible", lower, upper)
                continue
            d_lower = sum(lower)
            d_upper = sum(upper)
            reason = None
            if d_upper <= Fraction(19, 15):
                reason = "core"
            elif d_lower >= 3:
                reason = "far"
            elif lower[0] >= Fraction(9, 16):
                reason = "central"
            elif d_lower > 1:
                u = arb_exact(d_lower)
                f_lower = u * (-4 / u).exp() - (u - 1) ** 4 / u**3
                scalar_upper = upper[0] * (4 - d_lower) / 12
                if f_lower > arb_exact(scalar_upper):
                    reason = "scalar"
            if reason is not None:
                main_counts[reason] += 1
                update_transcript(main_digest, reason, lower, upper)
                continue

            center = tuple(
                (lower[index] + upper[index]) / 2 for index in range(4)
            )
            weighted_center = sum(
                main_weights[index] * center[index]
                for index in range(4)
            )
            if weighted_center > 4:
                scale = Fraction(4, 1) / weighted_center
                center = tuple(value * scale for value in center)
            if any(center[index] == 0 for index in range(1, 4)):
                raise AssertionError("unexpected confluent main center")
            center_value = arb_divided_difference(arb_f3, center)
            radii = [
                max(
                    center[index] - lower[index],
                    upper[index] - center[index],
                )
                for index in range(4)
            ]
            error = sum(
                main_weights[index] * radii[index]
                for index in range(4)
            )
            if center_value > arb_exact(error):
                main_counts["lipschitz"] += 1
                main_max_depth = max(main_max_depth, depth)
                update_transcript(main_digest, "L", lower, upper)
                continue
            if depth >= 40:
                raise AssertionError("n=4 main subdivision did not close")
            split_index = max(
                range(4),
                key=lambda index: main_weights[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            main_stack.append((tuple(lower_child), upper, depth + 1))
            main_stack.append((lower, tuple(upper_child), depth + 1))

        assert main_calls == 77401
        assert main_counts == {
            "core": 6190,
            "far": 0,
            "central": 23,
            "scalar": 7599,
            "lipschitz": 24889,
        }
        assert main_calls == 2 * sum(main_counts.values()) - 1
        assert main_max_depth == 28
        assert main_digest.hexdigest() == (
            "66dfba78573895b27c4ba6dd7616b68ed96a377f0c767107d30a5f08de8fbf94"
        )
    finally:
        ctx.prec = previous_precision

    return {
        "arb_precision_bits": 160,
        "extended_convex_core": {
            "largest_knot_upper": "19/15",
            "endpoint_exp_threshold": _fraction_record(core_exp_threshold),
            "endpoint_exp_upper": _fraction_record(core_exp_upper),
            "exact_margin": _fraction_record(core_endpoint_margin),
        },
        "central_affine_minorant": {
            "knot_interval": "[9/16,37/16]",
            "minorant": "q(u)>=(1-u)/6",
            "lower_endpoint_margin": _fraction_record(
                central_endpoint_margin
            ),
            "bernstein_minimum": _fraction_record(
                Fraction(11480818817, 25344)
            ),
        },
        "q2_linear_upper_bound": {
            "bound": "g_4''(u)/2<=u/4 for u>=0",
            "bernstein_minimum": _fraction_record(
                Fraction(90499, 640640)
            ),
        },
        "scalar_pruning_monotonicity": {
            "h_4_4_minus_h_4_1_scaled_margin": _fraction_record(
                h_endpoint_margin
            ),
            "f_4_4_minus_f_4_3_lower_margin": _fraction_record(
                f_four_minus_f_three_margin
            ),
        },
        "auxiliary_derivative_bound": {
            "bound": "abs((f_4''/2)')<3/5",
            "certified_leaf_intervals": derivative_leaves,
            "maximum_bisection_depth": derivative_max_depth,
        },
        "secondary_constrained_convexity": {
            "claim": "[b,c,d]f_4>=0 when 0<=b<=c<=d and b+c+d<=4",
            "parameterization": "(b,c,d)=(b,b+s,b+s+t), 3b+2s+t<=4",
            "initial_parameter_box": "[0,4/3]x[0,2]x[0,4]",
            "cumulative_coordinate_weights": [3, 2, 1],
            "box_preprocessing": (
                "tighten each upper endpoint from the other lower "
                "endpoints and the weighted budget"
            ),
            "evaluation_point_rule": (
                "coordinate midpoint, radially scaled to weighted budget "
                "4 when necessary"
            ),
            "split_rule": (
                "bisect the first coordinate maximizing weight*width; "
                "push lower child then upper child on a LIFO stack"
            ),
            "total_branch_calls": secondary_calls,
            "analytic_core_terminal_boxes": secondary_core_leaves,
            "lipschitz_terminal_boxes": secondary_lipschitz_leaves,
            "maximum_bisection_depth": secondary_max_depth,
            "terminal_transcript_sha256": secondary_digest.hexdigest(),
        },
        "main_four_positive_face": {
            "parameterization": (
                "(a,b,c,d)=(a,a+r,a+r+s,a+r+s+t), "
                "4a+3r+2s+t<=4"
            ),
            "initial_parameter_box": "[0,1]x[0,4/3]x[0,2]x[0,4]",
            "cumulative_coordinate_weights": [4, 3, 2, 1],
            "box_preprocessing": (
                "tighten each upper endpoint from the other lower "
                "endpoints and the weighted budget"
            ),
            "evaluation_point_rule": (
                "coordinate midpoint, radially scaled to weighted budget "
                "4 when necessary"
            ),
            "split_rule": (
                "bisect the first coordinate maximizing weight*width; "
                "push lower child then upper child on a LIFO stack"
            ),
            "total_branch_calls": main_calls,
            "terminal_box_counts": main_counts,
            "maximum_bisection_depth": main_max_depth,
            "terminal_transcript_sha256": main_digest.hexdigest(),
        },
    }


def verify_n5_four_positive_face() -> dict[str, object]:
    """Certify the complete four-positive coordinate face for ``n=5``.

    After two zero knots have been removed, the target is a third divided
    difference on four ordered knots.  The proof combines the existing
    all-``n`` sparse core and far cap with an affine minorant and one
    exhaustive directed real-ball subdivision of the residual compact set.
    """
    from flint import arb, ctx, fmpq
    from sympy import Rational, diff, exp, simplify, symbols

    lam = symbols("lambda", positive=True)
    u = symbols("u", positive=True)

    exponential_piece = exp(-5 / u)
    f3_lower = u**3 * exponential_piece
    f3_upper = f3_lower - (u - 1) ** 5 / u**2
    q_lower_u = exponential_piece * (
        1 + 5 / u + Rational(25, 2) / u**2 + Rational(125, 6) / u**3
    )
    q_upper_u = q_lower_u - 1 + 5 / u**4 - 4 / u**5
    assert simplify(diff(f3_lower, u, 3) / 6 - q_lower_u) == 0
    assert simplify(diff(f3_upper, u, 3) / 6 - q_upper_u) == 0

    poisson_three = exp(-lam) * (
        1 + lam + lam**2 / 2 + lam**3 / 6
    )
    q_upper_lambda = (
        poisson_three - 1 + lam**4 / 125 - 4 * lam**5 / 3125
    )
    q_upper_derivative = (
        exp(-lam) * lam**5 / 30
        - 4 * lam**5 * (5 - lam) / 3125
    )
    assert simplify(q_upper_u.subs(u, 5 / lam) - q_upper_lambda) == 0
    assert simplify(
        -lam**2 / 5 * diff(q_upper_lambda, lam) - q_upper_derivative
    ) == 0

    # On u<=1, q'(u)=exp(-lambda)lambda^5/30.  It decreases for
    # lambda>=5, and the exact Taylor lower bound makes its value at five
    # strictly smaller than one.
    q_lower_derivative = exp(-lam) * lam**5 / 30
    assert simplify(diff(q_lower_u, u) - q_lower_derivative.subs(lam, 5 / u)) == 0
    assert simplify(diff(q_lower_derivative, lam)) == (
        exp(-lam) * lam**4 * (5 - lam) / 30
    )
    exp_five_degree_six_lower = _exp_taylor_lower(Fraction(5), 6)
    assert exp_five_degree_six_lower == Fraction(16289, 144)
    assert exp_five_degree_six_lower > Fraction(625, 6)

    previous_precision = ctx.prec
    ctx.prec = 160

    def arb_exact(value: Fraction) -> arb:
        value = Fraction(value)
        return arb(fmpq(value.numerator, value.denominator))

    def arb_interval(lo: Fraction, hi: Fraction) -> arb:
        lo_q = fmpq(lo.numerator, lo.denominator)
        hi_q = fmpq(hi.numerator, hi.denominator)
        return arb((lo_q + hi_q) / 2, (hi_q - lo_q) / 2)

    def update_interval_digest(digest, label, lo, hi):
        digest.update(
            (
                f"{label}|{lo.numerator}/{lo.denominator}|"
                f"{hi.numerator}/{hi.denominator}\n"
            ).encode("ascii")
        )

    scalar_digest = hashlib.sha256()
    derivative_stack = [(Fraction(1), Fraction(5), 0)]
    derivative_leaves = 0
    derivative_max_depth = 0
    central_records = []

    try:
        while derivative_stack:
            lo, hi, depth = derivative_stack.pop()
            lam_ball = arb_interval(lo, hi)
            value = (
                (-lam_ball).exp() * lam_ball**5 / 30
                - 4 * lam_ball**5 * (5 - lam_ball) / 3125
            )
            if value < 1 and value > -1:
                derivative_leaves += 1
                derivative_max_depth = max(derivative_max_depth, depth)
                update_interval_digest(
                    scalar_digest, "derivative", lo, hi
                )
                continue
            if depth >= 30:
                raise AssertionError(
                    "n=5 four-positive derivative bound did not close"
                )
            midpoint = (lo + hi) / 2
            derivative_stack.append((lo, midpoint, depth + 1))
            derivative_stack.append((midpoint, hi, depth + 1))

        assert derivative_leaves == 17
        assert derivative_max_depth == 5

        def central_gap(value: arb, upper_piece: bool) -> arb:
            lambda_value = 5 / value
            result = (-lambda_value).exp() * (
                1
                + lambda_value
                + lambda_value**2 / 2
                + lambda_value**3 / 6
            )
            if upper_piece:
                result += -1 + 5 / value**4 - 4 / value**5
            return result - arb_exact(Fraction(2, 25)) * (
                arb_exact(Fraction(5, 4)) - value
            )

        for label, start, stop, upper_piece in (
            (
                "central-lower",
                Fraction(13, 20),
                Fraction(1),
                False,
            ),
            (
                "central-upper",
                Fraction(1),
                Fraction(61, 20),
                True,
            ),
        ):
            stack = [(start, stop, 0)]
            leaves = 0
            maximum_depth = 0
            while stack:
                lo, hi, depth = stack.pop()
                if central_gap(arb_interval(lo, hi), upper_piece) > 0:
                    leaves += 1
                    maximum_depth = max(maximum_depth, depth)
                    update_interval_digest(scalar_digest, label, lo, hi)
                    continue
                if depth >= 30:
                    raise AssertionError(
                        "n=5 four-positive affine minorant did not close"
                    )
                midpoint = (lo + hi) / 2
                stack.append((lo, midpoint, depth + 1))
                stack.append((midpoint, hi, depth + 1))
            central_records.append(
                {
                    "interval": f"[{start},{stop}]",
                    "certified_leaf_intervals": leaves,
                    "maximum_bisection_depth": maximum_depth,
                }
            )

        assert [
            record["certified_leaf_intervals"]
            for record in central_records
        ] == [7, 256]
        assert [
            record["maximum_bisection_depth"]
            for record in central_records
        ] == [5, 10]

        def arb_f3(value: Fraction) -> arb:
            if value == 0:
                return arb(0)
            point = arb_exact(value)
            result = point**3 * (-5 / point).exp()
            if value > 1:
                result -= (point - 1) ** 5 / point**2
            return result

        def arb_divided_difference(parameters):
            knots = []
            knot = Fraction()
            for parameter in parameters:
                knot += parameter
                knots.append(knot)
            values = [arb_f3(knot_value) for knot_value in knots]
            for order in range(1, len(knots)):
                values = [
                    (values[index + 1] - values[index])
                    / arb_exact(knots[index + order] - knots[index])
                    for index in range(len(values) - 1)
                ]
            return values[0]

        def tightened_upper(lower, upper, weights):
            result = list(upper)
            for index, weight in enumerate(weights):
                available = 5 - sum(
                    weights[other] * lower[other]
                    for other in range(len(weights))
                    if other != index
                )
                result[index] = min(result[index], available / weight)
            return tuple(result)

        def update_box_digest(digest, reason, lower, upper):
            digest.update(reason.encode("ascii") + b"|")
            for endpoint in (lower, upper):
                for value in endpoint:
                    digest.update(
                        f"{value.numerator}/{value.denominator},".encode(
                            "ascii"
                        )
                    )
            digest.update(b"\n")

        weights = (4, 3, 2, 1)
        stack = [
            (
                (Fraction(0),) * 4,
                (
                    Fraction(5, 4),
                    Fraction(5, 3),
                    Fraction(5, 2),
                    Fraction(5),
                ),
                0,
            )
        ]
        counts = {
            "sparse_core": 0,
            "far_cap": 0,
            "central_minorant": 0,
            "direct_lipschitz": 0,
        }
        calls = 0
        maximum_depth = 0
        box_digest = hashlib.sha256()
        while stack:
            lower, upper, depth = stack.pop()
            calls += 1
            upper = tightened_upper(lower, upper, weights)
            if any(
                upper[index] < lower[index] for index in range(4)
            ):
                update_box_digest(
                    box_digest, "infeasible", lower, upper
                )
                continue

            largest_lower = sum(lower)
            largest_upper = sum(upper)
            reason = None
            if largest_upper <= Fraction(25, 18):
                reason = "sparse_core"
            elif largest_lower >= 4:
                reason = "far_cap"
            elif lower[0] >= Fraction(13, 20):
                reason = "central_minorant"
            if reason is not None:
                counts[reason] += 1
                update_box_digest(box_digest, reason, lower, upper)
                continue

            center = tuple(
                (lower[index] + upper[index]) / 2
                for index in range(4)
            )
            weighted_center = sum(
                weights[index] * center[index] for index in range(4)
            )
            if weighted_center > 5:
                scale = Fraction(5) / weighted_center
                center = tuple(value * scale for value in center)
            if any(center[index] == 0 for index in range(1, 4)):
                raise AssertionError(
                    "unexpected confluent n=5 four-positive center"
                )

            center_value = arb_divided_difference(center)
            radii = [
                max(
                    center[index] - lower[index],
                    upper[index] - center[index],
                )
                for index in range(4)
            ]
            error = Fraction(1, 4) * sum(
                weights[index] * radii[index] for index in range(4)
            )
            if center_value > arb_exact(error):
                counts["direct_lipschitz"] += 1
                maximum_depth = max(maximum_depth, depth)
                update_box_digest(box_digest, "L", lower, upper)
                continue
            if depth >= 40:
                raise AssertionError(
                    "n=5 four-positive subdivision did not close"
                )
            split_index = max(
                range(4),
                key=lambda index: weights[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            stack.append((tuple(lower_child), upper, depth + 1))
            stack.append((lower, tuple(upper_child), depth + 1))

        assert calls == 81703
        assert counts == {
            "sparse_core": 1719,
            "far_cap": 350,
            "central_minorant": 82,
            "direct_lipschitz": 38701,
        }
        assert calls == 2 * sum(counts.values()) - 1
        assert maximum_depth == 28
        assert scalar_digest.hexdigest() == (
            "b5cc3d8cbda175722021249a795d275b0f29adbb6b76c658dada7a3903b2db2e"
        )
        assert box_digest.hexdigest() == (
            "f4bae2dcb30244286c97891596caffa5c8b1580f41aa80c600fba63c5af7d5a3"
        )
    finally:
        ctx.prec = previous_precision

    return {
        "arb_precision_bits": 160,
        "global_q_derivative_bound": {
            "claim": "abs(q'(u))<1 for 0<=u<=5",
            "analytic_lower_piece": "0<=u<=1",
            "upper_piece_certified_leaf_intervals": derivative_leaves,
            "upper_piece_maximum_bisection_depth": derivative_max_depth,
        },
        "central_affine_minorant": {
            "claim": "q(u)>=(2/25)(5/4-u) on [13/20,61/20]",
            "interval_records": central_records,
        },
        "scalar_terminal_transcript_sha256": scalar_digest.hexdigest(),
        "four_positive_face": {
            "parameterization": (
                "(a,b,c,d)=(a,a+r,a+r+s,a+r+s+t), "
                "4a+3r+2s+t<=5"
            ),
            "initial_parameter_box": "[0,5/4]x[0,5/3]x[0,5/2]x[0,5]",
            "cumulative_coordinate_weights": [4, 3, 2, 1],
            "box_preprocessing": (
                "tighten each upper endpoint from the other lower "
                "endpoints and the weighted budget"
            ),
            "evaluation_point_rule": (
                "coordinate midpoint, radially scaled to weighted budget "
                "5 when necessary"
            ),
            "split_rule": (
                "bisect the first coordinate maximizing weight*width; "
                "push lower child then upper child on a LIFO stack"
            ),
            "total_branch_calls": calls,
            "terminal_box_counts": counts,
            "maximum_bisection_depth": maximum_depth,
            "terminal_transcript_sha256": box_digest.hexdigest(),
        },
    }


def verify_n5_six_knot_theorem() -> dict[str, object]:
    """Certify the complete Dirichlet--Poissonization comparison at ``n=5``.

    The radial reduction leaves one zero knot and five ordered positive
    knots.  Hermite--Genocchi turns the remaining fourth divided difference
    into the expectation of a scalar function ``q``.  Exact scalar bounds,
    a recursive estimate from the already proved four-positive face, and a
    deterministic directed real-ball subdivision cover the constrained
    five-parameter region.
    """
    from flint import arb, ctx, fmpq
    from sympy import (
        Rational,
        diff,
        exp,
        factorial,
        integrate,
        simplify,
        symbols,
    )

    lam = symbols("lambda", positive=True)
    u = symbols("u", positive=True)

    exponential_piece = exp(-5 / u)
    f4_lower = u**4 * exponential_piece
    f4_upper = f4_lower - (u - 1) ** 5 / u
    q_lower_u = exponential_piece * sum(
        (5 / u) ** order / factorial(order) for order in range(5)
    )
    q_upper_u = q_lower_u - 1 + u**-5
    assert simplify(diff(f4_lower, u, 4) / 24 - q_lower_u) == 0
    assert simplify(diff(f4_upper, u, 4) / 24 - q_upper_u) == 0

    q_lower_derivative = exp(-lam) * lam**6 / 120
    q_upper_derivative = q_lower_derivative - lam**6 / 3125
    assert simplify(diff(q_lower_derivative, lam)) == (
        exp(-lam) * lam**5 * (6 - lam) / 120
    )
    assert simplify(
        diff(q_lower_u, u) - q_lower_derivative.subs(lam, 5 / u)
    ) == 0
    assert simplify(
        diff(q_upper_u, u) - q_upper_derivative.subs(lam, 5 / u)
    ) == 0

    # The positive core extends through u=19/16.  On [1,19/16], q is
    # decreasing because exp(lambda)>625/24; at the right endpoint its
    # positivity is the displayed exact upper bound for exp(80/19).
    core_lambda = Fraction(80, 19)
    core_poisson_polynomial = sum(
        core_lambda**order / math.factorial(order)
        for order in range(5)
    )
    core_exp_threshold = core_poisson_polynomial / (
        1 - (core_lambda / 5) ** 5
    )
    core_exp_upper = _exp_taylor_upper(core_lambda, 6)
    core_monotonicity_lower = _exp_taylor_lower(core_lambda, 4)
    core_endpoint_margin = core_exp_threshold - core_exp_upper
    assert core_poisson_polynomial == Fraction(5162241, 130321)
    assert core_exp_threshold == Fraction(32694193, 475841)
    assert core_exp_upper == Fraction(1819400743567, 26675014527)
    assert core_monotonicity_lower > Fraction(625, 24)
    assert core_endpoint_margin == Fraction(
        6372604003876864,
        12693065587542207,
    )
    assert core_endpoint_margin > 0

    # Scalar derivative bounds for the direct certificate.  On u<=1 the
    # derivative is positive, is maximized at lambda=6, and is below one.
    exp_six_lower = _exp_taylor_lower(Fraction(6), 11)
    assert exp_six_lower > Fraction(1944, 5)

    # On u>1 put r(lambda)=q'(5/lambda).  Its derivative has at most one
    # zero, from positive to negative, because
    # exp(-lambda)(6-lambda) is strictly decreasing on [1,5].  Hence the
    # minimum on any subinterval is at an endpoint.  Exact Taylor bounds at
    # lambda=5 and lambda=25/6 give the global 33/8 bound and the sharper
    # unit bound outside 1<u<6/5.
    r_derivative = lam**5 * (
        exp(-lam) * (6 - lam) / 120 - Rational(6, 3125)
    )
    assert simplify(diff(q_upper_derivative, lam) - r_derivative) == 0
    assert simplify(diff(exp(-lam) * (6 - lam), lam)) == (
        exp(-lam) * (lam - 7)
    )
    exp_one_upper = _exp_taylor_upper(Fraction(1), 1)
    exp_five_lower = _exp_taylor_lower(Fraction(5), 8)
    exp_five_upper = _exp_taylor_upper(Fraction(5), 9)
    transition_lambda = Fraction(25, 6)
    transition_exp_threshold = transition_lambda**6 / (
        120 * (transition_lambda**6 / 3125 - 1)
    )
    transition_exp_upper = _exp_taylor_upper(transition_lambda, 8)
    assert exp_one_upper < Fraction(625, 24)
    assert exp_five_lower > Fraction(3125, 24)
    assert exp_five_upper < Fraction(3125, 21)
    assert transition_exp_threshold == Fraction(48828125, 755256)
    assert transition_exp_upper < transition_exp_threshold

    # Check the four beta-tail moment formulas used to sharpen the
    # Lipschitz radii.  For a split with h high knots, their aggregate
    # Dirichlet weight S is Beta(h,5-h).  The formulas are the individual
    # low- and high-coordinate moments on {S>1-R}.
    s, radius = symbols("s radius", real=True)

    def symbolic_tail_moments(high_count):
        low_count = 5 - high_count
        density = (
            factorial(4)
            / (factorial(high_count - 1) * factorial(low_count - 1))
            * s ** (high_count - 1)
            * (1 - s) ** (low_count - 1)
        )
        low_moment = integrate(
            (1 - s) * density / low_count,
            (s, 1 - radius, 1),
        )
        high_moment = integrate(
            s * density / high_count,
            (s, 1 - radius, 1),
        )
        return simplify(low_moment), simplify(high_moment)

    symbolic_moments = {
        1: (
            radius**5 / 5,
            radius**4 * (5 - 4 * radius) / 5,
        ),
        2: (
            radius**4 - Rational(4, 5) * radius**5,
            2 * radius**3
            - 3 * radius**4
            + Rational(6, 5) * radius**5,
        ),
        3: (
            2 * radius**3
            - 3 * radius**4
            + Rational(6, 5) * radius**5,
            2 * radius**2
            - 4 * radius**3
            + 3 * radius**4
            - Rational(4, 5) * radius**5,
        ),
        4: (
            2 * radius**2
            - 4 * radius**3
            + 3 * radius**4
            - Rational(4, 5) * radius**5,
            radius
            - 2 * radius**2
            + 2 * radius**3
            - radius**4
            + radius**5 / 5,
        ),
    }
    for high_count, expected in symbolic_moments.items():
        actual = symbolic_tail_moments(high_count)
        assert all(
            simplify(actual[index] - expected[index]) == 0
            for index in range(2)
        )

    # Knot insertion for the recursive pruning rule.  If f3=F/u,
    # A=[a,b,c,d]f3, and B=[b,c,d,e]f3, then the target has the sign of
    # e*B-a*A.  Section 13 proves A,B>=0.  Also A<=1 because the
    # Hermite--Genocchi integrand for f3 is pointwise at most one.
    x0, x1, x2, x3, x4 = symbols("x0:5", positive=True)
    z0, z1, z2, z3, z4 = symbols("z0:5", real=True)

    def symbolic_dd(nodes, values):
        current = list(values)
        for order in range(1, len(nodes)):
            current = [
                (current[index + 1] - current[index])
                / (nodes[index + order] - nodes[index])
                for index in range(len(current) - 1)
            ]
        return current[0]

    nodes = [x0, x1, x2, x3, x4]
    values = [z0, z1, z2, z3, z4]
    target = symbolic_dd(
        nodes,
        [node * value for node, value in zip(nodes, values)],
    )
    left_face = symbolic_dd(nodes[:4], values[:4])
    right_face = symbolic_dd(nodes[1:], values[1:])
    assert simplify(
        (x4 - x0) * target - (x4 * right_face - x0 * left_face)
    ) == 0
    face_correction = 5 / u**4 - 4 / u**5
    assert simplify(
        diff(face_correction, u) - 20 * (1 - u) / u**6
    ) == 0
    assert face_correction.subs(u, 1) == 1

    previous_precision = ctx.prec
    ctx.prec = 160

    def arb_exact(value: Fraction) -> arb:
        value = Fraction(value)
        return arb(fmpq(value.numerator, value.denominator))

    def arb_interval(lo: Fraction, hi: Fraction) -> arb:
        lo_q = fmpq(lo.numerator, lo.denominator)
        hi_q = fmpq(hi.numerator, hi.denominator)
        return arb((lo_q + hi_q) / 2, (hi_q - lo_q) / 2)

    def update_interval_digest(digest, label, lo, hi):
        digest.update(
            (
                f"{label}|{lo.numerator}/{lo.denominator}|"
                f"{hi.numerator}/{hi.denominator}\n"
            ).encode("ascii")
        )

    scalar_digest = hashlib.sha256()
    central_records = []

    try:
        def central_gap(value: arb, upper_piece: bool) -> arb:
            lambda_value = 5 / value
            result = (-lambda_value).exp() * sum(
                lambda_value**order / math.factorial(order)
                for order in range(5)
            )
            if upper_piece:
                result += -1 + value**-5
            return result - arb_exact(Fraction(13, 50)) * (1 - value)

        for label, start, stop, upper_piece in (
            ("central-lower", Fraction(5, 8), Fraction(1), False),
            ("central-upper", Fraction(1), Fraction(5, 2), True),
        ):
            stack = [(start, stop, 0)]
            leaves = 0
            maximum_depth = 0
            while stack:
                lo, hi, depth = stack.pop()
                if central_gap(arb_interval(lo, hi), upper_piece) > 0:
                    leaves += 1
                    maximum_depth = max(maximum_depth, depth)
                    update_interval_digest(scalar_digest, label, lo, hi)
                    continue
                if depth >= 30:
                    raise AssertionError(
                        "n=5 complete affine minorant did not close"
                    )
                midpoint = (lo + hi) / 2
                stack.append((lo, midpoint, depth + 1))
                stack.append((midpoint, hi, depth + 1))
            central_records.append(
                {
                    "interval": f"[{start},{stop}]",
                    "certified_leaf_intervals": leaves,
                    "maximum_bisection_depth": maximum_depth,
                }
            )

        assert [
            record["certified_leaf_intervals"]
            for record in central_records
        ] == [9, 520]
        assert [
            record["maximum_bisection_depth"]
            for record in central_records
        ] == [7, 14]

        def arb_f3(value: Fraction) -> arb:
            if value == 0:
                return arb(0)
            point = arb_exact(value)
            result = point**3 * (-5 / point).exp()
            if value > 1:
                result -= (point - 1) ** 5 / point**2
            return result

        def arb_f4(value: Fraction) -> arb:
            if value == 0:
                return arb(0)
            point = arb_exact(value)
            result = point**4 * (-5 / point).exp()
            if value > 1:
                result -= (point - 1) ** 5 / point
            return result

        def arb_divided_difference(knots, function):
            values = [function(knot) for knot in knots]
            for order in range(1, len(knots)):
                values = [
                    (values[index + 1] - values[index])
                    / arb_exact(knots[index + order] - knots[index])
                    for index in range(len(values) - 1)
                ]
            return values[0]

        def tightened_upper(lower, upper, weights):
            result = list(upper)
            for index, weight in enumerate(weights):
                available = 5 - sum(
                    weights[other] * lower[other]
                    for other in range(len(weights))
                    if other != index
                )
                result[index] = min(result[index], available / weight)
            return tuple(result)

        def update_box_digest(digest, reason, lower, upper):
            digest.update(reason.encode("ascii") + b"|")
            for endpoint in (lower, upper):
                for value in endpoint:
                    digest.update(
                        f"{value.numerator}/{value.denominator},".encode(
                            "ascii"
                        )
                    )
            digest.update(b"\n")

        one_fifth = arb_exact(Fraction(1, 5))
        derivative_excess = arb_exact(Fraction(25, 8))
        four_fifths = arb_exact(Fraction(4, 5))
        six_fifths = arb_exact(Fraction(6, 5))

        def tail_moments(high_count: int, radius_value: Fraction):
            """Directed versions of the four checked beta-tail formulas."""
            value = arb_exact(radius_value)
            if high_count == 1:
                return (
                    value**5 / 5,
                    value**4 * (5 - 4 * value) / 5,
                )
            if high_count == 2:
                return (
                    value**4 - four_fifths * value**5,
                    2 * value**3
                    - 3 * value**4
                    + six_fifths * value**5,
                )
            if high_count == 3:
                return (
                    2 * value**3
                    - 3 * value**4
                    + six_fifths * value**5,
                    2 * value**2
                    - 4 * value**3
                    + 3 * value**4
                    - four_fifths * value**5,
                )
            if high_count == 4:
                return (
                    2 * value**2
                    - 4 * value**3
                    + 3 * value**4
                    - four_fifths * value**5,
                    value
                    - 2 * value**2
                    + 2 * value**3
                    - value**4
                    + value**5 / 5,
                )
            raise ValueError("high_count must lie in {1,2,3,4}")

        weights = (5, 4, 3, 2, 1)
        stack = [
            (
                (Fraction(0),) * 5,
                (
                    Fraction(1),
                    Fraction(5, 4),
                    Fraction(5, 3),
                    Fraction(5, 2),
                    Fraction(5),
                ),
                0,
            )
        ]
        counts = {
            "convex_core": 0,
            "central_minorant": 0,
            "recursive_global_A": 0,
            "recursive_local_A": 0,
            "direct_lipschitz": 0,
            "infeasible": 0,
        }
        calls = 0
        maximum_depth = 0
        box_digest = hashlib.sha256()
        while stack:
            lower, upper, depth = stack.pop()
            calls += 1
            upper = tightened_upper(lower, upper, weights)
            if any(upper[index] < lower[index] for index in range(5)):
                counts["infeasible"] += 1
                maximum_depth = max(maximum_depth, depth)
                update_box_digest(box_digest, "infeasible", lower, upper)
                continue

            largest_upper = sum(upper)
            if largest_upper <= Fraction(19, 16):
                counts["convex_core"] += 1
                maximum_depth = max(maximum_depth, depth)
                update_box_digest(box_digest, "core", lower, upper)
                continue
            if lower[0] >= Fraction(5, 8):
                counts["central_minorant"] += 1
                maximum_depth = max(maximum_depth, depth)
                update_box_digest(box_digest, "central", lower, upper)
                continue

            center = tuple(
                (lower[index] + upper[index]) / 2
                for index in range(5)
            )
            weighted_center = sum(
                weights[index] * center[index] for index in range(5)
            )
            if weighted_center > 5:
                weighted_lower = sum(
                    weights[index] * lower[index]
                    for index in range(5)
                )
                scale = (
                    Fraction(5) - weighted_lower
                ) / (weighted_center - weighted_lower)
                center = tuple(
                    lower[index]
                    + scale * (center[index] - lower[index])
                    for index in range(5)
                )
            if any(center[index] == 0 for index in range(1, 5)):
                raise AssertionError(
                    "unexpected confluent n=5 complete center"
                )

            radii = [
                max(
                    center[index] - lower[index],
                    upper[index] - center[index],
                )
                for index in range(5)
            ]
            knots = []
            knot = Fraction()
            for parameter in center:
                knot += parameter
                knots.append(knot)

            # First use the recursive four-positive-face estimate.  B is
            # bounded below at the center with the global |q3'|<1 radius.
            # Either A<=1 closes the box immediately, or a similarly
            # directed local upper bound for A supplies the refinement.
            right_center = arb_divided_difference(knots[1:], arb_f3)
            right_error = (
                arb_exact(radii[0])
                + arb_exact(radii[1])
                + arb_exact(Fraction(3, 4) * radii[2])
                + arb_exact(Fraction(1, 2) * radii[3])
                + arb_exact(Fraction(1, 4) * radii[4])
            )
            right_lower = right_center - right_error
            largest_lower = sum(lower)
            smallest_upper = upper[0]
            recursive_reason = None
            if right_lower > 0:
                if (
                    arb_exact(largest_lower) * right_lower
                    > arb_exact(smallest_upper)
                ):
                    recursive_reason = "recursive-global"
                else:
                    left_center = arb_divided_difference(
                        knots[:4], arb_f3
                    )
                    left_error = (
                        arb_exact(radii[0])
                        + arb_exact(Fraction(3, 4) * radii[1])
                        + arb_exact(Fraction(1, 2) * radii[2])
                        + arb_exact(Fraction(1, 4) * radii[3])
                    )
                    left_upper = left_center + left_error
                    if not left_upper < 1:
                        left_upper = arb(1)
                    if (
                        arb_exact(largest_lower) * right_lower
                        > arb_exact(smallest_upper) * left_upper
                    ):
                        recursive_reason = "recursive-local"
            if recursive_reason is not None:
                key = (
                    "recursive_global_A"
                    if recursive_reason == "recursive-global"
                    else "recursive_local_A"
                )
                counts[key] += 1
                maximum_depth = max(maximum_depth, depth)
                update_box_digest(
                    box_digest, recursive_reason, lower, upper
                )
                continue

            center_value = arb_divided_difference(knots, arb_f4)

            # For any feasible point of the box and the segment joining it
            # to the evaluation center, |q'|<1 except when 1<U<6/5 and is
            # globally below 33/8.  Bound E[W_i 1{1<U<6/5}] from both sides
            # using every admissible low/high knot split.  The excess over
            # the unit derivative bound is 33/8-1=25/8.
            knot_lower = []
            knot_upper = []
            lower_sum = Fraction()
            upper_sum = Fraction()
            for index in range(5):
                lower_sum += lower[index]
                upper_sum += upper[index]
                knot_lower.append(lower_sum)
                knot_upper.append(upper_sum)

            event_moments = [one_fifth] * 5
            largest_knot_upper = min(knot_upper[-1], Fraction(5))
            if largest_knot_upper > 1:
                for low_count in range(1, 5):
                    low_cap = knot_upper[low_count - 1]
                    if low_cap >= 1:
                        continue
                    radius_value = (
                        (largest_knot_upper - 1)
                        / (largest_knot_upper - low_cap)
                    )
                    low_moment, high_moment = tail_moments(
                        5 - low_count,
                        radius_value,
                    )
                    for index in range(low_count):
                        if low_moment < event_moments[index]:
                            event_moments[index] = low_moment
                    for index in range(low_count, 5):
                        if high_moment < event_moments[index]:
                            event_moments[index] = high_moment

            for low_count in range(1, 5):
                high_floor = knot_lower[low_count]
                if high_floor <= Fraction(6, 5):
                    continue
                radius_value = 1 - Fraction(6, 5) / high_floor
                low_tail, high_tail = tail_moments(
                    5 - low_count,
                    radius_value,
                )
                low_moment = one_fifth - low_tail
                high_moment = one_fifth - high_tail
                for index in range(low_count):
                    if low_moment < event_moments[index]:
                        event_moments[index] = low_moment
                for index in range(low_count, 5):
                    if high_moment < event_moments[index]:
                        event_moments[index] = high_moment

            knot_lipschitz = [
                one_fifth + derivative_excess * moment
                for moment in event_moments
            ]
            parameter_lipschitz = [
                sum(knot_lipschitz[index:]) for index in range(5)
            ]
            error = sum(
                parameter_lipschitz[index] * arb_exact(radii[index])
                for index in range(5)
            )
            if center_value > error:
                counts["direct_lipschitz"] += 1
                maximum_depth = max(maximum_depth, depth)
                update_box_digest(box_digest, "direct", lower, upper)
                continue

            if depth >= 50:
                raise AssertionError(
                    "n=5 complete subdivision did not close"
                )
            split_index = max(
                range(5),
                key=lambda index: weights[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            stack.append((tuple(lower_child), upper, depth + 1))
            stack.append((lower, tuple(upper_child), depth + 1))

        assert calls == 1669573
        assert counts == {
            "convex_core": 43171,
            "central_minorant": 6,
            "recursive_global_A": 105885,
            "recursive_local_A": 558019,
            "direct_lipschitz": 127706,
            "infeasible": 0,
        }
        assert calls == 2 * sum(counts.values()) - 1
        assert maximum_depth == 38
        assert scalar_digest.hexdigest() == (
            "bdf4f9223239168671e7eb3517d210b90f270b6394c668e6b21703065234fe97"
        )
        assert box_digest.hexdigest() == (
            "c07fa042526bb48d7a79cf67c80297fcbcd59481ac8b3560645bafd80e23a5f2"
        )
    finally:
        ctx.prec = previous_precision

    return {
        "arb_precision_bits": 160,
        "extended_convex_core": {
            "largest_knot_upper": "19/16",
            "endpoint_exp_threshold": _fraction_record(
                core_exp_threshold
            ),
            "endpoint_exp_upper": _fraction_record(core_exp_upper),
            "exact_margin": _fraction_record(core_endpoint_margin),
        },
        "central_affine_minorant": {
            "claim": "q(u)>=(13/50)(1-u) on [5/8,5/2]",
            "interval_records": central_records,
        },
        "derivative_bounds": {
            "global": (
                "q is globally 33/8-Lipschitz; abs(q'(u))<33/8 "
                "on each smooth piece"
            ),
            "outside_transition": (
                "abs(q'(u))<1 on the lower piece 0<u<1 and on "
                "the upper piece 6/5<=u<=5"
            ),
            "transition_interval": "1<u<6/5",
            "exp_5_upper_margin": _fraction_record(
                Fraction(3125, 21) - exp_five_upper
            ),
            "transition_endpoint_margin": _fraction_record(
                transition_exp_threshold - transition_exp_upper
            ),
        },
        "beta_tail_moment_formulas": "symbolically checked for h=1,2,3,4",
        "scalar_terminal_transcript_sha256": scalar_digest.hexdigest(),
        "five_positive_boundary": {
            "parameterization": (
                "(a,b,c,d,e)=(a,a+r,a+r+s,a+r+s+t,a+r+s+t+v), "
                "5a+4r+3s+2t+v<=5"
            ),
            "initial_parameter_box": (
                "[0,1]x[0,5/4]x[0,5/3]x[0,5/2]x[0,5]"
            ),
            "cumulative_coordinate_weights": [5, 4, 3, 2, 1],
            "evaluation_point_rule": (
                "coordinate midpoint, moved from the lower box corner "
                "along their connecting segment to weighted budget 5 "
                "when the midpoint is infeasible"
            ),
            "recursive_pruning": (
                "[a,b,c,d,e](u*f3)=(e*[b,c,d,e]f3-"
                "a*[a,b,c,d]f3)/(e-a), using the proved four-positive "
                "face, A<=1, and abs(q3')<1"
            ),
            "direct_lipschitz_pruning": (
                "unit derivative bound outside 1<U<6/5, global 33/8 "
                "bound, and directed Beta-tail coordinate moments"
            ),
            "total_branch_calls": calls,
            "terminal_box_counts": counts,
            "maximum_bisection_depth": maximum_depth,
            "terminal_transcript_sha256": box_digest.hexdigest(),
        },
    }


def _middle_knot_d_coefficient(lam: Fraction, order: int) -> Fraction:
    """Coefficient of ``z**order`` in the logarithm (15ak)."""
    lam = Fraction(lam)
    if order < 1 or lam <= 0:
        raise ValueError("require order>=1 and lambda>0")
    return (
        lam ** (order + 1) / (order + 1)
        - lam**order / order
        + (2 * lam / (lam + 2)) ** order / order
    )


def _quadratic_log_coefficient(
    linear: Fraction,
    quadratic: Fraction,
    order: int,
) -> Fraction:
    """Coefficient of ``z**order`` in ``-log(1-A*z-B*z**2)``."""
    linear = Fraction(linear)
    quadratic = Fraction(quadratic)
    if order < 1:
        raise ValueError("require order>=1")
    result = Fraction()
    # In (A*z+B*z**2)^power, choosing ``twos`` quadratic factors
    # gives total degree ``power+twos``.
    for power in range((order + 1) // 2, order + 1):
        twos = order - power
        if 0 <= twos <= power:
            result += (
                Fraction(math.comb(power, twos), power)
                * linear ** (power - twos)
                * quadratic**twos
            )
    return result


def verify_three_positive_middle_knot_region() -> None:
    """Check the all-``n`` three-positive middle-knot proof."""
    from sympy import (
        Poly,
        Rational,
        exp,
        factor,
        log,
        powsimp,
        simplify,
        symbols,
    )

    n = symbols("n", integer=True, positive=True)
    lam = symbols("lambda", positive=True)

    poisson_two = exp(-lam) * (1 + lam + lam**2 / 2)
    binomial_two = (1 - lam / n) ** (n - 2) * (
        1
        + (1 - 2 / n) * lam
        + (n - 1) * (n - 2) * lam**2 / (2 * n**2)
    )
    curvature_ratio = poisson_two / binomial_two
    p_poly = (
        (n - 1) * (n - 2) * lam**3
        + 2 * n * (n - 2) * lam**2
        + 4 * n * lam
        - 2 * n * (3 * n - 2)
    )
    q_poly = (
        (n - 1) * (n - 2) * lam**2
        + 2 * n * (n - 2) * lam
        + 2 * n**2
    )
    expected_log_derivative = (
        -lam**2 * p_poly
        / ((lam - n) * (lam**2 + 2 * lam + 2) * q_poly)
    )
    assert powsimp(
        factor(log(curvature_ratio).diff(lam) - expected_log_derivative),
        force=True,
    ) == 0
    assert simplify(
        p_poly.diff(lam)
        - (
            3 * (n - 1) * (n - 2) * lam**2
            + 4 * n * (n - 2) * lam
            + 4 * n
        )
    ) == 0

    h_formula = n / lam * (
        exp(-lam) * (lam + 2)
        - (1 - lam / n) ** (n - 1)
        * (lam + 2 - 2 * lam / n)
    )
    scaled_derivative = (
        exp(-lam) * (lam + 2)
        - (1 - lam / n) ** (n - 1)
        * (lam + 2 - 2 * lam / n)
    )
    assert simplify(h_formula - n * scaled_derivative / lam) == 0

    # lambda=1: compare d_1(z) with -log(1-z/6).
    lambda_one_differences = []
    for order in range(1, 4):
        difference = (
            _middle_knot_d_coefficient(Fraction(1), order)
            - Fraction(1, 6) ** order / order
        )
        lambda_one_differences.append(difference)
    assert lambda_one_differences == [
        Fraction(0),
        Fraction(1, 24),
        Fraction(1, 72),
    ]
    assert Fraction(1, 24) > Fraction(1, 240) + Fraction(1, 1000)

    # lambda=2: every coefficient after the first exceeds that of
    # -log(1-z) by this positive expression.
    order_symbol = symbols("r", integer=True, positive=True)
    d_two_symbolic = (
        2 ** (order_symbol + 1) / (order_symbol + 1)
        - 2**order_symbol / order_symbol
        + 1 / order_symbol
    )
    assert simplify(
        d_two_symbolic
        - 1 / order_symbol
        - 2**order_symbol
        * (order_symbol - 1)
        / (order_symbol * (order_symbol + 1))
    ) == 0

    # lambda=3: exact finite polynomial and positive tail comparison.
    lambda_three_differences = []
    for order in range(1, 9):
        difference = (
            Fraction(16, 5) ** order / order
            - _middle_knot_d_coefficient(Fraction(3), order)
        )
        lambda_three_differences.append(difference)
    assert lambda_three_differences == [
        Fraction(1, 2),
        Fraction(-1, 10),
        Fraction(-271, 300),
        Fraction(-1327, 500),
        Fraction(-7861, 1250),
        Fraction(-1636591, 131250),
        Fraction(-15185543, 875000),
        Fraction(4360661, 625000),
    ]
    finite_three_polynomial_at_quarter = sum(
        coefficient * Fraction(1, 4) ** order
        for order, coefficient in enumerate(
            lambda_three_differences[:7],
            start=1,
        )
    )
    assert finite_three_polynomial_at_quarter == Fraction(
        516239213,
        6144000000,
    )
    assert finite_three_polynomial_at_quarter > 0
    assert (
        Fraction(1, 15) * Fraction(16, 15) ** 8
        > Fraction(3, 9 * 10)
    )

    # lambda=4: coefficientwise upper bound from order three onward.
    normalized_four_difference = (
        (Rational(4, 3)) ** order_symbol
        - (3 * order_symbol - 1) / (order_symbol + 1)
        - Rational(1, 3) ** order_symbol
    )
    assert normalized_four_difference.subs(order_symbol, 3) == Rational(1, 3)
    assert Rational(4, 3) ** 4 > 3 + Rational(1, 81)

    # Rational exponential bounds used to order the four endpoint values.
    exp_three_lower = _exp_taylor_lower(Fraction(3), 4)
    e_lower = _exp_taylor_lower(Fraction(1), 4)
    assert exp_three_lower == Fraction(131, 8)
    assert exp_three_lower > 16
    assert e_lower == Fraction(65, 24)
    assert e_lower > Fraction(8, 3)

    # Algebra behind the lambda reflection r -> r/(r-2).
    r = symbols("r", positive=True)
    b = n / r
    reflected_lambda = simplify(n / (n - 2 * b))
    assert simplify(reflected_lambda - r / (r - 2)) == 0

    # Quantitative curvature inequality (15ap).  Its Poisson and binomial
    # pieces factor exactly as stated in the written proof.
    e_poisson = (lam**3 - 3 * lam - 6) * exp(-lam)
    e_r_poly = (
        (n - 1) * (n - 2) * lam**3
        + 3 * (n - 2) * lam**2
        - 3 * n * (n - 4) * lam
        - 6 * n**2
    )
    assert factor(
        2 * (lam - 2) * (1 + lam + lam**2 / 2)
        - (lam + 2)
        - (lam**3 - 3 * lam - 6)
    ) == 0
    t_symbol = symbols("t", positive=True)
    binomial_two_bracket = (
        t_symbol**2
        + lam * t_symbol
        + (n - 1) * lam**2 / (2 * n)
    )
    binomial_derivative_bracket = t_symbol * (
        lam + 2 - 2 * lam / n
    )
    assert simplify(
        (
            2 * (lam - 2) * binomial_two_bracket
            - binomial_derivative_bracket
        ).subs(t_symbol, 1 - lam / n)
        - e_r_poly / n**2
    ) == 0

    e_binomial = (1 - lam / n) ** (n - 2) * e_r_poly / n**2
    e_ratio = e_poisson / e_binomial
    e_w_poly = (
        lam**5 * n**2
        - 3 * lam**5 * n
        + 2 * lam**5
        - 2 * lam**4 * n**2
        + 9 * lam**4 * n
        - 10 * lam**4
        - 6 * lam**3 * n**2
        + 12 * lam**3 * n
        + 12 * lam**3
        + 3 * lam**2 * n**2
        - 45 * lam**2 * n
        + 6 * lam**2
        + 39 * lam * n**2
        - 15 * lam * n
        + 6 * lam
        - 27 * n**2
        + 18 * n
    )
    e_expected_log_derivative = (
        -lam**2 * e_w_poly
        / ((lam - n) * (lam**3 - 3 * lam - 6) * e_r_poly)
    )
    assert powsimp(
        factor(log(e_ratio).diff(lam) - e_expected_log_derivative),
        force=True,
    ) == 0

    x, m = symbols("x m", nonnegative=True)
    e_r_shifted = Poly(
        e_r_poly.subs({lam: x + 3, n: x + 3 + m}),
        x,
        m,
    )
    e_w_shifted = Poly(
        e_w_poly.subs({lam: x + 3, n: x + 3 + m}),
        x,
        m,
    )
    assert min(e_r_shifted.coeffs()) > 0
    assert min(e_w_shifted.coeffs()) > 0

    # The lambda=3 binomial piece is exactly the expression compared in
    # (15as), after cancelling the common factor 12.
    assert factor(e_r_poly.subs(lam, 3)) == 6 * n * (2 * n - 3)

    # At lambda=3, the logarithmic ratio is positive term by term after
    # t=3/n.  This identity is the coefficient in (15as).
    assert simplify(
        (2 - order_symbol) / (order_symbol * (order_symbol + 1))
        - 1 / (order_symbol * 2**order_symbol)
        + (
            (order_symbol - 2)
            / (order_symbol * (order_symbol + 1))
            + 1 / (order_symbol * 2**order_symbol)
        )
    ) == 0

    # Refined lambda=3/2 lower bound.  The first two logarithmic
    # coefficients agree exactly; all later differences are positive by the
    # two-root comparison in the written proof.
    lambda_three_halves_differences = [
        _middle_knot_d_coefficient(Fraction(3, 2), order)
        - _quadratic_log_coefficient(
            Fraction(27, 56),
            Fraction(225, 896),
            order,
        )
        for order in range(1, 9)
    ]
    assert lambda_three_halves_differences[:2] == [Fraction(0), Fraction(0)]
    assert all(value > 0 for value in lambda_three_halves_differences[2:])
    assert 3879 < 69**2  # both quadratic roots have modulus below 6/7
    assert Fraction(7, 4) ** 4 * Fraction(2, 2 * 5) > 1
    assert Fraction(7, 4) > 1

    # Refined lambda=3 upper bound.  The exact finite polynomial through
    # order nine is positive on z in [0,1/4] by its degree-eight Bernstein
    # coefficients; the remaining coefficients are positive analytically.
    refined_three_differences = [
        _quadratic_log_coefficient(
            Fraction(27, 10),
            Fraction(7, 4),
            order,
        )
        - _middle_knot_d_coefficient(Fraction(3), order)
        for order in range(1, 11)
    ]
    assert refined_three_differences[0] == 0
    assert refined_three_differences[1] == Fraction(7, 40)
    assert all(value < 0 for value in refined_three_differences[2:5])
    assert all(value > 0 for value in refined_three_differences[5:])

    # C(z)=sum_{r=2}^9 c_r z^(r-2), with z=x/4.  Convert its power
    # coefficients directly to degree-eight Bernstein coefficients.
    power_coefficients = [
        refined_three_differences[order - 1] / 4 ** (order - 2)
        for order in range(2, 10)
    ]
    bernstein_degree = 8
    bernstein_coefficients = []
    for index in range(bernstein_degree + 1):
        coefficient = sum(
            power_coefficients[power_order]
            * Fraction(
                math.comb(index, power_order),
                math.comb(bernstein_degree, power_order),
            )
            for power_order in range(min(index, 7) + 1)
        )
        bernstein_coefficients.append(coefficient)
    assert min(bernstein_coefficients) == Fraction(
        264755763361,
        68812800000000,
    )
    assert min(bernstein_coefficients) > 0
    assert 37**2 < 1429 < 39**2  # 16/5<p<33/10 and |q|<3/5
    tail_order = 11
    tail_margin = (
        Fraction(16, 5) ** tail_order
        - Fraction(3, 5) ** tail_order
        - 2 * 3**tail_order
        - Fraction(6, 5) ** tail_order
    )
    assert tail_margin > 0
    # If A_r is the left-minus-right tail margin, then
    # A_(r+1)-3*A_r is the displayed positive expression in the proof.
    next_tail_margin = (
        Fraction(16, 5) ** (tail_order + 1)
        - Fraction(3, 5) ** (tail_order + 1)
        - 2 * 3 ** (tail_order + 1)
        - Fraction(6, 5) ** (tail_order + 1)
    )
    assert next_tail_margin - 3 * tail_margin == (
        Fraction(1, 5) * Fraction(16, 5) ** tail_order
        + Fraction(12, 5) * Fraction(3, 5) ** tail_order
        + Fraction(9, 5) * Fraction(6, 5) ** tail_order
    )

    # Translate the refined logarithmic bounds into h_n(3/2)>h_n(3).
    endpoint_ratio = (
        (Fraction(9, 2) + Fraction(35, 12) * Fraction(1, 4))
        / (Fraction(9, 8) + Fraction(75, 128) * Fraction(1, 4))
    )
    assert endpoint_ratio == Fraction(8032, 1953)
    assert (
        Fraction(35, 12) * Fraction(9, 8)
        - Fraction(9, 2) * Fraction(75, 128)
        == Fraction(165, 256)
    )
    exp_three_degree_five_lower = _exp_taylor_lower(Fraction(3), 5)
    assert exp_three_degree_five_lower == Fraction(92, 5)
    assert exp_three_degree_five_lower > endpoint_ratio**2


def verify_three_positive_full_face() -> dict[str, object]:
    """Check the residual argument completing every three-positive face."""
    from flint import arb, ctx, fmpq
    from sympy import (
        Poly,
        Rational,
        diff,
        exp,
        expand,
        factor,
        fraction,
        log,
        simplify,
        symbols,
    )

    lam = symbols("lambda", positive=True)
    n = symbols("n", integer=True, positive=True)
    z = symbols("z", nonnegative=True)
    x = symbols("x", nonnegative=True)
    a, b, c = symbols("a b c", nonnegative=True)
    f_a, f_b, f_c = symbols("f_a f_b f_c")

    def bernstein_coefficients(poly, variable, lo, hi):
        """Exact Bernstein coefficients after mapping ``[lo,hi]`` to ``[0,1]``."""
        transformed = Poly(
            expand(poly.subs(variable, lo + (hi - lo) * x)),
            x,
        )
        degree = transformed.degree()
        power = [transformed.nth(index) for index in range(degree + 1)]
        return [
            factor(
                sum(
                    power[order]
                    * Rational(
                        math.comb(index, order),
                        math.comb(degree, order),
                    )
                    for order in range(index + 1)
                )
            )
            for index in range(degree + 1)
        ]

    # Multiplication by the coordinate function gives (15aw) exactly.
    first_ab = (b * f_b - a * f_a) / (b - a)
    first_bc = (c * f_c - b * f_b) / (c - b)
    second_g = (first_bc - first_ab) / (c - a)
    second_f_reduction = (
        c * (f_c - f_b) / (c - b)
        - a * (f_b - f_a) / (b - a)
    ) / (c - a)
    assert simplify(second_g - second_f_reduction) == 0

    # The derivative phi_n=f_n' and the exact one-maximum calculation.
    phi = exp(-lam) * (1 + lam) - (1 - lam / n) ** (n - 1) * (
        1 + lam - lam / n
    )
    t_symbol = symbols("t", positive=True)
    binomial_derivative = (
        -(n - 1)
        / n
        * t_symbol ** (n - 2)
        * (1 + (n - 1) * lam / n)
        + t_symbol ** (n - 1) * (n - 1) / n
    )
    assert factor(
        (
            binomial_derivative
            + (n - 1) * lam / n * t_symbol ** (n - 2)
        ).subs(lam, n * (1 - t_symbol))
    ) == 0
    curvature_log = (
        -lam
        - log(1 - 1 / n)
        - (n - 2) * log(1 - lam / n)
    )
    assert simplify(
        diff(curvature_log, lam) - (lam - 2) / (n - lam)
    ) == 0

    # For n>=13, this lower logarithmic bound puts the maximum of phi_n
    # before n/3.  The remaining n=5,...,12 cases use the degree-seven
    # Taylor lower bound for exp(3).
    curvature_log_lower = (
        -z
        - z**2 / (2 * (1 - z))
        + (1 / z - 2)
        * (-3 * z - 9 * z**2 / (2 * (1 - 3 * z)))
        + 3
    )
    assert simplify(
        curvature_log_lower
        - z
        * (15 * z**2 - 14 * z + 1)
        / (2 * (z - 1) * (3 * z - 1))
    ) == 0
    assert 15 * Fraction(1, 13) ** 2 - 14 * Fraction(1, 13) + 1 == Fraction(
        2,
        169,
    )
    exp_three_degree_seven = _exp_taylor_lower(Fraction(3), 7)
    assert exp_three_degree_seven == Fraction(5557, 280)
    finite_curvature_margins = [
        exp_three_degree_seven
        * Fraction(sample_n - 1, sample_n)
        * Fraction(sample_n - 3, sample_n) ** (sample_n - 2)
        - 1
        for sample_n in range(5, 13)
    ]
    assert min(finite_curvature_margins) == Fraction(353, 21875)

    # The n=4 endpoint comparison used when its maximum lies just to the
    # right of n/3.
    exp_five_halves_lower = _exp_taylor_lower(Fraction(5, 2), 4)
    assert exp_five_halves_lower > Fraction(256, 27)
    exp_half_upper = _exp_taylor_upper(Fraction(1, 2), 2)
    n4_endpoint_margin = (
        Fraction(4)
        - Fraction(7, 2) * exp_half_upper
        + Fraction(413, 4096) * _exp_taylor_lower(Fraction(3), 5)
    )
    assert exp_half_upper == Fraction(277, 168)
    assert n4_endpoint_margin == Fraction(1297, 15360)

    # Analytic upper bound for n*phi_n when n>=22.  The exact logarithmic
    # bound is encoded by the rational factor below.
    p_symbol = symbols("p", nonnegative=True)
    logarithmic_upper = (
        -p_symbol
        + (n - 1) * p_symbol**2 / (2 * (1 - p_symbol))
        + p_symbol / (1 + (n - 1) * p_symbol)
    )
    expected_logarithmic_upper = (
        (n - 1)
        * p_symbol**2
        * (n * p_symbol - 1 + p_symbol)
        / (
            2
            * (1 - p_symbol)
            * (1 + (n - 1) * p_symbol)
        )
    )
    assert simplify(logarithmic_upper - expected_logarithmic_upper) == 0
    upper_factor = (
        (1 - z)
        * (1 + lam)
        * lam**2
        * (lam - 1 + lam * z)
        / (2 * (1 - lam * z) * (1 + lam - lam * z))
    )
    assert simplify(
        upper_factor.subs(z, 1 / n)
        - n
        * (1 + lam)
        * expected_logarithmic_upper.subs(p_symbol, lam / n)
    ) == 0
    upper_factor_derivative_numerator = factor(diff(
        (1 - z)
        * (lam - 1 + lam * z)
        / ((1 - lam * z) * (1 + lam - lam * z)),
        z,
    )).as_numer_denom()[0]
    expected_upper_derivative_numerator = (
        lam**3 * z**2
        - 2 * lam**3 * z
        + lam**3
        + lam**2 * z**2
        + lam**2
        - 2 * lam * z
        - lam
        + 1
    )
    assert factor(
        upper_factor_derivative_numerator
        - expected_upper_derivative_numerator
    ) == 0
    shifted_upper_derivative = Poly(
        expand(expected_upper_derivative_numerator.subs(lam, 1 + x)),
        x,
    )
    # Every coefficient in x is a positive quadratic on 0<=z<=1/22.
    assert shifted_upper_derivative.all_coeffs() == [
        z**2 - 2 * z + 1,
        4 * z**2 - 6 * z + 4,
        5 * z**2 - 8 * z + 4,
        2 * z**2 - 4 * z + 2,
    ]
    for coefficient in shifted_upper_derivative.all_coeffs():
        assert coefficient.subs(z, Fraction(1, 22)) > 0
        assert diff(coefficient, z).subs(z, Fraction(1, 22)) < 0

    upper_factor_22 = factor(upper_factor.subs(z, Rational(1, 22)))
    exp_taylor_seven = sum(lam**order / math.factorial(order) for order in range(8))
    upper_22_gap = factor(Rational(3, 5) * exp_taylor_seven - upper_factor_22)
    upper_22_numerator, upper_22_denominator = fraction(upper_22_gap)
    assert simplify(
        upper_22_denominator
        - 8400 * (lam - 22) * (21 * lam + 22)
    ) == 0
    upper_22_positive_polynomial = -upper_22_numerator
    upper_22_breaks = [
        Rational(1),
        Rational(15, 8),
        Rational(11, 4),
        Rational(29, 8),
        Rational(9, 2),
    ]
    upper_22_bernstein = []
    for lo, hi in zip(upper_22_breaks, upper_22_breaks[1:]):
        upper_22_bernstein.extend(
            bernstein_coefficients(
                upper_22_positive_polynomial,
                lam,
                lo,
                hi,
            )
        )
    assert min(upper_22_bernstein) == Rational(
        6341019244847,
        25165824,
    )
    assert min(upper_22_bernstein) > 0

    # For lambda>=9/2, the p<=1/4 envelope is decreasing and below 3/5.
    envelope_log_derivative = factor(
        -1
        + 1 / (1 + lam)
        + 2 / lam
        + 1 / (lam - Rational(3, 4))
        - 1 / (lam + Rational(3, 4))
    )
    envelope_numerator = fraction(envelope_log_derivative)[0]
    shifted_envelope = Poly(
        expand((-envelope_numerator).subs(lam, x + Rational(9, 2))),
        x,
    )
    assert min(shifted_envelope.all_coeffs()) > 0
    exp_nine_halves_lower = _exp_taylor_lower(Fraction(9, 2), 9)
    assert exp_nine_halves_lower == Fraction(202948427, 2293760)
    assert exp_nine_halves_lower > Fraction(2475, 28)

    # In the p>=1/4 branch, n=22 is the worst analytic endpoint.
    exp_eleven_halves_lower = _exp_taylor_lower(Fraction(11, 2), 10)
    assert exp_eleven_halves_lower == Fraction(886288933661, 3715891200)
    assert exp_eleven_halves_lower > Fraction(715, 3)
    assert factor(
        5 * n * (n + 4) - 4 * (n + 1) * (n + 5)
    ) == n**2 - 4 * n - 20

    # Rigorous real-ball verification of max phi_n < 3/(5n) for the finite
    # prefix n=4,...,21.  The analytic argument above handles n>=22.
    previous_precision = ctx.prec
    ctx.prec = 128

    def arb_interval(lo: Fraction, hi: Fraction) -> arb:
        lo_q = fmpq(lo.numerator, lo.denominator)
        hi_q = fmpq(hi.numerator, hi.denominator)
        return arb((lo_q + hi_q) / 2, (hi_q - lo_q) / 2)

    finite_prefix_records = []
    try:
        for sample_n in range(4, 22):
            target = arb(fmpq(3, 5 * sample_n))
            stack = [(Fraction(1), Fraction(sample_n), 0)]
            leaves = 0
            maximum_depth = 0
            while stack:
                lo, hi, depth = stack.pop()
                lam_ball = arb_interval(lo, hi)
                value = (-lam_ball).exp() * (1 + lam_ball) - (
                    1 - lam_ball / sample_n
                ) ** (sample_n - 1) * (
                    1 + lam_ball - lam_ball / sample_n
                )
                if value < target:
                    leaves += 1
                    maximum_depth = max(maximum_depth, depth)
                    continue
                if depth >= 24:
                    raise AssertionError(
                        "finite f-prime upper-bound subdivision did not close"
                    )
                midpoint = (lo + hi) / 2
                stack.append((lo, midpoint, depth + 1))
                stack.append((midpoint, hi, depth + 1))
            finite_prefix_records.append(
                {
                    "n": sample_n,
                    "certified_leaf_intervals": leaves,
                    "maximum_bisection_depth": maximum_depth,
                }
            )
    finally:
        ctx.prec = previous_precision

    assert sum(
        record["certified_leaf_intervals"] for record in finite_prefix_records
    ) == 2006
    assert max(
        record["maximum_bisection_depth"] for record in finite_prefix_records
    ) == 10

    # Lower bound n*phi_n(n/lambda)>=(2lambda-3)/5.  For n>=9, retain only
    # the first positive logarithmic coefficient and use the n=9 endpoint.
    c_one = lam**2 * (lam - 1) / (2 * (lam + 1))
    lower_nine_rational = factor(
        exp(-lam)
        * (1 + lam)
        * c_one
        / (1 + c_one / 9)
    )
    lower_nine_without_exp = factor(lower_nine_rational / exp(-lam))
    lower_nine_derivative_gap = factor(
        Rational(5, 3)
        - (diff(lower_nine_without_exp, lam) - lower_nine_without_exp)
    )
    lower_nine_gap_numerator, lower_nine_gap_denominator = fraction(
        lower_nine_derivative_gap
    )
    assert lower_nine_gap_denominator.subs(lam, Rational(3, 2)) > 0
    lower_nine_breaks = [
        Rational(3, 2),
        Rational(7, 4),
        Rational(2),
        Rational(13, 6),
        Rational(20, 9),
    ]
    lower_nine_bernstein = []
    for lo, hi in zip(lower_nine_breaks, lower_nine_breaks[1:]):
        lower_nine_bernstein.extend(
            bernstein_coefficients(
                lower_nine_gap_numerator,
                lam,
                lo,
                hi,
            )
        )
    assert min(lower_nine_bernstein) == Rational(26185, 189)
    exp_three_degree_five = _exp_taylor_lower(Fraction(3), 5)
    assert exp_three_degree_five > Fraction(625, 36)
    exp_two_ninths_upper = _exp_taylor_upper(Fraction(2, 9), 1)
    assert exp_two_ninths_upper == Fraction(281, 225)
    assert exp_two_ninths_upper < Fraction(5, 4)
    e_upper = _exp_taylor_upper(Fraction(1), 6)
    assert e_upper < Fraction(11, 4)
    assert Fraction(11, 4) ** 2 * Fraction(5, 4) == Fraction(605, 64)
    lower_nine_endpoint = Rational(20, 9)
    lower_nine_endpoint_threshold = factor(
        lower_nine_without_exp.subs(lam, lower_nine_endpoint)
        / ((2 * lower_nine_endpoint - 3) / 5)
    )
    assert lower_nine_endpoint_threshold == Rational(2871000, 303433)
    assert Fraction(605, 64) < Fraction(2871000, 303433)

    # For n=4,...,8, three positive logarithmic coefficients suffice.
    finite_lower_bernstein_minima = []
    finite_lower_endpoint_margins = []
    for sample_n in range(4, 9):
        d_lower = sum(
            (
                lam ** (order + 1) / (order + 1)
                - lam**order / order
                + (lam / (lam + 1)) ** order / order
            )
            / sample_n**order
            for order in range(1, 4)
        )
        lower_rational = factor(
            sample_n * (1 + lam) * d_lower / (1 + d_lower)
        )
        derivative_gap = factor(
            Rational(2, 5)
            * sum(lam**order / math.factorial(order) for order in range(6))
            - (diff(lower_rational, lam) - lower_rational)
        )
        derivative_numerator, derivative_denominator = fraction(derivative_gap)
        assert derivative_denominator.subs(lam, Rational(3, 2)) > 0
        endpoint = Rational(2) + Rational(2, sample_n)
        coefficients = bernstein_coefficients(
            derivative_numerator,
            lam,
            Rational(3, 2),
            endpoint,
        )
        assert min(coefficients) > 0
        finite_lower_bernstein_minima.append(min(coefficients))

        endpoint_threshold = factor(
            5 * lower_rational.subs(lam, endpoint) / (2 * endpoint - 3)
        )
        taylor_order = 4 if sample_n == 4 else 3
        endpoint_exp_upper = _exp_taylor_upper(
            Fraction(int(endpoint.p), int(endpoint.q)),
            taylor_order,
        )
        endpoint_margin = Rational(
            endpoint_threshold.p,
            endpoint_threshold.q,
        ) - Rational(
            endpoint_exp_upper.numerator,
            endpoint_exp_upper.denominator,
        )
        assert endpoint_margin > 0
        finite_lower_endpoint_margins.append(endpoint_margin)

    # Final residual-region constant matching.
    assert simplify(
        Rational(3, 5) / n * (Rational(2, 3) * n - n / lam)
        - (2 * lam - 3) / (5 * lam)
    ) == 0
    convex_core = n**2 / (2 * (n + 1))
    assert simplify(n / convex_core - (2 + 2 / n)) == 0

    return {
        "arb_precision_bits": 128,
        "finite_fprime_upper_bound": finite_prefix_records,
        "finite_curvature_margin_minimum": _fraction_record(
            Fraction(353, 21875)
        ),
        "n4_endpoint_margin": _fraction_record(Fraction(1297, 15360)),
        "upper_n22_bernstein_minimum": _fraction_record(
            Fraction(6341019244847, 25165824)
        ),
        "lower_n9_bernstein_minimum": _fraction_record(
            Fraction(26185, 189)
        ),
        "finite_lower_bernstein_minimum": _fraction_record(
            Fraction(
                int(min(finite_lower_bernstein_minima).p),
                int(min(finite_lower_bernstein_minima).q),
            )
        ),
        "finite_lower_endpoint_margin_minimum": _fraction_record(
            Fraction(
                int(min(finite_lower_endpoint_margins).p),
                int(min(finite_lower_endpoint_margins).q),
            )
        ),
    }


def two_level_profile_regression(
    n: int,
    k: int,
    a: Fraction,
    b: Fraction,
) -> dict[str, object]:
    """Exact regression for one nontrivial two-level profile.

    The vector contains ``k`` copies of ``a`` and ``ell=n+1-k`` copies of
    ``b``.  This routine checks the exact radial-path and endpoint algebra,
    then certifies the resulting equal-block endpoint.  It is a regression
    check, not a substitute for the all-parameter argument.
    """
    if n < 1 or not 1 <= k <= n:
        raise ValueError("require n>=1 and 1<=k<=n")
    a = Fraction(a)
    b = Fraction(b)
    ell = n + 1 - k
    m = n + 1
    if not 0 < a < 1 < b:
        raise ValueError("require 0<a<1<b")
    coefficient_sum = k * a + ell * b
    if coefficient_sum > n:
        raise ValueError("the coefficient sum must be at most n")

    t_star = 1 / (1 - a)
    b_star = 1 + t_star * (b - 1)
    q = (1 - a) / (b - a)
    endpoint_sum = ell * b_star
    endpoint_sum_slack = Fraction(n) - endpoint_sum
    if endpoint_sum_slack < 0:
        raise AssertionError("the radial endpoint left the sum constraint")
    if b_star != 1 / q:
        raise AssertionError("the exact two-level endpoint identity failed")

    endpoint_lambda = Fraction(n) / b_star
    endpoint_check = equal_block_check(n, ell, endpoint_lambda)
    return {
        "n": n,
        "low_multiplicity": k,
        "high_multiplicity": ell,
        "low_coefficient": _fraction_record(a),
        "high_coefficient": _fraction_record(b),
        "coefficient_sum": _fraction_record(coefficient_sum),
        "coefficient_sum_slack": _fraction_record(Fraction(n) - coefficient_sum),
        "radial_path": {
            "t_star": _fraction_record(t_star),
            "coefficient_sum_slope": _fraction_record(coefficient_sum - m),
            "fixed_beta_threshold_q": _fraction_record(q),
            "endpoint_high_coefficient": _fraction_record(b_star),
            "endpoint_coefficient_sum": _fraction_record(endpoint_sum),
            "endpoint_sum_slack": _fraction_record(endpoint_sum_slack),
        },
        "equal_block_endpoint": {
            "positive_coordinates": ell,
            "common_positive_coordinate": _fraction_record(1 / q),
            "lambda": _fraction_record(endpoint_lambda),
            "certified_lower_margin": endpoint_check["certified_lower_margin"],
        },
        "conclusion": (
            "The radial lemma makes the poissonized tail nonincreasing "
            "toward the certified equal-block endpoint while the simplex "
            "tail P(Beta(ell,k)>q) stays fixed."
        ),
    }


def saffine_localization_obstruction() -> dict[str, object]:
    r"""Return an exact obstruction to the over-broad localization route.

    Put ``n=10`` and give ``Y`` the density proportional to

        (3/4 + 75y/388)^9,       0 <= y <= 97/75.

    This is a ``1/10``-concave probability law because its density to the
    power ``1/9`` is affine. Its mean is strictly below ``10/11``, the mean
    bound inherited from the coefficient-sum constraint. Direct integration
    gives

        P(Y>1) > 7/15,

    while, for ``S ~ Gamma(11,1)`` independent of ``Y``,

        P(SY>10) = C exp(-750/97) < 7/15,

    Thus checking every s-affine localization extremizer cannot establish
    the desired implication, even after retaining its necessary mean
    constraint. The actual Dirichlet-average class is narrower and is not
    refuted by this example.
    """
    n = 10
    support_upper = Fraction(97, 75)
    density_normalizer = Fraction(95984119, 196608000)
    mean = Fraction(729166363, 816359775)
    tail = Fraction(
        3108309643939756140704768,
        6633646218308706152889893,
    )
    comparison_level = Fraction(7, 15)
    smoothed_coefficient = Fraction(
        5573507995079350591862317513,
        5265884111440931688376513,
    )
    exp_lower, exp_upper = exact_exp_neg_bounds(Fraction(750, 97))
    smoothed_lower = smoothed_coefficient * exp_lower
    smoothed_upper = smoothed_coefficient * exp_upper
    mean_limit = Fraction(10, 11)

    if not mean < mean_limit:
        raise AssertionError("the exact s-affine mean constraint failed")
    if not tail > comparison_level:
        raise AssertionError("the exact s-affine upper-tail sign failed")
    if not smoothed_upper < comparison_level:
        raise AssertionError("the smoothed upper-tail sign failed")
    if not smoothed_upper < tail:
        raise AssertionError("the localization obstruction was not certified")

    return {
        "n": n,
        "law": {
            "support": ["0", str(support_upper)],
            "unnormalized_density": "(3/4 + 75*y/388)^9",
            "density_integral": _fraction_record(density_normalizer),
            "s_concavity_parameter": "1/10",
            "mean": _fraction_record(mean),
            "mean_constraint": "E[Y] < 10/11",
        },
        "tail_probability_P_Y_gt_1": _fraction_record(tail),
        "comparison_tail_probability": _fraction_record(comparison_level),
        "smoothed_probability": {
            "identity": "P(S*Y>10)=C*exp(-750/97)",
            "coefficient": _fraction_record(smoothed_coefficient),
            "exp_argument": "750/97",
            "exp_lower_bound": _large_fraction_record(exp_lower),
            "exp_upper_bound": _large_fraction_record(exp_upper),
            "probability_lower_bound": _large_fraction_record(smoothed_lower),
            "probability_upper_bound": _large_fraction_record(smoothed_upper),
        },
        "exact_comparisons": {
            "mean": "< 10/11",
            "P_Y_gt_1": "> 7/15",
            "P_S_Y_gt_10": "< 7/15",
            "conclusion": (
                "The generic s-affine localization class is too broad; "
                "this is not a counterexample for a Dirichlet average."
            ),
        },
    }


def build_certificate() -> dict[str, object]:
    verify_saffine_symbolic_identities()
    verify_radial_symbolic_identities()
    verify_two_level_symbolic_identities()
    verify_two_positive_knot_theorem()
    verify_n3_four_knot_theorem()
    n4_five_knot_record = verify_n4_five_knot_theorem()
    n5_four_positive_record = verify_n5_four_positive_face()
    n5_six_knot_record = verify_n5_six_knot_theorem()
    verify_three_positive_convex_core()
    sparse_convex_core_record = verify_sparse_convex_core()
    verify_three_positive_far_cap()
    four_positive_far_cap_record = verify_four_positive_far_cap()
    verify_three_positive_middle_knot_region()
    full_three_positive_record = verify_three_positive_full_face()
    checks = []
    # Include both the boundary lambda=k and strict lambda>k cases.  These
    # finite checks guard the beta/binomial and gamma/Poisson translations;
    # they are not a substitute for the cited all-parameter theorem.
    for n in (2, 3, 5, 10, 20):
        for k in range(1, n + 1):
            checks.append(equal_block_check(n, k, Fraction(k)))
            if k < n:
                checks.append(equal_block_check(n, k, Fraction(k + n, 2)))

    two_level_checks = []
    # Each pair contains an active-boundary case and a strict-sum case.
    for n, k, a, b in (
        (2, 2, Fraction(1, 4), Fraction(3, 2)),
        (2, 2, Fraction(1, 4), Fraction(5, 4)),
        (3, 2, Fraction(1, 4), Fraction(5, 4)),
        (3, 2, Fraction(1, 4), Fraction(9, 8)),
        (5, 3, Fraction(1, 3), Fraction(4, 3)),
        (5, 3, Fraction(1, 3), Fraction(5, 4)),
        (10, 5, Fraction(2, 5), Fraction(4, 3)),
        (10, 5, Fraction(2, 5), Fraction(5, 4)),
        (20, 15, Fraction(1, 2), Fraction(25, 12)),
        (20, 15, Fraction(1, 2), Fraction(2, 1)),
    ):
        two_level_checks.append(two_level_profile_regression(n, k, a, b))

    return {
        "schema_version": 16,
        "status": (
            "Research certificate: exact reductions and an obstruction, "
            "not an all-sample-size coverage certificate."
        ),
        "divided_difference_target": (
            "For y_i>=0 with sum_i y_i<=n, prove "
            "[y_0,...,y_n]H_n>=0, where "
            "H_n(y)=y^n exp(-n/y)-(y-1)_+^n."
        ),
        "equal_block_profiles": {
            "analytic_basis": (
                "Anderson--Samuels (1967), with continuity at lambda=k"
            ),
            "regression_checks": checks,
        },
        "radial_zero_knot_reduction": {
            "analytic_basis": (
                "The Laplace-transform convolution identity, the resulting "
                "mode-below-mean bound for weighted exponential sums, and "
                "the density-derivative argument in "
                "DIRICHLET-POISSONIZATION.md"
            ),
            "scope": (
                "Any strictly positive coefficient profile with sum_i "
                "y_i<=n can be moved to a zero-knot boundary without "
                "increasing the Dirichlet--Poissonization gap."
            ),
            "symbolic_identity_check": "passed",
        },
        "two_level_profiles": {
            "analytic_basis": (
                "The radial zero-knot reduction in "
                "DIRICHLET-POISSONIZATION.md and Anderson--Samuels (1967) "
                "at the resulting equal-block endpoint"
            ),
            "symbolic_identity_check": "passed",
            "exact_rational_regression_checks": two_level_checks,
        },
        "n2_three_knot_theorem": {
            "analytic_basis": (
                "The radial zero-knot reduction followed by the all-n "
                "two-positive-knot theorem"
            ),
            "scope": (
                "Every nonnegative three-coordinate profile with "
                "sum_i y_i<=2"
            ),
            "symbolic_identity_and_constant_check": "passed",
        },
        "n3_four_knot_theorem": {
            "analytic_basis": (
                "The radial zero-knot reduction followed by the constrained "
                "three-knot secant proof for H_3(u)/u"
            ),
            "scope": (
                "Every nonnegative four-coordinate profile with "
                "sum_i y_i<=3"
            ),
            "symbolic_identity_and_constant_check": "passed",
            "exact_taylor_bounds": {
                "exp_15_over_7_upper": _fraction_record(
                    Fraction(4605295343, 540244208)
                ),
                "exp_10_over_3_lower": _fraction_record(
                    Fraction(3781751, 137781)
                ),
                "exp_10_over_3_upper": _fraction_record(
                    Fraction(65696017, 2342277)
                ),
                "e_lower": _fraction_record(Fraction(685, 252)),
                "e_upper": _fraction_record(Fraction(31967, 11760)),
            },
        },
        "n4_five_knot_theorem": {
            "analytic_basis": (
                "The radial zero-knot reduction followed by exact scalar "
                "convex-core, affine-minorant, and comparison bounds, an "
                "auxiliary constrained-convexity certificate, and a "
                "deterministic 160-bit Arb branch certificate for the "
                "residual ordered four-knot face"
            ),
            "scope": (
                "Every nonnegative five-coordinate profile with "
                "sum_i y_i<=4"
            ),
            "symbolic_exact_and_interval_checks": "passed",
            "proof_record": n4_five_knot_record,
        },
        "n5_four_positive_face": {
            "analytic_basis": (
                "The all-n sparse core and four-positive far cap, a "
                "directed affine-minorant certificate, and a deterministic "
                "160-bit Arb branch certificate for the residual ordered "
                "four-knot face"
            ),
            "scope": (
                "Every n=5 profile having two zero coordinates and four "
                "ordered nonnegative coordinates with total sum at most 5"
            ),
            "symbolic_exact_and_interval_checks": "passed",
            "proof_record": n5_four_positive_record,
        },
        "n5_six_knot_theorem": {
            "analytic_basis": (
                "The radial zero-knot reduction, the complete n=5 "
                "four-positive face, exact scalar core and derivative "
                "bounds, a recursive knot-insertion estimate, and a "
                "deterministic directed 160-bit Arb branch certificate "
                "for the residual ordered five-positive boundary"
            ),
            "scope": (
                "Every nonnegative six-coordinate profile with "
                "sum_i y_i<=5"
            ),
            "symbolic_exact_and_interval_checks": "passed",
            "proof_record": n5_six_knot_record,
        },
        "three_positive_convex_core_all_n": {
            "analytic_basis": (
                "The divided-difference multiplication identity, the exact "
                "k=2 Poisson-minus-binomial CDF formula for g_n''/2, and "
                "Anderson--Samuels (1967)"
            ),
            "scope": (
                "Every n>=3 profile having n-2 zero coefficients and three "
                "ordered nonzero coefficients a<=b<=c with sum at most n "
                "and c<=n^2/(2(n+1))"
            ),
            "symbolic_identity_check": "passed",
        },
        "sparse_convex_core_all_n": {
            "analytic_basis": (
                "The repeated divided-difference multiplication identity, "
                "the arbitrary-order derivative-CDF identities proved by "
                "induction, and Anderson--Samuels (1967)"
            ),
            "scope": (
                "For every n>=2 and 1<=k<=n-1, every profile having n-k "
                "zero coefficients and k+1 ordered nonnegative coefficients "
                "x_0<=...<=x_k with sum at most n and "
                "x_k<=n^2/(k(n+1)); in particular, every n>=4 "
                "four-positive profile in the core d<=n^2/(3(n+1))"
            ),
            "symbolic_identity_and_threshold_check": "passed",
            "corroboration_record": sparse_convex_core_record,
        },
        "three_positive_far_cap_all_n": {
            "analytic_basis": (
                "The divided-difference multiplication identity, an exact "
                "monotonicity proof for g_n on [0,n], and exact endpoint "
                "bounds comparing the two secant slopes"
            ),
            "scope": (
                "Every n>=4 profile having n-2 zero coefficients and three "
                "ordered nonnegative coefficients a<=b<=c with sum at most "
                "n and c>=n-1"
            ),
            "symbolic_identity_and_constant_check": "passed",
            "exact_rational_bounds": {
                "exp_4_degree_9_lower": _fraction_record(
                    Fraction(153527, 2835)
                ),
                "exp_3_degree_4_lower": _fraction_record(
                    Fraction(131, 8)
                ),
                "n4_g4_at_3_lower": _fraction_record(
                    Fraction(17, 36)
                ),
                "n4_25_exp_neg_4_upper": _fraction_record(
                    Fraction(25, 54)
                ),
            },
        },
        "four_positive_far_cap_all_n": {
            "analytic_basis": (
                "The multiplication identity expressing the third divided "
                "difference through two weighted second divided "
                "differences, the complete three-positive machinery for "
                "g_n, and exact endpoint estimates"
            ),
            "scope": (
                "Every n>=4 profile having n-3 zero coefficients and four "
                "ordered nonnegative coefficients a<=b<=c<=d with sum at "
                "most n and d>=n-1"
            ),
            "symbolic_identity_and_constant_check": "passed",
            "exact_rational_bounds": four_positive_far_cap_record,
        },
        "three_positive_middle_knot_all_n": {
            "analytic_basis": (
                "The divided-difference multiplication identity, an exact "
                "one-crossing curvature proof, exact logarithmic-series "
                "endpoint comparisons, a quantitative curvature inequality, "
                "and the endpoint tangent-remainder argument"
            ),
            "scope": (
                "Every n>=4 profile having n-2 zero coefficients and three "
                "ordered nonnegative coefficients a<=b<=c with sum at most "
                "n and b<=n/3"
            ),
            "symbolic_identity_and_constant_check": "passed",
            "exact_rational_bounds": {
                "lambda_3_finite_polynomial_at_z_1_over_4": (
                    _fraction_record(Fraction(516239213, 6144000000))
                ),
                "exp_3_degree_4_lower": _fraction_record(
                    Fraction(131, 8)
                ),
                "e_degree_4_lower": _fraction_record(
                    Fraction(65, 24)
                ),
                "refined_lambda_3_bernstein_minimum": _fraction_record(
                    Fraction(264755763361, 68812800000000)
                ),
                "exp_3_degree_5_lower": _fraction_record(
                    Fraction(92, 5)
                ),
            },
        },
        "three_positive_full_face_all_n": {
            "analytic_basis": (
                "The divided-difference multiplication identity, the "
                "one-maximum derivative lemma, exact global upper and "
                "residual-region lower bounds for f_n', and a rigorous "
                "real-ball check of the finite n=4,...,21 prefix"
            ),
            "scope": (
                "Every n>=4 profile having n-2 zero coefficients and three "
                "ordered nonnegative coefficients a<=b<=c with sum at most n"
            ),
            "symbolic_identity_and_constant_check": "passed",
            "proof_record": full_three_positive_record,
        },
        "two_positive_knots_all_n": {
            "analytic_basis": (
                "A divided-difference multiplication identity and the "
                "one-variable monotonicity proof for H_n(u)/u^(n-1)"
            ),
            "scope": (
                "Every n>=2 profile having at most two nonzero "
                "coefficients and sum_i y_i<=n"
            ),
            "symbolic_identity_and_constant_check": "passed",
        },
        "localization_obstruction": saffine_localization_obstruction(),
        "exponential_series_pairs": EXACT_EXP_PAIRS,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args(argv)
    certificate = build_certificate()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    obstruction = certificate["localization_obstruction"]
    print("Dirichlet--Poissonization reduction checks: PASS")
    print(
        "equal-block rational regression cases:",
        len(certificate["equal_block_profiles"]["regression_checks"]),
    )
    print(
        "radial zero-knot reduction symbolic identities:",
        certificate["radial_zero_knot_reduction"]["symbolic_identity_check"],
    )
    print(
        "all-two-level exact rational regressions:",
        len(
            certificate["two_level_profiles"][
                "exact_rational_regression_checks"
            ]
        ),
    )
    print(
        "all-n two-positive-knot comparison:",
        certificate["two_positive_knots_all_n"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "complete n=2 three-knot comparison:",
        certificate["n2_three_knot_theorem"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "complete n=3 four-knot comparison:",
        certificate["n3_four_knot_theorem"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "complete n=4 five-knot comparison:",
        certificate["n4_five_knot_theorem"][
            "symbolic_exact_and_interval_checks"
        ],
    )
    print(
        "complete n=5 four-positive coordinate face:",
        certificate["n5_four_positive_face"][
            "symbolic_exact_and_interval_checks"
        ],
    )
    print(
        "complete n=5 six-knot comparison:",
        certificate["n5_six_knot_theorem"][
            "symbolic_exact_and_interval_checks"
        ],
    )
    print(
        "all-n three-positive convex core:",
        certificate["three_positive_convex_core_all_n"][
            "symbolic_identity_check"
        ],
    )
    print(
        "all-n sparse convex core:",
        certificate["sparse_convex_core_all_n"][
            "symbolic_identity_and_threshold_check"
        ],
    )
    print(
        "all-n three-positive far cap:",
        certificate["three_positive_far_cap_all_n"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "all-n four-positive far cap:",
        certificate["four_positive_far_cap_all_n"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "all-n three-positive middle-knot region:",
        certificate["three_positive_middle_knot_all_n"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "all-n complete three-positive face:",
        certificate["three_positive_full_face_all_n"][
            "symbolic_identity_and_constant_check"
        ],
    )
    print(
        "mean-constrained s-affine obstruction: P(Y>1)=%s; P(SY>10)<%s"
        % (
            obstruction["tail_probability_P_Y_gt_1"]["decimal"],
            obstruction["smoothed_probability"]["probability_upper_bound"][
                "decimal"
            ],
        )
    )
    print("first open complete simplex dimension: n=6")
    print("general Dirichlet-average inequality: OPEN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
