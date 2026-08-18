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
        "schema_version": 13,
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
    print("general Dirichlet-average inequality: OPEN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
