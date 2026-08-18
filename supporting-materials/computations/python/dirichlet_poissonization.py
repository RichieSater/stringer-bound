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
        "schema_version": 7,
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
