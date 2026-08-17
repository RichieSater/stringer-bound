"""Exact scalar calibrations for all-n Poisson Stringer bounds.

At the conventional levels ``alpha < exp(-1)``, let
``lambda_j(alpha)`` be the usual Poisson confidence limits and let

    Q_n(kappa) = P(V_(i:n) <= min(1, kappa*lambda_(i-1)/n), all i),

for uniform order statistics.  The corrected simultaneous-band lemma shows
that ``kappa * SB_P`` has coverage at least ``Q_n(kappa)``. It also applies
to an anchored path that fixes the ordinary zero-taint factor and multiplies
only the error increments by ``eta``. This module certifies adjacent dyadic
brackets for the smallest multiplier on each path whose corrected band event
has probability at least ``1-alpha``. On either path, capping the effective
calibrated factors at one preserves the same event and gives a pointwise
no-larger report; sample-level output encloses both cap conventions.

Floating-point arithmetic only locates a candidate dyadic cell.  Every
reported conclusion uses:

* exact rational brackets for each ``lambda_j``;
* exact rational evaluation of Bolshev's recursion; and
* opposite endpoint choices to enclose ``Q_n(kappa)`` monotonically.

The fixed certificate covers representative audit sample sizes. The
analytic theorem also handles other confidence levels by imposing the
additional terminal-factor lower bound recorded in the theory note; this
script is deliberately limited to the conventional regime where that lower
bound is one. A separate analytic corollary shows that multiplier two is
valid uniformly in the sample size when ``alpha <= 1-2/e``. A refinement
retains a finite prefix of the exact Poisson crossing probabilities and
bounds the infinite remainder geometrically, producing smaller uniform
choices at conventional levels. The all-sample-size theorem is recorded in
``theory/POISSON-BAND-CALIBRATION.md``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from mpmath import mp

from poisson_band_certificate import bolshev_probability
from stringer import (
    EXACT_EXP_PAIRS,
    exact_exp_neg_bounds,
    exact_poisson_cdf_bounds,
    exact_poisson_lambda_brackets,
)


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


CASES = {
    "0.10": (25, 50, 100, 200),
    "0.05": (25, 50, 100, 200),
    "0.01": (25, 50, 100, 200),
}
ELEMENTARY_UNIFORM_MULTIPLIERS = {
    "0.10": Fraction(44, 25),
    "0.05": Fraction(83, 50),
    "0.01": Fraction(38, 25),
}
REFINED_UNIFORM_MULTIPLIERS = {
    "0.10": (Fraction(762953, 500000), 200),
    "0.05": (Fraction(287227, 200000), 200),
    "0.01": (Fraction(1320081, 1000000), 200),
}
FACTOR_BITS = 64
KAPPA_BITS = 28
UNIFORM_PROBABILITY_BITS = 256
NUMERICAL_DPS = 80
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "poisson-band-calibration-certificate.json"
)
DECIMAL_POLICY = (
    "Fields named decimal are display approximations to the exact numerator "
    "and denominator. The fields kappa_upper.valid_decimal_ceiling_12 and "
    "eta_upper.valid_decimal_ceiling_12 are rounded upward and are "
    "themselves valid decimal multipliers."
)


def _mp_fraction(value: Fraction):
    value = Fraction(value)
    return mp.mpf(value.numerator) / value.denominator


def bolshev_probability_numeric(boundaries):
    """High-precision locator version of Bolshev's recursion.

    This function is never used for a certificate sign.  Its output only
    identifies the dyadic cell that is subsequently checked exactly.
    """
    boundaries = tuple(mp.mpf(value) for value in boundaries)
    q = [mp.mpf(1)]
    for m in range(1, len(boundaries) + 1):
        excluded = sum(
            math.comb(m, i)
            * (1 - boundaries[m - i]) ** i
            * q[m - i]
            for i in range(1, m + 1)
        )
        q.append(1 - excluded)
    return q[-1]


def _numeric_event_probability(n, lambdas, kappa):
    return bolshev_probability_numeric(
        min(mp.mpf(1), kappa * lambdas[j] / n)
        for j in range(n)
    )


def _numeric_zero_anchor_event_probability(n, lambdas, eta):
    p0 = min(mp.mpf(1), lambdas[0] / n)
    boundaries = [p0]
    boundaries.extend(
        min(
            mp.mpf(1),
            p0 + eta * (min(mp.mpf(1), lambdas[j] / n) - p0),
        )
        for j in range(1, n)
    )
    return bolshev_probability_numeric(boundaries)


def locate_multiplier(n: int, alpha: Fraction, brackets):
    """Numerically locate the minimal scalar calibration.

    ``brackets`` must contain rigorous lambda brackets through ``n-1``.
    The result has no proof status until :func:`certify_case` performs the
    exact endpoint checks.
    """
    if n < 1:
        raise ValueError("n must be positive")
    alpha = Fraction(alpha)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    if len(brackets) < n:
        raise ValueError("insufficient Poisson-limit brackets")

    with mp.workdps(NUMERICAL_DPS):
        lambdas = [
            (_mp_fraction(lower) + _mp_fraction(upper)) / 2
            for lower, upper in brackets[:n]
        ]
        target = 1 - _mp_fraction(alpha)
        if _numeric_event_probability(n, lambdas, 1) >= target:
            return mp.mpf(1)

        lower, upper = mp.mpf(1), mp.mpf(2)
        while _numeric_event_probability(n, lambdas, upper) < target:
            upper *= 2
        for _ in range(80):
            midpoint = (lower + upper) / 2
            if _numeric_event_probability(n, lambdas, midpoint) >= target:
                upper = midpoint
            else:
                lower = midpoint
        return +upper


def locate_zero_anchor_multiplier(n: int, alpha: Fraction, brackets):
    """Numerically locate the minimal zero-taint-preserving multiplier."""
    if n < 1:
        raise ValueError("n must be positive")
    alpha = Fraction(alpha)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    if len(brackets) < n:
        raise ValueError("insufficient Poisson-limit brackets")

    with mp.workdps(NUMERICAL_DPS):
        lambdas = [
            (_mp_fraction(lower) + _mp_fraction(upper)) / 2
            for lower, upper in brackets[:n]
        ]
        target = 1 - _mp_fraction(alpha)
        if _numeric_zero_anchor_event_probability(
                n, lambdas, mp.mpf(1)) >= target:
            return mp.mpf(1)

        lower, upper = mp.mpf(1), mp.mpf(2)
        while _numeric_zero_anchor_event_probability(
                n, lambdas, upper) < target:
            upper *= 2
        for _ in range(80):
            midpoint = (lower + upper) / 2
            if _numeric_zero_anchor_event_probability(
                    n, lambdas, midpoint) >= target:
                upper = midpoint
            else:
                lower = midpoint
        return +upper


def _exact_boundaries(n: int, brackets, kappa: Fraction,
                      lambda_endpoint: str):
    position = 0 if lambda_endpoint == "lower" else 1
    return tuple(
        min(Fraction(1), kappa * brackets[j][position] / n)
        for j in range(n)
    )


def _exact_zero_anchor_boundaries(
        n: int, brackets, eta: Fraction, event_endpoint: str):
    """Enclose the anchored event boundaries monotonically.

    For boundaries after the first, the anchored factor is increasing in
    ``lambda_j`` and decreasing in ``lambda_0`` because ``eta>=1``. Thus a
    lower event enclosure pairs lower ``lambda_j`` endpoints with the upper
    ``lambda_0`` endpoint, and an upper enclosure reverses those choices.
    """
    eta = Fraction(eta)
    if eta < 1:
        raise ValueError("eta must be at least one")
    if event_endpoint not in {"lower", "upper"}:
        raise ValueError("event_endpoint must be lower or upper")

    p0_lower = min(Fraction(1), brackets[0][0] / n)
    p0_upper = min(Fraction(1), brackets[0][1] / n)
    if event_endpoint == "lower":
        boundaries = [p0_lower]
        p0_for_later = p0_upper
        position = 0
    else:
        boundaries = [p0_upper]
        p0_for_later = p0_lower
        position = 1

    for j in range(1, n):
        pj = min(Fraction(1), brackets[j][position] / n)
        value = eta * pj - (eta - 1) * p0_for_later
        boundaries.append(max(Fraction(0), min(Fraction(1), value)))
    return tuple(boundaries)


def _decimal(value, digits=16):
    value = Fraction(value)
    with localcontext() as context:
        context.prec = digits + 30
        decimal = Decimal(value.numerator) / Decimal(value.denominator)
        return format(decimal, f".{digits}f")


def _decimal_ceiling(value, digits=12):
    """Round a nonnegative rational upward to a fixed decimal precision."""
    value = Fraction(value)
    if value < 0:
        raise ValueError("upward decimal helper requires a nonnegative value")
    scale = 10 ** digits
    numerator = (
        value.numerator * scale + value.denominator - 1
    ) // value.denominator
    return f"{numerator // scale}.{numerator % scale:0{digits}d}"


def _fraction_record(value, digits=16):
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": _decimal(value, digits),
    }


def _factor_bracket_records(brackets):
    """Serialize every exact Poisson-limit bracket used by a certificate."""
    denominator = 1 << FACTOR_BITS
    return [
        {
            "j": j,
            "dyadic_denominator": str(denominator),
            "lower_numerator": str(lower * denominator),
            "upper_numerator": str(upper * denominator),
            "lower_decimal": _decimal(lower),
            "upper_decimal": _decimal(upper),
            "width": str(upper - lower),
            "endpoint_signs": {
                "cdf_at_lower": "> alpha",
                "cdf_at_upper": "< alpha",
            },
        }
        for j, (lower, upper) in enumerate(brackets)
    ]


def _conventional_tail_level(alpha: Fraction):
    """Return whether ``alpha < exp(-1)`` without a float comparison."""
    alpha = Fraction(alpha)
    exp_lower, _ = exact_exp_neg_bounds(Fraction(1))
    return alpha < exp_lower


def certify_uniform_multiplier(alpha: Fraction, kappa: Fraction):
    """Certify the sample-size-uniform analytic multiplier condition.

    For ``kappa=p/q``, the condition

        alpha**(kappa-1) <= 1-kappa*exp(1-kappa)

    is proved without floating point by enclosing the exponential above and
    comparing the positive quantities after raising both sides to ``q``.
    """
    alpha = Fraction(alpha)
    kappa = Fraction(kappa)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    if kappa <= 1:
        raise ValueError("uniform kappa must exceed one")
    exp_one_lower, _ = exact_exp_neg_bounds(Fraction(1))
    if alpha > exp_one_lower:
        raise ValueError("uniform certificate requires alpha <= exp(-1)")

    _, exp_upper = exact_exp_neg_bounds(kappa - 1)
    rho_upper = kappa * exp_upper
    rhs_lower = 1 - rho_upper
    if rhs_lower <= 0:
        raise ArithmeticError("uniform multiplier right side is not positive")

    denominator = kappa.denominator
    exponent_numerator = kappa.numerator - denominator
    powered_lhs = alpha ** exponent_numerator
    powered_rhs_lower = rhs_lower ** denominator
    margin = powered_rhs_lower - powered_lhs
    if margin < 0:
        raise ArithmeticError("uniform multiplier condition is not certified")
    return {
        "alpha": str(alpha),
        "nominal_confidence": _decimal(1 - alpha),
        "kappa": _fraction_record(kappa, 12),
        "exp_upper_at_kappa_minus_one": _fraction_record(exp_upper),
        "rho_upper": _fraction_record(rho_upper),
        "one_minus_rho_lower": _fraction_record(rhs_lower),
        "powered_condition": (
            f"alpha^{exponent_numerator} <= "
            f"(one_minus_rho_lower)^{denominator}"
        ),
        "powered_lhs": _fraction_record(powered_lhs),
        "powered_rhs_lower": _fraction_record(powered_rhs_lower),
        "powered_margin_lower": _fraction_record(margin, 60),
        "conclusion": (
            "The displayed kappa is a rigorous valid full-scale multiplier "
            "for every positive sample size under the analytic geometric-"
            "union-bound criterion."
        ),
    }


def _dyadic_ceiling(value: Fraction, bits: int):
    """Round a nonnegative rational upward to the ``2**-bits`` grid."""
    value = Fraction(value)
    if value < 0:
        raise ValueError("dyadic ceiling requires a nonnegative value")
    denominator = 1 << bits
    numerator = (
        value.numerator * denominator + value.denominator - 1
    ) // value.denominator
    return Fraction(numerator, denominator)


def certify_refined_uniform_multiplier(
        alpha: Fraction, kappa: Fraction, prefix_terms: int, brackets):
    """Certify a finite-prefix plus geometric-tail uniform multiplier.

    The analytic corollary bounds the all-``n`` crossing probability by

        sum(P(Pois(kappa*lambda_j) <= j), j < J)
        + alpha**kappa * rho**J / (1-rho),

    where ``rho=kappa*exp(1-kappa)``.  Each finite-prefix term is bounded
    above by substituting the rigorous lower endpoint for ``lambda_j`` and
    then rounding the rational exponential enclosure upward to a common
    dyadic grid.  The infinite tail is bounded on the same grid.
    """
    alpha = Fraction(alpha)
    kappa = Fraction(kappa)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    if kappa <= 1:
        raise ValueError("refined uniform kappa must exceed one")
    if prefix_terms < 1:
        raise ValueError("prefix_terms must be positive")
    if len(brackets) < prefix_terms:
        raise ValueError("insufficient Poisson-limit brackets")
    exp_one_lower, _ = exact_exp_neg_bounds(Fraction(1))
    if alpha > exp_one_lower:
        raise ValueError("uniform certificate requires alpha <= exp(-1)")

    term_upper_bounds = []
    lambda_lower_numerators = []
    lambda_denominator = 1 << FACTOR_BITS
    for j in range(prefix_terms):
        lambda_lower = Fraction(brackets[j][0])
        scaled_lambda_lower = lambda_lower * lambda_denominator
        if scaled_lambda_lower.denominator != 1:
            raise ValueError(
                "Poisson-limit bracket is not on the configured dyadic grid")
        lambda_lower_numerators.append(str(scaled_lambda_lower.numerator))
        _, probability_upper = exact_poisson_cdf_bounds(
            kappa * lambda_lower, j, EXACT_EXP_PAIRS)
        term_upper_bounds.append(_dyadic_ceiling(
            probability_upper, UNIFORM_PROBABILITY_BITS))
    finite_prefix_upper = sum(term_upper_bounds, Fraction(0))

    _, exp_rho_upper = exact_exp_neg_bounds(
        kappa - 1, EXACT_EXP_PAIRS)
    rho_upper = _dyadic_ceiling(
        kappa * exp_rho_upper, UNIFORM_PROBABILITY_BITS)
    if rho_upper >= 1:
        raise ArithmeticError("geometric-tail ratio is not below one")

    # lambda_0=-log(alpha), and the stored lower endpoint is below lambda_0.
    # Hence alpha*exp(-(kappa-1)*lambda_0) is bounded above by the following
    # exact rational enclosure of alpha**kappa.
    _, exp_alpha_upper = exact_exp_neg_bounds(
        (kappa - 1) * Fraction(brackets[0][0]), EXACT_EXP_PAIRS)
    alpha_kappa_upper = _dyadic_ceiling(
        alpha * exp_alpha_upper, UNIFORM_PROBABILITY_BITS)
    tail_upper = _dyadic_ceiling(
        alpha_kappa_upper * rho_upper ** prefix_terms / (1 - rho_upper),
        UNIFORM_PROBABILITY_BITS,
    )
    total_upper = finite_prefix_upper + tail_upper
    margin_lower = alpha - total_upper
    if margin_lower <= 0:
        raise ArithmeticError(
            "refined uniform multiplier condition is not certified")

    probability_denominator = 1 << UNIFORM_PROBABILITY_BITS
    return {
        "alpha": str(alpha),
        "nominal_confidence": _decimal(1 - alpha),
        "kappa": _fraction_record(kappa, 12),
        "prefix_terms": prefix_terms,
        "lambda_lower_dyadic_bits": FACTOR_BITS,
        "lambda_lower_numerators": lambda_lower_numerators,
        "probability_upper_dyadic_bits": UNIFORM_PROBABILITY_BITS,
        "finite_prefix_term_upper_numerators": [
            str(value * probability_denominator)
            for value in term_upper_bounds
        ],
        "finite_prefix_upper": _fraction_record(
            finite_prefix_upper, 60),
        "rho_upper": _fraction_record(rho_upper, 60),
        "alpha_to_kappa_upper": _fraction_record(
            alpha_kappa_upper, 60),
        "geometric_tail_upper": _fraction_record(tail_upper, 60),
        "total_crossing_probability_upper": _fraction_record(
            total_upper, 60),
        "margin_below_alpha_lower": _fraction_record(
            margin_lower, 60),
        "conclusion": (
            "The displayed kappa is a rigorous valid full-scale multiplier "
            "for every positive sample size under the finite-prefix plus "
            "geometric-tail crossing bound."
        ),
    }


def _parse_taints(text: str | None, n: int):
    """Parse comma-separated decimal taints as exact rationals.

    Zero taints may be omitted, but no more than ``n`` entries may be
    supplied. This mirrors the zero-heavy command convention used by the
    Gaffke safeguard.
    """
    if text is None or not text.strip():
        return ()
    pieces = [piece.strip() for piece in text.split(",")]
    if any(not piece for piece in pieces):
        raise ValueError("taints must be nonempty comma-separated decimals")
    try:
        values = tuple(Fraction(piece) for piece in pieces)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"invalid taint: {error}") from error
    if len(values) > n:
        raise ValueError("more taints than sample units")
    if any(not 0 <= value <= 1 for value in values):
        raise ValueError("every taint must lie in [0,1]")
    return values


def exact_calibrated_report(
        n: int, brackets, kappa: Fraction, taints,
        zero_anchor_eta: Fraction | None = None):
    """Return rational upper enclosures for calibrated audit results.

    The Stringer expression is written by summation by parts, so every
    Poisson factor has a nonnegative coefficient. Substituting the upper
    endpoint of each exact factor bracket therefore encloses ordinary
    Poisson Stringer from above. Multiplication by the certified ``kappa``
    and capping at one gives a directly reportable rational upper bound.
    The pointwise no-larger variant that caps each calibrated factor at one
    is enclosed separately.
    If ``zero_anchor_eta`` is supplied, the function also encloses the
    affine calibration that fixes the ordinary zero-taint factor and its
    calibrated-factor-capped variant.
    """
    raw = tuple(Fraction(value) for value in taints)
    if len(raw) > n:
        raise ValueError("more taints than sample units")
    if any(not 0 <= value <= 1 for value in raw):
        raise ValueError("every taint must lie in [0,1]")
    ordered = sorted((value for value in raw if value), reverse=True)
    maximum_j = len(ordered)
    if len(brackets) <= maximum_j:
        raise ValueError("insufficient Poisson-limit brackets for sample")

    factors_upper = [
        brackets[j][1] / n
        for j in range(maximum_j + 1)
    ]

    def from_factor_upper_endpoints(factors):
        if not ordered:
            return factors[0]
        value = factors[0] * (1 - ordered[0])
        value += sum(
            factors[j] * (ordered[j - 1] - ordered[j])
            for j in range(1, maximum_j)
        )
        value += factors[maximum_j] * ordered[-1]
        return value

    ordinary_upper = from_factor_upper_endpoints(factors_upper)

    calibrated_upper = min(Fraction(1), Fraction(kappa) * ordinary_upper)
    factorwise_calibrated_factors_upper = [
        min(Fraction(1), Fraction(kappa) * factor)
        for factor in factors_upper
    ]
    factorwise_calibrated_upper = min(
        calibrated_upper,
        from_factor_upper_endpoints(factorwise_calibrated_factors_upper),
    )
    report = {
        "factor_convention": "untruncated_poisson_factors_then_final_cap",
        "reported_calibration_variants": [
            "final_cap",
            "calibrated_factorwise_cap",
        ],
        "sample_size": n,
        "nonzero_taint_count": maximum_j,
        "zero_taint_count": n - maximum_j,
        "sorted_nonzero_taints": [
            _fraction_record(value, 16) for value in ordered
        ],
        "ordinary_poisson_stringer_upper": _fraction_record(
            ordinary_upper, 16),
        "calibrated_poisson_upper": _fraction_record(
            calibrated_upper, 16),
        "factorwise_capped_calibrated_poisson_upper": _fraction_record(
            factorwise_calibrated_upper, 16),
        "cap_at_one_active": calibrated_upper == 1,
        "factorwise_capped_upper_is_lower": (
            factorwise_calibrated_upper < calibrated_upper),
        "semantics": (
            "The ordinary component substitutes exact upper endpoints for "
            "all untruncated Poisson factors in the nonnegative "
            "summation-by-parts formula. The calibrated value multiplies "
            "that rational upper enclosure by certified kappa_upper and "
            "caps only the final result at one. The factorwise-capped value "
            "instead caps every calibrated factor at one; it has the same "
            "coverage guarantee and is pointwise no larger."
        ),
    }
    if zero_anchor_eta is not None:
        eta = Fraction(zero_anchor_eta)
        if eta < 1:
            raise ValueError("zero-anchor eta must be at least one")
        p0_lower = brackets[0][0] / n
        p0_upper = brackets[0][1] / n
        if not ordered:
            anchored_upper = p0_upper
        else:
            coefficient_p0 = 1 - eta * ordered[0]
            anchored_upper = coefficient_p0 * (
                p0_upper if coefficient_p0 >= 0 else p0_lower)
            anchored_upper += eta * sum(
                factors_upper[j] * (ordered[j - 1] - ordered[j])
                for j in range(1, maximum_j)
            )
            anchored_upper += (
                eta * factors_upper[maximum_j] * ordered[-1])
        anchored_upper = min(Fraction(1), anchored_upper)
        anchored_factors_upper = [p0_upper]
        anchored_factors_upper.extend(
            min(
                Fraction(1),
                eta * factors_upper[j] - (eta - 1) * p0_lower,
            )
            for j in range(1, maximum_j + 1)
        )
        factorwise_anchored_upper = min(
            anchored_upper,
            from_factor_upper_endpoints(anchored_factors_upper),
        )
        report.update({
            "zero_anchor_eta_used": _fraction_record(eta, 12),
            "zero_anchor_calibrated_poisson_upper": _fraction_record(
                anchored_upper, 16),
            "zero_anchor_factorwise_capped_poisson_upper": _fraction_record(
                factorwise_anchored_upper, 16),
            "zero_anchor_cap_at_one_active": anchored_upper == 1,
            "zero_anchor_factorwise_capped_upper_is_lower": (
                factorwise_anchored_upper < anchored_upper),
            "zero_anchor_semantics": (
                "The anchored value fixes the ordinary untruncated Poisson "
                "zero-taint factor and multiplies only the error increments "
                "by certified eta_upper. Endpoint selection follows each "
                "factor's exact coefficient sign, and only the final result "
                "is capped at one. The factorwise-capped anchored value caps "
                "each effective calibrated factor at one; it has the same "
                "coverage guarantee and is pointwise no larger. Its "
                "reported enclosure is intersected with the final-cap "
                "enclosure to preserve that ordering despite independent "
                "factor brackets."
            ),
            "calibration_selection_warning": (
                "The full-scale or zero-anchor path must be selected before "
                "observing the sample. Taking the smaller result across "
                "paths is not certified. Within a preselected path, the "
                "calibrated-factor-capped variant is separately proved and "
                "may be used."
            ),
        })
    return report


def certify_case(n: int, alpha: Fraction, brackets):
    """Return an exact adjacent-dyadic bracket for the band multiplier.

    This certificate routine is for ``alpha < exp(-1)``, where the ordinary
    terminal Poisson factor already exceeds one and the permitted scalar
    domain starts at one.
    """
    alpha = Fraction(alpha)
    if not _conventional_tail_level(alpha):
        raise ValueError("certificate routine requires alpha < exp(-1)")
    nominal = 1 - alpha
    approximate = locate_multiplier(n, alpha, brackets)
    denominator = 1 << KAPPA_BITS
    with mp.workdps(NUMERICAL_DPS):
        lower_numerator = int(mp.floor(approximate * denominator))
    lower_numerator = max(denominator, lower_numerator)

    # The upper lambda endpoints give an upper bound for Q at a proposed
    # lower kappa.  Move left if numerical location landed just to the right
    # of the exact dyadic crossing.
    while True:
        lower = Fraction(lower_numerator, denominator)
        probability_upper = bolshev_probability(
            _exact_boundaries(n, brackets, lower, "upper"))
        if probability_upper < nominal or lower == 1:
            break
        lower_numerator -= 1

    upper_numerator = lower_numerator + 1
    # The lower lambda endpoints give a lower bound for Q at a proposed
    # upper kappa.  Factor brackets are far narrower than the kappa grid, so
    # the first adjacent point should certify; retain a defensive loop and
    # reject a nonadjacent result below rather than hiding it.
    while True:
        upper = Fraction(upper_numerator, denominator)
        probability_lower = bolshev_probability(
            _exact_boundaries(n, brackets, upper, "lower"))
        if probability_lower >= nominal:
            break
        upper_numerator += 1

    if upper_numerator != lower_numerator + 1:
        raise ArithmeticError(
            "factor uncertainty spans more than one kappa grid cell")
    if lower == 1 and probability_upper >= nominal:
        # Collapse the bracket to the unadjusted rule only when the lower
        # factor endpoints already certify it.  An upper-endpoint check at
        # one is not enough: in a near-frontier case the factor enclosure
        # could straddle the target, in which case the adjacent upper kappa
        # remains the certified choice.
        probability_at_one_lower = bolshev_probability(
            _exact_boundaries(n, brackets, lower, "lower"))
        if probability_at_one_lower >= nominal:
            upper = lower
            probability_lower = probability_at_one_lower

    no_error_lower = upper * brackets[0][0] / n
    no_error_upper = upper * brackets[0][1] / n
    kappa_upper_record = _fraction_record(upper, 12)
    kappa_upper_record["valid_decimal_ceiling_12"] = _decimal_ceiling(
        upper, 12)
    return {
        "n": n,
        "kappa_lower": _fraction_record(lower, 12),
        "kappa_upper": kappa_upper_record,
        "kappa_bracket_width": str(upper - lower),
        "event_probability_upper_at_kappa_lower":
            _fraction_record(probability_upper),
        "event_probability_lower_at_kappa_upper":
            _fraction_record(probability_lower),
        "nominal_confidence": _decimal(nominal),
        "certified_inflation_percent_upper":
            _decimal(100 * (upper - 1), 10),
        "adjusted_zero_taint_factor_lower":
            _fraction_record(no_error_lower, 12),
        "adjusted_zero_taint_factor_upper":
            _fraction_record(no_error_upper, 12),
        "conclusion": (
            "The unadjusted simultaneous-band criterion is rigorously "
            "certified, so kappa_upper equals one."
            if upper == lower else
            "The minimal scalar multiplier certified by the corrected "
            "simultaneous-band criterion lies in the displayed adjacent "
            "dyadic interval; kappa_upper is a rigorous valid choice."
        ),
    }


def certify_zero_anchor_case(n: int, alpha: Fraction, brackets):
    """Certify an adjacent-dyadic bracket for the anchored multiplier."""
    alpha = Fraction(alpha)
    if not _conventional_tail_level(alpha):
        raise ValueError("certificate routine requires alpha < exp(-1)")
    nominal = 1 - alpha
    approximate = locate_zero_anchor_multiplier(n, alpha, brackets)
    denominator = 1 << KAPPA_BITS
    with mp.workdps(NUMERICAL_DPS):
        lower_numerator = int(mp.floor(approximate * denominator))
    lower_numerator = max(denominator, lower_numerator)

    while True:
        lower = Fraction(lower_numerator, denominator)
        probability_upper = bolshev_probability(
            _exact_zero_anchor_boundaries(
                n, brackets, lower, "upper"))
        if probability_upper < nominal or lower == 1:
            break
        lower_numerator -= 1

    upper_numerator = lower_numerator + 1
    while True:
        upper = Fraction(upper_numerator, denominator)
        probability_lower = bolshev_probability(
            _exact_zero_anchor_boundaries(
                n, brackets, upper, "lower"))
        if probability_lower >= nominal:
            break
        upper_numerator += 1

    if upper_numerator != lower_numerator + 1:
        raise ArithmeticError(
            "factor uncertainty spans more than one eta grid cell")
    if lower == 1 and probability_upper >= nominal:
        probability_at_one_lower = bolshev_probability(
            _exact_zero_anchor_boundaries(
                n, brackets, lower, "lower"))
        if probability_at_one_lower >= nominal:
            upper = lower
            probability_lower = probability_at_one_lower

    eta_upper_record = _fraction_record(upper, 12)
    eta_upper_record["valid_decimal_ceiling_12"] = _decimal_ceiling(
        upper, 12)
    no_error_lower = brackets[0][0] / n
    no_error_upper = brackets[0][1] / n
    return {
        "n": n,
        "eta_lower": _fraction_record(lower, 12),
        "eta_upper": eta_upper_record,
        "eta_bracket_width": str(upper - lower),
        "event_probability_upper_at_eta_lower":
            _fraction_record(probability_upper),
        "event_probability_lower_at_eta_upper":
            _fraction_record(probability_lower),
        "nominal_confidence": _decimal(nominal),
        "ordinary_zero_taint_factor_lower":
            _fraction_record(no_error_lower, 12),
        "ordinary_zero_taint_factor_upper":
            _fraction_record(no_error_upper, 12),
        "conclusion": (
            "The unadjusted anchored simultaneous-band criterion is "
            "rigorously certified, so eta_upper equals one."
            if upper == lower else
            "The minimal zero-taint-preserving multiplier under the "
            "anchored simultaneous-band criterion lies in the displayed "
            "adjacent dyadic interval; eta_upper is a rigorous valid choice."
        ),
    }


def build_certificate():
    levels = []
    refined_uniform_cases = []
    for alpha_text, sample_sizes in CASES.items():
        alpha = Fraction(alpha_text)
        maximum_n = max(sample_sizes)
        brackets = exact_poisson_lambda_brackets(
            alpha_text, maximum_n - 1, FACTOR_BITS)
        factor_records = _factor_bracket_records(brackets)
        records = [
            certify_case(n, alpha, brackets)
            for n in sample_sizes
        ]
        zero_anchor_records = [
            certify_zero_anchor_case(n, alpha, brackets)
            for n in sample_sizes
        ]
        refined_kappa, prefix_terms = REFINED_UNIFORM_MULTIPLIERS[alpha_text]
        refined_uniform_cases.append(certify_refined_uniform_multiplier(
            alpha, refined_kappa, prefix_terms, brackets))
        levels.append({
            "alpha": alpha_text,
            "nominal_confidence": _decimal(1 - alpha),
            "sample_sizes": list(sample_sizes),
            "cases": records,
            "zero_anchor_cases": zero_anchor_records,
            "poisson_lambda_brackets": factor_records,
        })

    elementary_uniform_cases = [
        certify_uniform_multiplier(Fraction(alpha_text), kappa)
        for alpha_text, kappa in ELEMENTARY_UNIFORM_MULTIPLIERS.items()
    ]
    return {
        "schema_version": 4,
        "claim": (
            "Both the full-scale and zero-taint-preserving Poisson Stringer "
            "calibrations have corrected simultaneous-band probability at "
            "least 1-alpha when their certified upper multipliers are used. "
            "The displayed elementary and finite-prefix full-scale "
            "multipliers are valid uniformly over every positive sample "
            "size."
        ),
        "scope": (
            "The all-n validity theorem is analytic. These representative "
            "calculations certify the smallest multiplier on each of two "
            "particular one-parameter sufficient-event paths to a 2^-28 "
            "bracket; they do not assert that ordinary Stringer fails or "
            "that either path gives the globally shortest valid confidence "
            "bound. The uniform full-scale cases are sufficient analytic "
            "choices for every positive sample size, not optimality claims."
        ),
        "decimal_policy": DECIMAL_POLICY,
        "arithmetic": (
            "Poisson limits have exact adjacent 64-bit dyadic brackets with "
            "rational exponential endpoint signs; every endpoint used is "
            "stored below. Bolshev probabilities are evaluated with "
            "Fraction. Floating point only proposes multiplier cells, whose "
            "two sides are then certified with monotone opposite endpoint "
            "substitutions. The elementary uniform cases use exact rational "
            "exponential upper bounds and exact powered comparisons. The "
            "refined cases retain one hundred Poisson crossing terms using "
            "exact dyadic upper bounds and enclose the infinite remainder "
            "by an exact geometric tail."
        ),
        "factor_bits": FACTOR_BITS,
        "exponential_series_pairs": EXACT_EXP_PAIRS,
        "kappa_bits": KAPPA_BITS,
        "uniform_probability_bits": UNIFORM_PROBABILITY_BITS,
        "elementary_uniform_full_scale_cases": elementary_uniform_cases,
        "refined_uniform_full_scale_cases": refined_uniform_cases,
        "levels": levels,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    parser.add_argument("--n", type=int)
    parser.add_argument("--alpha")
    parser.add_argument(
        "--taints",
        help="comma-separated decimal taints; zero taints may be omitted",
    )
    arguments = parser.parse_args(argv)

    if (arguments.n is None) != (arguments.alpha is None):
        parser.error("--n and --alpha must be supplied together")
    if arguments.taints is not None and arguments.n is None:
        parser.error("--taints requires --n and --alpha")

    if arguments.n is not None:
        if arguments.n < 1:
            parser.error("--n must be positive")
        try:
            alpha = Fraction(arguments.alpha)
        except (ValueError, ZeroDivisionError) as error:
            parser.error(f"invalid --alpha: {error}")
        if not 0 < alpha < Fraction(1):
            parser.error("--alpha must lie strictly between zero and one")
        if not _conventional_tail_level(alpha):
            parser.error(
                "exact CLI calibration currently requires alpha < exp(-1)")
        try:
            taints = _parse_taints(arguments.taints, arguments.n)
        except ValueError as error:
            parser.error(str(error))
        maximum_j = max(
            arguments.n - 1,
            sum(value != 0 for value in taints),
        )
        brackets = exact_poisson_lambda_brackets(
            arguments.alpha, maximum_j, FACTOR_BITS)
        case = certify_case(arguments.n, alpha, brackets)
        zero_anchor_case = certify_zero_anchor_case(
            arguments.n, alpha, brackets)
        payload = {
            "schema_version": 4,
            "alpha": str(alpha),
            "factor_bits": FACTOR_BITS,
            "kappa_bits": KAPPA_BITS,
            "decimal_policy": DECIMAL_POLICY,
            "case": case,
            "zero_anchor_case": zero_anchor_case,
            "poisson_lambda_brackets": _factor_bracket_records(brackets),
        }
        for alpha_text, uniform_kappa in ELEMENTARY_UNIFORM_MULTIPLIERS.items():
            if Fraction(alpha_text) == alpha:
                payload["elementary_uniform_full_scale_case"] = (
                    certify_uniform_multiplier(alpha, uniform_kappa))
                refined_kappa, prefix_terms = (
                    REFINED_UNIFORM_MULTIPLIERS[alpha_text])
                uniform_brackets = brackets
                if len(uniform_brackets) < prefix_terms:
                    uniform_brackets = exact_poisson_lambda_brackets(
                        alpha_text, prefix_terms - 1, FACTOR_BITS)
                payload["refined_uniform_full_scale_case"] = (
                    certify_refined_uniform_multiplier(
                        alpha, refined_kappa, prefix_terms,
                        uniform_brackets))
                break
        if arguments.taints is not None:
            kappa_record = case["kappa_upper"]
            kappa = Fraction(
                int(kappa_record["numerator"]),
                int(kappa_record["denominator"]),
            )
            eta_record = zero_anchor_case["eta_upper"]
            eta = Fraction(
                int(eta_record["numerator"]),
                int(eta_record["denominator"]),
            )
            payload["report"] = exact_calibrated_report(
                arguments.n, brackets, kappa, taints, eta)
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if arguments.out is None:
            print(rendered, end="")
        else:
            arguments.out.parent.mkdir(parents=True, exist_ok=True)
            arguments.out.write_text(rendered, encoding="utf-8")
        return 0

    certificate = build_certificate()
    output = arguments.out or DEFAULT_OUTPUT
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for case in certificate["elementary_uniform_full_scale_cases"]:
        print(
            "alpha=%s: elementary sample-size-uniform kappa %s; "
            "powered margin %s"
            % (
                case["alpha"],
                case["kappa"]["decimal"],
                case["powered_margin_lower"]["decimal"],
            )
        )
    for case in certificate["refined_uniform_full_scale_cases"]:
        print(
            "alpha=%s: refined sample-size-uniform kappa %s; "
            "crossing margin %s"
            % (
                case["alpha"],
                case["kappa"]["decimal"],
                case["margin_below_alpha_lower"]["decimal"],
            )
        )
    for level in certificate["levels"]:
        for case, anchor_case in zip(
                level["cases"], level["zero_anchor_cases"]):
            print(
                "alpha=%s n=%d: exact kappa bracket [%s/%s, %s/%s]; "
                "valid kappa %s; valid zero-anchor eta %s"
                % (
                    level["alpha"],
                    case["n"],
                    case["kappa_lower"]["numerator"],
                    case["kappa_lower"]["denominator"],
                    case["kappa_upper"]["numerator"],
                    case["kappa_upper"]["denominator"],
                    case["kappa_upper"]["valid_decimal_ceiling_12"],
                    anchor_case["eta_upper"]["valid_decimal_ceiling_12"],
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
