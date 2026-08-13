"""Exact scalar calibrations for an all-n Poisson Stringer bound.

At the conventional levels ``alpha < exp(-1)``, let
``lambda_j(alpha)`` be the usual Poisson confidence limits and let

    Q_n(kappa) = P(V_(i:n) <= min(1, kappa*lambda_(i-1)/n), all i),

for uniform order statistics.  The corrected simultaneous-band lemma shows
that ``kappa * SB_P`` has coverage at least ``Q_n(kappa)``.  This module
certifies adjacent dyadic brackets for the smallest ``kappa >= 1`` at which
``Q_n(kappa) >= 1-alpha``.

Floating-point arithmetic only locates a candidate dyadic cell.  Every
reported conclusion uses:

* exact rational brackets for each ``lambda_j``;
* exact rational evaluation of Bolshev's recursion; and
* opposite endpoint choices to enclose ``Q_n(kappa)`` monotonically.

The fixed certificate covers representative audit sample sizes. The
analytic theorem also handles other confidence levels by imposing the
additional terminal-factor lower bound recorded in the theory note; this
script is deliberately limited to the conventional regime where that lower
bound is one. The all-sample-size theorem is recorded in
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
    exact_poisson_lambda_brackets,
)


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


CASES = {
    "0.10": (25, 50, 100, 200),
    "0.05": (25, 50, 100, 200),
    "0.01": (25, 50, 100, 200),
}
FACTOR_BITS = 64
KAPPA_BITS = 28
NUMERICAL_DPS = 80
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "poisson-band-calibration-certificate.json"
)
DECIMAL_POLICY = (
    "Fields named decimal are display approximations to the exact numerator "
    "and denominator. The field kappa_upper.valid_decimal_ceiling_12 is "
    "rounded upward and is itself a valid decimal multiplier."
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


def _exact_boundaries(n: int, brackets, kappa: Fraction,
                      lambda_endpoint: str):
    position = 0 if lambda_endpoint == "lower" else 1
    return tuple(
        min(Fraction(1), kappa * brackets[j][position] / n)
        for j in range(n)
    )


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


def exact_calibrated_report(n: int, brackets, kappa: Fraction, taints):
    """Return a rational upper enclosure of the calibrated audit result.

    The Stringer expression is written by summation by parts, so every
    Poisson factor has a nonnegative coefficient. Substituting the upper
    endpoint of each exact factor bracket therefore encloses ordinary
    Poisson Stringer from above. Multiplication by the certified ``kappa``
    and capping at one gives a directly reportable rational upper bound.
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
    if not ordered:
        ordinary_upper = factors_upper[0]
    else:
        ordinary_upper = factors_upper[0] * (1 - ordered[0])
        ordinary_upper += sum(
            factors_upper[j] * (ordered[j - 1] - ordered[j])
            for j in range(1, maximum_j)
        )
        ordinary_upper += factors_upper[maximum_j] * ordered[-1]

    calibrated_upper = min(Fraction(1), Fraction(kappa) * ordinary_upper)
    return {
        "factor_convention": "untruncated_poisson_factors_then_final_cap",
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
        "cap_at_one_active": calibrated_upper == 1,
        "semantics": (
            "The ordinary component substitutes exact upper endpoints for "
            "all untruncated Poisson factors in the nonnegative "
            "summation-by-parts formula. The calibrated value multiplies "
            "that rational upper enclosure by certified kappa_upper and "
            "caps only the final result at one."
        ),
    }


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


def build_certificate():
    levels = []
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
        levels.append({
            "alpha": alpha_text,
            "nominal_confidence": _decimal(1 - alpha),
            "sample_sizes": list(sample_sizes),
            "cases": records,
            "poisson_lambda_brackets": factor_records,
        })

    return {
        "schema_version": 1,
        "claim": (
            "Scaling every ordinary Poisson Stringer factor by the certified "
            "kappa_upper gives a corrected simultaneous-band event of "
            "probability at least 1-alpha at each listed case."
        ),
        "scope": (
            "The all-n validity theorem is analytic. These representative "
            "calculations certify the smallest scalar inflation under this "
            "particular sufficient-event criterion to a 2^-28 bracket; they "
            "do not assert that ordinary Stringer fails or that this is the "
            "globally shortest valid confidence bound."
        ),
        "decimal_policy": DECIMAL_POLICY,
        "arithmetic": (
            "Poisson limits have exact adjacent 64-bit dyadic brackets with "
            "rational exponential endpoint signs; every endpoint used is "
            "stored below. Bolshev probabilities are evaluated with "
            "Fraction. Floating point only proposes a kappa cell, whose two "
            "sides are then certified with opposite lambda endpoints."
        ),
        "factor_bits": FACTOR_BITS,
        "exponential_series_pairs": EXACT_EXP_PAIRS,
        "kappa_bits": KAPPA_BITS,
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
        payload = {
            "schema_version": 1,
            "alpha": str(alpha),
            "factor_bits": FACTOR_BITS,
            "kappa_bits": KAPPA_BITS,
            "decimal_policy": DECIMAL_POLICY,
            "case": case,
            "poisson_lambda_brackets": _factor_bracket_records(brackets),
        }
        if arguments.taints is not None:
            kappa_record = case["kappa_upper"]
            kappa = Fraction(
                int(kappa_record["numerator"]),
                int(kappa_record["denominator"]),
            )
            payload["report"] = exact_calibrated_report(
                arguments.n, brackets, kappa, taints)
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
    for level in certificate["levels"]:
        for case in level["cases"]:
            print(
                "alpha=%s n=%d: exact kappa bracket [%s/%s, %s/%s]; "
                "valid decimal choice %s"
                % (
                    level["alpha"],
                    case["n"],
                    case["kappa_lower"]["numerator"],
                    case["kappa_lower"]["denominator"],
                    case["kappa_upper"]["numerator"],
                    case["kappa_upper"]["denominator"],
                    case["kappa_upper"]["valid_decimal_ceiling_12"],
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
