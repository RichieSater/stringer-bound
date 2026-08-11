"""Exact simultaneous-band certificates for Poisson Stringer factors.

For Poisson limits ``lambda_j`` satisfying

    P(Poisson(lambda_j) <= j) = alpha,

the event

    V_(i:n) <= min(1, lambda_(i-1) / n),  i=1,...,n,

for uniform order statistics is sufficient for Stringer coverage.  This
module proves that the event itself has probability at least ``1-alpha`` in
the conventional-level ranges recorded in ``LEVELS``.

Every numerical-looking step has an exact certificate:

* each ``lambda_j`` is enclosed by adjacent dyadic rationals;
* the Poisson-CDF signs at those endpoints use rational alternating-series
  bounds for ``exp(-x)`` (see :mod:`stringer`);
* Bolshev's order-statistic recursion is evaluated with ``Fraction``; and
* the final comparison with ``1-alpha`` is a rational comparison.

The numerical root locator affects speed only.  This script also certifies
that the sufficient-event probability is below nominal at the next sample
size.  That latter fact is a limitation of this proof route, not evidence of
Stringer undercoverage.
"""

from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from stringer import (
    EXACT_EXP_PAIRS,
    exact_poisson_cdf_bounds,
    exact_poisson_lambda_brackets,
)


LEVELS = {
    "0.10": 8,
    "0.05": 11,
    "0.01": 20,
}
FACTOR_BITS = 80
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "poisson-simultaneous-band-certificate.json"
)


def bolshev_probability(boundaries):
    """Exact ``P(V_(i:m) <= a_i, i=1,...,m)`` for uniform order stats.

    ``boundaries`` is the nondecreasing sequence ``a_1,...,a_m``.  The
    recursion is Bolshev's formula as stated in Shorack and Wellner (1986,
    Section 9.1).
    """
    boundaries = tuple(Fraction(value) for value in boundaries)
    if any(not 0 <= value <= 1 for value in boundaries):
        raise ValueError("boundaries must lie in [0,1]")
    if any(left > right for left, right in zip(
            boundaries, boundaries[1:])):
        raise ValueError("boundaries must be nondecreasing")

    q = [Fraction(1)]
    for m in range(1, len(boundaries) + 1):
        excluded = sum(
            math.comb(m, i)
            * (1 - boundaries[m - i]) ** i
            * q[m - i]
            for i in range(1, m + 1)
        )
        q.append(1 - excluded)
        if not 0 <= q[-1] <= 1:
            raise ArithmeticError("Bolshev recursion left [0,1]")
    return q[-1]


def _event_boundaries(n, brackets, endpoint):
    position = 0 if endpoint == "lower" else 1
    return tuple(min(Fraction(1), brackets[j][position] / n)
                 for j in range(n))


def _decimal(value, digits=16):
    value = Fraction(value)
    with localcontext() as context:
        context.prec = digits + 20
        decimal = Decimal(value.numerator) / Decimal(value.denominator)
        return format(decimal, f".{digits}f")


def _fraction_record(value, digits=16):
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": _decimal(value, digits),
    }


def build_certificate():
    levels = []
    denominator = 1 << FACTOR_BITS
    for alpha_text, frontier in LEVELS.items():
        alpha = Fraction(alpha_text)
        nominal = 1 - alpha
        # The first unresolved size uses lambda_0,...,lambda_frontier.
        brackets = exact_poisson_lambda_brackets(
            alpha_text, frontier, FACTOR_BITS, EXACT_EXP_PAIRS)

        factor_records = []
        for j, (lower, upper) in enumerate(brackets):
            cdf_lower_at_lower, _ = exact_poisson_cdf_bounds(
                lower, j, EXACT_EXP_PAIRS)
            _, cdf_upper_at_upper = exact_poisson_cdf_bounds(
                upper, j, EXACT_EXP_PAIRS)
            if not cdf_lower_at_lower > alpha:
                raise AssertionError("lower Poisson endpoint sign failed")
            if not cdf_upper_at_upper < alpha:
                raise AssertionError("upper Poisson endpoint sign failed")
            factor_records.append({
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
            })

        certified_sizes = []
        frontier_probability = None
        for n in range(1, frontier + 1):
            boundaries = _event_boundaries(n, brackets, "lower")
            probability = bolshev_probability(boundaries)
            margin = probability - nominal
            if margin < 0:
                raise AssertionError(
                    f"simultaneous-band certificate failed at alpha={alpha_text}, n={n}")
            certified_sizes.append({
                "n": n,
                "active_boundaries": sum(value < 1 for value in boundaries),
                "probability_lower_bound": _fraction_record(probability),
                "margin_over_nominal": _fraction_record(margin),
            })
            frontier_probability = probability

        next_n = frontier + 1
        next_lower = bolshev_probability(
            _event_boundaries(next_n, brackets, "lower"))
        next_upper = bolshev_probability(
            _event_boundaries(next_n, brackets, "upper"))
        if not next_upper < nominal:
            raise AssertionError("next-size band-event limitation not certified")

        levels.append({
            "alpha": alpha_text,
            "nominal_confidence": _decimal(nominal),
            "certified_n_range": [1, frontier],
            "frontier_probability_lower_bound":
                _fraction_record(frontier_probability),
            "frontier_margin_over_nominal":
                certified_sizes[-1]["margin_over_nominal"],
            "sizes": certified_sizes,
            "next_size_limitation": {
                "n": next_n,
                "event_probability_lower_bound": _fraction_record(next_lower),
                "event_probability_upper_bound": _fraction_record(next_upper),
                "conclusion": (
                    "The sufficient simultaneous-band event has probability "
                    "below nominal. This does not imply Stringer undercoverage."
                ),
            },
            "poisson_lambda_brackets": factor_records,
        })

    return {
        "schema_version": 1,
        "claim": (
            "For Poisson Stringer factors, the corrected simultaneous-band "
            "event proves distribution-free coverage for every n<=8 at 90%, "
            "every n<=11 at 95%, and every n<=20 at 99% nominal confidence."
        ),
        "arithmetic": (
            "Poisson limits have adjacent 80-bit dyadic brackets. Endpoint "
            "signs use exact rational alternating-series enclosures of "
            "exp(-x); uniform-order-statistic probabilities use exact "
            "Fraction evaluation of Bolshev's recursion."
        ),
        "factor_bits": FACTOR_BITS,
        "exponential_series_pairs": EXACT_EXP_PAIRS,
        "levels": levels,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args(argv)
    certificate = build_certificate()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for level in certificate["levels"]:
        print(
            "alpha=%s: certified n=1..%d; frontier band probability=%s; "
            "next n=%d has band probability below nominal"
            % (
                level["alpha"],
                level["certified_n_range"][1],
                level["frontier_probability_lower_bound"]["decimal"],
                level["next_size_limitation"]["n"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
