"""Exact arithmetic for the SRSWOR conditioning bridge.

The proof is in theory/SRSWOR-CONDITIONING.md.  This module only regenerates
the finite arithmetic reported there and exposes small reusable functions.
It uses O(n) time and O(1) working state apart from integer numerator sizes.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def distinct_probability(population_size: int, sample_size: int) -> Fraction:
    """Return (N)_n/N^n exactly."""
    if population_size < 1:
        raise ValueError("population_size must be positive")
    if sample_size < 1 or sample_size > population_size:
        raise ValueError("require 1 <= sample_size <= population_size")
    result = Fraction(1, 1)
    for k in range(sample_size):
        result *= Fraction(population_size - k, population_size)
    return result


def adjusted_tail(alpha: Fraction, population_size: int, sample_size: int) -> Fraction:
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0,1)")
    return alpha * distinct_probability(population_size, sample_size)


def minimum_population_size(sample_size: int, required_collision_factor: Fraction) -> int:
    """Smallest N>=n for which c_(N,n) is at least the required factor."""
    if not 0 < required_collision_factor < 1:
        raise ValueError("required_collision_factor must lie in (0,1)")
    population_size = sample_size
    while distinct_probability(population_size, sample_size) < required_collision_factor:
        population_size += 1
    return population_size


def decimal_string(value: Fraction, digits: int = 18) -> str:
    with localcontext() as ctx:
        ctx.prec = digits + 8
        number = Decimal(value.numerator) / Decimal(value.denominator)
        return f"{number:.{digits}f}"


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal_18": decimal_string(value),
    }


def build_certificate() -> dict[str, object]:
    thresholds = []
    requirements = (
        ("95pct_iid_to_90pct_srswor", Fraction(1, 2)),
        ("99pct_iid_to_90pct_srswor", Fraction(1, 10)),
        ("99pct_iid_to_95pct_srswor", Fraction(1, 5)),
    )
    for n in range(4, 8):
        row: dict[str, object] = {"sample_size": n}
        for label, required in requirements:
            N = minimum_population_size(n, required)
            current = distinct_probability(N, n)
            previous = distinct_probability(N - 1, n) if N > n else None
            row[label] = {
                "required_factor": _fraction_record(required),
                "minimum_population_size": N,
                "factor_at_minimum": _fraction_record(current),
                "factor_at_previous_population_size": (
                    _fraction_record(previous) if previous is not None else None
                ),
            }
        thresholds.append(row)

    examples = []
    for N, n in ((1_000, 25), (10_000, 25), (10_000, 100),
                 (100_000, 100), (100_000, 200), (1_000_000, 200)):
        c = distinct_probability(N, n)
        examples.append({
            "population_size": N,
            "sample_size": n,
            "collision_free_probability_decimal_18": decimal_string(c),
            "adjusted_tail_at_alpha_0.05_decimal_18": decimal_string(Fraction(1, 20) * c),
        })

    return {
        "schema_version": 1,
        "theorem": "SRSWOR noncoverage <= eta / ((N)_n/N^n)",
        "arithmetic_status": "exact for threshold rows; rounded display only for illustrative large-N rows",
        "fixed_level_population_thresholds": thresholds,
        "illustrative_collision_factors": examples,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    parser.add_argument("--population-size", type=int)
    parser.add_argument("--sample-size", type=int)
    parser.add_argument("--alpha", default="0.05")
    args = parser.parse_args()

    if args.population_size is not None or args.sample_size is not None:
        if args.population_size is None or args.sample_size is None:
            parser.error("--population-size and --sample-size must be supplied together")
        alpha = Fraction(args.alpha)
        c = distinct_probability(args.population_size, args.sample_size)
        payload: dict[str, object] = {
            "population_size": args.population_size,
            "sample_size": args.sample_size,
            "collision_free_probability": _fraction_record(c),
            "adjusted_tail": _fraction_record(alpha * c),
        }
    else:
        payload = build_certificate()

    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
