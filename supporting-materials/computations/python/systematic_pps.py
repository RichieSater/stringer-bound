"""Exact finite-population calculations for one-start systematic PPS/MUS.

The mathematical proofs are in ``theory/SYSTEMATIC-PPS.md``.  This module
keeps the design calculation separate from the i.i.d. coverage code:

* item book amounts and taints are represented by exact fractions;
* all random-start phases are enumerated exactly;
* a small binary census locates the first unit-book-value counterexample at
  95% confidence; and
* the practical-size periodic example is certified with the repository's
  exact binomial- and Poisson-factor enclosures.

Endpoints of item intervals have design probability zero in the continuous
model.  Each returned phase is therefore represented by an open interval and
evaluated at its exact rational midpoint.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from stringer import (
    EXACT_FACTOR_BITS,
    exact_binomial_factor_brackets,
    exact_poisson_cdf_bounds,
    exact_poisson_lambda_brackets,
)


@dataclass(frozen=True)
class SystematicPhase:
    """One positive-probability random-start phase."""

    start_lower: Fraction
    start_upper: Fraction
    probability: Fraction
    item_indices: tuple[int, ...]
    taints: tuple[Fraction, ...]


def _fractions(values: Iterable[Fraction | int | str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def _validated_population(
    weights: Sequence[Fraction | int | str],
    taints: Sequence[Fraction | int | str],
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...], Fraction]:
    ws = _fractions(weights)
    ts = _fractions(taints)
    if not ws or len(ws) != len(ts):
        raise ValueError("weights and taints must have the same positive length")
    if any(weight <= 0 for weight in ws):
        raise ValueError("all book weights must be positive")
    if any(not 0 <= taint <= 1 for taint in ts):
        raise ValueError("taints must lie in [0,1]")
    return ws, ts, sum(ws, Fraction(0))


def finite_population_target(
    weights: Sequence[Fraction | int | str],
    taints: Sequence[Fraction | int | str],
) -> Fraction:
    """Return total overstatement divided by total recorded amount."""
    ws, ts, total = _validated_population(weights, taints)
    numerator = sum(
        (weight * taint for weight, taint in zip(ws, ts)), Fraction(0)
    )
    return numerator / total


def systematic_phases(
    weights: Sequence[Fraction | int | str],
    taints: Sequence[Fraction | int | str],
    sample_size: int,
) -> tuple[SystematicPhase, ...]:
    """Enumerate the exact phases of normalized one-start systematic PPS.

    On a book-value line of length ``W``, the design draws
    ``R, R+W/n, ..., R+(n-1)W/n`` for ``R`` uniform on ``[0,W/n)``.  The
    item containing a point supplies its taint.  Breakpoints are the item
    boundaries reduced modulo the sampling interval.
    """
    ws, ts, total = _validated_population(weights, taints)
    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    interval = total / sample_size

    cumulative: list[Fraction] = []
    running = Fraction(0)
    for weight in ws:
        running += weight
        cumulative.append(running)

    breakpoints = {Fraction(0), interval}
    for boundary in cumulative[:-1]:
        residue = boundary % interval
        if residue:
            breakpoints.add(residue)
    ordered = sorted(breakpoints)

    def item_at(position: Fraction) -> int:
        for index, upper in enumerate(cumulative):
            if position < upper:
                return index
        raise ArithmeticError("systematic point fell outside the population line")

    phases = []
    for lower, upper in zip(ordered, ordered[1:]):
        if lower == upper:
            continue
        start = (lower + upper) / 2
        indices = tuple(
            item_at(start + draw * interval) for draw in range(sample_size)
        )
        phases.append(
            SystematicPhase(
                start_lower=lower,
                start_upper=upper,
                probability=(upper - lower) / interval,
                item_indices=indices,
                taints=tuple(ts[index] for index in indices),
            )
        )
    if sum((phase.probability for phase in phases), Fraction(0)) != 1:
        raise ArithmeticError("random-start phases do not partition the interval")
    return tuple(phases)


def phase_sample_mean(phase: SystematicPhase) -> Fraction:
    return sum(phase.taints, Fraction(0)) / len(phase.taints)


def systematic_mean_expectation(phases: Sequence[SystematicPhase]) -> Fraction:
    return sum(
        (phase.probability * phase_sample_mean(phase) for phase in phases),
        Fraction(0),
    )


def markov_design_bound(
    sample_taints: Sequence[Fraction | int | str], alpha: Fraction | str
) -> Fraction:
    """The design-valid randomization bound ``1-alpha+alpha*mean(taints)``."""
    ts = _fractions(sample_taints)
    a = Fraction(alpha)
    if not ts:
        raise ValueError("sample_taints must be nonempty")
    if any(not 0 <= taint <= 1 for taint in ts):
        raise ValueError("taints must lie in [0,1]")
    if not 0 < a < 1:
        raise ValueError("alpha must lie in (0,1)")
    return 1 - a + a * sum(ts, Fraction(0)) / len(ts)


def completion_bound(
    weights: Sequence[Fraction | int | str],
    taints: Sequence[Fraction | int | str],
    sampled_item_indices: Sequence[int],
) -> Fraction:
    """Deterministic bound using audited identities and their book weights.

    Every distinct sampled item contributes its known weighted taint.  Every
    unobserved item is conservatively assigned taint one.  Repeated hits on
    the same item are counted once.
    """
    ws, ts, _ = _validated_population(weights, taints)
    selected = set(sampled_item_indices)
    if any(index < 0 or index >= len(ws) for index in selected):
        raise ValueError("sampled item index is outside the population")
    return completion_bound_from_observed(
        ws, {index: ts[index] for index in selected}
    )


def completion_bound_from_observed(
    weights: Sequence[Fraction | int | str],
    observed_taints: Mapping[int, Fraction | int | str],
) -> Fraction:
    """Completion bound using only book weights and observed item taints.

    This is the operational form of :func:`completion_bound`: unobserved
    taints are not inputs because the bound assigns them their maximum value
    one.  Mapping keys are zero-based population item indices.
    """
    ws = _fractions(weights)
    if not ws or any(weight <= 0 for weight in ws):
        raise ValueError("all book weights must be positive")
    total = sum(ws, Fraction(0))
    numerator = total
    for index, raw_taint in observed_taints.items():
        if index < 0 or index >= len(ws):
            raise ValueError("observed item index is outside the population")
        taint = Fraction(raw_taint)
        if not 0 <= taint <= 1:
            raise ValueError("observed taints must lie in [0,1]")
        numerator -= ws[index] * (1 - taint)
    return numerator / total


def hybrid_design_bound(
    weights: Sequence[Fraction | int | str],
    taints: Sequence[Fraction | int | str],
    phase: SystematicPhase,
    alpha: Fraction | str,
) -> Fraction:
    """Minimum of the randomization and deterministic completion bounds."""
    return min(
        markov_design_bound(phase.taints, alpha),
        completion_bound(weights, taints, phase.item_indices),
    )


def periodic_population(
    sample_size: int,
    high_fraction: Fraction | str,
    low_taint: Fraction | str = Fraction(0),
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    """Return the two-item-per-cell periodic counterexample family.

    In each of ``sample_size`` systematic cells, a leading item of book
    weight ``q/n`` has taint one and the remaining item of weight
    ``(1-q)/n`` has taint ``low_taint``.
    """
    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    q = Fraction(high_fraction)
    y = Fraction(low_taint)
    if not 0 < q < 1:
        raise ValueError("high_fraction must lie in (0,1)")
    if not 0 <= y <= 1:
        raise ValueError("low_taint must lie in [0,1]")
    weights = []
    taints = []
    for _ in range(sample_size):
        weights.extend((q / sample_size, (1 - q) / sample_size))
        taints.extend((Fraction(1), y))
    return tuple(weights), tuple(taints)


def _decimal_string(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits + 12
        decimal = Decimal(value.numerator) / Decimal(value.denominator)
        return f"{decimal:.{digits}f}"


def _fraction_record(value: Fraction, digits: int = 18) -> dict[str, str]:
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        f"decimal_{digits}": _decimal_string(value, digits),
    }


def _interval_record(
    lower: Fraction, upper: Fraction, digits: int = 18
) -> dict[str, object]:
    return {
        "lower": _fraction_record(lower, digits),
        "upper": _fraction_record(upper, digits),
    }


def _binary_phase_counts(
    mask: int, total_units: int, sample_size: int
) -> tuple[int, ...]:
    interval = total_units // sample_size
    return tuple(
        sum((mask >> (start + draw * interval)) & 1 for draw in range(sample_size))
        for start in range(interval)
    )


def exhaustive_binary_census(
    alpha: str = "0.05", maximum_total_units: int = 9
) -> dict[str, object]:
    """Exhaust all binary unit-item populations through a finite size.

    For every factorization ``N=n*m``, the population consists of ``N``
    ordered unit-book-value items and the random start is uniform over the
    ``m`` integer positions in a systematic interval.  Because a binary
    sample with ``j`` ones has Stringer bound ``p_n(j)``, exact factor
    brackets decide every phase.
    """
    a = Fraction(alpha)
    if not 0 < a < 1:
        raise ValueError("alpha must lie in (0,1)")
    if maximum_total_units < 2:
        raise ValueError("maximum_total_units must be at least two")

    first_failure: dict[str, object] | None = None
    total_population_design_pairs = 0
    total_phase_comparisons = 0
    by_size = []
    for total_units in range(2, maximum_total_units + 1):
        designs = []
        for sample_size in range(2, total_units + 1):
            if total_units % sample_size:
                continue
            interval = total_units // sample_size
            factors = exact_binomial_factor_brackets(sample_size, alpha)
            violating = 0
            smallest_coverage = Fraction(1)
            first_mask: int | None = None
            first_failed_starts: tuple[int, ...] = ()
            for mask in range(1 << total_units):
                target = Fraction(mask.bit_count(), total_units)
                failed_starts = []
                for start, count in enumerate(
                    _binary_phase_counts(mask, total_units, sample_size)
                ):
                    lower, upper = factors[count]
                    if upper < target:
                        failed_starts.append(start)
                    elif lower < target <= upper:
                        raise ArithmeticError(
                            "factor bracket is too wide to decide a census phase"
                        )
                coverage = 1 - Fraction(len(failed_starts), interval)
                smallest_coverage = min(smallest_coverage, coverage)
                if coverage < 1 - a:
                    violating += 1
                    if first_mask is None:
                        first_mask = mask
                        first_failed_starts = tuple(failed_starts)

            population_count = 1 << total_units
            total_population_design_pairs += population_count
            total_phase_comparisons += population_count * interval
            design_record: dict[str, object] = {
                "sample_size": sample_size,
                "systematic_interval_units": interval,
                "binary_populations_checked": population_count,
                "phase_comparisons_checked": population_count * interval,
                "violating_population_count": violating,
                "smallest_coverage": _fraction_record(smallest_coverage),
            }
            if first_mask is not None:
                bits = tuple((first_mask >> index) & 1 for index in range(total_units))
                counts = _binary_phase_counts(first_mask, total_units, sample_size)
                target = Fraction(sum(bits), total_units)
                witness = {
                    "ordered_taints": list(bits),
                    "target": _fraction_record(target),
                    "sample_one_counts_by_start": list(counts),
                    "failed_zero_based_starts": list(first_failed_starts),
                    "coverage": _fraction_record(
                        1 - Fraction(len(first_failed_starts), interval)
                    ),
                    "all_zero_factor_bracket": _interval_record(*factors[0]),
                }
                design_record["first_witness"] = witness
                if first_failure is None:
                    first_failure = {
                        "total_units": total_units,
                        "sample_size": sample_size,
                        "systematic_interval_units": interval,
                        **witness,
                    }
            designs.append(design_record)
        by_size.append({"total_units": total_units, "designs": designs})

    return {
        "scope": (
            "all ordered binary-taint populations of unit-book-value items, "
            "for every integer systematic design N=n*m in the stated range"
        ),
        "alpha": alpha,
        "maximum_total_units": maximum_total_units,
        "population_design_pairs_checked": total_population_design_pairs,
        "phase_comparisons_checked": total_phase_comparisons,
        "by_total_units": by_size,
        "first_failure": first_failure,
    }


def build_certificate() -> dict[str, object]:
    sample_size = 100
    q = Fraction(1, 10)
    weights, taints = periodic_population(sample_size, q)
    phases = systematic_phases(weights, taints, sample_size)
    target = finite_population_target(weights, taints)
    if target != q or systematic_mean_expectation(phases) != target:
        raise ArithmeticError("periodic population failed its design identity")

    low_phase = next(phase for phase in phases if not any(phase.taints))
    high_phase = next(phase for phase in phases if all(t == 1 for t in phase.taints))
    if low_phase.probability != 1 - q or high_phase.probability != q:
        raise ArithmeticError("periodic population did not have the expected phases")

    # The single exact inequality P(Pois(10)<=10)>0.1 also implies
    # P(Pois(10)<=100)>alpha for every alpha below, hence lambda_100>10.
    high_cdf_lower, high_cdf_upper = exact_poisson_cdf_bounds(Fraction(10), 10)
    if high_cdf_lower <= Fraction(1, 10):
        raise ArithmeticError("Poisson high-phase coverage witness failed")

    levels = []
    for alpha in ("0.10", "0.05", "0.01"):
        a = Fraction(alpha)
        binomial_lower, binomial_upper = exact_binomial_factor_brackets(
            sample_size, alpha
        )[0]
        lambda_lower, lambda_upper = exact_poisson_lambda_brackets(alpha, 0)[0]
        poisson_lower = lambda_lower / sample_size
        poisson_upper = lambda_upper / sample_size
        if not binomial_upper < target or not poisson_upper < target:
            raise ArithmeticError("all-zero ordinary bound did not fail")
        if not low_phase.probability > a:
            raise ArithmeticError("noncoverage probability is not above alpha")

        markov_low = markov_design_bound(low_phase.taints, a)
        completion_low = completion_bound(weights, taints, low_phase.item_indices)
        hybrid_low = min(markov_low, completion_low)
        if hybrid_low != target:
            raise ArithmeticError("hybrid bound should be exact on the low phase")

        levels.append(
            {
                "alpha": alpha,
                "nominal_confidence": _fraction_record(1 - a),
                "ordinary_exact_coverage_both_factor_conventions": (
                    _fraction_record(q)
                ),
                "ordinary_coverage_shortfall": _fraction_record(1 - a - q),
                "ordinary_noncoverage_excess_over_alpha": _fraction_record(1 - q - a),
                "all_zero_phase": {
                    "probability": _fraction_record(low_phase.probability),
                    "binomial_stringer": _interval_record(
                        binomial_lower, binomial_upper
                    ),
                    "poisson_stringer": _interval_record(
                        poisson_lower, poisson_upper
                    ),
                    "markov_design_bound": _fraction_record(markov_low),
                    "completion_bound": _fraction_record(completion_low),
                    "hybrid_design_bound": _fraction_record(hybrid_low),
                    "safe_binomial_output": _fraction_record(target),
                    "safe_poisson_output": _fraction_record(target),
                    "safe_binomial_uplift": _interval_record(
                        target - binomial_upper, target - binomial_lower
                    ),
                    "safe_poisson_uplift": _interval_record(
                        target - poisson_upper, target - poisson_lower
                    ),
                },
            }
        )

    return {
        "schema_version": 1,
        "design": (
            "one uniform random start on one systematic PPS interval; "
            "n equally spaced book-unit points"
        ),
        "estimand": "sum_i(weight_i*taint_i)/sum_i(weight_i)",
        "arithmetic_status": (
            "exact rational phase enumeration and exact 80-bit dyadic factor "
            "enclosures; decimals are display only"
        ),
        "small_binary_census_at_95pct": exhaustive_binary_census(),
        "periodic_100_draw_example": {
            "sample_size": sample_size,
            "total_recorded_amount_integer_version": 1000,
            "systematic_interval_integer_version": 10,
            "item_pattern_per_cell_integer_version": [
                {"recorded_amount": 1, "taint": 1},
                {"recorded_amount": 9, "taint": 0},
            ],
            "number_of_items": len(weights),
            "target": _fraction_record(target),
            "phases": [
                {
                    "sample_taint": int(phase.taints[0]),
                    "sample_taint_repetitions": sample_size,
                    "probability": _fraction_record(phase.probability),
                    "distinct_items_sampled": len(set(phase.item_indices)),
                    "sampled_book_weight": _fraction_record(
                        sum(
                            (
                                weights[index]
                                for index in set(phase.item_indices)
                            ),
                            Fraction(0),
                        )
                    ),
                }
                for phase in phases
            ],
            "exact_design_unbiasedness": {
                "expected_sample_mean": _fraction_record(
                    systematic_mean_expectation(phases)
                ),
                "finite_population_target": _fraction_record(target),
            },
            "poisson_high_phase_witness": {
                "statement": (
                    "P(Poisson(10)<=10)>0.10, so lambda_100(alpha)>10 "
                    "for alpha in {0.10,0.05,0.01}"
                ),
                "cdf_bracket": _interval_record(high_cdf_lower, high_cdf_upper),
            },
            "levels": levels,
        },
        "factor_bracket_bits": EXACT_FACTOR_BITS,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
