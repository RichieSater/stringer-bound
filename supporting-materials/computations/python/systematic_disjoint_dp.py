"""Exact disjoint-phase finite-frame bounds by pseudo-polynomial DP.

For a phase-partition sampling design, the population target is the
probability-weighted average of the phase means.  If phase ``j`` is observed
with mean ``y``, exact finite-frame inversion therefore reduces to finding a
minimum-probability coalition that contains ``j`` and has probability
strictly above ``alpha``.  This module solves that minimum-excess subset-sum
problem exactly after clearing rational denominators.

The same dynamic program handles the all-zero observation from independent
random starts.  If the coalition has phase probability ``rho``, the exact
lower-tail probability for ``r`` starts is ``rho**r``.  All comparisons,
back-pointers, bounds, and certificate records use exact integer or rational
arithmetic.  The scaled probability denominator controls the
pseudo-polynomial complexity, so an explicit limit fails closed before a
partial answer can be reported.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path
from typing import Iterable, Mapping, Sequence


class DPComplexityLimitError(RuntimeError):
    """Raised before a DP whose exact state space exceeds its declared limit."""


@dataclass(frozen=True)
class ExactCoalitionDPResult:
    """Exact minimum-mass coalition and the DP evidence used to obtain it."""

    probabilities: tuple[Fraction, ...]
    alpha: Fraction
    independent_starts: int
    required_indices: tuple[int, ...]
    scaled_total_mass: int
    scaled_phase_masses: tuple[int, ...]
    strict_scaled_threshold: int
    required_scaled_mass: int
    added_scaled_mass: int
    coalition: tuple[int, ...]
    coalition_scaled_mass: int
    coalition_probability: Fraction
    product_lower_tail_probability: Fraction
    reachable_sum_count: int
    transitions_considered: int


@dataclass(frozen=True)
class DisjointPhaseBoundResult:
    """One-start disjoint-phase bound for an arbitrary observed phase mean."""

    observed_phase_index: int
    observed_mean: Fraction
    bound: Fraction
    coalition: ExactCoalitionDPResult


@dataclass(frozen=True)
class MultipleStartAllZeroBoundResult:
    """All-zero disjoint-phase bound for independent random starts."""

    observed_phase_indices: tuple[int, ...]
    independent_starts: int
    bound: Fraction
    coalition: ExactCoalitionDPResult


@dataclass(frozen=True)
class DisjointSystematicFrameBoundResult:
    """Specialized DP result after verifying an actual systematic-PPS frame."""

    sample_size: int
    observed_phase_index: int
    observed_atom_index: int
    observed_mean: Fraction
    raw_phase_count: int
    phase_atom_count: int
    bound: Fraction
    coalition: ExactCoalitionDPResult


@dataclass(frozen=True)
class EqualPhaseFormulaResult:
    """Closed-form specialization for equiprobable disjoint phases."""

    phase_count: int
    alpha: Fraction
    observed_mean: Fraction
    required_low_phase_count: int
    lower_tail_probability: Fraction
    bound: Fraction


@dataclass(frozen=True)
class EqualPhaseMultipleStartResult:
    """Closed-form all-zero specialization for equal phases and many starts."""

    phase_count: int
    alpha: Fraction
    independent_starts: int
    distinct_observed_phase_count: int
    required_zero_phase_count: int
    lower_tail_probability: Fraction
    bound: Fraction


def _fractions(values: Iterable[Fraction | int | str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def _validated_probabilities(
    probabilities: Sequence[Fraction | int | str],
) -> tuple[Fraction, ...]:
    values = _fractions(probabilities)
    if not values or any(value <= 0 for value in values):
        raise ValueError("phase probabilities must be positive")
    if sum(values, Fraction(0)) != 1:
        raise ValueError("phase probabilities must sum exactly to one")
    return values


def _scaled_probability_masses(
    probabilities: tuple[Fraction, ...], *, max_scaled_total: int
) -> tuple[tuple[int, ...], int]:
    if max_scaled_total < 1:
        raise ValueError("max_scaled_total must be positive")
    denominator = 1
    for probability in probabilities:
        denominator = lcm(denominator, probability.denominator)
    masses = tuple(
        probability.numerator * (denominator // probability.denominator)
        for probability in probabilities
    )
    common = 0
    for mass in masses:
        common = gcd(common, mass)
    if common > 1:
        masses = tuple(mass // common for mass in masses)
        denominator //= common
    if sum(masses) != denominator:
        raise ArithmeticError("scaled phase masses do not sum to their denominator")
    if denominator > max_scaled_total:
        raise DPComplexityLimitError(
            "exact disjoint-phase DP requires scaled total mass "
            f"{denominator}, exceeding max_scaled_total={max_scaled_total}; "
            "no bound was certified"
        )
    return masses, denominator


def _strict_scaled_threshold(
    scaled_total: int, alpha: Fraction, independent_starts: int
) -> int:
    """Smallest integer c with (c/scaled_total)**starts > alpha."""
    if independent_starts < 1:
        raise ValueError("independent_starts must be positive")
    low = 0
    high = scaled_total
    while low < high:
        middle = (low + high) // 2
        left = pow(middle, independent_starts) * alpha.denominator
        right = alpha.numerator * pow(scaled_total, independent_starts)
        if left > right:
            high = middle
        else:
            low = middle + 1
    if not Fraction(low, scaled_total) ** independent_starts > alpha:
        raise ArithmeticError("strict product-tail threshold failed its inequality")
    if low and Fraction(low - 1, scaled_total) ** independent_starts > alpha:
        raise ArithmeticError("strict product-tail threshold was not minimal")
    return low


def minimum_strict_tail_coalition(
    probabilities: Sequence[Fraction | int | str],
    alpha: Fraction | int | str,
    required_indices: Sequence[int],
    *,
    independent_starts: int = 1,
    max_scaled_total: int = 2_000_000,
) -> ExactCoalitionDPResult:
    """Return an exact minimum-mass required coalition with rho**r > alpha.

    Rational phase probabilities are scaled to positive integer masses.  A
    one-dimensional 0--1 subset-sum DP enumerates every reachable added mass
    from phases that are not required by the observation.  The first
    reachable mass crossing the strict product-tail threshold is optimal.
    """
    values = _validated_probabilities(probabilities)
    a = Fraction(alpha)
    if not 0 < a < 1:
        raise ValueError("alpha must lie in (0,1)")
    if independent_starts < 1:
        raise ValueError("independent_starts must be positive")
    required = tuple(sorted(set(int(index) for index in required_indices)))
    if not required:
        raise ValueError("at least one observed phase must be required")
    if required[0] < 0 or required[-1] >= len(values):
        raise ValueError("required phase index is outside the phase list")

    masses, scaled_total = _scaled_probability_masses(
        values, max_scaled_total=max_scaled_total
    )
    threshold = _strict_scaled_threshold(scaled_total, a, independent_starts)
    required_mass = sum(masses[index] for index in required)
    optional = tuple(index for index in range(len(values)) if index not in required)
    optional_total = scaled_total - required_mass

    reachable = bytearray(optional_total + 1)
    reachable[0] = 1
    parent_item = [-1] * (optional_total + 1)
    parent_sum = [-1] * (optional_total + 1)
    current_max = 0
    transitions = 0
    for index in optional:
        mass = masses[index]
        for subtotal in range(current_max, -1, -1):
            transitions += 1
            if not reachable[subtotal]:
                continue
            new_total = subtotal + mass
            if not reachable[new_total]:
                reachable[new_total] = 1
                parent_item[new_total] = index
                parent_sum[new_total] = subtotal
        current_max += mass

    needed = max(0, threshold - required_mass)
    added_mass = next(
        (mass for mass in range(needed, optional_total + 1) if reachable[mass]),
        None,
    )
    if added_mass is None:
        raise ArithmeticError("all phases failed to cross a threshold below one")

    chosen = list(required)
    remainder = added_mass
    while remainder:
        index = parent_item[remainder]
        previous = parent_sum[remainder]
        if index < 0 or previous < 0 or previous >= remainder:
            raise ArithmeticError("subset-sum back-pointer chain is invalid")
        chosen.append(index)
        remainder = previous
    coalition = tuple(sorted(chosen))
    coalition_mass = sum(masses[index] for index in coalition)
    probability = sum((values[index] for index in coalition), Fraction(0))
    if coalition_mass != required_mass + added_mass:
        raise ArithmeticError("reconstructed coalition changed the DP mass")
    if probability != Fraction(coalition_mass, scaled_total):
        raise ArithmeticError("scaled and rational coalition masses disagree")
    product_tail = probability**independent_starts
    if not product_tail > a:
        raise ArithmeticError("DP coalition does not cross the strict tail threshold")
    if coalition_mass < threshold:
        raise ArithmeticError("DP coalition lies below the integer threshold")
    for mass in range(needed, added_mass):
        if reachable[mass]:
            raise ArithmeticError("DP did not return the minimum reachable crossing")

    return ExactCoalitionDPResult(
        probabilities=values,
        alpha=a,
        independent_starts=independent_starts,
        required_indices=required,
        scaled_total_mass=scaled_total,
        scaled_phase_masses=masses,
        strict_scaled_threshold=threshold,
        required_scaled_mass=required_mass,
        added_scaled_mass=added_mass,
        coalition=coalition,
        coalition_scaled_mass=coalition_mass,
        coalition_probability=probability,
        product_lower_tail_probability=product_tail,
        reachable_sum_count=sum(reachable),
        transitions_considered=transitions,
    )


def exact_disjoint_phase_bound(
    probabilities: Sequence[Fraction | int | str],
    observed_phase_index: int,
    observed_mean: Fraction | int | str,
    alpha: Fraction | int | str,
    *,
    max_scaled_total: int = 2_000_000,
) -> DisjointPhaseBoundResult:
    """Compute the one-start closed-form bound using the exact coalition DP."""
    y = Fraction(observed_mean)
    if not 0 <= y <= 1:
        raise ValueError("observed_mean must lie in [0,1]")
    coalition = minimum_strict_tail_coalition(
        probabilities,
        alpha,
        (observed_phase_index,),
        max_scaled_total=max_scaled_total,
    )
    bound = 1 - (1 - y) * coalition.coalition_probability
    if not 0 <= bound <= 1:
        raise ArithmeticError("disjoint-phase formula left the unit interval")
    return DisjointPhaseBoundResult(
        observed_phase_index=observed_phase_index,
        observed_mean=y,
        bound=bound,
        coalition=coalition,
    )


def exact_multiple_start_all_zero_bound(
    probabilities: Sequence[Fraction | int | str],
    observed_phase_indices: Sequence[int],
    independent_starts: int,
    alpha: Fraction | int | str,
    *,
    max_scaled_total: int = 2_000_000,
) -> MultipleStartAllZeroBoundResult:
    """Compute the exact all-zero product-design bound for disjoint phases."""
    observed = tuple(int(index) for index in observed_phase_indices)
    if len(observed) != independent_starts:
        raise ValueError(
            "observed_phase_indices must contain one phase for every start"
        )
    required = tuple(sorted(set(observed)))
    coalition = minimum_strict_tail_coalition(
        probabilities,
        alpha,
        required,
        independent_starts=independent_starts,
        max_scaled_total=max_scaled_total,
    )
    bound = 1 - coalition.coalition_probability
    return MultipleStartAllZeroBoundResult(
        observed_phase_indices=observed,
        independent_starts=independent_starts,
        bound=bound,
        coalition=coalition,
    )


def exact_disjoint_systematic_frame_bound(
    weights: Sequence[Fraction | int | str],
    sample_size: int,
    observed_phase_index: int,
    observed_taints: Mapping[int, Fraction | int | str],
    alpha: Fraction | int | str,
    *,
    max_scaled_total: int = 2_000_000,
) -> DisjointSystematicFrameBoundResult:
    """Verify phase disjointness and apply the specialized DP to a frame.

    The generic frame enumerator is imported lazily to keep this module's
    probability-mass algorithm independent.  No formula is returned unless
    every population item belongs to exactly one distinct phase atom.
    """
    from systematic_minimax import systematic_frame_atoms

    ws = _fractions(weights)
    if not ws or any(weight <= 0 for weight in ws):
        raise ValueError("all book weights must be positive")
    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    atoms, raw_to_atom, raw_items = systematic_frame_atoms(ws, sample_size)
    if observed_phase_index < 0 or observed_phase_index >= len(raw_items):
        raise ValueError("observed_phase_index is outside the phase enumeration")

    owners: list[int | None] = [None] * len(ws)
    for atom_index, atom in enumerate(atoms):
        for item_index, multiplicity in atom.item_multiplicities:
            if multiplicity < 1:
                raise ArithmeticError("phase atom contains a nonpositive hit count")
            owner = owners[item_index]
            if owner is not None and owner != atom_index:
                raise ValueError(
                    "specialized DP requires pairwise disjoint phase-atom supports"
                )
            owners[item_index] = atom_index
    if any(owner is None for owner in owners):
        raise ValueError("specialized DP requires phase supports to cover the frame")

    observed = {int(index): Fraction(value) for index, value in observed_taints.items()}
    required_items = set(raw_items[observed_phase_index])
    if set(observed) != required_items:
        raise ValueError(
            "observed_taints keys must equal the distinct sampled item identities"
        )
    if any(not 0 <= value <= 1 for value in observed.values()):
        raise ValueError("observed taints must lie in [0,1]")
    observed_sequence = tuple(
        observed[index] for index in raw_items[observed_phase_index]
    )
    observed_mean = sum(observed_sequence, Fraction(0)) / sample_size
    atom_index = raw_to_atom[observed_phase_index]
    specialized = exact_disjoint_phase_bound(
        tuple(atom.probability for atom in atoms),
        atom_index,
        observed_mean,
        alpha,
        max_scaled_total=max_scaled_total,
    )
    return DisjointSystematicFrameBoundResult(
        sample_size=sample_size,
        observed_phase_index=observed_phase_index,
        observed_atom_index=atom_index,
        observed_mean=observed_mean,
        raw_phase_count=len(raw_items),
        phase_atom_count=len(atoms),
        bound=specialized.bound,
        coalition=specialized.coalition,
    )


def equal_disjoint_phase_formula(
    phase_count: int,
    observed_mean: Fraction | int | str,
    alpha: Fraction | int | str,
) -> EqualPhaseFormulaResult:
    """Evaluate 1-(floor(alpha*m)+1)(1-y)/m exactly."""
    if phase_count < 1:
        raise ValueError("phase_count must be positive")
    y = Fraction(observed_mean)
    a = Fraction(alpha)
    if not 0 <= y <= 1:
        raise ValueError("observed_mean must lie in [0,1]")
    if not 0 < a < 1:
        raise ValueError("alpha must lie in (0,1)")
    count = (a.numerator * phase_count) // a.denominator + 1
    if not 1 <= count <= phase_count:
        raise ArithmeticError("equal-phase crossing count is outside its range")
    probability = Fraction(count, phase_count)
    if not probability > a:
        raise ArithmeticError("equal-phase formula did not cross alpha strictly")
    if count > 1 and Fraction(count - 1, phase_count) > a:
        raise ArithmeticError("equal-phase crossing count was not minimal")
    return EqualPhaseFormulaResult(
        phase_count=phase_count,
        alpha=a,
        observed_mean=y,
        required_low_phase_count=count,
        lower_tail_probability=probability,
        bound=1 - probability * (1 - y),
    )


def equal_multiple_start_all_zero_formula(
    phase_count: int,
    independent_starts: int,
    distinct_observed_phase_count: int,
    alpha: Fraction | int | str,
) -> EqualPhaseMultipleStartResult:
    """Evaluate the exact equal-phase all-zero multiple-start formula."""
    if phase_count < 1 or independent_starts < 1:
        raise ValueError("phase_count and independent_starts must be positive")
    if not 1 <= distinct_observed_phase_count <= min(
        phase_count, independent_starts
    ):
        raise ValueError(
            "distinct_observed_phase_count must be between one and "
            "min(phase_count, independent_starts)"
        )
    a = Fraction(alpha)
    if not 0 < a < 1:
        raise ValueError("alpha must lie in (0,1)")
    count = distinct_observed_phase_count
    while Fraction(count, phase_count) ** independent_starts <= a:
        count += 1
    if count > phase_count:
        raise ArithmeticError("all phases failed to cross alpha below one")
    tail = Fraction(count, phase_count) ** independent_starts
    if count > distinct_observed_phase_count:
        previous = Fraction(count - 1, phase_count) ** independent_starts
        if previous > a:
            raise ArithmeticError("multiple-start phase count was not minimal")
    return EqualPhaseMultipleStartResult(
        phase_count=phase_count,
        alpha=a,
        independent_starts=independent_starts,
        distinct_observed_phase_count=distinct_observed_phase_count,
        required_zero_phase_count=count,
        lower_tail_probability=tail,
        bound=1 - Fraction(count, phase_count),
    )


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


def _coalition_record(result: ExactCoalitionDPResult) -> dict[str, object]:
    return {
        "independent_starts": result.independent_starts,
        "required_phase_indices_zero_based": list(result.required_indices),
        "scaled_total_mass": result.scaled_total_mass,
        "scaled_phase_masses": list(result.scaled_phase_masses),
        "strict_scaled_threshold": result.strict_scaled_threshold,
        "required_scaled_mass": result.required_scaled_mass,
        "added_scaled_mass": result.added_scaled_mass,
        "minimizing_coalition_zero_based": list(result.coalition),
        "coalition_scaled_mass": result.coalition_scaled_mass,
        "coalition_probability": _fraction_record(result.coalition_probability),
        "product_lower_tail_probability": _fraction_record(
            result.product_lower_tail_probability
        ),
        "reachable_sum_count": result.reachable_sum_count,
        "transitions_considered": result.transitions_considered,
    }


def solve_payload(payload: Mapping[str, object]) -> dict[str, object]:
    """Solve either documented one-start or all-zero multiple-start JSON."""
    if "weights" in payload:
        weights = payload["weights"]
        observed = payload.get("observed_taints")
        if not isinstance(weights, list):
            raise ValueError("weights must be a JSON list")
        if not isinstance(observed, dict):
            raise ValueError("observed_taints must be a JSON object")
        result = exact_disjoint_systematic_frame_bound(
            weights,  # type: ignore[arg-type]
            int(payload["sample_size"]),
            int(payload["observed_phase_index"]),
            {int(index): value for index, value in observed.items()},
            str(payload.get("alpha")),
            max_scaled_total=int(payload.get("max_scaled_total", 2_000_000)),
        )
        return {
            "schema_version": 1,
            "mode": "verified phase-partitioned systematic-PPS frame",
            "bound": _fraction_record(result.bound),
            "observed_mean": _fraction_record(result.observed_mean),
            "observed_phase_index": result.observed_phase_index,
            "observed_atom_index": result.observed_atom_index,
            "raw_phase_count": result.raw_phase_count,
            "phase_atom_count": result.phase_atom_count,
            "coalition_dp": _coalition_record(result.coalition),
        }
    probabilities = payload.get("probabilities")
    if not isinstance(probabilities, list):
        raise ValueError("probabilities must be a JSON list")
    alpha = str(payload.get("alpha"))
    max_scaled_total = int(payload.get("max_scaled_total", 2_000_000))
    if "observed_mean" in payload:
        result = exact_disjoint_phase_bound(
            probabilities,  # type: ignore[arg-type]
            int(payload["observed_phase_index"]),
            str(payload["observed_mean"]),
            alpha,
            max_scaled_total=max_scaled_total,
        )
        return {
            "schema_version": 1,
            "mode": "one-start disjoint-phase bound",
            "bound": _fraction_record(result.bound),
            "observed_mean": _fraction_record(result.observed_mean),
            "coalition_dp": _coalition_record(result.coalition),
        }
    observed = payload.get("observed_phase_indices")
    if not isinstance(observed, list):
        raise ValueError(
            "multiple-start input requires observed_phase_indices as a JSON list"
        )
    result = exact_multiple_start_all_zero_bound(
        probabilities,  # type: ignore[arg-type]
        [int(index) for index in observed],
        int(payload["independent_starts"]),
        alpha,
        max_scaled_total=max_scaled_total,
    )
    return {
        "schema_version": 1,
        "mode": "independent-start disjoint-phase all-zero bound",
        "bound": _fraction_record(result.bound),
        "observed_phase_indices_zero_based": list(result.observed_phase_indices),
        "coalition_dp": _coalition_record(result.coalition),
    }


def build_certificate() -> dict[str, object]:
    """Build deterministic exact examples for every theorem and algorithm."""
    alpha = Fraction(1, 20)
    equal = equal_disjoint_phase_formula(20, 0, alpha)
    equal_dp = exact_disjoint_phase_bound((Fraction(1, 20),) * 20, 0, 0, alpha)
    if equal.bound != Fraction(9, 10) or equal_dp.bound != equal.bound:
        raise ArithmeticError("equal-phase formula and DP did not agree at 9/10")
    frame_weights = (1,) * 2000
    frame_observed = {20 * index: 0 for index in range(100)}
    equal_frame = exact_disjoint_systematic_frame_bound(
        frame_weights, 100, 0, frame_observed, alpha
    )
    if equal_frame.bound != equal.bound or equal_frame.phase_atom_count != 20:
        raise ArithmeticError("verified 2,000-item frame did not match the formula")

    unequal_probabilities = (
        Fraction(1, 10),
        Fraction(1, 5),
        Fraction(3, 10),
        Fraction(2, 5),
    )
    unequal = exact_disjoint_phase_bound(
        unequal_probabilities, 0, Fraction(1, 3), Fraction(7, 20)
    )
    if unequal.bound != Fraction(11, 15):
        raise ArithmeticError("unequal-phase DP fixture did not equal 11/15")
    if unequal.coalition.coalition != (0, 2):
        raise ArithmeticError("unequal-phase DP chose an unexpected coalition")

    # The reduction maps a subset-sum instance (a_i,T) to phase weights
    # (1,2a_1,...,2a_k), an all-zero observed first phase, and alpha=2T/W.
    subset_yes = (3, 5, 6)
    target_yes = 8
    weights_yes = (1,) + tuple(2 * value for value in subset_yes)
    total_yes = sum(weights_yes)
    reduction_yes = exact_disjoint_phase_bound(
        tuple(Fraction(weight, total_yes) for weight in weights_yes),
        0,
        0,
        Fraction(2 * target_yes, total_yes),
    )
    if reduction_yes.coalition.coalition_probability != Fraction(17, 29):
        raise ArithmeticError("yes-instance reduction did not attain target sum")

    subset_no = (3, 5, 9)
    target_no = 7
    weights_no = (1,) + tuple(2 * value for value in subset_no)
    total_no = sum(weights_no)
    reduction_no = exact_disjoint_phase_bound(
        tuple(Fraction(weight, total_no) for weight in weights_no),
        0,
        0,
        Fraction(2 * target_no, total_no),
    )
    if reduction_no.coalition.coalition_probability != Fraction(17, 35):
        raise ArithmeticError("no-instance reduction did not attain next larger sum")

    multiple_records = []
    expected_bounds = (Fraction(9, 10), Fraction(3, 4), Fraction(3, 5))
    for starts, expected in enumerate(expected_bounds, start=1):
        formula = equal_multiple_start_all_zero_formula(20, starts, starts, alpha)
        dp = exact_multiple_start_all_zero_bound(
            (Fraction(1, 20),) * 20,
            tuple(range(starts)),
            starts,
            alpha,
        )
        if formula.bound != expected or dp.bound != expected:
            raise ArithmeticError("multiple-start formula and DP did not agree")
        multiple_records.append(
            {
                "independent_starts": starts,
                "distinct_observed_phases": starts,
                "required_zero_phases": formula.required_zero_phase_count,
                "exact_bound": _fraction_record(formula.bound),
                "exact_product_lower_tail_probability": _fraction_record(
                    formula.lower_tail_probability
                ),
                "dp": _coalition_record(dp.coalition),
            }
        )

    return {
        "schema_version": 1,
        "method": (
            "exact denominator-cleared subset-sum dynamic programming for "
            "disjoint phase coalitions"
        ),
        "arithmetic_status": (
            "all phase masses, strict product-tail comparisons, DP states, "
            "back-pointers, and bound calculations are exact"
        ),
        "equal_20_phase_one_start": {
            "phase_count": equal.phase_count,
            "alpha": _fraction_record(equal.alpha),
            "observed_mean": _fraction_record(equal.observed_mean),
            "required_low_phases": equal.required_low_phase_count,
            "exact_bound": _fraction_record(equal.bound),
            "verified_systematic_frame": {
                "population_items": len(frame_weights),
                "sample_size": equal_frame.sample_size,
                "raw_phase_count": equal_frame.raw_phase_count,
                "phase_atom_count": equal_frame.phase_atom_count,
                "observed_phase_zero_based": equal_frame.observed_phase_index,
                "observed_atom_zero_based": equal_frame.observed_atom_index,
                "exact_bound": _fraction_record(equal_frame.bound),
            },
            "dp": _coalition_record(equal_dp.coalition),
        },
        "unequal_phase_fixture": {
            "phase_probabilities": [
                _fraction_record(value) for value in unequal_probabilities
            ],
            "observed_phase_zero_based": 0,
            "observed_mean": _fraction_record(Fraction(1, 3)),
            "alpha": _fraction_record(Fraction(7, 20)),
            "exact_bound": _fraction_record(unequal.bound),
            "dp": _coalition_record(unequal.coalition),
        },
        "subset_sum_reduction_fixtures": {
            "yes_instance": {
                "integers": list(subset_yes),
                "target": target_yes,
                "phase_weights": list(weights_yes),
                "alpha": _fraction_record(Fraction(2 * target_yes, total_yes)),
                "minimum_coalition_probability": _fraction_record(
                    reduction_yes.coalition.coalition_probability
                ),
                "exact_bound": _fraction_record(reduction_yes.bound),
                "decodes_exact_target": True,
            },
            "no_instance": {
                "integers": list(subset_no),
                "target": target_no,
                "phase_weights": list(weights_no),
                "alpha": _fraction_record(Fraction(2 * target_no, total_no)),
                "minimum_coalition_probability": _fraction_record(
                    reduction_no.coalition.coalition_probability
                ),
                "exact_bound": _fraction_record(reduction_no.bound),
                "decodes_exact_target": False,
            },
        },
        "equal_20_phase_independent_starts_all_zero": multiple_records,
        "complexity": {
            "time": "O(KD)",
            "memory": "O(D)",
            "parameters": (
                "K is the number of phase atoms and D is their least common "
                "scaled total mass after exact denominator clearing"
            ),
            "failure_mode": (
                "the implementation raises DPComplexityLimitError before "
                "allocation when D exceeds max_scaled_total"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        help="solve documented disjoint-phase JSON instead of the certificate",
    )
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    payload = (
        solve_payload(json.loads(args.input.read_text()))
        if args.input
        else build_certificate()
    )
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
