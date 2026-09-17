"""Exact least-favourable completion bounds for a fixed systematic-PPS frame.

For a known ordered frame, one random start has finitely many phase atoms.
After a phase is observed, the audited item identities and taints fix some of
the population taints.  This module computes the largest population target
among all compatible completions for which the observed grid mean is not in
the bottom ``alpha`` probability of its own exact randomization distribution.

The resulting bound is design-valid.  It is also no larger than either the
Markov randomization bound or the deterministic completion bound.  The
optimization is a finite union of rational linear programs: the code
enumerates every inclusion-minimal phase coalition having probability
strictly above ``alpha`` and solves every surviving LP with SymPy's exact
rational simplex.  It checks primal feasibility, dual feasibility, and exact
equality of their objectives.  It never certifies a floating-point MILP
answer.

The enumeration can be exponential in the number of distinct phase atoms.
Explicit search limits therefore fail closed with ``ComplexityLimitError``.
Equal-weight items having identical phase-incidence columns are aggregated
exactly, which makes the 2,000-item benchmark in the accompanying certificate
small despite its practical sample size.

When distinct phase atoms have disjoint item supports, the LP union reduces
to a minimum-excess subset sum.  ``systematic_disjoint_dp.py`` implements the
corresponding exact pseudo-polynomial algorithm and the independent-start
all-zero specialization; this module remains the general overlapping-frame
solver and an independent cross-check of the disjoint formulas.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from sympy import Rational
from sympy.solvers.simplex import InfeasibleLPError, linprog

from stringer import (
    EXACT_FACTOR_BITS,
    exact_binomial_factor_brackets,
    exact_poisson_lambda_brackets,
)
from systematic_pps import (
    completion_bound_from_observed,
    markov_design_bound,
    systematic_phases,
)


class ComplexityLimitError(RuntimeError):
    """Raised before an incomplete coalition search could be reported."""


@dataclass(frozen=True)
class FrameAtom:
    """One distinct grid-incidence atom in the random-start distribution."""

    probability: Fraction
    item_multiplicities: tuple[tuple[int, int], ...]
    source_phase_indices: tuple[int, ...]


@dataclass(frozen=True)
class ExactMinimaxResult:
    """Exact output and a least-favourable completion certificate."""

    alpha: Fraction
    observed_phase_index: int
    observed_atom_index: int
    observed_mean: Fraction
    bound: Fraction
    markov_bound: Fraction
    completion_bound: Fraction
    raw_phase_count: int
    phase_atom_count: int
    minimal_coalition_count: int
    coalition_lps_evaluated: int
    infeasible_coalition_count: int
    coalition_search_nodes: int
    maximizing_coalition: tuple[int, ...] | None
    maximizing_coalition_probability: Fraction | None
    witness_taints: tuple[Fraction, ...] | None
    witness_lower_tail_probability: Fraction | None

    @property
    def prior_hybrid_bound(self) -> Fraction:
        return min(self.markov_bound, self.completion_bound)


@dataclass(frozen=True)
class _VariableGroup:
    """Unknown equal-weight items with one common phase-incidence column."""

    weight: Fraction
    signature: tuple[int, ...]
    item_indices: tuple[int, ...]


def _fractions(values: Iterable[Fraction | int | str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def _validated_weights(
    weights: Sequence[Fraction | int | str],
) -> tuple[tuple[Fraction, ...], Fraction]:
    ws = _fractions(weights)
    if not ws or any(weight <= 0 for weight in ws):
        raise ValueError("all book weights must be positive")
    return ws, sum(ws, Fraction(0))


def systematic_frame_atoms(
    weights: Sequence[Fraction | int | str], sample_size: int
) -> tuple[tuple[FrameAtom, ...], tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """Return exact phase atoms, a raw-to-atom map, and raw item sequences.

    Raw phase intervals with the same item multiplicities define the same
    grid-mean functional and are combined into one atom.  Their probabilities
    are added exactly.  The raw item sequences retain the observed sample
    identities, including duplicate hits.
    """
    ws, _ = _validated_weights(weights)
    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    raw = systematic_phases(ws, (Fraction(0),) * len(ws), sample_size)

    atom_by_signature: dict[tuple[tuple[int, int], ...], int] = {}
    atom_probabilities: list[Fraction] = []
    atom_sources: list[list[int]] = []
    raw_to_atom = []
    raw_items = []
    for raw_index, phase in enumerate(raw):
        signature = tuple(sorted(Counter(phase.item_indices).items()))
        atom_index = atom_by_signature.get(signature)
        if atom_index is None:
            atom_index = len(atom_probabilities)
            atom_by_signature[signature] = atom_index
            atom_probabilities.append(Fraction(0))
            atom_sources.append([])
        atom_probabilities[atom_index] += phase.probability
        atom_sources[atom_index].append(raw_index)
        raw_to_atom.append(atom_index)
        raw_items.append(phase.item_indices)

    atoms = tuple(
        FrameAtom(
            probability=atom_probabilities[index],
            item_multiplicities=signature,
            source_phase_indices=tuple(atom_sources[index]),
        )
        for signature, index in sorted(
            atom_by_signature.items(), key=lambda item: item[1]
        )
    )
    if sum((atom.probability for atom in atoms), Fraction(0)) != 1:
        raise ArithmeticError("frame atoms do not have total probability one")
    return atoms, tuple(raw_to_atom), tuple(raw_items)


def _minimal_probability_coalitions(
    probabilities: Sequence[Fraction],
    alpha: Fraction,
    *,
    max_coalitions: int,
    max_search_nodes: int,
) -> tuple[tuple[tuple[int, ...], ...], int]:
    """Enumerate every inclusion-minimal coalition of mass above alpha."""
    if max_coalitions < 1 or max_search_nodes < 1:
        raise ValueError("complexity limits must be positive")
    if not probabilities or any(probability <= 0 for probability in probabilities):
        raise ValueError("phase probabilities must be positive")
    if sum(probabilities, Fraction(0)) != 1:
        raise ValueError("phase probabilities must sum to one")

    # Large atoms first normally makes a crossing, and therefore a prune,
    # occur early.  Returned coalitions are sorted in the original atom order.
    order = tuple(
        sorted(range(len(probabilities)), key=lambda i: (-probabilities[i], i))
    )
    suffix = [Fraction(0)] * (len(order) + 1)
    for position in range(len(order) - 1, -1, -1):
        suffix[position] = suffix[position + 1] + probabilities[order[position]]

    coalitions: list[tuple[int, ...]] = []
    selected: list[int] = []
    nodes = 0

    def visit(position: int, total: Fraction, smallest: Fraction | None) -> None:
        nonlocal nodes
        nodes += 1
        if nodes > max_search_nodes:
            raise ComplexityLimitError(
                "minimal-coalition search exceeded max_search_nodes="
                f"{max_search_nodes}; no bound was certified"
            )
        if total > alpha:
            if smallest is None or total - smallest > alpha:
                return
            coalitions.append(tuple(sorted(selected)))
            if len(coalitions) > max_coalitions:
                raise ComplexityLimitError(
                    "minimal-coalition search exceeded max_coalitions="
                    f"{max_coalitions}; no bound was certified"
                )
            return
        if position == len(order) or total + suffix[position] <= alpha:
            return

        atom = order[position]
        probability = probabilities[atom]
        selected.append(atom)
        visit(
            position + 1,
            total + probability,
            probability if smallest is None else min(smallest, probability),
        )
        selected.pop()
        visit(position + 1, total, smallest)

    visit(0, Fraction(0), None)
    coalitions.sort()
    return tuple(coalitions), nodes


def _to_sympy(value: Fraction | int) -> Rational:
    value = Fraction(value)
    return Rational(value.numerator, value.denominator)


def _to_fraction(value: object) -> Fraction:
    if isinstance(value, (int, Fraction)):
        return Fraction(value)
    numerator, denominator = value.as_numer_denom()  # type: ignore[attr-defined]
    return Fraction(int(numerator), int(denominator))


def _variable_groups(
    weights: tuple[Fraction, ...],
    atoms: tuple[FrameAtom, ...],
    observed_indices: set[int],
) -> tuple[_VariableGroup, ...]:
    multiplicities = [dict(atom.item_multiplicities) for atom in atoms]
    grouped: dict[tuple[Fraction, tuple[int, ...]], list[int]] = {}
    for item_index, weight in enumerate(weights):
        if item_index in observed_indices:
            continue
        signature = tuple(
            atom_multiplicities.get(item_index, 0)
            for atom_multiplicities in multiplicities
        )
        grouped.setdefault((weight, signature), []).append(item_index)
    return tuple(
        _VariableGroup(weight, signature, tuple(indices))
        for (weight, signature), indices in sorted(
            grouped.items(), key=lambda entry: entry[1][0]
        )
    )


def _solve_coalition_lp(
    *,
    weights: tuple[Fraction, ...],
    total_weight: Fraction,
    atoms: tuple[FrameAtom, ...],
    coalition: tuple[int, ...],
    sample_size: int,
    observed_taints: Mapping[int, Fraction],
    observed_mean: Fraction,
    groups: tuple[_VariableGroup, ...],
) -> tuple[Fraction, tuple[Fraction, ...]] | None:
    """Solve one coalition LP and return target plus a full exact completion."""
    observed_weighted_total = sum(
        (weights[index] * taint for index, taint in observed_taints.items()),
        Fraction(0),
    )
    atom_multiplicities = [dict(atom.item_multiplicities) for atom in atoms]

    rows: list[tuple[int, ...]] = []
    right_sides: list[Fraction] = []
    for atom_index in coalition:
        counts = atom_multiplicities[atom_index]
        known = sum(
            (counts.get(index, 0) * taint for index, taint in observed_taints.items()),
            Fraction(0),
        )
        right_side = sample_size * observed_mean - known
        if right_side < 0:
            return None
        rows.append(tuple(group.signature[atom_index] for group in groups))
        right_sides.append(right_side)

    # Nonnegative rows with zero right side force every incident group to zero.
    forced_zero = {
        group_index
        for row, right_side in zip(rows, right_sides)
        if right_side == 0
        for group_index, coefficient in enumerate(row)
        if coefficient > 0
    }
    active_groups = tuple(
        index for index in range(len(groups)) if index not in forced_zero
    )

    reduced_rows: list[list[Rational]] = []
    reduced_right_sides: list[Rational] = []
    for row, right_side in zip(rows, right_sides):
        reduced = [row[index] for index in active_groups]
        maximum = sum(
            (
                Fraction(coefficient * len(groups[group_index].item_indices), 1)
                for coefficient, group_index in zip(reduced, active_groups)
            ),
            Fraction(0),
        )
        if maximum <= right_side:
            continue
        if not any(reduced):
            return None
        reduced_rows.append([_to_sympy(coefficient) for coefficient in reduced])
        reduced_right_sides.append(_to_sympy(right_side))

    group_sums = [Fraction(0)] * len(groups)
    if active_groups:
        if not reduced_rows:
            for group_index in active_groups:
                group_sums[group_index] = Fraction(
                    len(groups[group_index].item_indices), 1
                )
        else:
            objective_coefficients = [
                groups[group_index].weight for group_index in active_groups
            ]
            objective = [
                -_to_sympy(coefficient) for coefficient in objective_coefficients
            ]
            bounds = [
                (Rational(0), Rational(len(groups[group_index].item_indices)))
                for group_index in active_groups
            ]
            try:
                primal_minimum, solution = linprog(
                    objective,
                    reduced_rows,
                    reduced_right_sides,
                    bounds=bounds,
                )
            except InfeasibleLPError:
                return None
            primal = tuple(_to_fraction(value) for value in solution)
            for row, right_side in zip(reduced_rows, reduced_right_sides):
                exact_row = tuple(_to_fraction(value) for value in row)
                if sum(
                    (coefficient * value for coefficient, value in zip(exact_row, primal)),
                    Fraction(0),
                ) > _to_fraction(right_side):
                    raise ArithmeticError("rational simplex returned an infeasible primal")
            for value, (_, upper) in zip(primal, bounds):
                if not 0 <= value <= _to_fraction(upper):
                    raise ArithmeticError("rational simplex violated a primal bound")

            # An independently checked dual witness turns the simplex output
            # into an exact optimality certificate for this LP.  For
            # max c*x subject to A*x<=b and 0<=x<=u, the dual is
            # min b*y+u*z subject to A^T*y+z>=c, y,z>=0.
            row_count = len(reduced_rows)
            variable_count = len(active_groups)
            dual_objective = list(reduced_right_sides) + [
                upper for _, upper in bounds
            ]
            dual_rows = []
            dual_right_sides = []
            for variable_index, coefficient in enumerate(objective_coefficients):
                dual_rows.append(
                    [-row[variable_index] for row in reduced_rows]
                    + [
                        Rational(-1) if index == variable_index else Rational(0)
                        for index in range(variable_count)
                    ]
                )
                dual_right_sides.append(-_to_sympy(coefficient))
            try:
                dual_minimum, dual_solution = linprog(
                    dual_objective,
                    dual_rows,
                    dual_right_sides,
                    bounds=(0, None),
                )
            except InfeasibleLPError as error:
                raise ArithmeticError("dual LP was unexpectedly infeasible") from error
            dual = tuple(_to_fraction(value) for value in dual_solution)
            if any(value < 0 for value in dual):
                raise ArithmeticError("rational simplex returned a negative dual value")
            dual_y = dual[:row_count]
            dual_z = dual[row_count:]
            for variable_index, coefficient in enumerate(objective_coefficients):
                left = sum(
                    (
                        _to_fraction(reduced_rows[row_index][variable_index])
                        * dual_y[row_index]
                        for row_index in range(row_count)
                    ),
                    Fraction(0),
                ) + dual_z[variable_index]
                if left < coefficient:
                    raise ArithmeticError("rational simplex returned an infeasible dual")
            primal_value = sum(
                (
                    coefficient * value
                    for coefficient, value in zip(objective_coefficients, primal)
                ),
                Fraction(0),
            )
            dual_value = sum(
                (
                    _to_fraction(coefficient) * value
                    for coefficient, value in zip(dual_objective, dual)
                ),
                Fraction(0),
            )
            if primal_value != dual_value:
                raise ArithmeticError("primal and dual exact objectives do not agree")
            if -_to_fraction(primal_minimum) != primal_value:
                raise ArithmeticError("reported primal objective is inconsistent")
            if _to_fraction(dual_minimum) != dual_value:
                raise ArithmeticError("reported dual objective is inconsistent")

            for group_index, value in zip(active_groups, primal):
                group_sums[group_index] = value

    completion: list[Fraction | None] = [None] * len(weights)
    for item_index, taint in observed_taints.items():
        completion[item_index] = taint
    for group, group_sum in zip(groups, group_sums):
        if not 0 <= group_sum <= len(group.item_indices):
            raise ArithmeticError("exact LP returned a group sum outside its bounds")
        remaining = group_sum
        for item_index in group.item_indices:
            taint = min(Fraction(1), remaining)
            completion[item_index] = taint
            remaining -= taint
        if remaining:
            raise ArithmeticError("group-sum expansion did not terminate exactly")
    if any(value is None for value in completion):
        raise ArithmeticError("least-favourable completion omitted an item")
    completed = tuple(value for value in completion if value is not None)

    for atom_index in coalition:
        grid_total = sum(
            (
                multiplicity * completed[item_index]
                for item_index, multiplicity in atoms[
                    atom_index
                ].item_multiplicities
            ),
            Fraction(0),
        )
        if grid_total > sample_size * observed_mean:
            raise ArithmeticError("exact LP witness violates a coalition constraint")

    target = sum(
        (weight * taint for weight, taint in zip(weights, completed)),
        Fraction(0),
    ) / total_weight
    objective_target = (
        observed_weighted_total
        + sum(
            (
                group.weight * group_sum
                for group, group_sum in zip(groups, group_sums)
            ),
            Fraction(0),
        )
    ) / total_weight
    if target != objective_target:
        raise ArithmeticError("expanded completion changed the exact LP objective")
    return target, completed


def exact_minimax_frame_bound(
    weights: Sequence[Fraction | int | str],
    sample_size: int,
    observed_phase_index: int,
    observed_taints: Mapping[int, Fraction | int | str],
    alpha: Fraction | str,
    *,
    max_coalitions: int = 250_000,
    max_search_nodes: int = 2_000_000,
) -> ExactMinimaxResult:
    """Compute the exact finite-frame least-favourable completion bound.

    ``observed_phase_index`` indexes the raw phase tuple returned by
    :func:`systematic_pps.systematic_phases`.  ``observed_taints`` must give
    one value for every distinct item selected in that phase and no others.
    All inputs are converted to :class:`fractions.Fraction`.

    The confidence set uses the exact lower-tail rank

    ``P_u{grid_mean <= observed_mean} > alpha``.

    Its supremum target is returned, with the convention zero if the inverted
    confidence set is empty.  That convention can occur only on an outcome
    the exact test rejects for every compatible completion.
    """
    ws, total_weight = _validated_weights(weights)
    a = Fraction(alpha)
    if not 0 < a < 1:
        raise ValueError("alpha must lie in (0,1)")
    atoms, raw_to_atom, raw_items = systematic_frame_atoms(ws, sample_size)
    if observed_phase_index < 0 or observed_phase_index >= len(raw_items):
        raise ValueError("observed_phase_index is outside the phase enumeration")

    observed = {index: Fraction(value) for index, value in observed_taints.items()}
    required = set(raw_items[observed_phase_index])
    if set(observed) != required:
        raise ValueError(
            "observed_taints keys must equal the distinct sampled item identities"
        )
    if any(index < 0 or index >= len(ws) for index in observed):
        raise ValueError("observed item index is outside the population")
    if any(not 0 <= taint <= 1 for taint in observed.values()):
        raise ValueError("observed taints must lie in [0,1]")

    observed_sequence = tuple(observed[index] for index in raw_items[observed_phase_index])
    observed_mean = sum(observed_sequence, Fraction(0)) / sample_size
    markov = markov_design_bound(observed_sequence, a)
    completion = completion_bound_from_observed(ws, observed)

    probabilities = tuple(atom.probability for atom in atoms)
    coalitions, nodes = _minimal_probability_coalitions(
        probabilities,
        a,
        max_coalitions=max_coalitions,
        max_search_nodes=max_search_nodes,
    )
    groups = _variable_groups(ws, atoms, set(observed))

    best_target = Fraction(0)
    best_coalition: tuple[int, ...] | None = None
    best_completion: tuple[Fraction, ...] | None = None
    infeasible = 0
    for coalition in coalitions:
        solved = _solve_coalition_lp(
            weights=ws,
            total_weight=total_weight,
            atoms=atoms,
            coalition=coalition,
            sample_size=sample_size,
            observed_taints=observed,
            observed_mean=observed_mean,
            groups=groups,
        )
        if solved is None:
            infeasible += 1
            continue
        target, witness = solved
        if (
            best_coalition is None
            or target > best_target
            or (target == best_target and coalition < best_coalition)
        ):
            best_target = target
            best_coalition = coalition
            best_completion = witness

    witness_probability: Fraction | None = None
    coalition_probability: Fraction | None = None
    if best_completion is not None and best_coalition is not None:
        coalition_probability = sum(
            (probabilities[index] for index in best_coalition), Fraction(0)
        )
        witness_probability = sum(
            (
                atom.probability
                for atom in atoms
                if sum(
                    (
                        multiplicity * best_completion[item_index]
                        for item_index, multiplicity in atom.item_multiplicities
                    ),
                    Fraction(0),
                )
                <= sample_size * observed_mean
            ),
            Fraction(0),
        )
        if not coalition_probability > a or not witness_probability > a:
            raise ArithmeticError("least-favourable witness lacks tail mass above alpha")
        for index, taint in observed.items():
            if best_completion[index] != taint:
                raise ArithmeticError("least-favourable witness changed an observation")

    # These are theorem-level dominance checks and also strong implementation
    # invariants.  They use exact comparisons only.
    if best_target > markov or best_target > completion:
        raise ArithmeticError("finite-frame inversion failed its dominance theorem")

    return ExactMinimaxResult(
        alpha=a,
        observed_phase_index=observed_phase_index,
        observed_atom_index=raw_to_atom[observed_phase_index],
        observed_mean=observed_mean,
        bound=best_target,
        markov_bound=markov,
        completion_bound=completion,
        raw_phase_count=len(raw_items),
        phase_atom_count=len(atoms),
        minimal_coalition_count=len(coalitions),
        coalition_lps_evaluated=len(coalitions),
        infeasible_coalition_count=infeasible,
        coalition_search_nodes=nodes,
        maximizing_coalition=best_coalition,
        maximizing_coalition_probability=coalition_probability,
        witness_taints=best_completion,
        witness_lower_tail_probability=witness_probability,
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


def _interval_record(
    lower: Fraction, upper: Fraction, digits: int = 18
) -> dict[str, object]:
    return {
        "lower": _fraction_record(lower, digits),
        "upper": _fraction_record(upper, digits),
    }


def result_record(
    result: ExactMinimaxResult, *, include_witness: bool = True
) -> dict[str, object]:
    """Serialize a solver result without converting exact values to floats."""
    payload: dict[str, object] = {
        "alpha": _fraction_record(result.alpha),
        "observed_phase_index": result.observed_phase_index,
        "observed_atom_index": result.observed_atom_index,
        "observed_mean": _fraction_record(result.observed_mean),
        "exact_minimax_bound": _fraction_record(result.bound),
        "markov_bound": _fraction_record(result.markov_bound),
        "completion_bound": _fraction_record(result.completion_bound),
        "prior_hybrid_bound": _fraction_record(result.prior_hybrid_bound),
        "raw_phase_count": result.raw_phase_count,
        "phase_atom_count": result.phase_atom_count,
        "minimal_coalition_count": result.minimal_coalition_count,
        "coalition_lps_evaluated": result.coalition_lps_evaluated,
        "infeasible_coalition_count": result.infeasible_coalition_count,
        "coalition_search_nodes": result.coalition_search_nodes,
        "maximizing_coalition": (
            list(result.maximizing_coalition)
            if result.maximizing_coalition is not None
            else None
        ),
        "maximizing_coalition_probability": (
            _fraction_record(result.maximizing_coalition_probability)
            if result.maximizing_coalition_probability is not None
            else None
        ),
        "witness_lower_tail_probability": (
            _fraction_record(result.witness_lower_tail_probability)
            if result.witness_lower_tail_probability is not None
            else None
        ),
    }
    if include_witness:
        payload["witness_taints"] = (
            [
                f"{value.numerator}/{value.denominator}"
                for value in result.witness_taints
            ]
            if result.witness_taints is not None
            else None
        )
    return payload


def solve_payload(
    payload: Mapping[str, object], *, include_witness: bool = True
) -> dict[str, object]:
    """Solve the documented JSON input schema and return exact JSON data."""
    required = {
        "weights",
        "sample_size",
        "observed_phase_index",
        "observed_taints",
        "alpha",
    }
    missing = required - set(payload)
    if missing:
        raise ValueError(f"solver input is missing keys: {sorted(missing)}")
    raw_weights = payload["weights"]
    raw_observed = payload["observed_taints"]
    if not isinstance(raw_weights, list):
        raise ValueError("weights must be a JSON list")
    if not isinstance(raw_observed, dict):
        raise ValueError("observed_taints must be a JSON object")
    observed = {int(index): value for index, value in raw_observed.items()}
    max_coalitions = int(payload.get("max_coalitions", 250_000))
    max_search_nodes = int(payload.get("max_search_nodes", 2_000_000))
    result = exact_minimax_frame_bound(
        raw_weights,  # type: ignore[arg-type]
        int(payload["sample_size"]),
        int(payload["observed_phase_index"]),
        observed,  # type: ignore[arg-type]
        str(payload["alpha"]),
        max_coalitions=max_coalitions,
        max_search_nodes=max_search_nodes,
    )
    return {
        "schema_version": 1,
        "method": "exact finite-frame least-favourable completion inversion",
        "input": {
            "weights": [str(value) for value in raw_weights],
            "sample_size": int(payload["sample_size"]),
            "observed_phase_index": int(payload["observed_phase_index"]),
            "observed_taints": {
                str(index): str(value) for index, value in sorted(observed.items())
            },
            "alpha": str(payload["alpha"]),
            "max_coalitions": max_coalitions,
            "max_search_nodes": max_search_nodes,
        },
        "result": result_record(result, include_witness=include_witness),
    }


def build_certificate() -> dict[str, object]:
    """Build the deterministic exact certificate committed with the paper."""
    alpha = Fraction(1, 20)
    sample_size = 100
    phase_count = 20
    item_count = sample_size * phase_count
    weights = (Fraction(1),) * item_count
    atoms, _, raw_items = systematic_frame_atoms(weights, sample_size)
    observed_phase = 0
    observed = {index: Fraction(0) for index in set(raw_items[observed_phase])}
    benchmark = exact_minimax_frame_bound(
        weights, sample_size, observed_phase, observed, alpha
    )
    if benchmark.bound != Fraction(9, 10):
        raise ArithmeticError("2,000-item minimax benchmark did not equal 9/10")
    if benchmark.prior_hybrid_bound != Fraction(19, 20):
        raise ArithmeticError("2,000-item prior hybrid did not equal 19/20")
    if benchmark.minimal_coalition_count != 190:
        raise ArithmeticError("uniform 20-phase frame should have 190 coalitions")
    if benchmark.witness_taints is None:
        raise ArithmeticError("2,000-item benchmark lacks a witness completion")

    binomial_lower, binomial_upper = exact_binomial_factor_brackets(
        sample_size, "0.05"
    )[0]
    lambda_lower, lambda_upper = exact_poisson_lambda_brackets("0.05", 0)[0]
    poisson_lower = lambda_lower / sample_size
    poisson_upper = lambda_upper / sample_size
    if not binomial_upper < benchmark.bound or not poisson_upper < benchmark.bound:
        raise ArithmeticError("ordinary all-zero bound unexpectedly controls overlay")

    # A nonbinary four-item fixture exercises an interior rational LP endpoint.
    small_weights = (1, 1, 1, 1)
    small_atoms, _, small_raw_items = systematic_frame_atoms(small_weights, 1)
    small_observed = {small_raw_items[0][0]: Fraction(1, 2)}
    small = exact_minimax_frame_bound(
        small_weights, 1, 0, small_observed, Fraction(1, 4)
    )
    if small.bound != Fraction(3, 4):
        raise ArithmeticError("four-item rational LP fixture did not equal 3/4")

    witness_zero_count = sum(taint == 0 for taint in benchmark.witness_taints)
    witness_one_count = sum(taint == 1 for taint in benchmark.witness_taints)
    if (witness_zero_count, witness_one_count) != (200, 1800):
        raise ArithmeticError("benchmark witness has an unexpected taint profile")

    return {
        "schema_version": 1,
        "method": (
            "exact inversion of the lower-tail randomization rank of the "
            "systematic grid mean, maximized over compatible taint completions"
        ),
        "arithmetic_status": (
            "phase probabilities, coalitions, LP coefficients, primal-dual "
            "simplex solutions, witnesses, and comparisons are exact rationals; "
            "no floating-point optimizer is used"
        ),
        "scope": (
            "one-start systematic PPS on a fixed known rational frame with "
            "rational observed taints in [0,1]"
        ),
        "algorithm": {
            "phase_step": "combine identical grid-incidence atoms",
            "coalition_step": (
                "enumerate every inclusion-minimal atom coalition with "
                "probability strictly above alpha"
            ),
            "optimization_step": (
                "solve each compatible least-favourable-completion LP by "
                "exact rational simplex, check matching feasible primal and "
                "dual objectives, and take the exact maximum"
            ),
            "complexity_note": (
                "coalition enumeration is exponential in the worst case and "
                "fails closed at explicit limits"
            ),
        },
        "practical_95pct_all_zero_benchmark": {
            "population_items": item_count,
            "unit_book_weight_per_item": 1,
            "sample_size": sample_size,
            "integer_random_start_positions": phase_count,
            "items_per_phase": sample_size,
            "observed_phase_zero_based": observed_phase,
            "distinct_observed_items": len(observed),
            "observed_taint": 0,
            "alpha": _fraction_record(alpha),
            "raw_phase_count": benchmark.raw_phase_count,
            "distinct_phase_atom_count": benchmark.phase_atom_count,
            "phase_probabilities": [
                _fraction_record(atom.probability) for atom in atoms
            ],
            "minimal_coalitions_enumerated": benchmark.minimal_coalition_count,
            "coalition_lps_evaluated": benchmark.coalition_lps_evaluated,
            "infeasible_coalitions": benchmark.infeasible_coalition_count,
            "coalition_search_nodes": benchmark.coalition_search_nodes,
            "observed_grid_mean": _fraction_record(benchmark.observed_mean),
            "markov_bound": _fraction_record(benchmark.markov_bound),
            "completion_bound": _fraction_record(benchmark.completion_bound),
            "prior_hybrid_bound": _fraction_record(benchmark.prior_hybrid_bound),
            "exact_minimax_bound": _fraction_record(benchmark.bound),
            "improvement_over_prior_hybrid": _fraction_record(
                benchmark.prior_hybrid_bound - benchmark.bound
            ),
            "ordinary_binomial_all_zero": _interval_record(
                binomial_lower, binomial_upper
            ),
            "ordinary_poisson_all_zero": _interval_record(
                poisson_lower, poisson_upper
            ),
            "safe_binomial_overlay": _fraction_record(benchmark.bound),
            "safe_poisson_overlay": _fraction_record(benchmark.bound),
            "least_favourable_witness": {
                "maximizing_atom_coalition_zero_based": list(
                    benchmark.maximizing_coalition or ()
                ),
                "coalition_probability": _fraction_record(
                    benchmark.maximizing_coalition_probability or Fraction(0)
                ),
                "exact_lower_tail_probability": _fraction_record(
                    benchmark.witness_lower_tail_probability or Fraction(0)
                ),
                "zero_taints": witness_zero_count,
                "one_taints": witness_one_count,
                "target": _fraction_record(benchmark.bound),
            },
        },
        "four_item_nonbinary_fixture": {
            "weights": list(small_weights),
            "sample_size": 1,
            "phase_probabilities": [
                _fraction_record(atom.probability) for atom in small_atoms
            ],
            "alpha": _fraction_record(Fraction(1, 4)),
            "observed_taint": _fraction_record(Fraction(1, 2)),
            "exact_minimax_bound": _fraction_record(small.bound),
            "markov_bound": _fraction_record(small.markov_bound),
            "completion_bound": _fraction_record(small.completion_bound),
            "least_favourable_lower_tail_probability": _fraction_record(
                small.witness_lower_tail_probability or Fraction(0)
            ),
        },
        "factor_bracket_bits": EXACT_FACTOR_BITS,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        help="solve a frame JSON instead of regenerating the paper certificate",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument(
        "--omit-witness",
        action="store_true",
        help="omit the complete taint witness from generic solver JSON",
    )
    args = parser.parse_args()
    if args.input:
        payload = solve_payload(
            json.loads(args.input.read_text()),
            include_witness=not args.omit_witness,
        )
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
