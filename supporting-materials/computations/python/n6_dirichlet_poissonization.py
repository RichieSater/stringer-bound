"""Directed certificate for the complete ``n=6`` Poissonization comparison.

The final six-positive boundary is partitioned into a fixed collection of
exact-rational root boxes.  Each root is independent, so the expensive
directed-real-ball subdivision can be regenerated in parallel without making
the mathematical output depend on process scheduling.
"""

from __future__ import annotations

import hashlib
import math
import multiprocessing as mp
import os
from fractions import Fraction

from flint import arb, ctx, fmpq


ARB_PRECISION_BITS = 160
SIX_BOUNDARY_BASE_ROOT_COUNT = 1024
SIX_BOUNDARY_CORNER_REFINEMENT = 64
SIX_BOUNDARY_CORNER_REFINEMENT_LEVELS = 2
SIX_BOUNDARY_ROOT_COUNT = (
    SIX_BOUNDARY_BASE_ROOT_COUNT
    - SIX_BOUNDARY_CORNER_REFINEMENT_LEVELS
    + SIX_BOUNDARY_CORNER_REFINEMENT_LEVELS
    * SIX_BOUNDARY_CORNER_REFINEMENT
)


def _fraction_record(value: Fraction) -> dict[str, str]:
    value = Fraction(value)
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": f"{float(value):.16g}",
    }


def _exp_taylor_lower(value: Fraction, degree: int) -> Fraction:
    value = Fraction(value)
    return sum(
        (
            value**order / math.factorial(order)
            for order in range(degree + 1)
        ),
        Fraction(),
    )


def _arb_exact(value: Fraction) -> arb:
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def _tightened_upper(
    lower: tuple[Fraction, ...],
    upper: tuple[Fraction, ...],
    weights: tuple[int, ...],
) -> tuple[Fraction, ...]:
    result = list(upper)
    for index, weight in enumerate(weights):
        available = 6 - sum(
            weights[other] * lower[other]
            for other in range(len(weights))
            if other != index
        )
        result[index] = min(result[index], available / weight)
    return tuple(result)


def _cumulative_knots(
    parameters: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    knots = []
    value = Fraction()
    for parameter in parameters:
        value += parameter
        knots.append(value)
    return tuple(knots)


def _feasible_center(
    lower: tuple[Fraction, ...],
    upper: tuple[Fraction, ...],
    weights: tuple[int, ...],
) -> tuple[Fraction, ...]:
    center = tuple(
        (lower[index] + upper[index]) / 2
        for index in range(len(weights))
    )
    weighted_center = sum(
        weights[index] * center[index] for index in range(len(weights))
    )
    if weighted_center > 6:
        weighted_lower = sum(
            weights[index] * lower[index]
            for index in range(len(weights))
        )
        scale = (Fraction(6) - weighted_lower) / (
            weighted_center - weighted_lower
        )
        center = tuple(
            lower[index] + scale * (center[index] - lower[index])
            for index in range(len(weights))
        )
    return center


def _f3(value: Fraction) -> arb:
    if value == 0:
        return arb(0)
    point = _arb_exact(value)
    result = point**3 * (-6 / point).exp()
    if value > 1:
        result -= (point - 1) ** 6 / point**3
    return result


def _f4(value: Fraction) -> arb:
    if value == 0:
        return arb(0)
    point = _arb_exact(value)
    result = point**4 * (-6 / point).exp()
    if value > 1:
        result -= (point - 1) ** 6 / point**2
    return result


def _divided_difference(
    knots: tuple[Fraction, ...],
    function,
) -> arb:
    values = [function(knot) for knot in knots]
    for order in range(1, len(knots)):
        values = [
            (values[index + 1] - values[index])
            / _arb_exact(knots[index + order] - knots[index])
            for index in range(len(values) - 1)
        ]
    return values[0]


def _tail_coefficients(
    total: int,
) -> dict[
    int,
    tuple[
        tuple[tuple[Fraction, int], ...],
        tuple[tuple[Fraction, int], ...],
    ],
]:
    """Return exact individual-coordinate Beta-tail moment polynomials."""
    result = {}
    for high_count in range(1, total):
        low_count = total - high_count
        constant = Fraction(
            math.factorial(total - 1),
            math.factorial(high_count - 1)
            * math.factorial(low_count - 1),
        )
        low_terms = tuple(
            (
                constant
                * Fraction(
                    (-1) ** order * math.comb(high_count - 1, order),
                    low_count * (low_count + order + 1),
                ),
                low_count + order + 1,
            )
            for order in range(high_count)
        )
        high_terms = tuple(
            (
                constant
                * Fraction(
                    (-1) ** order * math.comb(high_count, order),
                    high_count * (low_count + order),
                ),
                low_count + order,
            )
            for order in range(high_count + 1)
        )
        result[high_count] = (low_terms, high_terms)
    return result


_TAIL_COEFFICIENTS_5 = _tail_coefficients(5)
_WORKER_ARB_CONSTANTS = None


def _worker_arb_constants():
    """Create process-local 160-bit constants after setting Arb precision."""
    global _WORKER_ARB_CONSTANTS
    if _WORKER_ARB_CONSTANTS is None:
        coefficient_balls = {
            high_count: (
                tuple(
                    (_arb_exact(coefficient), power)
                    for coefficient, power in low_terms
                ),
                tuple(
                    (_arb_exact(coefficient), power)
                    for coefficient, power in high_terms
                ),
            )
            for high_count, (low_terms, high_terms) in (
                _TAIL_COEFFICIENTS_5.items()
            )
        }
        _WORKER_ARB_CONSTANTS = {
            "tail_coefficients": coefficient_balls,
            "base": _arb_exact(Fraction(1, 5)),
            "envelope": tuple(
                _arb_exact(value)
                for value in (
                    Fraction(1, 100),
                    Fraction(1, 5),
                    Fraction(3, 10),
                    Fraction(13, 50),
                    Fraction(41, 100),
                )
            ),
            "q4_cap": _arb_exact(Fraction(149, 500)),
        }
    return _WORKER_ARB_CONSTANTS


def _tail_moments_arb(high_count: int, radius: Fraction) -> tuple[arb, arb]:
    constants = _worker_arb_constants()
    radius_ball = _arb_exact(radius)

    def evaluate(terms):
        return sum(
            (
                coefficient * radius_ball**power
                for coefficient, power in terms
            ),
            arb(0),
        )

    low_terms, high_terms = constants["tail_coefficients"][high_count]
    return evaluate(low_terms), evaluate(high_terms)


def _event_coordinate_moments(
    knot_lower: tuple[Fraction, ...],
    knot_upper: tuple[Fraction, ...],
    event_lower: Fraction,
    event_upper: Fraction,
) -> list[arb]:
    """Upper-bound ``E[W_i 1{event_lower<=U<=event_upper}]``.

    The bounds are the minimum of the trivial Dirichlet moment and every
    admissible low/high-knot split from either side of the event interval.
    """
    total = len(knot_lower)
    if total != 5:
        raise AssertionError("the q4 radius uses five Dirichlet coordinates")
    base = _worker_arb_constants()["base"]
    moments = [base] * total
    if knot_upper[-1] < event_lower or knot_lower[0] > event_upper:
        return [arb(0)] * total

    largest_upper = min(knot_upper[-1], Fraction(6))
    if largest_upper > event_lower:
        for low_count in range(1, total):
            low_cap = knot_upper[low_count - 1]
            if low_cap >= event_lower:
                continue
            radius = (largest_upper - event_lower) / (
                largest_upper - low_cap
            )
            low_moment, high_moment = _tail_moments_arb(
                total - low_count,
                radius,
            )
            for index in range(low_count):
                if low_moment < moments[index]:
                    moments[index] = low_moment
            for index in range(low_count, total):
                if high_moment < moments[index]:
                    moments[index] = high_moment

    for low_count in range(1, total):
        high_floor = knot_lower[low_count]
        if high_floor <= event_upper:
            continue
        radius = 1 - event_upper / high_floor
        low_tail, high_tail = _tail_moments_arb(
            total - low_count,
            radius,
        )
        low_moment = base - low_tail
        high_moment = base - high_tail
        for index in range(low_count):
            if low_moment < moments[index]:
                moments[index] = low_moment
        for index in range(low_count, total):
            if high_moment < moments[index]:
                moments[index] = high_moment
    return moments


def _q4_knot_lipschitz_bounds(
    knot_lower: tuple[Fraction, ...],
    knot_upper: tuple[Fraction, ...],
) -> list[arb]:
    """Return the localized five-knot derivative bounds for ``E q4(U)``."""
    outer = _event_coordinate_moments(
        knot_lower,
        knot_upper,
        Fraction(2, 5),
        Fraction(19, 10),
    )
    middle = _event_coordinate_moments(
        knot_lower,
        knot_upper,
        Fraction(11, 20),
        Fraction(8, 5),
    )
    lower_inner = _event_coordinate_moments(
        knot_lower,
        knot_upper,
        Fraction(7, 10),
        Fraction(1),
    )
    upper_inner = _event_coordinate_moments(
        knot_lower,
        knot_upper,
        Fraction(1),
        Fraction(3, 2),
    )
    base, outer_step, middle_step, lower_step, upper_step = (
        _worker_arb_constants()["envelope"]
    )
    return [
        base
        + outer_step * outer[index]
        + middle_step * middle[index]
        + lower_step * lower_inner[index]
        + upper_step * upper_inner[index]
        for index in range(5)
    ]


def _update_box_digest(
    digest,
    reason: str,
    lower: tuple[Fraction, ...],
    upper: tuple[Fraction, ...],
) -> None:
    digest.update(reason.encode("ascii") + b"|")
    for endpoint in (lower, upper):
        for value in endpoint:
            digest.update(
                f"{value.numerator}/{value.denominator},".encode("ascii")
            )
    digest.update(b"\n")


def _six_boundary_worker(task) -> dict[str, object]:
    """Certify one exact root box of the six-positive boundary."""
    root_index, root_lower, root_upper, root_depth = task
    ctx.prec = ARB_PRECISION_BITS
    _worker_arb_constants()
    weights = (6, 5, 4, 3, 2, 1)
    stack = [(root_lower, root_upper, root_depth)]
    counts = {
        "scalar_core": 0,
        "central_minorant": 0,
        "recursive_global_A": 0,
        "recursive_local_A": 0,
        "infeasible": 0,
    }
    calls = 0
    maximum_depth = 0
    digest = hashlib.sha256()
    q4_cap = _worker_arb_constants()["q4_cap"]

    while stack:
        lower, upper, depth = stack.pop()
        calls += 1
        upper = _tightened_upper(lower, upper, weights)
        if any(upper[index] < lower[index] for index in range(6)):
            counts["infeasible"] += 1
            maximum_depth = max(maximum_depth, depth)
            _update_box_digest(digest, "inf", lower, upper)
            continue
        if sum(upper) <= Fraction(37, 32):
            counts["scalar_core"] += 1
            maximum_depth = max(maximum_depth, depth)
            _update_box_digest(digest, "core", lower, upper)
            continue
        if lower[0] >= Fraction(27, 40):
            counts["central_minorant"] += 1
            maximum_depth = max(maximum_depth, depth)
            _update_box_digest(digest, "central", lower, upper)
            continue

        center = _feasible_center(lower, upper, weights)
        if any(center[index] == 0 for index in range(1, 6)):
            raise AssertionError("unexpected confluent n=6 boundary center")
        radii = [
            max(
                center[index] - lower[index],
                upper[index] - center[index],
            )
            for index in range(6)
        ]
        knots = _cumulative_knots(center)
        knot_lower = _cumulative_knots(lower)
        knot_upper = _cumulative_knots(upper)

        right_knot_bounds = _q4_knot_lipschitz_bounds(
            knot_lower[1:],
            knot_upper[1:],
        )
        right_parameter_bounds = [
            sum(right_knot_bounds[index:], arb(0)) for index in range(5)
        ]
        right_center = _divided_difference(knots[1:], _f4)
        right_error = (
            right_parameter_bounds[0] * _arb_exact(radii[0])
            + right_parameter_bounds[0] * _arb_exact(radii[1])
            + right_parameter_bounds[1] * _arb_exact(radii[2])
            + right_parameter_bounds[2] * _arb_exact(radii[3])
            + right_parameter_bounds[3] * _arb_exact(radii[4])
            + right_parameter_bounds[4] * _arb_exact(radii[5])
        )
        right_lower = right_center - right_error

        reason = None
        if right_lower > 0:
            if (
                _arb_exact(sum(lower)) * right_lower
                > _arb_exact(Fraction(149, 500) * upper[0])
            ):
                reason = "recursive_global_A"
            else:
                left_knot_bounds = _q4_knot_lipschitz_bounds(
                    knot_lower[:5],
                    knot_upper[:5],
                )
                left_parameter_bounds = [
                    sum(left_knot_bounds[index:], arb(0))
                    for index in range(5)
                ]
                left_center = _divided_difference(knots[:5], _f4)
                left_error = sum(
                    (
                        left_parameter_bounds[index]
                        * _arb_exact(radii[index])
                        for index in range(5)
                    ),
                    arb(0),
                )
                left_upper = left_center + left_error
                if not left_upper < q4_cap:
                    left_upper = q4_cap
                if (
                    _arb_exact(sum(lower)) * right_lower
                    > _arb_exact(upper[0]) * left_upper
                ):
                    reason = "recursive_local_A"

        if reason is not None:
            counts[reason] += 1
            maximum_depth = max(maximum_depth, depth)
            digest_reason = {
                "recursive_global_A": "Rg",
                "recursive_local_A": "Rl",
            }[reason]
            _update_box_digest(digest, digest_reason, lower, upper)
            continue
        if depth >= 70:
            raise AssertionError(
                "n=6 six-positive subdivision did not close at "
                f"root {root_index}: {lower!r}, {upper!r}"
            )

        sensitivity = [
            right_parameter_bounds[0],
            right_parameter_bounds[0],
            right_parameter_bounds[1],
            right_parameter_bounds[2],
            right_parameter_bounds[3],
            right_parameter_bounds[4],
        ]
        scores = [
            sensitivity[index] * _arb_exact(upper[index] - lower[index])
            for index in range(6)
        ]
        split_index = 0
        for index in range(1, 6):
            if scores[index] > scores[split_index]:
                split_index = index
        midpoint = (lower[split_index] + upper[split_index]) / 2
        lower_child = list(lower)
        lower_child[split_index] = midpoint
        upper_child = list(upper)
        upper_child[split_index] = midpoint
        stack.append((tuple(lower_child), upper, depth + 1))
        stack.append((lower, tuple(upper_child), depth + 1))

    return {
        "root_index": root_index,
        "total_branch_calls": calls,
        "terminal_box_counts": counts,
        "maximum_bisection_depth": maximum_depth,
        "terminal_transcript_sha256": digest.hexdigest(),
    }


def _six_boundary_roots() -> list[
    tuple[
        int,
        tuple[Fraction, ...],
        tuple[Fraction, ...],
        int,
    ]
]:
    weights = (6, 5, 4, 3, 2, 1)
    boxes = [
        (
            (Fraction(),) * 6,
            (
                Fraction(1),
                Fraction(6, 5),
                Fraction(3, 2),
                Fraction(2),
                Fraction(3),
                Fraction(6),
            ),
            0,
        )
    ]
    while len(boxes) < SIX_BOUNDARY_BASE_ROOT_COUNT:
        lower, upper, depth = boxes.pop(0)
        upper = _tightened_upper(lower, upper, weights)
        split_index = max(
            range(6),
            key=lambda index: weights[index]
            * (upper[index] - lower[index]),
        )
        midpoint = (lower[split_index] + upper[split_index]) / 2
        lower_child = list(lower)
        lower_child[split_index] = midpoint
        upper_child = list(upper)
        upper_child[split_index] = midpoint
        boxes.append((lower, tuple(upper_child), depth + 1))
        boxes.append((tuple(lower_child), upper, depth + 1))

    # The all-zero corner is much more expensive than the other base roots.
    # Refine it twice by the same exact, budget-aware 64-way rule so
    # parallel regeneration does not end with a single long-running worker.
    # This is a fixed proof partition, not runtime-dependent work stealing.
    corner = boxes[0]
    fixed_siblings = []
    for _ in range(SIX_BOUNDARY_CORNER_REFINEMENT_LEVELS):
        refined = [corner]
        while len(refined) < SIX_BOUNDARY_CORNER_REFINEMENT:
            lower, upper, depth = refined.pop(0)
            upper = _tightened_upper(lower, upper, weights)
            split_index = max(
                range(6),
                key=lambda index: weights[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            refined.append((lower, tuple(upper_child), depth + 1))
            refined.append((tuple(lower_child), upper, depth + 1))
        corner = refined[0]
        fixed_siblings.extend(refined[1:])
    boxes = [corner] + fixed_siblings + boxes[1:]
    assert len(boxes) == SIX_BOUNDARY_ROOT_COUNT
    return [
        (index, lower, upper, depth)
        for index, (lower, upper, depth) in enumerate(boxes)
    ]


def _root_partition_digest(tasks) -> str:
    digest = hashlib.sha256()
    for index, lower, upper, depth in tasks:
        digest.update(f"{index}|{depth}|".encode("ascii"))
        for endpoint in (lower, upper):
            for value in endpoint:
                digest.update(
                    f"{value.numerator}/{value.denominator},".encode("ascii")
                )
        digest.update(b"\n")
    return digest.hexdigest()


def _verify_six_positive_boundary() -> dict[str, object]:
    tasks = _six_boundary_roots()
    worker_text = os.environ.get("STRINGER_N6_WORKERS")
    if worker_text is None:
        worker_count = min(9, os.cpu_count() or 1)
    else:
        worker_count = int(worker_text)
        if worker_count < 1:
            raise ValueError("STRINGER_N6_WORKERS must be positive")

    if worker_count == 1:
        results = [_six_boundary_worker(task) for task in tasks]
    else:
        spawn_context = mp.get_context("spawn")
        with spawn_context.Pool(processes=worker_count) as pool:
            results = list(
                pool.imap_unordered(
                    _six_boundary_worker,
                    tasks,
                    chunksize=1,
                )
            )
    results.sort(key=lambda result: result["root_index"])

    counts = {
        key: sum(
            result["terminal_box_counts"][key] for result in results
        )
        for key in results[0]["terminal_box_counts"]
    }
    calls = sum(result["total_branch_calls"] for result in results)
    maximum_depth = max(
        result["maximum_bisection_depth"] for result in results
    )
    transcript_digest = hashlib.sha256()
    for result in results:
        transcript_digest.update(
            (
                f"{result['root_index']}|"
                f"{result['terminal_transcript_sha256']}\n"
            ).encode("ascii")
        )

    # Filled from the first complete fixed-root regeneration.  Keeping these
    # values explicit turns accidental changes to any branch decision into a
    # failing proof-regression test.
    expected_calls = 46317164
    expected_counts = {
        "scalar_core": 809393,
        "central_minorant": 348,
        "recursive_global_A": 10442544,
        "recursive_local_A": 11906872,
        "infeasible": 0,
    }
    expected_depth = 45
    expected_root_digest = (
        "2c8df24aedb29c4b796687d8f79ebeed0bc85b6c6784682da480deef727fbf36"
    )
    expected_transcript_digest = (
        "0c6ff61e6a587c8e00975053ad00848c008ee530c6bfcab961175397ee801431"
    )
    assert calls == expected_calls
    assert counts == expected_counts
    assert maximum_depth == expected_depth
    assert _root_partition_digest(tasks) == expected_root_digest
    assert transcript_digest.hexdigest() == expected_transcript_digest
    assert calls == 2 * sum(counts.values()) - len(tasks)

    return {
        "parameterization": (
            "(a,b,c,d,e,f)=(a,a+r,a+r+s,a+r+s+t,"
            "a+r+s+t+v,a+r+s+t+v+w), "
            "6a+5r+4s+3t+2v+w<=6"
        ),
        "cumulative_coordinate_weights": [6, 5, 4, 3, 2, 1],
        "recursive_pruning": (
            "knot insertion from the proved five-positive face, "
            "A<=149/500, and a localized step envelope for abs(q4')"
        ),
        "localized_q4_derivative_envelope": (
            "1/20 off [2/5,19/10], 1/4 on [2/5,19/10] "
            "off [11/20,8/5], 11/20 on [11/20,8/5] off the "
            "two inner intervals, 81/100 on [7/10,1], and "
            "24/25 on [1,3/2]"
        ),
        "initial_root_boxes": len(tasks),
        "root_partition_sha256": _root_partition_digest(tasks),
        "total_branch_calls": calls,
        "terminal_box_counts": counts,
        "maximum_bisection_depth": maximum_depth,
        "terminal_transcript_sha256": transcript_digest.hexdigest(),
    }


def verify_n6_seven_knot_theorem() -> dict[str, object]:
    """Certify the complete constrained comparison at ``n=6``.

    The radial reduction leaves at most six positive knots.  Three nested
    certificates prove the four-, five-, and six-positive faces in that
    order.  Every compact numerical sign is a directed 160-bit Arb
    inclusion; all subdivision endpoints and analytic constants are exact
    rationals.
    """
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
    exponential_piece = exp(-6 / u)

    f3_lower = u**3 * exponential_piece
    f3_upper = f3_lower - (u - 1) ** 6 / u**3
    f4_lower = u**4 * exponential_piece
    f4_upper = f4_lower - (u - 1) ** 6 / u**2
    f5_lower = u**5 * exponential_piece
    f5_upper = f5_lower - (u - 1) ** 6 / u

    q3_lower_u = exponential_piece * sum(
        (6 / u) ** order / factorial(order) for order in range(4)
    )
    q3_upper_u = (
        q3_lower_u - 1 + 15 / u**4 - 24 / u**5 + 10 / u**6
    )
    q4_lower_u = exponential_piece * sum(
        (6 / u) ** order / factorial(order) for order in range(5)
    )
    q4_upper_u = q4_lower_u - 1 + 6 / u**5 - 5 / u**6
    q5_lower_u = exponential_piece * sum(
        (6 / u) ** order / factorial(order) for order in range(6)
    )
    q5_upper_u = q5_lower_u - 1 + u**-6
    assert simplify(diff(f3_lower, u, 3) / 6 - q3_lower_u) == 0
    assert simplify(diff(f3_upper, u, 3) / 6 - q3_upper_u) == 0
    assert simplify(diff(f4_lower, u, 4) / 24 - q4_lower_u) == 0
    assert simplify(diff(f4_upper, u, 4) / 24 - q4_upper_u) == 0
    assert simplify(diff(f5_lower, u, 5) / 120 - q5_lower_u) == 0
    assert simplify(diff(f5_upper, u, 5) / 120 - q5_upper_u) == 0

    q3_lower_derivative = exp(-lam) * lam**5 / 36
    q3_upper_derivative = (
        q3_lower_derivative
        - 5 * lam**5 * (6 - lam) ** 2 / 23328
    )
    q4_lower_derivative = exp(-lam) * lam**6 / 144
    q4_upper_derivative = (
        q4_lower_derivative + 5 * lam**6 * (lam - 6) / 46656
    )
    q5_lower_derivative = exp(-lam) * lam**7 / 720
    # The correction denominator is 6^6 because d(u^-6)/du=-6u^-7.
    q5_upper_derivative = q5_lower_derivative - lam**7 / 46656
    for lower_u, upper_u, lower_lam, upper_lam in (
        (
            q3_lower_u,
            q3_upper_u,
            q3_lower_derivative,
            q3_upper_derivative,
        ),
        (
            q4_lower_u,
            q4_upper_u,
            q4_lower_derivative,
            q4_upper_derivative,
        ),
        (
            q5_lower_u,
            q5_upper_u,
            q5_lower_derivative,
            q5_upper_derivative,
        ),
    ):
        assert simplify(
            diff(lower_u, u) - lower_lam.subs(lam, 6 / u)
        ) == 0
        assert simplify(
            diff(upper_u, u) - upper_lam.subs(lam, 6 / u)
        ) == 0
    assert simplify(diff(q3_lower_derivative, lam)) == (
        exp(-lam) * lam**4 * (5 - lam) / 36
    )
    assert simplify(diff(q4_lower_derivative, lam)) == (
        exp(-lam) * lam**5 * (6 - lam) / 144
    )

    # Exact Taylor comparisons used on the analytic lower pieces of q3' and
    # q4'.  The four q4 values prove the successive step-envelope bounds.
    exp_six_degree_11 = _exp_taylor_lower(Fraction(6), 11)
    exp_six_degree_13 = _exp_taylor_lower(Fraction(6), 13)
    exp_sixty_sevenths_degree_14 = _exp_taylor_lower(
        Fraction(60, 7),
        14,
    )
    exp_one_twenty_elevenths_degree_14 = _exp_taylor_lower(
        Fraction(120, 11),
        14,
    )
    exp_fifteen_degree_15 = _exp_taylor_lower(Fraction(15), 15)
    assert exp_six_degree_11 == Fraction(760997, 1925)
    assert exp_six_degree_13 == Fraction(10059173, 25025)
    assert exp_sixty_sevenths_degree_14 == Fraction(
        24353880812705672483,
        4752309071452943,
    )
    assert exp_one_twenty_elevenths_degree_14 == Fraction(
        125192034360679668826127,
        2660907083917769687,
    )
    assert exp_fifteen_degree_15 == Fraction(
        13324966405463,
        7175168,
    )
    assert exp_six_degree_11 > Fraction(4320, 11)
    assert exp_six_degree_11 > Fraction(57500, 149)
    assert exp_six_degree_13 > 400
    assert exp_sixty_sevenths_degree_14 > (
        Fraction(20, 11) * Fraction(60, 7) ** 6 / 144
    )
    assert exp_one_twenty_elevenths_degree_14 > (
        4 * Fraction(120, 11) ** 6 / 144
    )
    assert exp_fifteen_degree_15 > (
        20 * Fraction(15) ** 6 / 144
    )

    # Symbolically rederive every individual-coordinate Beta(h,5-h) moment
    # polynomial used by the localized final-face radius.
    s, radius = symbols("s radius", real=True)
    for high_count, (low_terms, high_terms) in (
        _TAIL_COEFFICIENTS_5.items()
    ):
        low_count = 5 - high_count
        density = (
            factorial(4)
            / (factorial(high_count - 1) * factorial(low_count - 1))
            * s ** (high_count - 1)
            * (1 - s) ** (low_count - 1)
        )
        symbolic_low = integrate(
            (1 - s) * density / low_count,
            (s, 1 - radius, 1),
        )
        symbolic_high = integrate(
            s * density / high_count,
            (s, 1 - radius, 1),
        )
        finite_low = sum(
            Rational(coefficient.numerator, coefficient.denominator)
            * radius**power
            for coefficient, power in low_terms
        )
        finite_high = sum(
            Rational(coefficient.numerator, coefficient.denominator)
            * radius**power
            for coefficient, power in high_terms
        )
        assert simplify(symbolic_low - finite_low) == 0
        assert simplify(symbolic_high - finite_high) == 0

    # Knot insertion is the recursive bridge between consecutive faces.
    nodes = symbols("x0:6", positive=True)
    values = symbols("z0:6", real=True)

    def symbolic_dd(input_nodes, input_values):
        current = list(input_values)
        for order in range(1, len(input_nodes)):
            current = [
                (current[index + 1] - current[index])
                / (input_nodes[index + order] - input_nodes[index])
                for index in range(len(current) - 1)
            ]
        return current[0]

    for length in (5, 6):
        selected_nodes = nodes[:length]
        selected_values = values[:length]
        target = symbolic_dd(
            selected_nodes,
            [
                node * value
                for node, value in zip(selected_nodes, selected_values)
            ],
        )
        left = symbolic_dd(selected_nodes[:-1], selected_values[:-1])
        right = symbolic_dd(selected_nodes[1:], selected_values[1:])
        assert simplify(
            (selected_nodes[-1] - selected_nodes[0]) * target
            - (
                selected_nodes[-1] * right
                - selected_nodes[0] * left
            )
        ) == 0

    previous_precision = ctx.prec
    ctx.prec = ARB_PRECISION_BITS

    def arb_interval(lower: Fraction, upper: Fraction) -> arb:
        lower_q = fmpq(lower.numerator, lower.denominator)
        upper_q = fmpq(upper.numerator, upper.denominator)
        return arb((lower_q + upper_q) / 2, (upper_q - lower_q) / 2)

    def update_interval_digest(digest, label, lower, upper):
        digest.update(
            (
                f"{label}|{lower.numerator}/{lower.denominator}|"
                f"{upper.numerator}/{upper.denominator}\n"
            ).encode("ascii")
        )

    def arb_q(value: arb, order: int, upper_piece: bool) -> arb:
        lambda_value = 6 / value
        result = (-lambda_value).exp() * sum(
            lambda_value**index / math.factorial(index)
            for index in range(order + 1)
        )
        if upper_piece:
            if order == 3:
                result += (
                    -1
                    + 15 / value**4
                    - 24 / value**5
                    + 10 / value**6
                )
            elif order == 4:
                result += -1 + 6 / value**5 - 5 / value**6
            elif order == 5:
                result += -1 + value**-6
            else:
                raise ValueError("unsupported scalar order")
        return result

    def arb_q_derivative(lambda_value: arb, order: int) -> arb:
        if order == 3:
            return (
                (-lambda_value).exp() * lambda_value**5 / 36
                - 5
                * lambda_value**5
                * (6 - lambda_value) ** 2
                / 23328
            )
        if order == 4:
            return (
                (-lambda_value).exp() * lambda_value**6 / 144
                + 5 * lambda_value**6 * (lambda_value - 6) / 46656
            )
        raise ValueError("unsupported derivative order")

    scalar_digest = hashlib.sha256()
    scalar_records: dict[str, dict[str, object]] = {}

    def certify_scalar(
        label: str,
        start: Fraction,
        stop: Fraction,
        evaluator,
        maximum_depth: int = 40,
    ) -> None:
        stack = [(Fraction(start), Fraction(stop), 0)]
        leaves = 0
        depth_seen = 0
        while stack:
            lower, upper, depth = stack.pop()
            if evaluator(arb_interval(lower, upper)) > 0:
                leaves += 1
                depth_seen = max(depth_seen, depth)
                update_interval_digest(
                    scalar_digest,
                    label,
                    lower,
                    upper,
                )
                continue
            if depth >= maximum_depth:
                raise AssertionError(
                    f"n=6 scalar certificate did not close: {label}"
                )
            midpoint = (lower + upper) / 2
            stack.append((lower, midpoint, depth + 1))
            stack.append((midpoint, upper, depth + 1))
        scalar_records[label] = {
            "interval": f"[{start},{stop}]",
            "certified_leaf_intervals": leaves,
            "maximum_bisection_depth": depth_seen,
        }

    def certify_abs_upper_derivative(
        label: str,
        start: Fraction,
        stop: Fraction,
        bound: Fraction,
    ) -> None:
        certify_scalar(
            f"{label}-upper",
            start,
            stop,
            lambda value: _arb_exact(bound) - arb_q_derivative(value, 4),
        )
        certify_scalar(
            f"{label}-lower",
            start,
            stop,
            lambda value: _arb_exact(bound) + arb_q_derivative(value, 4),
        )

    branch_records: dict[str, dict[str, object]] = {}
    try:
        # Scalar estimates for the four-positive face.
        certify_scalar(
            "q3-core",
            Fraction(1),
            Fraction(31, 16),
            lambda value: arb_q(value, 3, True),
        )
        certify_scalar(
            "q3-central-lower",
            Fraction(3, 4),
            Fraction(1),
            lambda value: arb_q(value, 3, False)
            - _arb_exact(Fraction(1, 20))
            * (_arb_exact(Fraction(3, 2)) - value),
        )
        certify_scalar(
            "q3-central-upper",
            Fraction(1),
            Fraction(15, 4),
            lambda value: arb_q(value, 3, True)
            - _arb_exact(Fraction(1, 20))
            * (_arb_exact(Fraction(3, 2)) - value),
        )
        certify_scalar(
            "q3-derivative-upper",
            Fraction(1),
            Fraction(6),
            lambda value: _arb_exact(Fraction(11, 20))
            - arb_q_derivative(value, 3),
        )
        certify_scalar(
            "q3-derivative-lower",
            Fraction(1),
            Fraction(6),
            lambda value: _arb_exact(Fraction(11, 20))
            + arb_q_derivative(value, 3),
        )
        certify_scalar(
            "q3-upper-cap",
            Fraction(1),
            Fraction(6),
            lambda value: _arb_exact(Fraction(41, 200))
            - arb_q(value, 3, True),
        )

        # Scalar estimates for the five-positive face and final recursion.
        certify_scalar(
            "q4-core",
            Fraction(1),
            Fraction(23, 16),
            lambda value: arb_q(value, 4, True),
        )
        certify_scalar(
            "q4-central-lower",
            Fraction(7, 10),
            Fraction(1),
            lambda value: arb_q(value, 4, False)
            - _arb_exact(Fraction(13, 100))
            * (_arb_exact(Fraction(6, 5)) - value),
        )
        certify_scalar(
            "q4-central-upper",
            Fraction(1),
            Fraction(16, 5),
            lambda value: arb_q(value, 4, True)
            - _arb_exact(Fraction(13, 100))
            * (_arb_exact(Fraction(6, 5)) - value),
        )
        certify_scalar(
            "q4-derivative-upper",
            Fraction(1),
            Fraction(6),
            lambda value: _arb_exact(Fraction(24, 25))
            - arb_q_derivative(value, 4),
        )
        certify_scalar(
            "q4-derivative-lower",
            Fraction(1),
            Fraction(6),
            lambda value: _arb_exact(Fraction(24, 25))
            + arb_q_derivative(value, 4),
        )
        certify_scalar(
            "q4-upper-cap",
            Fraction(1),
            Fraction(6),
            lambda value: _arb_exact(Fraction(149, 500))
            - arb_q(value, 4, True),
        )
        certify_abs_upper_derivative(
            "q4-local-outer",
            Fraction(1),
            Fraction(60, 19),
            Fraction(1, 20),
        )
        certify_abs_upper_derivative(
            "q4-local-middle",
            Fraction(60, 19),
            Fraction(15, 4),
            Fraction(1, 4),
        )
        certify_abs_upper_derivative(
            "q4-local-shoulder",
            Fraction(15, 4),
            Fraction(4),
            Fraction(11, 20),
        )

        # Scalar estimates for the final six-positive boundary.
        certify_scalar(
            "q5-core",
            Fraction(1),
            Fraction(37, 32),
            lambda value: arb_q(value, 5, True),
        )
        certify_scalar(
            "q5-central-lower",
            Fraction(27, 40),
            Fraction(1),
            lambda value: arb_q(value, 5, False)
            - _arb_exact(Fraction(3, 8)) * (1 - value),
        )
        certify_scalar(
            "q5-central-upper",
            Fraction(1),
            Fraction(21, 8),
            lambda value: arb_q(value, 5, True)
            - _arb_exact(Fraction(3, 8)) * (1 - value),
        )

        expected_scalar = {
            "q3-core": (444, 18),
            "q3-central-lower": (5, 4),
            "q3-central-upper": (529, 12),
            "q3-derivative-upper": (26, 8),
            "q3-derivative-lower": (25, 6),
            "q3-upper-cap": (4834, 18),
            "q4-core": (64, 9),
            "q4-central-lower": (7, 5),
            "q4-central-upper": (1151, 15),
            "q4-derivative-upper": (12, 6),
            "q4-derivative-lower": (125, 11),
            "q4-upper-cap": (1467, 18),
            "q4-local-outer-upper": (102, 9),
            "q4-local-outer-lower": (13, 6),
            "q4-local-middle-upper": (3, 2),
            "q4-local-middle-lower": (9, 5),
            "q4-local-shoulder-upper": (1, 0),
            "q4-local-shoulder-lower": (2, 1),
            "q5-core": (17, 9),
            "q5-central-lower": (12, 9),
            "q5-central-upper": (164, 11),
        }
        assert {
            key: (
                value["certified_leaf_intervals"],
                value["maximum_bisection_depth"],
            )
            for key, value in scalar_records.items()
        } == expected_scalar
        assert scalar_digest.hexdigest() == (
            "45911d92da792ce6dd0ce418ed617c0530f033ac98f80a4a5d0ece0c8d909db1"
        )

        # Four-positive face: [a,b,c,d]f3.
        weights4 = (4, 3, 2, 1)
        stack4 = [
            (
                (Fraction(),) * 4,
                (
                    Fraction(3, 2),
                    Fraction(2),
                    Fraction(3),
                    Fraction(6),
                ),
                0,
            )
        ]
        counts4 = {
            "scalar_core": 0,
            "far_cap": 0,
            "central_minorant": 0,
            "direct_lipschitz": 0,
            "infeasible": 0,
        }
        calls4 = 0
        depth4 = 0
        digest4 = hashlib.sha256()
        while stack4:
            lower, upper, depth = stack4.pop()
            calls4 += 1
            upper = _tightened_upper(lower, upper, weights4)
            if any(upper[index] < lower[index] for index in range(4)):
                counts4["infeasible"] += 1
                depth4 = max(depth4, depth)
                _update_box_digest(digest4, "infeasible", lower, upper)
                continue
            if sum(upper) <= Fraction(31, 16):
                counts4["scalar_core"] += 1
                depth4 = max(depth4, depth)
                _update_box_digest(digest4, "core", lower, upper)
                continue
            if sum(lower) >= 5:
                counts4["far_cap"] += 1
                depth4 = max(depth4, depth)
                _update_box_digest(digest4, "far", lower, upper)
                continue
            if lower[0] >= Fraction(3, 4):
                counts4["central_minorant"] += 1
                depth4 = max(depth4, depth)
                _update_box_digest(digest4, "central", lower, upper)
                continue

            center = _feasible_center(lower, upper, weights4)
            if any(center[index] == 0 for index in range(1, 4)):
                raise AssertionError(
                    "unexpected confluent n=6 four-positive center"
                )
            radii = [
                max(
                    center[index] - lower[index],
                    upper[index] - center[index],
                )
                for index in range(4)
            ]
            center_value = _divided_difference(
                _cumulative_knots(center),
                _f3,
            )
            error = _arb_exact(Fraction(11, 80)) * sum(
                _arb_exact(weights4[index] * radii[index])
                for index in range(4)
            )
            if center_value > error:
                counts4["direct_lipschitz"] += 1
                depth4 = max(depth4, depth)
                _update_box_digest(digest4, "direct", lower, upper)
                continue
            if depth >= 45:
                raise AssertionError(
                    "n=6 four-positive subdivision did not close"
                )
            split_index = max(
                range(4),
                key=lambda index: weights4[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            stack4.append((tuple(lower_child), upper, depth + 1))
            stack4.append((lower, tuple(upper_child), depth + 1))

        assert calls4 == 36203
        assert counts4 == {
            "scalar_core": 334,
            "far_cap": 152,
            "central_minorant": 1,
            "direct_lipschitz": 17615,
            "infeasible": 0,
        }
        assert depth4 == 24
        assert calls4 == 2 * sum(counts4.values()) - 1
        assert digest4.hexdigest() == (
            "796d71cff3aebbbaf38f6b63814f8a25d03eaa2a0baf701d9b62e23689ec888b"
        )
        branch_records["four_positive_face"] = {
            "parameterization": (
                "(a,b,c,d)=(a,a+r,a+r+s,a+r+s+t), "
                "4a+3r+2s+t<=6"
            ),
            "cumulative_coordinate_weights": list(weights4),
            "total_branch_calls": calls4,
            "terminal_box_counts": counts4,
            "maximum_bisection_depth": depth4,
            "terminal_transcript_sha256": digest4.hexdigest(),
        }

        # Five-positive face: [a,b,c,d,e]f4, recursively from f3.
        weights5 = (5, 4, 3, 2, 1)
        stack5 = [
            (
                (Fraction(),) * 5,
                (
                    Fraction(6, 5),
                    Fraction(3, 2),
                    Fraction(2),
                    Fraction(3),
                    Fraction(6),
                ),
                0,
            )
        ]
        counts5 = {
            "scalar_core": 0,
            "central_minorant": 0,
            "recursive_global_A": 0,
            "recursive_local_A": 0,
            "direct_lipschitz": 0,
            "infeasible": 0,
        }
        calls5 = 0
        depth5 = 0
        digest5 = hashlib.sha256()
        q3_lipschitz = _arb_exact(Fraction(11, 20))
        q3_cap = _arb_exact(Fraction(41, 200))
        while stack5:
            lower, upper, depth = stack5.pop()
            calls5 += 1
            upper = _tightened_upper(lower, upper, weights5)
            if any(upper[index] < lower[index] for index in range(5)):
                counts5["infeasible"] += 1
                depth5 = max(depth5, depth)
                _update_box_digest(digest5, "infeasible", lower, upper)
                continue
            if sum(upper) <= Fraction(23, 16):
                counts5["scalar_core"] += 1
                depth5 = max(depth5, depth)
                _update_box_digest(digest5, "core", lower, upper)
                continue
            if lower[0] >= Fraction(7, 10):
                counts5["central_minorant"] += 1
                depth5 = max(depth5, depth)
                _update_box_digest(digest5, "central", lower, upper)
                continue

            center = _feasible_center(lower, upper, weights5)
            if any(center[index] == 0 for index in range(1, 5)):
                raise AssertionError(
                    "unexpected confluent n=6 five-positive center"
                )
            radii = [
                max(
                    center[index] - lower[index],
                    upper[index] - center[index],
                )
                for index in range(5)
            ]
            knots = _cumulative_knots(center)
            right_center = _divided_difference(knots[1:], _f3)
            right_error = q3_lipschitz * (
                _arb_exact(radii[0])
                + _arb_exact(radii[1])
                + _arb_exact(Fraction(3, 4) * radii[2])
                + _arb_exact(Fraction(1, 2) * radii[3])
                + _arb_exact(Fraction(1, 4) * radii[4])
            )
            right_lower = right_center - right_error
            reason = None
            if right_lower > 0:
                if (
                    _arb_exact(sum(lower)) * right_lower
                    > _arb_exact(upper[0]) * q3_cap
                ):
                    reason = "recursive_global_A"
                else:
                    left_center = _divided_difference(knots[:4], _f3)
                    left_error = q3_lipschitz * (
                        _arb_exact(radii[0])
                        + _arb_exact(Fraction(3, 4) * radii[1])
                        + _arb_exact(Fraction(1, 2) * radii[2])
                        + _arb_exact(Fraction(1, 4) * radii[3])
                    )
                    left_upper = left_center + left_error
                    if not left_upper < q3_cap:
                        left_upper = q3_cap
                    if (
                        _arb_exact(sum(lower)) * right_lower
                        > _arb_exact(upper[0]) * left_upper
                    ):
                        reason = "recursive_local_A"
            if reason is not None:
                counts5[reason] += 1
                depth5 = max(depth5, depth)
                _update_box_digest(digest5, reason, lower, upper)
                continue

            center_value = _divided_difference(knots, _f4)
            error = _arb_exact(Fraction(24, 125)) * sum(
                _arb_exact(weights5[index] * radii[index])
                for index in range(5)
            )
            if center_value > error:
                counts5["direct_lipschitz"] += 1
                depth5 = max(depth5, depth)
                _update_box_digest(digest5, "direct", lower, upper)
                continue
            if depth >= 50:
                raise AssertionError(
                    "n=6 five-positive subdivision did not close"
                )
            split_index = max(
                range(5),
                key=lambda index: weights5[index]
                * (upper[index] - lower[index]),
            )
            midpoint = (lower[split_index] + upper[split_index]) / 2
            lower_child = list(lower)
            lower_child[split_index] = midpoint
            upper_child = list(upper)
            upper_child[split_index] = midpoint
            stack5.append((tuple(lower_child), upper, depth + 1))
            stack5.append((lower, tuple(upper_child), depth + 1))

        assert calls5 == 725061
        assert counts5 == {
            "scalar_core": 21295,
            "central_minorant": 164,
            "recursive_global_A": 189416,
            "recursive_local_A": 144720,
            "direct_lipschitz": 6936,
            "infeasible": 0,
        }
        assert depth5 == 37
        assert calls5 == 2 * sum(counts5.values()) - 1
        assert digest5.hexdigest() == (
            "d3442eed3224f9c4dc1fb65624c1d375aa5d223670830948f7a13009ac9039b8"
        )
        branch_records["five_positive_face"] = {
            "parameterization": (
                "(a,b,c,d,e)=(a,a+r,a+r+s,a+r+s+t,a+r+s+t+v), "
                "5a+4r+3s+2t+v<=6"
            ),
            "cumulative_coordinate_weights": list(weights5),
            "recursive_pruning": (
                "knot insertion from the proved four-positive face, "
                "A<=41/200, and abs(q3')<11/20"
            ),
            "total_branch_calls": calls5,
            "terminal_box_counts": counts5,
            "maximum_bisection_depth": depth5,
            "terminal_transcript_sha256": digest5.hexdigest(),
        }

        branch_records["six_positive_boundary"] = (
            _verify_six_positive_boundary()
        )
    finally:
        ctx.prec = previous_precision

    return {
        "arb_precision_bits": ARB_PRECISION_BITS,
        "exact_exponential_lower_bounds": {
            "exp_6_degree_11": _fraction_record(exp_six_degree_11),
            "exp_6_degree_13": _fraction_record(exp_six_degree_13),
            "exp_60_over_7_degree_14": _fraction_record(
                exp_sixty_sevenths_degree_14
            ),
            "exp_120_over_11_degree_14": _fraction_record(
                exp_one_twenty_elevenths_degree_14
            ),
            "exp_15_degree_15": _fraction_record(exp_fifteen_degree_15),
        },
        "scalar_interval_records": scalar_records,
        "scalar_terminal_transcript_sha256": scalar_digest.hexdigest(),
        "derivative_bounds": {
            "q3_global": "abs(q3'(u))<11/20",
            "q4_global": "abs(q4'(u))<24/25",
            "q4_localized": (
                "the nested 1/20, 1/4, 11/20, 81/100, and 24/25 "
                "step envelope recorded in the six-positive proof"
            ),
        },
        "pointwise_upper_bounds": {
            "q3": "q3(u)<41/200",
            "q4": "q4(u)<149/500",
        },
        "beta_tail_moment_formulas": (
            "symbolically checked for Beta(h,5-h), h=1,...,4"
        ),
        **branch_records,
    }


if __name__ == "__main__":
    import json

    print(
        json.dumps(
            verify_n6_seven_knot_theorem(),
            indent=2,
            sort_keys=True,
        )
    )
