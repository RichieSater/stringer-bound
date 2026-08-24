"""Directed certificate for the complete ``n=7`` Poissonization comparison.

The proof closes the four-, five-, six-, and seven-positive boundary faces
in order.  All subdivision endpoints and analytic constants are exact
rationals; every compact transcendental sign is a directed 160-bit Arb
inclusion.  The expensive seven-positive face has a fixed root partition so
that independent shards can be regenerated without schedule-dependent proof
output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import pathlib
import sys
import time
from fractions import Fraction

from flint import arb, ctx
from n6_dirichlet_poissonization import (
    _arb_exact,
    _cumulative_knots,
    _divided_difference,
    _tail_coefficients,
    _update_box_digest,
)
from n7_dirichlet_poissonization_scalars import (
    ARB_PRECISION_BITS,
    FACE6_Q4_EVENTS,
    FACE7_Q5_EVENTS,
    verify_n7_scalar_bounds,
)

TOTAL = 7
FACE7_BASE_ROOT_COUNT = 1024
FACE7_PRIMARY_REFINEMENT = 4096
FACE7_PRIMARY_REFINED_BASE_INDICES = (*range(17), 32)
FACE7_EXTRA_REFINEMENT = 4096
FACE7_EXTRA_REFINEMENT_CALL_CUTOFF = 5_000_000
# Populated only after the complete primary partition has been regenerated.
# An entry is a root identifier in that primary partition whose children
# replace the parent in the final, schedule-independent proof partition.
FACE7_EXTRA_REFINED_PARENT_IDS: tuple[int, ...] = (
    0,
    1,
    2,
    3,
    4,
    5,
    8,
    9,
    16,
    17,
    24,
    32,
    33,
    40,
    48,
    64,
    65,
    80,
    128,
    129,
    132,
    136,
    144,
    160,
    192,
    526,
    770,
    772,
    16384,
    32768,
    32769,
    32770,
    32771,
    32772,
    32773,
    32774,
    32776,
    32784,
    32800,
    32832,
    32896,
    33024,
    73774,
    73838,
    74222,
)
FACE7_SHARD_COUNT = 256
FACE7_SHARD_SALT = "n7-dirichlet-poissonization-v1"
FACE7_WORKER_CALL_LIMIT = int(os.environ.get("N7_FACE7_CALL_LIMIT", "0"))

SCALAR_EXPECTATION = {
    "compact_terminal_intervals": 140197,
    "compact_maximum_bisection_depth": 25,
    "compact_transcript_sha256": (
        "87ca0938a8faa34e014dc792c0627b750a54179ddba0bd00e8ff7ad6af7c16f4"
    ),
}


# The corresponding helpers in the n=6 verifier deliberately hard-code a
# budget of six.  Reusing them here would omit the shell 6 < sum(y_i) <= 7.
def _tightened_upper(lower, upper, weights):
    result = list(upper)
    for index, weight in enumerate(weights):
        available = Fraction(TOTAL) - sum(
            weights[other] * lower[other]
            for other in range(len(weights))
            if other != index
        )
        result[index] = min(result[index], available / weight)
    return tuple(result)


def _feasible_center(lower, upper, weights):
    center = tuple((lower[index] + upper[index]) / 2 for index in range(len(weights)))
    weighted_center = sum(
        weights[index] * center[index] for index in range(len(weights))
    )
    if weighted_center > TOTAL:
        weighted_lower = sum(
            weights[index] * lower[index] for index in range(len(weights))
        )
        scale = (Fraction(TOTAL) - weighted_lower) / (weighted_center - weighted_lower)
        center = tuple(
            lower[index] + scale * (center[index] - lower[index])
            for index in range(len(weights))
        )
    return center


def _f(k, v):
    """Evaluate ``u^k(exp(-7/u) - (1-1/u)_+^7)`` as an Arb ball."""
    if v == 0:
        return arb(0)
    x = _arb_exact(v)
    out = x**k * (-7 / x).exp()
    if v > 1:
        out -= (x - 1) ** 7 / x ** (7 - k)
    return out


def _base_roots(dim, count=1024):
    """Build a deterministic longest-weighted-edge box partition."""
    weights = tuple(range(dim, 0, -1))
    boxes = [((Fraction(),) * dim, tuple(Fraction(7, j) for j in range(dim, 0, -1)), 0)]
    while len(boxes) < count:
        lower, upper, depth = boxes.pop(0)
        upper = _tightened_upper(lower, upper, weights)
        index = max(range(dim), key=lambda i: weights[i] * (upper[i] - lower[i]))
        midpoint = (lower[index] + upper[index]) / 2
        lower_child = list(lower)
        lower_child[index] = midpoint
        upper_child = list(upper)
        upper_child[index] = midpoint
        boxes.extend(
            (
                (lower, tuple(upper_child), depth + 1),
                (tuple(lower_child), upper, depth + 1),
            )
        )
    return boxes


def _refine(box, count, dim):
    """Apply the same exact partition rule inside one existing root box."""
    weights = tuple(range(dim, 0, -1))
    boxes = [box]
    while len(boxes) < count:
        lower, upper, depth = boxes.pop(0)
        upper = _tightened_upper(lower, upper, weights)
        index = max(range(dim), key=lambda i: weights[i] * (upper[i] - lower[i]))
        midpoint = (lower[index] + upper[index]) / 2
        lower_child = list(lower)
        lower_child[index] = midpoint
        upper_child = list(upper)
        upper_child[index] = midpoint
        boxes.extend(
            (
                (lower, tuple(upper_child), depth + 1),
                (tuple(lower_child), upper, depth + 1),
            )
        )
    return boxes


def face4_worker(task):
    """Certify one root of the four-positive coordinate face."""
    root_index, lower0, upper0, depth0 = task
    ctx.prec = ARB_PRECISION_BITS
    weights = (4, 3, 2, 1)
    stack = [(lower0, upper0, depth0)]
    counts = {x: 0 for x in ["core", "far", "central", "direct", "infeasible"]}
    calls = maximum_depth = 0
    digest = hashlib.sha256()
    while stack:
        lower, upper, depth = stack.pop()
        calls += 1
        upper = _tightened_upper(lower, upper, weights)
        reason = None
        if any(upper[i] < lower[i] for i in range(4)):
            reason = "infeasible"
        elif sum(upper) <= Fraction(9, 4):
            reason = "core"
        elif sum(lower) >= 6:
            reason = "far"
        elif lower[0] >= Fraction(7, 8):
            reason = "central"
        else:
            center = _feasible_center(lower, upper, weights)
            radii = [max(center[i] - lower[i], upper[i] - center[i]) for i in range(4)]
            value = _divided_difference(_cumulative_knots(center), lambda v: _f(3, v))
            error = _arb_exact(Fraction(39, 400)) * sum(
                (_arb_exact(weights[i] * radii[i]) for i in range(4)), arb(0)
            )
            if value > error:
                reason = "direct"
        if reason:
            counts[reason] += 1
            maximum_depth = max(maximum_depth, depth)
            _update_box_digest(digest, reason, lower, upper)
            continue
        if depth >= 50:
            raise RuntimeError(
                ("face4 fail", root_index, lower, upper, depth, value, error)
            )
        index = max(range(4), key=lambda i: weights[i] * (upper[i] - lower[i]))
        midpoint = (lower[index] + upper[index]) / 2
        lower_child = list(lower)
        lower_child[index] = midpoint
        upper_child = list(upper)
        upper_child[index] = midpoint
        stack.extend(
            (
                (tuple(lower_child), upper, depth + 1),
                (lower, tuple(upper_child), depth + 1),
            )
        )
    return {
        "i": root_index,
        "calls": calls,
        "counts": counts,
        "depth": maximum_depth,
        "digest": digest.hexdigest(),
    }


def face5_worker(task):
    """Certify one root of the five-positive coordinate face."""
    root_index, lower0, upper0, depth0 = task
    ctx.prec = ARB_PRECISION_BITS
    weights = (5, 4, 3, 2, 1)
    stack = [(lower0, upper0, depth0)]
    counts = {
        x: 0
        for x in ["core", "central", "rec_global", "rec_local", "direct", "infeasible"]
    }
    calls = maximum_depth = 0
    digest = hashlib.sha256()
    qlip = _arb_exact(Fraction(39, 100))
    qcap = _arb_exact(Fraction(81, 500))
    while stack:
        lower, upper, depth = stack.pop()
        calls += 1
        upper = _tightened_upper(lower, upper, weights)
        reason = None
        if any(upper[i] < lower[i] for i in range(5)):
            reason = "infeasible"
        elif sum(upper) <= Fraction(27, 16):
            reason = "core"
        elif lower[0] >= Fraction(4, 5):
            reason = "central"
        else:
            center = _feasible_center(lower, upper, weights)
            radii = [max(center[i] - lower[i], upper[i] - center[i]) for i in range(5)]
            knots = _cumulative_knots(center)
            right_center = _divided_difference(knots[1:], lambda v: _f(3, v))
            right_error = qlip * sum(
                (
                    _arb_exact(coefficient * radii[i])
                    for i, coefficient in enumerate(
                        (
                            Fraction(1),
                            Fraction(1),
                            Fraction(3, 4),
                            Fraction(1, 2),
                            Fraction(1, 4),
                        )
                    )
                ),
                arb(0),
            )
            right_lower = right_center - right_error
            if right_lower > 0:
                if _arb_exact(sum(lower)) * right_lower > _arb_exact(upper[0]) * qcap:
                    reason = "rec_global"
                else:
                    left_center = _divided_difference(knots[:4], lambda v: _f(3, v))
                    left_error = qlip * sum(
                        (
                            _arb_exact(coefficient * radii[i])
                            for i, coefficient in enumerate(
                                (
                                    Fraction(1),
                                    Fraction(3, 4),
                                    Fraction(1, 2),
                                    Fraction(1, 4),
                                )
                            )
                        ),
                        arb(0),
                    )
                    left_upper = left_center + left_error
                    if not left_upper < qcap:
                        left_upper = qcap
                    if (
                        _arb_exact(sum(lower)) * right_lower
                        > _arb_exact(upper[0]) * left_upper
                    ):
                        reason = "rec_local"
            if reason is None:
                value = _divided_difference(knots, lambda v: _f(4, v))
                error = _arb_exact(Fraction(16, 125)) * sum(
                    (_arb_exact(weights[i] * radii[i]) for i in range(5)), arb(0)
                )
                if value > error:
                    reason = "direct"
        if reason:
            counts[reason] += 1
            maximum_depth = max(maximum_depth, depth)
            _update_box_digest(digest, reason, lower, upper)
            continue
        if depth >= 55:
            raise RuntimeError(("face5 fail", root_index, lower, upper, depth))
        index = max(range(5), key=lambda i: weights[i] * (upper[i] - lower[i]))
        midpoint = (lower[index] + upper[index]) / 2
        lower_child = list(lower)
        lower_child[index] = midpoint
        upper_child = list(upper)
        upper_child[index] = midpoint
        stack.extend(
            (
                (tuple(lower_child), upper, depth + 1),
                (lower, tuple(upper_child), depth + 1),
            )
        )
    return {
        "i": root_index,
        "calls": calls,
        "counts": counts,
        "depth": maximum_depth,
        "digest": digest.hexdigest(),
    }


class _LocalEnvelope:
    """Directed knot sensitivities from scalar event envelopes.

    Each event moment is bounded by the trivial Dirichlet coordinate moment
    and by every admissible ordered low/high-knot split.  ``base_der`` is
    deliberately left undivided by the number of knots, which is conservative.
    """

    def __init__(self, total, events, base_der, cap):
        self.total = total
        self.tc = _tail_coefficients(total)
        self.events = events
        self.base = Fraction(1, total)
        self.base_der = base_der
        self.cap = cap
        self.cache = None

    def init(self):
        if self.cache is None:
            self.cache = {
                "tc": {
                    h: (
                        tuple((_arb_exact(c), p) for c, p in lo),
                        tuple((_arb_exact(c), p) for c, p in hi),
                    )
                    for h, (lo, hi) in self.tc.items()
                },
                "base": _arb_exact(self.base),
                "base_der": _arb_exact(self.base_der),
                "steps": tuple(_arb_exact(x[2]) for x in self.events),
                "cap": _arb_exact(self.cap),
            }
        return self.cache

    def tail(self, h, r):
        rb = _arb_exact(r)
        lo, hi = self.init()["tc"][h]
        return sum((c * rb**p for c, p in lo), arb(0)), sum(
            (c * rb**p for c, p in hi), arb(0)
        )

    def moments(self, kl, ku, el, eu, limit=7):
        """Bound ``E[W_i 1{el <= sum_j x_j W_j <= eu}]`` for every knot."""
        total = self.total
        base = self.init()["base"]
        ms = [base] * total
        if ku[-1] < el or kl[0] > eu:
            return [arb(0)] * total
        largest = min(ku[-1], Fraction(limit))
        if largest > el:
            for lc in range(1, total):
                cap = ku[lc - 1]
                if cap >= el:
                    continue
                lm, hm = self.tail(total - lc, (largest - el) / (largest - cap))
                for i in range(lc):
                    if lm < ms[i]:
                        ms[i] = lm
                for i in range(lc, total):
                    if hm < ms[i]:
                        ms[i] = hm
        for lc in range(1, total):
            hf = kl[lc]
            if hf <= eu:
                continue
            lt, ht = self.tail(total - lc, 1 - eu / hf)
            lm = base - lt
            hm = base - ht
            for i in range(lc):
                if lm < ms[i]:
                    ms[i] = lm
            for i in range(lc, total):
                if hm < ms[i]:
                    ms[i] = hm
        return ms

    def lips(self, kl, ku):
        """Return directed coordinatewise Lipschitz bounds for one face."""
        out = [self.init()["base_der"]] * self.total
        for step, (el, eu, _) in zip(self.init()["steps"], self.events):
            ms = self.moments(kl, ku, el, eu)
            for i in range(self.total):
                out[i] += step * ms[i]
        return out


ENV6 = _LocalEnvelope(
    5,
    FACE6_Q4_EVENTS,
    Fraction(1, 500),
    Fraction(223, 1000),
)
ENV7 = _LocalEnvelope(
    6,
    FACE7_Q5_EVENTS,
    Fraction(1, 600),
    Fraction(499, 1600),
)


class FaceCallLimitExceeded(RuntimeError):
    """Signal that an audit traversal has proved a root exceeds a call limit."""

    def __init__(self, root_index: int, calls: int):
        super().__init__(root_index, calls)
        self.root_index = root_index
        self.calls = calls


def recursive_worker(task, dim, env, core, central, scalar_k, maxdepth):
    """Certify one six- or seven-positive root by recursive knot insertion."""
    root_index, lower0, upper0, depth0 = task
    ctx.prec = ARB_PRECISION_BITS
    env.init()
    weights = tuple(range(dim, 0, -1))
    stack = [(lower0, upper0, depth0)]
    counts = {
        x: 0 for x in ["core", "central", "rec_global", "rec_local", "infeasible"]
    }
    calls = maximum_depth = 0
    digest = hashlib.sha256()
    cap = env.init()["cap"]
    while stack:
        lower, upper, depth = stack.pop()
        calls += 1
        if FACE7_WORKER_CALL_LIMIT and calls > FACE7_WORKER_CALL_LIMIT:
            raise FaceCallLimitExceeded(root_index, calls)
        upper = _tightened_upper(lower, upper, weights)
        reason = None
        sensitivities = None
        if any(upper[i] < lower[i] for i in range(dim)):
            reason = "infeasible"
        elif sum(upper) <= core:
            reason = "core"
        elif lower[0] >= central:
            reason = "central"
        else:
            center = _feasible_center(lower, upper, weights)
            radii = [
                max(center[i] - lower[i], upper[i] - center[i]) for i in range(dim)
            ]
            knots = _cumulative_knots(center)
            lower_knots = _cumulative_knots(lower)
            upper_knots = _cumulative_knots(upper)
            right_bounds = env.lips(lower_knots[1:], upper_knots[1:])
            right_partials = [sum(right_bounds[i:], arb(0)) for i in range(dim - 1)]
            right_center = _divided_difference(knots[1:], lambda v: _f(scalar_k, v))
            right_error = right_partials[0] * _arb_exact(radii[0]) + sum(
                (right_partials[i - 1] * _arb_exact(radii[i]) for i in range(1, dim)),
                arb(0),
            )
            right_lower = right_center - right_error
            if right_lower > 0:
                if _arb_exact(sum(lower)) * right_lower > _arb_exact(
                    env.cap * upper[0]
                ):
                    reason = "rec_global"
                else:
                    left_bounds = env.lips(
                        lower_knots[: dim - 1], upper_knots[: dim - 1]
                    )
                    left_partials = [
                        sum(left_bounds[i:], arb(0)) for i in range(dim - 1)
                    ]
                    left_center = _divided_difference(
                        knots[: dim - 1], lambda v: _f(scalar_k, v)
                    )
                    left_error = sum(
                        (
                            left_partials[i] * _arb_exact(radii[i])
                            for i in range(dim - 1)
                        ),
                        arb(0),
                    )
                    left_upper = left_center + left_error
                    if not left_upper < cap:
                        left_upper = cap
                    if (
                        _arb_exact(sum(lower)) * right_lower
                        > _arb_exact(upper[0]) * left_upper
                    ):
                        reason = "rec_local"
            sensitivities = (
                right_partials[0],
                right_partials[0],
                *right_partials[1:],
            )
        if reason:
            counts[reason] += 1
            maximum_depth = max(maximum_depth, depth)
            _update_box_digest(digest, reason, lower, upper)
            continue
        if depth >= maxdepth:
            raise RuntimeError(
                (f"face{dim} fail", root_index, lower, upper, depth, right_lower)
            )
        scores = [
            sensitivities[i] * _arb_exact(upper[i] - lower[i]) for i in range(dim)
        ]
        index = 0
        for i in range(1, dim):
            if scores[i] > scores[index]:
                index = i
        midpoint = (lower[index] + upper[index]) / 2
        lower_child = list(lower)
        lower_child[index] = midpoint
        upper_child = list(upper)
        upper_child[index] = midpoint
        stack.extend(
            (
                (tuple(lower_child), upper, depth + 1),
                (lower, tuple(upper_child), depth + 1),
            )
        )
    return {
        "i": root_index,
        "calls": calls,
        "counts": counts,
        "depth": maximum_depth,
        "digest": digest.hexdigest(),
    }


def face6_worker(task):
    return recursive_worker(task, 6, ENV6, Fraction(173, 128), Fraction(47, 64), 4, 75)


def face7_worker(task):
    return recursive_worker(task, 7, ENV7, Fraction(9, 8), Fraction(91, 128), 5, 90)


def _face7_primary_tasks():
    """Return the first-stage seven-positive partition with stable IDs."""
    base = _base_roots(7, FACE7_BASE_ROOT_COUNT)
    hard = set(FACE7_PRIMARY_REFINED_BASE_INDICES)
    boxes = []
    for base_index, box in enumerate(base):
        boxes.extend(
            _refine(box, FACE7_PRIMARY_REFINEMENT, 7) if base_index in hard else [box]
        )
    return [(root_index, *box) for root_index, box in enumerate(boxes)]


def _face_tasks(face):
    """Return the fixed rational roots for one ordered coordinate face."""
    if face == 4:
        return [(0, *_base_roots(4, 1)[0])]
    if face == 5:
        return [(0, *_base_roots(5, 1)[0])]
    if face == 6:
        base = _base_roots(6, 1024)
        boxes = _refine(base[0], 4096, 6) + base[1:]
        return [(i, *box) for i, box in enumerate(boxes)]
    if face == 7:
        primary = _face7_primary_tasks()
        extra = set(FACE7_EXTRA_REFINED_PARENT_IDS)
        if len(extra) != len(FACE7_EXTRA_REFINED_PARENT_IDS):
            raise AssertionError("duplicate extra-refinement parent")
        primary_ids = {task[0] for task in primary}
        if not extra <= primary_ids:
            raise AssertionError("unknown extra-refinement parent")
        next_root_index = len(primary)
        tasks = []
        for root_index, lower, upper, depth in primary:
            if root_index not in extra:
                tasks.append((root_index, lower, upper, depth))
                continue
            children = _refine((lower, upper, depth), FACE7_EXTRA_REFINEMENT, 7)
            tasks.extend(
                (next_root_index + child_index, *child)
                for child_index, child in enumerate(children)
            )
            next_root_index += len(children)
        return tasks
    raise ValueError(face)


def _face7_shard(root_index: int, shard_count: int = FACE7_SHARD_COUNT) -> int:
    """Assign a root to a deterministic hash shard.

    Hashing stable root identifiers spreads spatially adjacent hard boxes
    more evenly than a contiguous or modulo partition while leaving the
    mathematical partition and every branch decision unchanged.
    """
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    payload = f"{FACE7_SHARD_SALT}:{root_index}".encode("ascii")
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return value % shard_count


def _root_partition_digest(face_tasks) -> str:
    root_digest = hashlib.sha256()
    for root_index, lower, upper, initial_depth in face_tasks:
        root_digest.update(f"{root_index}|{initial_depth}|".encode("ascii"))
        for endpoint in (lower, upper):
            for value in endpoint:
                root_digest.update(
                    f"{value.numerator}/{value.denominator},".encode("ascii")
                )
        root_digest.update(b"\n")
    return root_digest.hexdigest()


def _aggregate(face_tasks, results):
    by_index = {result["i"]: result for result in results}
    assert len(by_index) == len(face_tasks)
    transcript_digest = hashlib.sha256()
    calls = 0
    counts = {}
    depth = 0
    for root_index, _lower, _upper, _initial_depth in face_tasks:
        result = by_index[root_index]
        calls += result["calls"]
        depth = max(depth, result["depth"])
        transcript_digest.update(f"{root_index}|{result['digest']}\n".encode("ascii"))
        for key, value in result["counts"].items():
            counts[key] = counts.get(key, 0) + value
    terminals = sum(counts.values())
    assert calls == 2 * terminals - len(face_tasks)
    return {
        "initial_root_boxes": len(face_tasks),
        "total_branch_calls": calls,
        "terminal_box_counts": counts,
        "terminal_boxes": terminals,
        "maximum_bisection_depth": depth,
        "root_partition_sha256": _root_partition_digest(face_tasks),
        "terminal_transcript_sha256": transcript_digest.hexdigest(),
    }


def run_face(
    face: int,
    *,
    workers: int,
    checkpoint: pathlib.Path | None = None,
    shard_index: int | None = None,
    shard_count: int | None = None,
    progress: bool = False,
) -> dict[str, object]:
    """Regenerate one complete face or a fixed hash shard of face seven."""
    if workers < 1:
        raise ValueError("workers must be positive")
    face_tasks = _face_tasks(face)
    if (shard_index is None) != (shard_count is None):
        raise ValueError("shard_index and shard_count must be supplied together")
    if shard_count is not None:
        if face != 7:
            raise ValueError("only the seven-positive face is sharded")
        if shard_count < 1 or not 0 <= shard_index < shard_count:
            raise ValueError("invalid shard")
        face_tasks = [
            task
            for task in face_tasks
            if _face7_shard(task[0], shard_count) == shard_index
        ]

    completed = {}
    if checkpoint is not None and checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            if line.strip():
                result = json.loads(line)
                completed[result["i"]] = result
    required_indices = {task[0] for task in face_tasks}
    completed = {
        index: result
        for index, result in completed.items()
        if index in required_indices
    }
    pending = [task for task in reversed(face_tasks) if task[0] not in completed]
    pending_count = len(pending)
    worker = {4: face4_worker, 5: face5_worker, 6: face6_worker, 7: face7_worker}[face]
    if workers == 1:
        iterator = map(worker, pending)
        pool = None
    else:
        pool = mp.get_context("spawn").Pool(workers)
        iterator = pool.imap_unordered(worker, pending, chunksize=1)
    started = time.monotonic()
    try:
        for completed_now, result in enumerate(iterator, 1):
            if checkpoint is not None:
                checkpoint.parent.mkdir(parents=True, exist_ok=True)
                with checkpoint.open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(result, sort_keys=True) + "\n")
            completed[result["i"]] = result
            if progress and (
                completed_now % 100 == 0
                or result["calls"] >= 100000
                or completed_now == pending_count
            ):
                print(
                    f"face={face} completed={completed_now}/{pending_count} "
                    f"root={result['i']} calls={result['calls']} "
                    f"elapsed_seconds={time.monotonic() - started:.1f}",
                    file=sys.stderr,
                    flush=True,
                )
    finally:
        if pool:
            pool.close()
            pool.join()
    record = _aggregate(face_tasks, list(completed.values()))
    if shard_count is not None:
        record["shard"] = {
            "index": shard_index,
            "count": shard_count,
            "assignment": "sha256(salt:root_index)[0:8] modulo shard_count",
            "salt": FACE7_SHARD_SALT,
        }
    return record


LOWER_FACE_EXPECTATIONS = {
    4: {
        "initial_root_boxes": 1,
        "total_branch_calls": 36057,
        "terminal_box_counts": {
            "central": 1,
            "core": 300,
            "direct": 17629,
            "far": 99,
            "infeasible": 0,
        },
        "terminal_boxes": 18029,
        "maximum_bisection_depth": 23,
        "root_partition_sha256": (
            "576e13355c9bb811ec6c63129886ff8c4fc076dcc3bf705bf74c4a6ab815e516"
        ),
        "terminal_transcript_sha256": (
            "cfd4b924faa5f5ec0e5f8adb9084ebcc5a26dcac742ecae40de2e956e8d784dc"
        ),
    },
    5: {
        "initial_root_boxes": 1,
        "total_branch_calls": 657565,
        "terminal_box_counts": {
            "central": 224,
            "core": 12309,
            "direct": 10093,
            "infeasible": 0,
            "rec_global": 203771,
            "rec_local": 102386,
        },
        "terminal_boxes": 328783,
        "maximum_bisection_depth": 35,
        "root_partition_sha256": (
            "b112af6f22e1770655801cf933de13276284f6c959afb817f1fa28353101b76b"
        ),
        "terminal_transcript_sha256": (
            "c0d540f02ae6e263e69c032b338fffc9f63a3641affdc140ebe9fdbbffff3804"
        ),
    },
    6: {
        "initial_root_boxes": 5119,
        "total_branch_calls": 17760001,
        "terminal_box_counts": {
            "central": 6015,
            "core": 431363,
            "infeasible": 0,
            "rec_global": 4361284,
            "rec_local": 4083898,
        },
        "terminal_boxes": 8882560,
        "maximum_bisection_depth": 44,
        "root_partition_sha256": (
            "0908c297bd6af9acc65a80f4b3fdf7ad78a633b01915e191f3bbe29cf60ba026"
        ),
        "terminal_transcript_sha256": (
            "ddea7fda9cf5a4a28fd13c4e9104bb0697dc44722b68b81e07a06b0f1f9c17a8"
        ),
    },
}

# These locks are filled from a complete regeneration of the final partition.
# Keeping the full-face record separate from the shard records makes both the
# mathematical aggregate and every CI regeneration unit independently
# checkable.
FACE7_EXPECTATION: dict[str, object] = {}
FACE7_SHARD_EXPECTATIONS: dict[int, dict[str, object]] = {}


def _require_face7_expectations() -> None:
    if not FACE7_EXPECTATION:
        raise RuntimeError("the seven-positive full-face lock is not populated")
    if set(FACE7_SHARD_EXPECTATIONS) != set(range(FACE7_SHARD_COUNT)):
        raise RuntimeError("the seven-positive shard locks are not populated")


def verify_face_record(
    face: int,
    record: dict[str, object],
    *,
    shard_index: int | None = None,
    shard_count: int | None = None,
) -> dict[str, object]:
    """Compare one regenerated face record with its committed proof lock."""
    if face in LOWER_FACE_EXPECTATIONS:
        if shard_index is not None or shard_count is not None:
            raise ValueError("lower faces are not sharded")
        expected = LOWER_FACE_EXPECTATIONS[face]
    elif face == 7:
        _require_face7_expectations()
        if shard_index is None and shard_count is None:
            expected = FACE7_EXPECTATION
        else:
            if shard_index is None or shard_count is None:
                raise ValueError("both shard arguments are required")
            if shard_count != FACE7_SHARD_COUNT:
                raise ValueError(
                    f"the locked partition uses {FACE7_SHARD_COUNT} shards"
                )
            expected = FACE7_SHARD_EXPECTATIONS[shard_index]
    else:
        raise ValueError(face)
    if record != expected:
        raise AssertionError(
            "regenerated face record differs from its committed proof lock"
        )
    return record


def verify_lower_faces(workers: int) -> dict[str, object]:
    """Regenerate and lock the four-, five-, and six-positive faces."""
    records = {}
    for face in (4, 5, 6):
        record = run_face(face, workers=workers)
        verify_face_record(face, record)
        records[f"{face}_positive_face"] = record
    return records


def verify_scalar_bounds() -> dict[str, object]:
    """Regenerate the scalar partitions and check their locked summary."""
    record = verify_n7_scalar_bounds()
    for key, expected in SCALAR_EXPECTATION.items():
        assert record[key] == expected
    return record


def build_certificate() -> dict[str, object]:
    """Build the compact, machine-readable complete ``n=7`` certificate."""
    _require_face7_expectations()
    primary_tasks = _face7_primary_tasks()
    return {
        "schema_version": 1,
        "status": (
            "Proves the complete Dirichlet--Poissonization comparison through "
            "n=7; this is not by itself a general-sample-size Stringer "
            "coverage theorem."
        ),
        "theorem": (
            "For every nonnegative eight-coordinate vector y with sum(y)<=7, "
            "P{sum_i y_i E_i>7} >= P{sum_i y_i D_i>1}, where the E_i are "
            "independent unit exponentials and D is uniform Dirichlet."
        ),
        "first_open_complete_simplex_dimension": 8,
        "arithmetic": {
            "arb_precision_bits": ARB_PRECISION_BITS,
            "transcendental_signs": "directed Arb inclusions",
            "partition_endpoints_and_budget_operations": "exact rationals",
            "beta_event_moments": "exact rational polynomials",
        },
        "scalar_bounds": verify_scalar_bounds(),
        "coordinate_faces": {
            **{
                f"{face}_positive_face": expectation
                for face, expectation in LOWER_FACE_EXPECTATIONS.items()
            },
            "7_positive_face": FACE7_EXPECTATION,
        },
        "seven_positive_partition": {
            "base_root_boxes": FACE7_BASE_ROOT_COUNT,
            "primary_refinement_size": FACE7_PRIMARY_REFINEMENT,
            "primary_refined_base_indices": list(FACE7_PRIMARY_REFINED_BASE_INDICES),
            "primary_root_boxes": len(primary_tasks),
            "primary_root_partition_sha256": _root_partition_digest(primary_tasks),
            "extra_refinement_size": FACE7_EXTRA_REFINEMENT,
            "extra_refinement_selection": (
                "replace every primary root whose first complete regeneration "
                "used strictly more than the recorded call cutoff; this is a "
                "load-balancing choice and not a sign-dependent proof rule"
            ),
            "extra_refinement_primary_call_cutoff": (
                FACE7_EXTRA_REFINEMENT_CALL_CUTOFF
            ),
            "extra_refined_primary_root_ids": list(FACE7_EXTRA_REFINED_PARENT_IDS),
            "final_root_boxes": FACE7_EXPECTATION["initial_root_boxes"],
            "final_root_partition_sha256": FACE7_EXPECTATION["root_partition_sha256"],
        },
        "seven_positive_shards": {
            "count": FACE7_SHARD_COUNT,
            "assignment": "sha256(salt:root_index)[0:8] modulo shard count",
            "salt": FACE7_SHARD_SALT,
            "records": [
                FACE7_SHARD_EXPECTATIONS[index] for index in range(FACE7_SHARD_COUNT)
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--face", type=int, choices=[4, 5, 6, 7])
    mode.add_argument("--scalar", action="store_true")
    mode.add_argument("--certificate", action="store_true")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--checkpoint", type=pathlib.Path)
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--shard-count", type=int)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--out", type=pathlib.Path)
    args = parser.parse_args()
    if args.certificate:
        if (
            any(
                value is not None
                for value in (args.checkpoint, args.shard_index, args.shard_count)
            )
            or args.progress
            or args.check
        ):
            parser.error("face-regeneration options do not apply to --certificate")
        record = build_certificate()
    elif args.scalar:
        if (
            any(
                value is not None
                for value in (args.checkpoint, args.shard_index, args.shard_count)
            )
            or args.progress
        ):
            parser.error("face-regeneration options do not apply to --scalar")
        record = verify_scalar_bounds()
    else:
        record = run_face(
            args.face,
            workers=args.workers,
            checkpoint=args.checkpoint,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            progress=args.progress,
        )
        if args.check:
            verify_face_record(
                args.face,
                record,
                shard_index=args.shard_index,
                shard_count=args.shard_count,
            )
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.out is None:
        print(rendered, end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
