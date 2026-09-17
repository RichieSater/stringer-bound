"""Rigorous binary64 enclosures for large tensor-Bernstein certificates.

The sign decision in this module is not an unverified floating-point test.
Integer power coefficients below ``2**53`` are represented exactly.  Every
rational Bernstein weight is compared exactly with its binary64 rounding, and
an exact rational forward-error bound is propagated through each multiply,
add, and midpoint de Casteljau subdivision.  A leaf closes only when its
binary64 lower center exceeds that rational bound strictly.

For an axis transform, the invariant is coefficientwise but is proved using
the l1 norm of the disjoint source slice.  Positive Bernstein weights do not
increase that norm.  The uniform ratio returned after each axis therefore
multiplies the original global l1 norm to give a valid (slightly coarser)
absolute error bound for every coefficient.  Midpoint subdivision preserves
the same invariant.  The code deliberately uses ``2^-52``, twice the usual
binary64 unit roundoff, throughout the standard ``gamma_n`` estimates.
"""

from __future__ import annotations

import gc
import hashlib
import math
import sys
from collections import Counter
from fractions import Fraction

import numpy as np
import sympy as sp


UNIT_ROUNDOFF_BOUND = Fraction(1, 2**52)


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _gamma(operation_count: int) -> Fraction:
    value = operation_count * UNIT_ROUNDOFF_BOUND
    if value >= 1:
        raise ArithmeticError("roundoff estimate left its valid range")
    return value / (1 - value)


def assert_binary64_environment() -> None:
    """Fail unless Python and NumPy expose round-to-nearest binary64."""
    info = np.finfo(np.float64)
    if (
        sys.float_info.radix != 2
        or sys.float_info.mant_dig != 53
        or info.nmant != 52
        or info.maxexp != 1024
        or info.minexp != -1022
    ):
        raise AssertionError("IEEE-754 binary64 is required")
    half_ulp = 2.0**-53
    if 1.0 + half_ulp != 1.0:
        raise AssertionError("round-to-nearest, ties-to-even is required")
    if 1.0 + 3.0 * half_ulp != 1.0 + 2.0**-51:
        raise AssertionError("round-to-nearest, ties-to-even is required")
    one = np.float64(1.0)
    if np.add(one, np.float64(half_ulp)) != one:
        raise AssertionError("NumPy binary64 rounding mode changed")


def _array_digest(array: np.ndarray) -> str:
    normalized = np.ascontiguousarray(array, dtype="<f8")
    return hashlib.sha256(normalized.view(np.uint8)).hexdigest()


def _bernstein64(
    poly: sp.Poly, target_degree: tuple[int, ...]
) -> tuple[np.ndarray, Fraction, int, int, dict[str, object]]:
    """Construct centers and a rational error ratio for one tensor."""
    denominator, integer_poly = poly.clear_denoms(convert=True)
    denominator = int(denominator)
    coefficients = [int(value) for value in integer_poly.coeffs()]
    if not coefficients:
        raise AssertionError("a strict sign certificate cannot use zero")
    if max(abs(value) for value in coefficients) >= 2**53:
        raise AssertionError("a power coefficient is not exact in binary64")
    l1_norm = sum(abs(value) for value in coefficients)
    if l1_norm >= 2**53:
        raise AssertionError("the power l1 norm is not exact in binary64")

    shape = tuple(degree + 1 for degree in target_degree)
    array = np.zeros(shape, dtype=np.float64)
    for powers, coefficient in integer_poly.terms():
        array[powers] = int(coefficient)

    error_ratio = Fraction(0)
    weight_count = 0
    largest_weight_relative_error = Fraction(0)
    with np.errstate(over="raise", under="raise", invalid="raise"):
        for axis, degree in enumerate(target_degree):
            view = np.moveaxis(array, axis, -1)
            output = np.empty_like(view)
            for index in range(degree + 1):
                accumulator = np.zeros(view.shape[:-1], dtype=np.float64)
                for power in range(index + 1):
                    weight = Fraction(
                        math.comb(index, power), math.comb(degree, power)
                    )
                    rounded_weight = np.float64(float(weight))
                    rounded_fraction = Fraction.from_float(
                        float(rounded_weight)
                    )
                    weight_error = abs(rounded_fraction - weight)
                    if (
                        weight_error
                        > UNIT_ROUNDOFF_BOUND * abs(rounded_fraction)
                    ):
                        raise AssertionError("Bernstein weight bound failed")
                    if not 0.0 <= rounded_weight <= 1.0:
                        raise AssertionError("Bernstein weight left [0,1]")
                    if (
                        rounded_weight
                        and abs(rounded_weight) < np.finfo(np.float64).tiny
                    ):
                        raise AssertionError("subnormal Bernstein weight")
                    relative_error = (
                        weight_error / abs(rounded_fraction)
                        if rounded_fraction
                        else Fraction(0)
                    )
                    largest_weight_relative_error = max(
                        largest_weight_relative_error, relative_error
                    )
                    weight_count += 1
                    product = np.multiply(view[..., power], rounded_weight)
                    accumulator = np.add(accumulator, product)
                output[..., index] = accumulator
            del view, array
            array = np.moveaxis(output, -1, axis)
            if not np.isfinite(array).all():
                raise AssertionError("nonfinite Bernstein center")
            error_ratio = (
                error_ratio
                + UNIT_ROUNDOFF_BOUND * (1 + error_ratio)
                + _gamma(2 * (degree + 1))
                * (1 + UNIT_ROUNDOFF_BOUND)
                * (1 + error_ratio)
            )
            gc.collect()

    return (
        array,
        error_ratio,
        l1_norm,
        denominator,
        {
            "weight_count": weight_count,
            "largest_weight_relative_error": _fraction_record(
                largest_weight_relative_error
            ),
        },
    )


def _split_pair(
    array: np.ndarray, axis: int
) -> tuple[np.ndarray, np.ndarray]:
    degree = array.shape[axis] - 1
    work = np.moveaxis(array, axis, 0).copy()
    left = np.empty_like(work)
    right = np.empty_like(work)
    left[0] = work[0]
    right[degree] = work[degree]
    with np.errstate(over="raise", under="raise", invalid="raise"):
        for level in range(1, degree + 1):
            work[: degree - level + 1] = np.multiply(
                np.add(
                    work[: degree - level + 1],
                    work[1 : degree - level + 2],
                ),
                np.float64(0.5),
            )
            left[level] = work[0]
            right[degree - level] = work[degree - level]
    del work
    left = np.moveaxis(left, 0, axis)
    right = np.moveaxis(right, 0, axis)
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise AssertionError("nonfinite subdivision center")
    return left, right


def _subdivided_error(
    error_ratio: Fraction, degree: int
) -> Fraction:
    # Scaling by 1/2 is exact for normal binary64 values.  The enclosing
    # transforms fail on underflow, so only the rounded addition contributes.
    for _ in range(degree):
        error_ratio += UNIT_ROUNDOFF_BOUND * (1 + error_ratio)
    return error_ratio


def _witness_trie(leaves) -> tuple[dict[str, object], dict[str, int]]:
    root: dict[str, object] = {}
    for path, closure, multiplier in leaves:
        if any(bit not in "01" for bit in path):
            raise AssertionError("nonbinary subdivision path")
        if closure not in {"T", "G4", "M4"}:
            raise AssertionError("unknown leaf closure")
        if closure == "M4":
            if multiplier is None or Fraction(*multiplier) < 0:
                raise AssertionError("invalid nonnegative multiplier")
        elif multiplier is not None:
            raise AssertionError("unexpected multiplier")
        node = root
        for bit in path:
            if "leaf" in node:
                raise AssertionError("a leaf path is a prefix")
            node = node.setdefault(bit, {})
        if node:
            raise AssertionError("duplicate path or path with descendants")
        node["leaf"] = (closure, multiplier)

    def check(node) -> tuple[int, int, int]:
        if "leaf" in node:
            if len(node) != 1:
                raise AssertionError("leaf has descendants")
            return 1, 1, 0
        if set(node) != {"0", "1"}:
            raise AssertionError("subdivision tree does not cover the cube")
        left_nodes, left_leaves, left_depth = check(node["0"])
        right_nodes, right_leaves, right_depth = check(node["1"])
        return (
            1 + left_nodes + right_nodes,
            left_leaves + right_leaves,
            1 + max(left_depth, right_depth),
        )

    nodes, leaf_count, maximum_depth = check(root)
    if leaf_count != len(leaves):
        raise AssertionError("leaf count mismatch")
    return root, {
        "nodes": nodes,
        "leaves": leaf_count,
        "maximum_depth": maximum_depth,
    }


def verify_binary64_bernstein_tree(
    target: sp.Poly,
    gap4: sp.Poly,
    leaves,
    axis_cycle: tuple[int, ...],
) -> dict[str, object]:
    """Certify a dyadic cover by exact forward-error inequalities.

    ``T`` closes a leaf by target positivity; ``G4`` closes it because the
    ordered-centroid condition would require ``gap4 >= 0``; and ``M4`` uses
    ``target - mu*gap4 > 0`` for a nonnegative integer ``mu``.
    """
    assert_binary64_environment()
    target_degree = tuple(map(int, target.degree_list()))
    if not axis_cycle or any(
        axis not in range(len(target_degree)) for axis in axis_cycle
    ):
        raise AssertionError("invalid subdivision-axis cycle")
    root, tree_record = _witness_trie(leaves)

    (
        target_array,
        target_error,
        target_l1,
        target_scale,
        target_weight_record,
    ) = _bernstein64(target, target_degree)
    (
        gap_array,
        gap_error,
        gap_l1,
        gap_scale,
        gap_weight_record,
    ) = _bernstein64(gap4, target_degree)
    if target_scale != gap_scale:
        raise AssertionError("target and gap need the same integer scale")

    target_digest = _array_digest(target_array)
    gap_digest = _array_digest(gap_array)
    target_range = {
        "smallest_center": _fraction_record(
            Fraction.from_float(float(target_array.min()))
        ),
        "largest_center": _fraction_record(
            Fraction.from_float(float(target_array.max()))
        ),
    }
    gap_range = {
        "smallest_center": _fraction_record(
            Fraction.from_float(float(gap_array.min()))
        ),
        "largest_center": _fraction_record(
            Fraction.from_float(float(gap_array.max()))
        ),
    }
    records = []

    def walk(
        node,
        local_target,
        local_gap,
        target_ratio: Fraction,
        gap_ratio: Fraction,
        path: str,
    ) -> None:
        if "leaf" in node:
            closure, multiplier_pair = node["leaf"]
            target_absolute_error = target_ratio * target_l1
            gap_absolute_error = gap_ratio * gap_l1
            multiplier = None
            if closure == "T":
                stored_lower = Fraction.from_float(
                    float(local_target.min())
                )
                error = target_absolute_error
            elif closure == "G4":
                stored_lower = -Fraction.from_float(
                    float(local_gap.max())
                )
                error = gap_absolute_error
            else:
                multiplier_fraction = Fraction(*multiplier_pair)
                if (
                    multiplier_fraction.denominator != 1
                    or multiplier_fraction.numerator >= 2**53
                ):
                    raise AssertionError(
                        "binary64 enclosure requires an exact integer multiplier"
                    )
                multiplier = multiplier_fraction.numerator
                # The nextafter endpoints enclose the exact operations on the
                # stored dyadic centers.  Producing the smallest subnormal from
                # zero is intentional here, so gradual underflow is allowed.
                with np.errstate(
                    over="raise", under="ignore", invalid="raise"
                ):
                    product_upper = np.nextafter(
                        np.multiply(np.float64(multiplier), local_gap),
                        np.float64(np.inf),
                    )
                    difference_lower = np.nextafter(
                        np.subtract(local_target, product_upper),
                        np.float64(-np.inf),
                    )
                if not np.isfinite(difference_lower).all():
                    raise AssertionError("nonfinite multiplier enclosure")
                stored_lower = Fraction.from_float(
                    float(difference_lower.min())
                )
                error = (
                    target_absolute_error
                    + multiplier * gap_absolute_error
                )
                del product_upper, difference_lower
            margin = stored_lower - error
            if margin <= 0:
                raise AssertionError(
                    ("uncertified Bernstein leaf", path, closure, margin)
                )
            records.append(
                {
                    "path": path,
                    "closure": closure,
                    "multiplier": multiplier,
                    "stored_lower_bound": _fraction_record(stored_lower),
                    "absolute_error_bound": _fraction_record(error),
                    "certified_margin": _fraction_record(margin),
                }
            )
            return

        axis = axis_cycle[len(path) % len(axis_cycle)]
        target_left, target_right = _split_pair(local_target, axis)
        gap_left, gap_right = _split_pair(local_gap, axis)
        degree = target_degree[axis]
        next_target_error = _subdivided_error(target_ratio, degree)
        next_gap_error = _subdivided_error(gap_ratio, degree)
        del local_target, local_gap
        gc.collect()
        walk(
            node["0"],
            target_left,
            gap_left,
            next_target_error,
            next_gap_error,
            path + "0",
        )
        del target_left, gap_left
        gc.collect()
        walk(
            node["1"],
            target_right,
            gap_right,
            next_target_error,
            next_gap_error,
            path + "1",
        )
        del target_right, gap_right
        gc.collect()

    walk(
        root,
        target_array,
        gap_array,
        target_error,
        gap_error,
        "",
    )
    margins = [
        Fraction(
            int(record["certified_margin"]["numerator"]),
            int(record["certified_margin"]["denominator"]),
        )
        for record in records
    ]
    closure_counts = Counter(record["closure"] for record in records)
    return {
        "method": (
            "IEEE-754 binary64 centers with exact rational forward-error "
            "bounds and outward-rounded multiplier differences"
        ),
        "unit_roundoff_bound": _fraction_record(UNIT_ROUNDOFF_BOUND),
        "axis_cycle": list(axis_cycle),
        "subdivision": tree_record,
        "closure_counts": dict(sorted(closure_counts.items())),
        "minimum_certified_margin": _fraction_record(min(margins)),
        "target": {
            "degree": list(target_degree),
            "power_terms": len(target.terms()),
            "integer_scale": target_scale,
            "power_l1_norm": target_l1,
            "initial_error_ratio": _fraction_record(target_error),
            "binary64_bernstein_sha256": target_digest,
            **target_range,
            **target_weight_record,
        },
        "gap4": {
            "native_degree": list(map(int, gap4.degree_list())),
            "degree_elevated_to": list(target_degree),
            "power_terms": len(gap4.terms()),
            "integer_scale": gap_scale,
            "power_l1_norm": gap_l1,
            "initial_error_ratio": _fraction_record(gap_error),
            "binary64_bernstein_sha256": gap_digest,
            **gap_range,
            **gap_weight_record,
        },
        "leaves": records,
    }
