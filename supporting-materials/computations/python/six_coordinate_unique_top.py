"""Exact certificate for the six-coordinate unique-top ``L=4`` family.

The normalized knots are

    (-1, -u, -v, -w, t*a, t),

with ``1 >= u >= v >= w > 0``, ``0 < a < 1``, and ``t > 0``.  Confluent
divided differences give the cap probability, its section centroid, five
adjacent centroid-gap polynomials ``G_i``, and a comparison polynomial ``T``.
On the ordered relative interior, ``G_i >= 0`` and

    T >= 0  <=>  p <= I_{gamma_0}(1,5) = 1-(1-gamma_0)^5.

The singular small-scale corner is resolved by scale and projective charts.
Every remaining polynomial sign is checked in an exact tensor-Bernstein
basis over a finite rational-box cover.  Floating point is not used by this
module.

The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative.
The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
from fractions import Fraction
from math import prod
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx

from five_coordinate_mixed_section import (
    _fraction_record,
    _power_digest,
    exact_bernstein_coefficients,
    exact_coefficient_statistics,
)
from six_coordinate_face32 import _confluent_divided_difference
from six_coordinate_unique_top_witness import (
    B_GE_V_FINAL_LEAVES,
    LOWER_LOW_Z_LEAVES,
    LOWER_MID_CORNER_LEAVES,
    MIDDLE_LOW_Z_LEAVES,
    Z_LE_A_FINAL_LEAVES,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates" / "six-coordinate-unique-top-certificate.json"
)


def _poly(expression, variables) -> sp.Poly:
    return sp.Poly(expression, *variables, domain=sp.QQ)


def _polynomial_record(poly: sp.Poly) -> dict[str, object]:
    return {
        "degree": list(map(int, poly.degree_list())),
        "power_terms": len(poly.terms()),
        "sha256_power_expansion": _power_digest(poly),
    }


def _coefficientwise_positive(poly: sp.Poly) -> None:
    if not poly.terms() or any(coefficient <= 0 for _, coefficient in poly.terms()):
        raise AssertionError("expected a coefficientwise-positive polynomial")


def _remove_positive_monomial(poly: sp.Poly) -> tuple[sp.Poly, tuple[int, ...]]:
    minima = tuple(
        min(powers[index] for powers, _ in poly.terms())
        for index in range(len(poly.gens))
    )
    monomial = _poly(
        sp.prod(poly.gens[index] ** minima[index] for index in range(len(minima))),
        poly.gens,
    )
    return sp.exquo(poly, monomial), minima


def _chart_polynomial(
    poly: sp.Poly,
    substitution: dict[sp.Symbol, sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Poly, tuple[int, ...]]:
    return _remove_positive_monomial(
        _poly(poly.as_expr().subs(substitution), variables)
    )


def _chart_family(
    target: sp.Poly,
    gaps: tuple[sp.Poly, ...],
    substitution: dict[sp.Symbol, sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> dict[str, object]:
    chart_target, target_monomial = _chart_polynomial(
        target, substitution, variables
    )
    gaps_with_monomials = tuple(
        _chart_polynomial(gap, substitution, variables) for gap in gaps
    )
    return {
        "variables": variables,
        "target": chart_target,
        "gaps": tuple(value[0] for value in gaps_with_monomials),
        "target_monomial": target_monomial,
        "gap_monomials": tuple(value[1] for value in gaps_with_monomials),
    }


def _bernstein_statistics(
    poly: sp.Poly,
    axes: tuple[int, ...] | None = None,
) -> dict[str, object]:
    coefficients = exact_bernstein_coefficients(poly, axes=axes)
    return exact_coefficient_statistics(
        coefficients, tuple(map(int, poly.degree_list()))
    )


def derive_unique_top_l4() -> dict[str, object]:
    """Derive the cap target and five ordered-centroid gap polynomials."""
    u, v, w, a, t = sp.symbols("u v w a t")
    variables = (u, v, w, a, t)
    specialized_knots = (-1, -u, -v, -w, t * a, t)
    abstract_knots = sp.symbols("Y0:6")
    specialization = dict(zip(abstract_knots, specialized_knots, strict=True))

    cap_abstract = _confluent_divided_difference(
        abstract_knots, (1,) * 6, range(4, 6), 5
    )
    density_abstract = 5 * _confluent_divided_difference(
        abstract_knots, (1,) * 6, range(4, 6), 4
    )
    cap = sp.cancel(cap_abstract.subs(specialization))
    density = sp.cancel(density_abstract.subs(specialization))
    centroids = tuple(
        sp.cancel(sp.diff(cap_abstract, knot).subs(specialization) / density)
        for knot in abstract_knots
    )

    gamma_numerator, gamma_denominator = (
        _poly(expression, variables)
        for expression in sp.fraction(centroids[0])
    )
    cap_numerator, cap_denominator = (
        _poly(expression, variables) for expression in sp.fraction(cap)
    )
    _coefficientwise_positive(gamma_denominator)
    _coefficientwise_positive(cap_denominator)

    raw_target = (
        gamma_denominator**5
        - (gamma_denominator - gamma_numerator) ** 5
    ) * cap_denominator - cap_numerator * gamma_denominator**5
    common = sp.gcd(raw_target, gamma_denominator**5 * cap_denominator)
    expected_common = _poly((t + 1) * (a * t + 1), variables)
    if common != expected_common:
        raise AssertionError("unique-top positive target factor changed")
    _coefficientwise_positive(common)
    reduced_target = sp.exquo(raw_target, common)
    reduced_target, initial_target_monomial = _remove_positive_monomial(
        reduced_target
    )
    target_content, reduced_target = reduced_target.primitive()
    if target_content <= 0:
        raise AssertionError("unique-top target acquired negative content")
    if (
        tuple(map(int, reduced_target.degree_list())),
        len(reduced_target.terms()),
        initial_target_monomial,
    ) != ((7, 7, 7, 23, 30), 87_047, (0, 0, 0, 0, 1)):
        raise AssertionError("unique-top target structure changed")

    U, V, W, A, q = sp.symbols("U V W A q")
    cube_variables = (U, V, W, A, q)
    ordered_substitution = {
        u: U,
        v: U * V,
        w: U * V * W,
        a: A,
        t: q,
    }
    cube_target, cube_target_monomial = _chart_polynomial(
        reduced_target, ordered_substitution, cube_variables
    )
    if (
        tuple(map(int, cube_target.degree_list())),
        len(cube_target.terms()),
        cube_target_monomial,
    ) != ((21, 14, 7, 23, 30), 87_047, (0, 0, 0, 0, 0)):
        raise AssertionError("unique-top cube target structure changed")

    # The raw adjacent-gap numerators are these nonnegative collision factors
    # times G_i.  Their denominators have strictly positive power coefficients.
    gap_cofactors = (
        _poly(q * (1 - U), cube_variables),
        _poly(U * q * (1 - V), cube_variables),
        _poly(U * V * q * (1 - W), cube_variables),
        _poly(1, cube_variables),
        _poly(1 - A, cube_variables),
    )
    expected_gaps = (
        ((4, 2, 1, 5, 6), 74),
        ((5, 3, 1, 5, 6), 74),
        ((5, 4, 2, 5, 6), 74),
        ((7, 5, 3, 7, 8), 218),
        ((9, 6, 3, 6, 10), 572),
    )
    gaps = []
    gap_denominators = []
    cofactor_statistics = []
    for index, cofactor in enumerate(gap_cofactors):
        gap = sp.cancel(centroids[index + 1] - centroids[index])
        numerator, denominator = sp.fraction(gap)
        numerator_poly = _poly(
            numerator.subs(ordered_substitution), cube_variables
        )
        denominator_poly = _poly(
            denominator.subs(ordered_substitution), cube_variables
        )
        _coefficientwise_positive(denominator_poly)
        gap_poly = sp.exquo(numerator_poly, cofactor)
        gap_content, gap_poly = gap_poly.primitive()
        if gap_content <= 0:
            raise AssertionError("unique-top gap acquired negative content")
        if (
            tuple(map(int, gap_poly.degree_list())),
            len(gap_poly.terms()),
        ) != expected_gaps[index]:
            raise AssertionError(f"unique-top gap {index} structure changed")
        if numerator_poly != gap_content * cofactor * gap_poly:
            raise AssertionError("unique-top gap factorization changed")
        gaps.append(gap_poly)
        gap_denominators.append(denominator_poly)
        cofactor_statistics.append(_bernstein_statistics(cofactor))

    rho = sp.symbols("rho")
    global_scale = _poly(
        (-gaps[4].as_expr()).subs(q, 3 + rho), (U, V, W, A, rho)
    )
    global_scale_statistics = _bernstein_statistics(
        global_scale, axes=(0, 1, 2, 3)
    )
    if (
        tuple(map(int, global_scale.degree_list())),
        len(global_scale.terms()),
        global_scale_statistics["positive_coefficients"],
        global_scale_statistics["zero_coefficients"],
        global_scale_statistics["negative_coefficients"],
    ) != ((9, 6, 3, 6, 10), 1_607, 10_279, 11_281, 0):
        raise AssertionError("unique-top global scale confinement changed")

    return {
        "variables": variables,
        "cube_variables": cube_variables,
        "cap": cap,
        "gamma0": centroids[0],
        "target": cube_target,
        "gaps": tuple(gaps),
        "gap_denominators": tuple(gap_denominators),
        "gap_cofactors": gap_cofactors,
        "cofactor_statistics": tuple(cofactor_statistics),
        "initial_target_monomial": initial_target_monomial,
        "cube_target_monomial": cube_target_monomial,
        "global_scale": global_scale,
        "global_scale_statistics": global_scale_statistics,
    }


def _confinement(
    gap: sp.Poly,
    substitution: dict[sp.Symbol, sp.Expr],
    variables: tuple[sp.Symbol, ...],
    bounded_axes: tuple[int, ...] | None,
    expected: tuple[tuple[int, ...], int, int, int],
) -> dict[str, object]:
    polynomial = _poly((-gap.as_expr()).subs(substitution), variables)
    statistics = _bernstein_statistics(polynomial, axes=bounded_axes)
    actual = (
        tuple(map(int, polynomial.degree_list())),
        len(polynomial.terms()),
        int(statistics["positive_coefficients"]),
        int(statistics["zero_coefficients"]),
    )
    if actual != expected or statistics["negative_coefficients"] != 0:
        raise AssertionError("unique-top scale confinement changed")
    return {"polynomial": polynomial, "statistics": statistics}


def derive_unique_top_l4_charts(derived: dict[str, object]) -> dict[str, object]:
    """Construct the scale and projective charts used by the proof."""
    U, V, W, A, q = derived["cube_variables"]
    target = derived["target"]
    gaps = derived["gaps"]
    Z, R, B, S, C, D, Q, z, rho = sp.symbols("Z R B S C D Q z rho")

    q_le_u = _chart_family(
        target, gaps, {q: U * Z}, (U, V, W, A, Z)
    )
    u_le_q = _chart_family(
        target, gaps, {U: q * Z}, (q, V, W, A, Z)
    )

    # Z <= A: successive consequences of the terminal gap turn the singular
    # corner into q=12*A*B^3*D and Z=A*B.
    z_le_a = _chart_family(
        u_le_q["target"], u_le_q["gaps"], {Z: A * B}, (q, V, W, A, B)
    )
    z_le_a_q_bound = _confinement(
        z_le_a["gaps"][4],
        {q: 2 * A + rho},
        (A, B, V, W, rho),
        (0, 1, 2, 3),
        ((9, 9, 6, 3, 3), 924, 9_380, 1_820),
    )
    z_le_a_r = _chart_family(
        z_le_a["target"], z_le_a["gaps"], {q: 2 * A * R},
        (A, B, V, W, R),
    )
    z_le_a_r_bound = _confinement(
        z_le_a_r["gaps"][4],
        {R: B + rho},
        (A, B, V, W, rho),
        None,
        ((8, 12, 6, 3, 3), 882, 12_852, 252),
    )
    z_le_a_s = _chart_family(
        z_le_a_r["target"], z_le_a_r["gaps"], {R: B * S},
        (A, B, V, W, S),
    )
    z_le_a_s_bound = _confinement(
        z_le_a_s["gaps"][4],
        {S: B + rho},
        (A, B, V, W, rho),
        None,
        ((8, 14, 6, 3, 3), 1_073, 14_868, 252),
    )
    z_le_a_c = _chart_family(
        z_le_a_s["target"], z_le_a_s["gaps"], {S: B * C},
        (A, B, V, W, C),
    )
    z_le_a_c_bound = _confinement(
        z_le_a_c["gaps"][4],
        {C: 6 * B + rho},
        (A, B, V, W, rho),
        None,
        ((8, 16, 6, 3, 3), 1_203, 16_884, 252),
    )
    z_le_a_final = _chart_family(
        z_le_a_c["target"], z_le_a_c["gaps"], {C: 6 * B * D},
        (A, B, V, W, D),
    )

    # A <= Z: compare B=A/Z with V and with V*W.  The three charts below
    # cover B>=V, V*W<=B<=V, and B<=V*W, respectively.
    a_le_z = _chart_family(
        u_le_q["target"], u_le_q["gaps"], {A: Z * B},
        (q, V, W, Z, B),
    )
    b_ge_v = _chart_family(
        a_le_z["target"], a_le_z["gaps"], {V: B * C},
        (q, W, Z, B, C),
    )
    b_ge_v_bound = _confinement(
        b_ge_v["gaps"][4],
        {q: 2 * Z + rho},
        (Z, W, B, C, rho),
        (0, 1, 2, 3),
        ((9, 3, 6, 6, 3), 882, 6_356, 1_484),
    )
    b_ge_v_final = _chart_family(
        b_ge_v["target"], b_ge_v["gaps"], {q: 2 * Z * R},
        (R, W, Z, B, C),
    )

    middle = _chart_family(
        a_le_z["target"], a_le_z["gaps"],
        {B: V * D, W: D * C}, (q, V, Z, D, C),
    )
    middle_low_bound = _confinement(
        middle["gaps"][4],
        {Z: z / 2, q: 2 * z + rho},
        (z, V, D, C, rho),
        (0, 1, 2, 3),
        ((9, 6, 6, 3, 3), 890, 5_816, 2_024),
    )
    middle_low = _chart_family(
        middle["target"], middle["gaps"],
        {Z: z / 2, q: 2 * z * R}, (R, V, z, D, C),
    )
    middle_scale = _chart_family(
        middle["target"], middle["gaps"], {q: 3 * Q},
        (Q, V, Z, D, C),
    )

    lower = _chart_family(
        a_le_z["target"], a_le_z["gaps"], {B: V * W * C},
        (q, V, W, Z, C),
    )
    lower_low_bound = _confinement(
        lower["gaps"][4],
        {Z: z / 3, q: 3 * z + rho},
        (z, V, W, C, rho),
        (0, 1, 2, 3),
        ((9, 6, 6, 6, 3), 890, 9_460, 4_260),
    )
    lower_low = _chart_family(
        lower["target"], lower["gaps"],
        {Z: z / 3, q: 3 * z * R}, (R, V, W, z, C),
    )
    lower_scale = _chart_family(
        lower["target"], lower["gaps"], {q: 3 * Q},
        (Q, V, W, Z, C),
    )

    expected_targets = {
        "q_le_u": ((14, 14, 7, 23, 30), 87_047),
        "z_le_a_final": ((36, 61, 14, 7, 14), 87_047),
        "b_ge_v_final": ((14, 7, 36, 23, 14), 87_047),
        "middle_low": ((14, 23, 36, 23, 7), 87_047),
        "middle_scale": ((14, 23, 23, 23, 7), 87_047),
        "lower_low": ((14, 23, 23, 36, 23), 87_047),
        "lower_scale": ((14, 23, 23, 23, 23), 87_047),
    }
    charts = {
        "q_le_u": q_le_u,
        "u_le_q": u_le_q,
        "z_le_a": z_le_a,
        "z_le_a_r": z_le_a_r,
        "z_le_a_s": z_le_a_s,
        "z_le_a_c": z_le_a_c,
        "z_le_a_final": z_le_a_final,
        "a_le_z": a_le_z,
        "b_ge_v": b_ge_v,
        "b_ge_v_final": b_ge_v_final,
        "middle": middle,
        "middle_low": middle_low,
        "middle_scale": middle_scale,
        "lower": lower,
        "lower_low": lower_low,
        "lower_scale": lower_scale,
    }
    for name, expected in expected_targets.items():
        chart_target = charts[name]["target"]
        actual = (
            tuple(map(int, chart_target.degree_list())),
            len(chart_target.terms()),
        )
        if actual != expected:
            raise AssertionError(f"unique-top chart {name} structure changed")

    return {
        **charts,
        "confinements": {
            "z_le_a_q": z_le_a_q_bound,
            "z_le_a_r": z_le_a_r_bound,
            "z_le_a_s": z_le_a_s_bound,
            "z_le_a_c": z_le_a_c_bound,
            "b_ge_v": b_ge_v_bound,
            "middle_low": middle_low_bound,
            "lower_low": lower_low_bound,
        },
    }


def _validate_leaf_partition(leaves, allowed: set[str]) -> dict[str, int]:
    trie: dict[str, object] = {}
    maximum_depth = 0
    for path, closure, multiplier in leaves:
        if not path or any(bit not in "01" for bit in path):
            raise AssertionError("invalid binary subdivision path")
        if closure not in allowed:
            raise AssertionError("invalid unique-top leaf closure")
        if closure.startswith("M"):
            if multiplier is None or Fraction(*multiplier) < 0:
                raise AssertionError("invalid unique-top multiplier")
        elif multiplier is not None:
            raise AssertionError("unexpected multiplier on a direct leaf")
        maximum_depth = max(maximum_depth, len(path))
        node = trie
        for bit in path:
            if "leaf" in node:
                raise AssertionError("subdivision witness is not prefix-free")
            node = node.setdefault(bit, {})
        if node:
            raise AssertionError("subdivision witness is not prefix-free")
        node["leaf"] = True

    def check_complete(node) -> tuple[int, int]:
        if "leaf" in node:
            if len(node) != 1:
                raise AssertionError("leaf has descendants")
            return 1, 1
        if set(node) != {"0", "1"}:
            raise AssertionError("subdivision tree does not cover the box")
        left_nodes, left_leaves = check_complete(node["0"])
        right_nodes, right_leaves = check_complete(node["1"])
        return 1 + left_nodes + right_nodes, left_leaves + right_leaves

    nodes, leaf_count = check_complete(trie)
    if leaf_count != len(leaves):
        raise AssertionError("leaf count mismatch")
    return {"nodes": nodes, "leaves": leaf_count, "maximum_depth": maximum_depth}


def _to_flint_polynomial(poly: sp.Poly, context):
    return context.from_dict(
        {powers: fmpq(str(coefficient)) for powers, coefficient in poly.terms()}
    )


def _flint_fraction_record(value) -> dict[str, str]:
    return {"numerator": str(int(value.p)), "denominator": str(int(value.q))}


def _dyadic_box(
    path: str,
    order: tuple[int, ...],
    base_box: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    intervals = list(base_box)
    for depth, bit in enumerate(path):
        axis = order[depth % len(order)]
        lower, upper = intervals[axis]
        middle = (lower + upper) / 2
        intervals[axis] = (middle, upper) if bit == "1" else (lower, middle)
    return tuple(intervals)


def _local_bernstein_sign_record(
    poly,
    variables: tuple[sp.Symbol, ...],
    box: tuple[tuple[Fraction, Fraction], ...],
    degree: tuple[int, ...],
) -> dict[str, object]:
    generators = poly.context().gens()
    substitutions = []
    for generator, (lower, upper) in zip(generators, box, strict=True):
        width = upper - lower
        substitutions.append(
            fmpq(lower.numerator, lower.denominator)
            + fmpq(width.numerator, width.denominator) * generator
        )
    local = poly.compose(*substitutions)

    count = len(variables)
    homogeneous_context = fmpq_mpoly_ctx.get(
        tuple(f"bz{index}" for index in range(count))
        + tuple(f"bh{index}" for index in range(count)),
        "lex",
    )
    homogeneous = homogeneous_context.from_dict(
        {
            tuple(powers)
            + tuple(degree[index] - powers[index] for index in range(count)):
            coefficient
            for powers, coefficient in local.terms()
        }
    )
    del local
    transformed = homogeneous.compose(
        *(tuple(generators) + tuple(1 + generator for generator in generators))
    )
    del homogeneous

    coefficients = transformed.coeffs()
    positive = sum(coefficient > 0 for coefficient in coefficients)
    negative = sum(coefficient < 0 for coefficient in coefficients)
    smallest = min(coefficients) if coefficients else None
    largest = max(coefficients) if coefficients else None
    total = prod(value + 1 for value in degree)
    del coefficients, transformed
    gc.collect()
    if negative:
        raise AssertionError(
            f"unique-top box has {negative} negative Bernstein coefficients"
        )
    if smallest is None or largest is None or smallest <= 0:
        raise AssertionError("unique-top box certificate is not strict")
    return {
        "degree": list(degree),
        "total_coefficients": total,
        "positive_coefficients": positive,
        "zero_coefficients": total - positive,
        "negative_coefficients": 0,
        "smallest_positive_scaled": _flint_fraction_record(smallest),
        "largest_scaled": _flint_fraction_record(largest),
    }


def _verify_boxes(
    name: str,
    chart: dict[str, object],
    boxes: tuple[
        tuple[
            str,
            str,
            tuple[int, int] | None,
            tuple[tuple[Fraction, Fraction], ...],
        ],
        ...,
    ],
    order: tuple[int, ...] | None = None,
    tree=None,
) -> dict[str, object]:
    variables = chart["variables"]
    target = chart["target"]
    gaps = chart["gaps"]
    context = fmpq_mpoly_ctx.get(tuple(map(str, variables)), "lex")
    target_flint = _to_flint_polynomial(target, context)
    gap_flint = tuple(_to_flint_polynomial(gap, context) for gap in gaps)
    target_degree = tuple(map(int, target.degree_list()))

    if tree is not None:
        tree_record = _validate_leaf_partition(tree, {"T", "G4", "M4"})
    else:
        tree_record = None

    closure_counts: dict[str, int] = {}
    total_coefficients = 0
    positive_coefficients = 0
    zero_coefficients = 0
    smallest_positive = None
    largest = None
    leaf_digest = hashlib.sha256()
    multiplier_records = []
    for path, closure, multiplier_pair, base_box in boxes:
        closure_counts[closure] = closure_counts.get(closure, 0) + 1
        if closure == "T":
            source = target_flint
            degree = target_degree
        elif closure == "G4":
            source = -gap_flint[4]
            degree = tuple(map(int, gaps[4].degree_list()))
        elif closure == "M4":
            multiplier = Fraction(*multiplier_pair)
            source = target_flint - fmpq(
                multiplier.numerator, multiplier.denominator
            ) * gap_flint[4]
            degree = target_degree
            multiplier_records.append(
                {
                    "path": path,
                    "multiplier": _fraction_record(multiplier),
                }
            )
        else:
            raise AssertionError("unknown unique-top closure")
        box = (
            _dyadic_box(path, order, base_box)
            if path
            else base_box
        )
        record = _local_bernstein_sign_record(source, variables, box, degree)
        leaf_digest.update(
            (
                f"{name}:{path}:{closure}:{multiplier_pair}:{box}:"
                f"{record['degree']}:{record['positive_coefficients']}:"
                f"{record['zero_coefficients']}:"
                f"{record['smallest_positive_scaled']}:"
                f"{record['largest_scaled']}\n"
            ).encode()
        )
        total_coefficients += int(record["total_coefficients"])
        positive_coefficients += int(record["positive_coefficients"])
        zero_coefficients += int(record["zero_coefficients"])
        local_smallest = fmpq(
            int(record["smallest_positive_scaled"]["numerator"]),
            int(record["smallest_positive_scaled"]["denominator"]),
        )
        local_largest = fmpq(
            int(record["largest_scaled"]["numerator"]),
            int(record["largest_scaled"]["denominator"]),
        )
        if smallest_positive is None or local_smallest < smallest_positive:
            smallest_positive = local_smallest
        if largest is None or local_largest > largest:
            largest = local_largest

    return {
        "box_count": len(boxes),
        "subdivision": (
            {"axis_cycle": list(order), **tree_record} if tree_record else None
        ),
        "closure_counts": closure_counts,
        "exact_multipliers": multiplier_records,
        "bernstein_signs": {
            "basis": (
                "local tensor-Bernstein basis; coefficients are recorded "
                "after positive binomial scaling"
            ),
            "total_coefficients": total_coefficients,
            "positive_coefficients": positive_coefficients,
            "zero_coefficients": zero_coefficients,
            "negative_coefficients": 0,
            "smallest_positive_scaled": _flint_fraction_record(
                smallest_positive
            ),
            "largest_scaled": _flint_fraction_record(largest),
            "sha256_leaf_sign_summaries": leaf_digest.hexdigest(),
        },
    }


def _unit_box(count: int):
    return tuple((Fraction(0), Fraction(1)) for _ in range(count))


def _tree_boxes(leaves, order, base_box):
    return tuple(
        (path, closure, multiplier, base_box)
        for path, closure, multiplier in leaves
    )


def _fixed_box(*intervals):
    return tuple((Fraction(lower), Fraction(upper)) for lower, upper in intervals)


def verify_unique_top_l4(derived: dict[str, object]) -> dict[str, object]:
    """Replay the complete exact rational-box cover for the ``L=4`` row."""
    charts = derive_unique_top_l4_charts(derived)
    unit = _unit_box(5)

    q_le_u = _verify_boxes(
        "q_le_u", charts["q_le_u"], (("", "T", None, unit),)
    )
    z_le_a = _verify_boxes(
        "z_le_a_final",
        charts["z_le_a_final"],
        _tree_boxes(Z_LE_A_FINAL_LEAVES, (4, 1, 0), unit),
        order=(4, 1, 0),
        tree=Z_LE_A_FINAL_LEAVES,
    )
    b_ge_v = _verify_boxes(
        "b_ge_v_final",
        charts["b_ge_v_final"],
        _tree_boxes(B_GE_V_FINAL_LEAVES, (2, 0, 1, 3, 4), unit),
        order=(2, 0, 1, 3, 4),
        tree=B_GE_V_FINAL_LEAVES,
    )
    middle_low = _verify_boxes(
        "middle_low",
        charts["middle_low"],
        _tree_boxes(MIDDLE_LOW_Z_LEAVES, (2, 0, 1, 3, 4), unit),
        order=(2, 0, 1, 3, 4),
        tree=MIDDLE_LOW_Z_LEAVES,
    )
    middle_high_boxes = (
        (
            "", "T", None,
            _fixed_box((0, 1), (0, Fraction(1, 2)),
                       (Fraction(1, 2), 1), (0, 1), (0, 1)),
        ),
        (
            "", "T", None,
            _fixed_box((0, 1), (Fraction(1, 2), 1),
                       (Fraction(1, 2), 1), (0, 1), (0, 1)),
        ),
    )
    middle_high = _verify_boxes(
        "middle_high", charts["middle_scale"], middle_high_boxes
    )
    lower_low = _verify_boxes(
        "lower_low",
        charts["lower_low"],
        _tree_boxes(LOWER_LOW_Z_LEAVES, (0, 3, 0, 1), unit),
        order=(0, 3, 0, 1),
        tree=LOWER_LOW_Z_LEAVES,
    )
    lower_large_boxes = (
        (
            "", "T", None,
            _fixed_box((0, 1), (0, Fraction(1, 2)), (0, 1),
                       (Fraction(1, 2), 1), (0, 1)),
        ),
        (
            "", "T", None,
            _fixed_box((0, 1), (Fraction(1, 2), 1), (0, 1),
                       (Fraction(1, 2), 1), (0, 1)),
        ),
        (
            "", "T", None,
            _fixed_box((0, 1), (Fraction(1, 2), 1), (0, 1),
                       (Fraction(1, 3), Fraction(1, 2)), (0, 1)),
        ),
        (
            "", "T", None,
            _fixed_box((0, Fraction(1, 2)), (0, Fraction(1, 2)), (0, 1),
                       (Fraction(1, 3), Fraction(1, 2)), (0, 1)),
        ),
    )
    lower_large = _verify_boxes(
        "lower_large", charts["lower_scale"], lower_large_boxes
    )
    lower_corner_base = _fixed_box(
        (Fraction(1, 2), 1),
        (0, Fraction(1, 2)),
        (0, 1),
        (Fraction(1, 3), Fraction(1, 2)),
        (0, 1),
    )
    lower_corner = _verify_boxes(
        "lower_mid_corner",
        charts["lower_scale"],
        _tree_boxes(
            LOWER_MID_CORNER_LEAVES,
            (3, 0, 1, 2, 4),
            lower_corner_base,
        ),
        order=(3, 0, 1, 2, 4),
        tree=LOWER_MID_CORNER_LEAVES,
    )
    return {
        "charts": charts,
        "q_le_u": q_le_u,
        "z_le_a": z_le_a,
        "b_ge_v": b_ge_v,
        "middle_low": middle_low,
        "middle_high": middle_high,
        "lower_low": lower_low,
        "lower_large": lower_large,
        "lower_mid_corner": lower_corner,
    }


def _chart_monomials(chart: dict[str, object]) -> dict[str, object]:
    return {
        "target": list(chart["target_monomial"]),
        "gaps": [list(value) for value in chart["gap_monomials"]],
    }


def unique_top_l4_certificate(
    derived: dict[str, object], verified: dict[str, object]
) -> dict[str, object]:
    charts = verified["charts"]
    confinement_records = {
        name: {
            "polynomial": _polynomial_record(value["polynomial"]),
            "bernstein_signs": value["statistics"],
        }
        for name, value in charts["confinements"].items()
    }
    return {
        "schema_version": 1,
        "mathematical_status": (
            "ordinary interpolation and projective reductions with exact "
            "rational polynomial certificates"
        ),
        "scope": {
            "hypothesis": (
                'The regular uniform six-coordinate section has normalized knots (-1,-u,-v,-w,t*a,t), ordered centroid, and the reduced gap inequalities of its chart. The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative.'
            ),
            "proved_unique_top_families": ["L=4"],
            "companion_unique_top_families": ["L=2", "L=3"],
            "conclusion": "p<=I_{gamma_0}(1,5)=1-(1-gamma_0)^5",
            "global_consequence": (
                'The companion certificates address all nine local first-coordinate families. The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.'
            ),
        },
        "normalized_shape": "(-1,-u,-v,-w,t*a,t)",
        "base_variables": [str(value) for value in derived["cube_variables"]],
        "base_parameterization": "u=U, v=U*V, w=U*V*W, a=A, t=q",
        "comparison_identity": (
            "I_{gamma_0}(1,5)-p is a positive rational factor times T"
        ),
        "gap_identity": (
            "gamma_{i+1}-gamma_i is a nonnegative collision factor times "
            "G_i over a coefficientwise-positive denominator"
        ),
        "projective_cover": (
            "Split q<=U and U<=q.  In the latter write U=q*Z and split "
            "Z<=A from A<=Z.  The first branch is resolved by nested scale "
            "bounds; the second by comparing B=A/Z with V and V*W."
        ),
        "scale_bounds": {
            "global": {
                "conclusion": "orderedness forces q<3",
                "polynomial": _polynomial_record(derived["global_scale"]),
                "bernstein_signs": derived["global_scale_statistics"],
            },
            **confinement_records,
        },
        "removed_positive_monomials": {
            "initial_target": list(derived["initial_target_monomial"]),
            "cube_target": list(derived["cube_target_monomial"]),
            **{
                name: _chart_monomials(charts[name])
                for name in (
                    "q_le_u", "u_le_q", "z_le_a", "z_le_a_r",
                    "z_le_a_s", "z_le_a_c", "z_le_a_final", "a_le_z",
                    "b_ge_v", "b_ge_v_final", "middle", "middle_low",
                    "middle_scale", "lower", "lower_low", "lower_scale",
                )
            },
        },
        "polynomials": {
            "target_T": _polynomial_record(derived["target"]),
            "centroid_gaps": [_polynomial_record(gap) for gap in derived["gaps"]],
            **{
                f"{name}_target": _polynomial_record(charts[name]["target"])
                for name in (
                    "q_le_u", "z_le_a_final", "b_ge_v_final",
                    "middle_low", "middle_scale", "lower_low", "lower_scale",
                )
            },
        },
        "positive_gap_factors": list(derived["cofactor_statistics"]),
        "bernstein_covers": {
            name: verified[name]
            for name in (
                "q_le_u", "z_le_a", "b_ge_v", "middle_low",
                "middle_high", "lower_low", "lower_large",
                "lower_mid_corner",
            )
        },
        "arithmetic": {
            "symbolic_derivation": "SymPy exact rational arithmetic",
            "bernstein_signs": "FLINT exact rational multihomogenization",
            "floating_point_role": "none",
        },
    }


def build_certificate() -> dict[str, object]:
    derived = derive_unique_top_l4()
    verified = verify_unique_top_l4(derived)
    return unique_top_l4_certificate(derived, verified)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        "six-coordinate unique-top L=4 family certified exactly; "
        "companion certificates complete all nine canonical families"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
