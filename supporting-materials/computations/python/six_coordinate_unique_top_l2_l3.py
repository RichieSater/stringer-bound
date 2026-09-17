"""Certificates for the remaining six-coordinate unique-top families.

For ``L=2`` the normalized knots are

    (-1, -u, t*a, t*b, t*c, t),

and for ``L=3`` they are

    (-1, -u, -v, t*a, t*b, t).

Ordered ratios put each family on a five-dimensional cube.  A sharp rational
minorant of ``I_g(1,5)`` reduces the comparison to a much smaller polynomial:

    5*g/(1 + (10974/4651)*g) <= 1-(1-g)^5,  0 <= g <= 1/6.

Projective charts resolve the simultaneous small-scale corners.  All but one
leaf cover are checked in exact rational tensor-Bernstein arithmetic.  The
large ``L=2`` lower chart is checked by binary64 centers with exact rational
forward-error bounds; no unbounded floating-point sign decision is trusted.

The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative.
The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp

from certified_binary64_bernstein import verify_binary64_bernstein_tree
from six_coordinate_face32 import _confluent_divided_difference
from six_coordinate_unique_top import (
    _bernstein_statistics,
    _chart_family,
    _chart_monomials,
    _coefficientwise_positive,
    _polynomial_record,
    _poly,
    _remove_positive_monomial,
    _tree_boxes,
    _unit_box,
    _verify_boxes,
)
from six_coordinate_unique_top_l2_l3_witness import (
    L2_C_LE_Z,
    L2_Q_LE_U,
    L2_Z_LE_C_LOWER,
    L2_Z_LE_C_MIDDLE,
    L2_Z_LE_C_UPPER,
    L3_B_LE_Z_LOWER,
    L3_B_LE_Z_MIDDLE,
    L3_B_LE_Z_UPPER,
    L3_Q_LE_U,
    L3_Z_LE_B_LOWER,
    L3_Z_LE_B_MIDDLE,
    L3_Z_LE_B_UPPER,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent
    / "certificates"
    / "six-coordinate-unique-top-l2-l3-certificate.json"
)
MINORANT_CONSTANT = sp.Rational(10974, 4651)


def _positive_constant(poly: sp.Poly) -> None:
    if poly.total_degree() != 0 or poly.LC() <= 0:
        raise AssertionError("expected a positive constant common factor")


def derive_beta_minorant() -> dict[str, object]:
    """Check the sharp linear-fractional minorant on ``[0,1/6]``."""
    g = sp.symbols("g")
    beta_tail = 1 - (1 - g) ** 5
    required_constant = sp.cancel((5 * g / beta_tail - 1) / g)
    expected = -(
        g**3 - 5 * g**2 + 10 * g - 10
    ) / (
        g**4 - 5 * g**3 + 10 * g**2 - 10 * g + 5
    )
    if sp.cancel(required_constant - expected) != 0:
        raise AssertionError("minorant constant function changed")
    derivative = sp.factor(sp.diff(required_constant, g))
    expected_derivative = (
        (g**2 - 5 * g + 5)
        * (g**4 - 5 * g**3 + 15 * g**2 - 20 * g + 10)
        / (g**4 - 5 * g**3 + 10 * g**2 - 10 * g + 5) ** 2
    )
    if sp.cancel(derivative - expected_derivative) != 0:
        raise AssertionError("minorant derivative factorization changed")
    endpoint = sp.cancel(required_constant.subs(g, sp.Rational(1, 6)))
    if endpoint != MINORANT_CONSTANT:
        raise AssertionError("minorant endpoint constant changed")
    return {
        "interval": ["0", "1/6"],
        "beta_tail": "1-(1-g)^5",
        "minorant": "5*g/(1+(10974/4651)*g)",
        "sharp_constant": {
            "numerator": "10974",
            "denominator": "4651",
        },
        "required_constant_function": sp.sstr(required_constant),
        "derivative_factorization": sp.sstr(derivative),
        "endpoint_equality": "g=1/6",
    }


def _family_data(level: int):
    if level == 2:
        u, a, b, c, t = sp.symbols("u a b c t")
        natural_variables = (u, a, b, c, t)
        knots = (-1, -u, t * a, t * b, t * c, t)
        positive_indices = range(2, 6)
        U, A, B, C, q = sp.symbols("U A B C q")
        cube_variables = (U, A, B, C, q)
        ordered_substitution = {
            u: U,
            a: A * B * C,
            b: B * C,
            c: C,
            t: q,
        }
        gap_cofactors = (
            q * (1 - U),
            sp.Integer(1),
            B * C * (1 - A),
            C * (1 - B),
            1 - C,
        )
        expected_target = (
            (7, 4, 4, 4, 14),
            (7, 4, 8, 12, 14),
            2_779,
        )
        expected_gaps = (
            ((6, 3, 6, 9, 10), 572),
            ((7, 3, 5, 7, 8), 218),
            ((5, 2, 4, 5, 6), 74),
            ((5, 1, 3, 5, 6), 74),
            ((5, 1, 2, 4, 6), 74),
        )
    elif level == 3:
        u, v, a, b, t = sp.symbols("u v a b t")
        natural_variables = (u, v, a, b, t)
        knots = (-1, -u, -v, t * a, t * b, t)
        positive_indices = range(3, 6)
        U, V, A, B, q = sp.symbols("U V A B q")
        cube_variables = (U, V, A, B, q)
        ordered_substitution = {
            u: U,
            v: U * V,
            a: A * B,
            b: B,
            t: q,
        }
        gap_cofactors = (
            q * (1 - U),
            U * q * (1 - V),
            sp.Integer(1),
            B * (1 - A),
            1 - B,
        )
        expected_target = (
            (5, 5, 6, 6, 15),
            (10, 5, 6, 12, 15),
            4_687,
        )
        expected_gaps = (
            ((6, 2, 4, 8, 9), 394),
            ((8, 4, 4, 8, 9), 394),
            ((8, 5, 5, 8, 9), 392),
            ((8, 4, 4, 8, 9), 394),
            ((8, 4, 2, 6, 9), 394),
        )
    else:  # pragma: no cover - internal caller fixes the two levels
        raise ValueError(level)
    return {
        "natural_variables": natural_variables,
        "knots": knots,
        "positive_indices": positive_indices,
        "cube_variables": cube_variables,
        "ordered_substitution": ordered_substitution,
        "gap_cofactors": tuple(
            _poly(value, cube_variables) for value in gap_cofactors
        ),
        "expected_target": expected_target,
        "expected_gaps": expected_gaps,
    }


def derive_unique_top_minorant(level: int) -> dict[str, object]:
    """Derive one cap comparison and its adjacent centroid gaps."""
    data = _family_data(level)
    natural_variables = data["natural_variables"]
    cube_variables = data["cube_variables"]
    abstract_knots = sp.symbols("Y0:6")
    specialization = dict(
        zip(abstract_knots, data["knots"], strict=True)
    )
    cap_abstract = _confluent_divided_difference(
        abstract_knots,
        (1,) * 6,
        data["positive_indices"],
        5,
    )
    density_abstract = 5 * _confluent_divided_difference(
        abstract_knots,
        (1,) * 6,
        data["positive_indices"],
        4,
    )
    cap = sp.cancel(cap_abstract.subs(specialization))
    density = sp.cancel(density_abstract.subs(specialization))
    centroids = tuple(
        sp.cancel(
            sp.diff(cap_abstract, knot).subs(specialization) / density
        )
        for knot in abstract_knots
    )

    gamma_numerator, gamma_denominator = (
        _poly(value, natural_variables)
        for value in sp.fraction(centroids[0])
    )
    cap_numerator, cap_denominator = (
        _poly(value, natural_variables) for value in sp.fraction(cap)
    )
    raw_target = (
        5 * gamma_numerator * cap_denominator
        - cap_numerator * gamma_denominator
        - MINORANT_CONSTANT * cap_numerator * gamma_numerator
    )
    comparison_denominator = (
        cap_denominator
        * (gamma_denominator + MINORANT_CONSTANT * gamma_numerator)
    )
    common = sp.gcd(raw_target, comparison_denominator)
    _positive_constant(common)
    reduced_target = sp.exquo(raw_target, common)
    reduced_target, initial_target_monomial = _remove_positive_monomial(
        reduced_target
    )
    target_content, reduced_target = reduced_target.primitive()
    if target_content != sp.Rational(1, 4651):
        raise AssertionError("minorant target content changed")
    cube_target, cube_target_monomial = _remove_positive_monomial(
        _poly(
            reduced_target.as_expr().subs(data["ordered_substitution"]),
            cube_variables,
        )
    )
    natural_degree, cube_degree, term_count = data["expected_target"]
    if (
        tuple(map(int, reduced_target.degree_list())) != natural_degree
        or tuple(map(int, cube_target.degree_list())) != cube_degree
        or len(reduced_target.terms()) != term_count
        or len(cube_target.terms()) != term_count
        or initial_target_monomial != (0, 0, 0, 0, 1)
        or cube_target_monomial != (0, 0, 0, 0, 0)
    ):
        raise AssertionError(f"L={level} minorant target changed")

    # These coefficientwise-positive cube denominators fix the orientations
    # of both the comparison identity and the five adjacent-gap identities.
    for denominator in (
        cap_denominator,
        gamma_denominator + MINORANT_CONSTANT * gamma_numerator,
    ):
        _coefficientwise_positive(
            _poly(
                denominator.as_expr().subs(data["ordered_substitution"]),
                cube_variables,
            )
        )

    gaps = []
    gap_denominators = []
    cofactor_statistics = []
    for index, cofactor in enumerate(data["gap_cofactors"]):
        gap = sp.cancel(centroids[index + 1] - centroids[index])
        numerator, denominator = sp.fraction(gap)
        numerator_poly = _poly(
            numerator.subs(data["ordered_substitution"]), cube_variables
        )
        denominator_poly = _poly(
            denominator.subs(data["ordered_substitution"]), cube_variables
        )
        _coefficientwise_positive(denominator_poly)
        quotient = sp.exquo(numerator_poly, cofactor)
        gap_content, gap_poly = quotient.primitive()
        if gap_content <= 0:
            raise AssertionError("adjacent-gap orientation changed")
        if numerator_poly != gap_content * cofactor * gap_poly:
            raise AssertionError("adjacent-gap factorization changed")
        if (
            tuple(map(int, gap_poly.degree_list())),
            len(gap_poly.terms()),
        ) != data["expected_gaps"][index]:
            raise AssertionError(f"L={level} gap {index} changed")
        gaps.append(gap_poly)
        gap_denominators.append(denominator_poly)
        cofactor_statistics.append(_bernstein_statistics(cofactor))

    return {
        "level": level,
        "natural_variables": natural_variables,
        "cube_variables": cube_variables,
        "cap": cap,
        "gamma0": centroids[0],
        "target": cube_target,
        "gaps": tuple(gaps),
        "gap_denominators": tuple(gap_denominators),
        "gap_cofactors": data["gap_cofactors"],
        "cofactor_statistics": tuple(cofactor_statistics),
        "initial_target_monomial": initial_target_monomial,
        "cube_target_monomial": cube_target_monomial,
    }


def _monomial_transform(poly, variables, exponent_map):
    coefficients = defaultdict(lambda: sp.S.Zero)
    for powers, coefficient in poly.terms():
        coefficients[exponent_map(powers)] += coefficient
    coefficients = {
        powers: coefficient
        for powers, coefficient in coefficients.items()
        if coefficient
    }
    minima = tuple(
        min(powers[index] for powers in coefficients)
        for index in range(len(variables))
    )
    reduced = {
        tuple(
            powers[index] - minima[index]
            for index in range(len(variables))
        ): coefficient
        for powers, coefficient in coefficients.items()
    }
    return sp.Poly.from_dict(reduced, *variables, domain=sp.QQ), minima


def _monomial_chart_family(chart, variables, exponent_map):
    target, target_monomial = _monomial_transform(
        chart["target"], variables, exponent_map
    )
    gap_pairs = tuple(
        _monomial_transform(gap, variables, exponent_map)
        for gap in chart["gaps"]
    )
    return {
        "variables": variables,
        "target": target,
        "gaps": tuple(pair[0] for pair in gap_pairs),
        "target_monomial": target_monomial,
        "gap_monomials": tuple(pair[1] for pair in gap_pairs),
    }


def _confinement_record(
    gap: sp.Poly,
    substitution: dict[sp.Symbol, sp.Expr],
    variables: tuple[sp.Symbol, ...],
    expected: tuple[tuple[int, ...], int, int, int],
) -> dict[str, object]:
    polynomial = _poly((-gap.as_expr()).subs(substitution), variables)
    statistics = _bernstein_statistics(polynomial, axes=(0, 1, 2, 3))
    actual = (
        tuple(map(int, polynomial.degree_list())),
        len(polynomial.terms()),
        int(statistics["positive_coefficients"]),
        int(statistics["zero_coefficients"]),
    )
    if actual != expected or statistics["negative_coefficients"] != 0:
        raise AssertionError("projective scale confinement changed")
    return {"polynomial": polynomial, "statistics": statistics}


def derive_l2_charts(derived: dict[str, object]) -> dict[str, object]:
    U, A, B, C, q = derived["cube_variables"]
    target, gaps = derived["target"], derived["gaps"]
    Z, D, R, E, F, S, rho = sp.symbols("Z D R E F S rho")

    q_le_u = _chart_family(
        target, gaps, {q: U * Z}, (U, A, B, C, Z)
    )
    u_le_q = _chart_family(
        target, gaps, {U: q * Z}, (q, A, B, C, Z)
    )
    z_le_c = _chart_family(
        u_le_q["target"], u_le_q["gaps"],
        {Z: C * D}, (q, A, B, C, D),
    )
    c_le_z = _chart_family(
        u_le_q["target"], u_le_q["gaps"],
        {C: Z * D}, (q, A, B, Z, D),
    )
    z_le_c_final = _chart_family(
        z_le_c["target"], z_le_c["gaps"],
        {q: C * R / 4}, (R, A, B, C, D),
    )
    c_le_z_final = _chart_family(
        c_le_z["target"], c_le_z["gaps"],
        {q: Z * R}, (R, A, B, Z, D),
    )

    lower = _monomial_chart_family(
        z_le_c_final,
        (R, A, B, C, E),
        lambda m: (m[0], m[1] + m[4], m[2] + m[4], m[3], m[4]),
    )
    middle = _monomial_chart_family(
        z_le_c_final,
        (R, F, B, C, E),
        lambda m: (m[0], m[1], m[2] + m[4], m[3], m[1] + m[4]),
    )
    upper = _monomial_chart_family(
        z_le_c_final,
        (R, A, D, C, E),
        lambda m: (m[0], m[1], m[2] + m[4], m[3], m[2]),
    )
    lower_scaled = _chart_family(
        lower["target"], lower["gaps"],
        {R: 4 * A**4 * B**3 * E**5 * S},
        (S, A, B, C, E),
    )
    middle_scaled = _chart_family(
        middle["target"], middle["gaps"],
        {R: 4 * B**3 * E**4 * S},
        (S, F, B, C, E),
    )
    upper_scaled = _chart_family(
        upper["target"], upper["gaps"],
        {R: 4 * D**3 * S},
        (S, A, D, C, E),
    )

    expected_targets = {
        "q_le_u": ((10, 4, 8, 12, 14), 2_779),
        "c_le_z_final": ((10, 4, 8, 19, 12), 2_779),
        "lower_scaled": ((10, 45, 37, 19, 53), 2_779),
        "middle_scaled": ((10, 4, 37, 19, 45), 2_779),
        "upper_scaled": ((10, 4, 37, 19, 8), 2_779),
    }
    charts = {
        "q_le_u": q_le_u,
        "u_le_q": u_le_q,
        "z_le_c": z_le_c,
        "c_le_z": c_le_z,
        "z_le_c_final": z_le_c_final,
        "c_le_z_final": c_le_z_final,
        "lower": lower,
        "middle": middle,
        "upper": upper,
        "lower_scaled": lower_scaled,
        "middle_scaled": middle_scaled,
        "upper_scaled": upper_scaled,
    }
    for name, expected in expected_targets.items():
        polynomial = charts[name]["target"]
        if (
            tuple(map(int, polynomial.degree_list())),
            len(polynomial.terms()),
        ) != expected:
            raise AssertionError(f"L=2 chart {name} changed")

    confinements = {
        "global": _confinement_record(
            gaps[4], {q: 1 + rho}, (U, A, B, C, rho),
            ((5, 1, 2, 4, 6), 195, 554, 706),
        ),
        "z_le_c_q": _confinement_record(
            z_le_c["gaps"][4], {q: C / 4 + rho},
            (A, B, C, D, rho), ((1, 2, 9, 5, 5), 185, 967, 1_193),
        ),
        "c_le_z_q": _confinement_record(
            c_le_z["gaps"][4], {q: Z + rho},
            (A, B, Z, D, rho), ((1, 2, 9, 4, 5), 181, 938, 862),
        ),
        "lower_r": _confinement_record(
            lower["gaps"][4], {R: 4 * A**4 * B**3 * E**5 + rho},
            (A, B, C, E, rho),
            ((24, 19, 8, 29, 5), 213, 401_670, 408_330),
        ),
        "middle_r": _confinement_record(
            middle["gaps"][4], {R: 4 * B**3 * E**4 + rho},
            (F, B, C, E, rho),
            ((1, 19, 8, 24, 5), 213, 28_377, 25_623),
        ),
        "upper_r": _confinement_record(
            upper["gaps"][4], {R: 4 * D**3 + rho},
            (A, D, C, E, rho),
            ((1, 19, 8, 2, 5), 211, 3_538, 2_942),
        ),
    }
    return {**charts, "confinements": confinements}


def derive_l3_charts(derived: dict[str, object]) -> dict[str, object]:
    U, V, A, B, q = derived["cube_variables"]
    target, gaps = derived["target"], derived["gaps"]
    Z, D, R, E, F, S, z, rho = sp.symbols("Z D R E F S z rho")

    q_le_u = _chart_family(
        target, gaps, {q: U * Z}, (U, V, A, B, Z)
    )
    u_le_q = _chart_family(
        target, gaps, {U: q * Z}, (q, V, A, B, Z)
    )
    z_le_b = _chart_family(
        u_le_q["target"], u_le_q["gaps"],
        {Z: B * D}, (q, V, A, B, D),
    )
    b_le_z = _chart_family(
        u_le_q["target"], u_le_q["gaps"],
        {B: Z * D}, (q, V, A, Z, D),
    )
    z_le_b_final = _chart_family(
        z_le_b["target"], z_le_b["gaps"],
        {q: B * R}, (R, V, A, B, D),
    )
    b_le_z_final = _chart_family(
        b_le_z["target"], b_le_z["gaps"],
        {q: 4 * Z * R}, (R, V, A, Z, D),
    )

    z_upper = _monomial_chart_family(
        z_le_b_final,
        (R, V, A, B, E),
        lambda m: (m[0], m[1], m[2] + m[4], m[3], m[4]),
    )
    z_middle = _monomial_chart_family(
        z_le_b_final,
        (R, F, D, B, E),
        lambda m: (m[0], m[1], m[2] + m[4], m[3], m[1] + m[2]),
    )
    z_lower = _monomial_chart_family(
        z_le_b_final,
        (R, V, D, B, E),
        lambda m: (m[0], m[1] + m[2], m[2] + m[4], m[3], m[2]),
    )
    b_upper = _monomial_chart_family(
        b_le_z_final,
        (R, A, D, Z, E),
        lambda m: (m[0], m[2] + m[1], m[4] + m[1], m[3], m[1]),
    )
    b_middle = _monomial_chart_family(
        b_le_z_final,
        (R, F, D, Z, E),
        lambda m: (m[0], m[2], m[1] + m[4], m[3], m[1] + m[2]),
    )
    b_lower = _monomial_chart_family(
        b_le_z_final,
        (R, V, A, Z, E),
        lambda m: (m[0], m[1] + m[4], m[2], m[3], m[4]),
    )
    z_upper_scaled = _chart_family(
        z_upper["target"], z_upper["gaps"],
        {R: 8 * A**3 * E**4 * S}, (S, V, A, B, E),
    )
    z_middle_scaled = _chart_family(
        z_middle["target"], z_middle["gaps"],
        {R: 4 * D**3 * S}, (S, F, D, B, E),
    )
    z_lower_scaled = _chart_family(
        z_lower["target"], z_lower["gaps"],
        {R: 4 * D**3 * S}, (S, V, D, B, E),
    )

    expected_targets = {
        "q_le_u": ((7, 5, 6, 12, 15), 4_687),
        "z_upper_scaled": ((7, 5, 29, 15, 35), 4_687),
        "z_middle_scaled": ((7, 5, 29, 15, 8), 4_687),
        "z_lower_scaled": ((7, 8, 29, 15, 6), 4_687),
        "b_upper": ((7, 8, 11, 15, 5), 4_687),
        "b_middle": ((7, 6, 11, 15, 8), 4_687),
        "b_lower": ((7, 11, 6, 15, 12), 4_687),
    }
    charts = {
        "q_le_u": q_le_u,
        "u_le_q": u_le_q,
        "z_le_b": z_le_b,
        "b_le_z": b_le_z,
        "z_le_b_final": z_le_b_final,
        "b_le_z_final": b_le_z_final,
        "z_upper": z_upper,
        "z_middle": z_middle,
        "z_lower": z_lower,
        "b_upper": b_upper,
        "b_middle": b_middle,
        "b_lower": b_lower,
        "z_upper_scaled": z_upper_scaled,
        "z_middle_scaled": z_middle_scaled,
        "z_lower_scaled": z_lower_scaled,
    }
    for name, expected in expected_targets.items():
        polynomial = charts[name]["target"]
        if (
            tuple(map(int, polynomial.degree_list())),
            len(polynomial.terms()),
        ) != expected:
            raise AssertionError(f"L=3 chart {name} changed")

    confinements = {
        "global": _confinement_record(
            gaps[4], {q: 2 + rho}, (U, V, A, B, rho),
            ((8, 4, 2, 6, 9), 1_125, 3_959, 5_491),
        ),
        "z_le_b_q": _confinement_record(
            z_le_b["gaps"][4], {q: B + rho},
            (V, A, B, D, rho), ((4, 2, 9, 8, 4), 791, 4_068, 2_682),
        ),
        "b_le_z_low_q": _confinement_record(
            b_le_z["gaps"][4], {Z: z / 2, q: 2 * z + rho},
            (V, A, z, D, rho), ((4, 2, 9, 6, 4), 748, 2_824, 2_426),
        ),
        "z_upper_r": _confinement_record(
            z_upper["gaps"][4], {R: 8 * A**3 * E**4 + rho},
            (V, A, B, E, rho),
            ((4, 19, 8, 24, 4), 1_024, 70_485, 42_015),
        ),
        "z_middle_r": _confinement_record(
            z_middle["gaps"][4], {R: 4 * D**3 + rho},
            (F, D, B, E, rho),
            ((4, 19, 8, 4, 4), 983, 14_920, 7_580),
        ),
        "z_lower_r": _confinement_record(
            z_lower["gaps"][4], {R: 4 * D**3 + rho},
            (V, D, B, E, rho),
            ((4, 19, 8, 2, 4), 983, 8_739, 4_761),
        ),
    }
    return {**charts, "confinements": confinements}


def verify_l2(derived: dict[str, object]) -> dict[str, object]:
    charts = derive_l2_charts(derived)
    unit = _unit_box(5)
    exact = {
        "q_le_u": _verify_boxes(
            "l2_q_le_u", charts["q_le_u"],
            _tree_boxes(L2_Q_LE_U, (4, 0, 1, 2, 3), unit),
            order=(4, 0, 1, 2, 3), tree=L2_Q_LE_U,
        ),
        "c_le_z": _verify_boxes(
            "l2_c_le_z", charts["c_le_z_final"],
            _tree_boxes(L2_C_LE_Z, (4, 0, 1, 2, 3), unit),
            order=(4, 0, 1, 2, 3), tree=L2_C_LE_Z,
        ),
        "z_le_c_middle": _verify_boxes(
            "l2_z_le_c_middle", charts["middle_scaled"],
            _tree_boxes(L2_Z_LE_C_MIDDLE, (2, 0, 4, 3, 1), unit),
            order=(2, 0, 4, 3, 1), tree=L2_Z_LE_C_MIDDLE,
        ),
        "z_le_c_upper": _verify_boxes(
            "l2_z_le_c_upper", charts["upper_scaled"],
            _tree_boxes(L2_Z_LE_C_UPPER, (0, 2, 0, 3, 4), unit),
            order=(0, 2, 0, 3, 4), tree=L2_Z_LE_C_UPPER,
        ),
    }
    lower = verify_binary64_bernstein_tree(
        charts["lower_scaled"]["target"],
        charts["lower_scaled"]["gaps"][4],
        L2_Z_LE_C_LOWER,
        (2, 0, 4, 3, 1),
    )
    return {"charts": charts, "exact_covers": exact, "lower_cover": lower}


def verify_l3(derived: dict[str, object]) -> dict[str, object]:
    charts = derive_l3_charts(derived)
    unit = _unit_box(5)
    specifications = {
        "q_le_u": ("q_le_u", L3_Q_LE_U, (4, 0, 1, 2, 3)),
        "z_le_b_upper": (
            "z_upper_scaled", L3_Z_LE_B_UPPER, (2, 0, 4, 3, 1)
        ),
        "z_le_b_middle": (
            "z_middle_scaled", L3_Z_LE_B_MIDDLE, (2, 0, 4, 3, 1)
        ),
        "z_le_b_lower": (
            "z_lower_scaled", L3_Z_LE_B_LOWER, (0, 1, 2, 3, 4)
        ),
        "b_le_z_upper": (
            "b_upper", L3_B_LE_Z_UPPER, (0, 1, 2, 3, 4)
        ),
        "b_le_z_middle": (
            "b_middle", L3_B_LE_Z_MIDDLE, (0, 1, 2, 3, 4)
        ),
        "b_le_z_lower": (
            "b_lower", L3_B_LE_Z_LOWER, (0, 1, 2, 3, 4)
        ),
    }
    covers = {}
    for name, (chart_name, leaves, order) in specifications.items():
        covers[name] = _verify_boxes(
            f"l3_{name}", charts[chart_name],
            _tree_boxes(leaves, order, unit),
            order=order, tree=leaves,
        )
    return {"charts": charts, "exact_covers": covers}


def _family_certificate(derived, verified):
    charts = verified["charts"]
    selected = (
        ("q_le_u", "c_le_z_final", "lower_scaled", "middle_scaled", "upper_scaled")
        if derived["level"] == 2
        else (
            "q_le_u", "z_upper_scaled", "z_middle_scaled", "z_lower_scaled",
            "b_upper", "b_middle", "b_lower",
        )
    )
    return {
        "base_variables": [str(value) for value in derived["cube_variables"]],
        "target": _polynomial_record(derived["target"]),
        "centroid_gaps": [
            _polynomial_record(gap) for gap in derived["gaps"]
        ],
        "positive_gap_factors": list(derived["cofactor_statistics"]),
        "selected_chart_targets": {
            name: _polynomial_record(charts[name]["target"])
            for name in selected
        },
        "selected_chart_monomials": {
            name: _chart_monomials(charts[name]) for name in selected
        },
        "scale_bounds": {
            name: {
                "polynomial": _polynomial_record(value["polynomial"]),
                "bernstein_signs": value["statistics"],
            }
            for name, value in charts["confinements"].items()
        },
        "exact_rational_covers": verified["exact_covers"],
        **(
            {"certified_binary64_cover": verified["lower_cover"]}
            if derived["level"] == 2
            else {}
        ),
    }


def build_certificate() -> dict[str, object]:
    minorant = derive_beta_minorant()
    l2 = derive_unique_top_minorant(2)
    l3 = derive_unique_top_minorant(3)
    l2_verified = verify_l2(l2)
    l3_verified = verify_l3(l3)
    return {
        "schema_version": 1,
        "mathematical_status": (
            "ordinary interpolation and projective reductions; exact rational "
            "Bernstein certificates except for one rigorously enclosed "
            "binary64 Bernstein cover"
        ),
        "scope": {
            "proved_unique_top_families": ["L=2", "L=3"],
            "companion_exact_family": "L=4",
            "conclusion": "p<=I_{gamma_0}(1,5)=1-(1-gamma_0)^5",
            "global_consequence": (
                'The companion certificates address all nine local first-coordinate families. The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative. The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.'
            ),
        },
        "comparison": {
            "stronger_conclusion": (
                "p<=5*gamma_0/(1+(10974/4651)*gamma_0)"
            ),
            "minorant": minorant,
            "reason_interval_applies": (
                "ordered barycentric centroid coordinates sum to one, hence "
                "0<=gamma_0<=1/6"
            ),
        },
        "families": {
            "L=2": _family_certificate(l2, l2_verified),
            "L=3": _family_certificate(l3, l3_verified),
        },
        "arithmetic": {
            "symbolic_derivation": "SymPy exact rational arithmetic",
            "exact_bernstein_signs": (
                "FLINT exact rational multihomogenization"
            ),
            "large_l2_lower_chart": (
                "NumPy binary64 tensor centers with exact rational forward-error "
                "bounds and outward-rounded dyadic enclosures"
            ),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        "six-coordinate unique-top L=2 and L=3 families certified; "
        "all nine canonical non-elementary families are proved"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
