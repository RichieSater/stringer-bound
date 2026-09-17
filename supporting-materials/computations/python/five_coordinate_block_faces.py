"""Exact first-coordinate certificates for three five-coordinate knot strata.

The normalized repeated-top shapes are (-1,-u,-u*s,t,t),
(-1,-u,t*v,t,t), and (-1,-u,t,t,t). The encoded reduced centroid-gap
constraints imply p <= I_{gamma_0}(1,4) on their relative interiors.
Exact order polynomials bound the scale by 3/2, 1, and 2/3, respectively.
Division and sparse rational systems prove the required Bernstein signs.

Limits are covered only when the reduced gap constraints persist. An
ordered centroid at a collision does not justify canceling a zero factor.
These first-coordinate inequalities do not prove the active-prefix bound
needed for the unresolved global monotone-weight cap theorem.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import prod
from pathlib import Path

import numpy as np
import sympy as sp
from flint import fmpq_mat

from five_coordinate_block_witness import (
    THREE_BELOW_TOP_PAIR_ACTIVE_ROWS_FLAT,
    THREE_BELOW_TOP_PAIR_MULTIPLIER_DEGREE,
    THREE_BELOW_TOP_PAIR_SUPPORT_FLAT,
    THREE_BELOW_TOP_PAIR_TARGET_DEGREE,
    TWO_BELOW_TOP_TRIPLE_ACTIVE_ROWS_FLAT,
    TWO_BELOW_TOP_TRIPLE_MULTIPLIER_DEGREE,
    TWO_BELOW_TOP_TRIPLE_SUPPORT_FLAT,
    TWO_BELOW_TOP_TRIPLE_TARGET_DEGREE,
)
from five_coordinate_mixed_section import (
    _fraction,
    _fraction_record,
    _from_fmpq,
    _index_digest,
    _poly,
    _power_digest,
    _product_entry,
    _to_fmpq,
    exact_bernstein_coefficients,
    exact_bernstein_statistics,
    exact_coefficient_statistics,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent
    / "certificates"
    / "five-coordinate-block-faces-certificate.json"
)


def _positive_polynomial(poly: sp.Poly) -> None:
    if not poly.terms() or any(coefficient <= 0 for _, coefficient in poly.terms()):
        raise AssertionError("expected a coefficientwise-positive polynomial")


def _nontrivial_factor(numerator, variables) -> sp.Poly:
    """Remove the positive monomial scale factor from a factored numerator."""
    unit, factors = sp.factor_list(numerator)
    candidates = []
    for factor, exponent in factors:
        polynomial = _poly(factor, variables)
        if len(polynomial.terms()) > 1:
            candidates.append((polynomial, exponent))
    if len(candidates) != 1 or candidates[0][1] != 1:
        raise AssertionError("comparison numerator factorization changed")
    return unit * candidates[0][0]


def _cap_data(positive_indices, substitution):
    z = sp.symbols("z0:5")
    cap = sum(
        z[index] ** 4
        / sp.prod(z[index] - z[j] for j in range(5) if j != index)
        for index in positive_indices
    )
    density = sum(sp.diff(cap, coordinate) for coordinate in z)
    specialized_cap = sp.factor(cap.subs(substitution))
    centroid = tuple(
        sp.factor(
            sp.diff(cap, coordinate).subs(substitution)
            / density.subs(substitution)
        )
        for coordinate in z
    )
    return specialized_cap, centroid


def _take_limit(cap, centroid, variable, value):
    return (
        sp.factor(sp.limit(cap, variable, value)),
        tuple(sp.factor(sp.limit(entry, variable, value)) for entry in centroid),
    )


def _derive_case(
    cap,
    centroid,
    substitution,
    variables,
    upper_index,
    lower_index,
    order_kind,
):
    cap = sp.factor(cap.subs(substitution))
    centroid = tuple(sp.factor(entry.subs(substitution)) for entry in centroid)
    if sp.factor(sum(centroid) - 1) != 0:
        raise AssertionError("section centroid no longer sums to one")

    margin = sp.factor(sp.cancel(1 - (1 - centroid[0]) ** 4 - cap))
    margin_numerator, margin_denominator = map(
        sp.factor, sp.fraction(margin)
    )
    target = _nontrivial_factor(margin_numerator, variables)
    if sp.factor(margin_numerator - variables[-1] * target.as_expr()) != 0:
        raise AssertionError("first-step margin factorization changed")
    _positive_polynomial(_poly(margin_denominator, variables))

    order_gap = sp.factor(
        sp.cancel(centroid[upper_index] - centroid[lower_index])
    )
    order_numerator, order_denominator = map(
        sp.factor, sp.fraction(order_gap)
    )
    _positive_polynomial(_poly(order_denominator, variables))
    if order_kind == "direct":
        order_factor = _poly(order_numerator, variables)
        order_identity = "adjacent centroid gap=G/Delta_order"
    elif order_kind == "v_minus_one":
        v = variables[1]
        quotient = sp.exquo(
            _poly(order_numerator, variables), _poly(v - 1, variables)
        )
        order_factor = -quotient
        order_identity = "adjacent centroid gap=(v-1)*R/Delta_order; G=-R"
    else:
        raise AssertionError("unknown order factorization")

    return {
        "variables": variables,
        "cap": cap,
        "centroid": centroid,
        "margin_denominator": margin_denominator,
        "target": target,
        "order_denominator": order_denominator,
        "order_factor": order_factor,
        "order_identity": order_identity,
    }


def derive_block_face_polynomials():
    """Derive the three repeated-top comparison problems from the cap formula."""
    u, v, w, t, s, tau, rho = sp.symbols("u v w t s tau rho")
    z = sp.symbols("z0:5")

    cap32, centroid32 = _cap_data(
        (3, 4),
        {z[0]: -1, z[1]: -u, z[2]: -v, z[3]: t * w, z[4]: t},
    )
    cap32, centroid32 = _take_limit(cap32, centroid32, w, 1)
    three_below_pair = _derive_case(
        cap32,
        centroid32,
        {v: u * s},
        (u, s, t),
        upper_index=3,
        lower_index=2,
        order_kind="direct",
    )
    if (
        tuple(three_below_pair["target"].degree_list()) != (12, 6, 15)
        or len(three_below_pair["target"].terms()) != 369
        or tuple(three_below_pair["order_factor"].degree_list()) != (4, 2, 5)
        or len(three_below_pair["order_factor"].terms()) != 20
    ):
        raise AssertionError("three-below/top-pair polynomials changed")

    cap23, centroid23 = _cap_data(
        (2, 3, 4),
        {z[0]: -1, z[1]: -u, z[2]: t * v, z[3]: t * w, z[4]: t},
    )
    cap23, centroid23 = _take_limit(cap23, centroid23, w, 1)
    two_below_pair = _derive_case(
        cap23,
        centroid23,
        {},
        (u, v, t),
        upper_index=3,
        lower_index=2,
        order_kind="v_minus_one",
    )
    if (
        tuple(two_below_pair["target"].degree_list()) != (11, 9, 18)
        or len(two_below_pair["target"].terms()) != 910
        or tuple(two_below_pair["order_factor"].degree_list()) != (4, 2, 5)
        or len(two_below_pair["order_factor"].terms()) != 20
    ):
        raise AssertionError("two-below/top-pair polynomials changed")

    cap23_triple, centroid23_triple = _take_limit(
        cap23, centroid23, v, 1
    )
    two_below_triple = _derive_case(
        cap23_triple,
        centroid23_triple,
        {},
        (u, t),
        upper_index=2,
        lower_index=1,
        order_kind="direct",
    )
    if (
        tuple(two_below_triple["target"].degree_list()) != (11, 14)
        or len(two_below_triple["target"].terms()) != 105
        or tuple(two_below_triple["order_factor"].degree_list()) != (3, 4)
        or len(two_below_triple["order_factor"].terms()) != 10
    ):
        raise AssertionError("two-below/top-triple polynomials changed")

    # Orderedness implies G>=0.  Exact mixed-basis expansions prove that this
    # is impossible at or beyond the displayed rational scale bounds.
    compactness = {}
    for name, case, bound in (
        ("three_below_top_pair", three_below_pair, sp.Rational(3, 2)),
        ("two_below_top_triple", two_below_triple, sp.Rational(2, 3)),
    ):
        base_variables = case["variables"][:-1]
        shifted = _poly(
            (-case["order_factor"].as_expr()).subs(t, bound + rho),
            (*base_variables, rho),
        )
        coefficients = exact_bernstein_coefficients(
            shifted, axes=tuple(range(len(base_variables)))
        )
        statistics = exact_coefficient_statistics(
            coefficients, tuple(map(int, shifted.degree_list()))
        )
        if not any(
            coefficients[(*index, 0)] > 0
            for index in np.ndindex(coefficients.shape[:-1])
        ):
            raise AssertionError("compactness boundary lost strict positivity")
        compactness[name] = {
            "bound": bound,
            "shifted_polynomial": shifted,
            "statistics": statistics,
        }

    # For the top-pair 2--3 crossing, G=-R and orderedness gives R<=0.
    # R=(t^2*u^3-u^4)+Q with Q coefficientwise positive, hence t^2<u<1.
    R = -two_below_pair["order_factor"]
    confinement = R - _poly(t**2 * u**3 - u**4, (u, v, t))
    if len(confinement.terms()) != 19:
        raise AssertionError("two-below/top-pair confinement changed")
    _positive_polynomial(confinement)
    two_below_pair["confinement"] = confinement

    # Map each exact ordered domain into a unit cube.
    scales = {
        "three_below_top_pair": sp.Rational(3, 2),
        "two_below_top_pair": sp.Rational(1),
        "two_below_top_triple": sp.Rational(2, 3),
    }
    cases = {
        "three_below_top_pair": three_below_pair,
        "two_below_top_pair": two_below_pair,
        "two_below_top_triple": two_below_triple,
    }
    for name, case in cases.items():
        cube_variables = (*case["variables"][:-1], tau)
        scale = scales[name]
        case["scale"] = scale
        case["cube_variables"] = cube_variables
        case["cube_target"] = _poly(
            case["target"].as_expr().subs(t, scale * tau), cube_variables
        )
        case["cube_order_factor"] = _poly(
            case["order_factor"].as_expr().subs(t, scale * tau),
            cube_variables,
        )

    return {
        "symbols": {"u": u, "v": v, "s": s, "t": t, "tau": tau},
        "cases": cases,
        "compactness": compactness,
    }


def _sparse_exact_multiplier(
    target: sp.Poly,
    order_factor: sp.Poly,
    multiplier_degree,
    support_flat,
    active_rows_flat,
):
    target_degree = tuple(map(int, target.degree_list()))
    order_degree = tuple(map(int, order_factor.degree_list()))
    multiplier_degree = tuple(multiplier_degree)
    if tuple(
        order_degree[index] + multiplier_degree[index]
        for index in range(len(target_degree))
    ) != target_degree:
        raise AssertionError("multiplier and order degrees do not reach target")

    target_coefficients = exact_bernstein_coefficients(target)
    order_coefficients = exact_bernstein_coefficients(order_factor)
    target_shape = tuple(value + 1 for value in target_degree)
    multiplier_shape = tuple(value + 1 for value in multiplier_degree)
    if len(support_flat) != len(active_rows_flat):
        raise AssertionError("sparse witness is not square")
    if len(set(support_flat)) != len(support_flat):
        raise AssertionError("multiplier support is not unique")
    if len(set(active_rows_flat)) != len(active_rows_flat):
        raise AssertionError("target rows are not unique")
    if not all(0 <= value < prod(multiplier_shape) for value in support_flat):
        raise AssertionError("multiplier support index is out of range")
    if not all(0 <= value < prod(target_shape) for value in active_rows_flat):
        raise AssertionError("target row index is out of range")

    multiplier_indices = tuple(
        tuple(map(int, np.unravel_index(value, multiplier_shape)))
        for value in support_flat
    )
    target_indices = tuple(
        tuple(map(int, np.unravel_index(value, target_shape)))
        for value in active_rows_flat
    )
    matrix = fmpq_mat(
        [
            [
                _to_fmpq(
                    _product_entry(
                        order_coefficients,
                        target_degree,
                        multiplier_degree,
                        target_index,
                        multiplier_index,
                    )
                )
                for multiplier_index in multiplier_indices
            ]
            for target_index in target_indices
        ]
    )
    rhs = fmpq_mat(
        [[_to_fmpq(target_coefficients[index])] for index in target_indices]
    )
    solution = matrix.solve(rhs)
    multiplier_values = tuple(
        _from_fmpq(solution[index, 0])
        for index in range(len(multiplier_indices))
    )
    if any(value <= 0 for value in multiplier_values):
        raise AssertionError("exact multiplier has a nonpositive support value")

    multiplier_coefficients = np.empty(multiplier_shape, dtype=object)
    multiplier_coefficients.fill(Fraction(0))
    for index, value in zip(
        multiplier_indices, multiplier_values, strict=True
    ):
        multiplier_coefficients[index] = value

    residual_coefficients = np.empty(target_shape, dtype=object)
    for target_index in np.ndindex(target_shape):
        residual_coefficients[target_index] = target_coefficients[
            target_index
        ] - sum(
            (
                value
                * _product_entry(
                    order_coefficients,
                    target_degree,
                    multiplier_degree,
                    target_index,
                    multiplier_index,
                )
                for multiplier_index, value in zip(
                    multiplier_indices, multiplier_values, strict=True
                )
            ),
            Fraction(0),
        )
    if any(
        residual_coefficients.flat[index] != 0
        for index in active_rows_flat
    ):
        raise AssertionError("active subsystem was not solved exactly")

    return {
        "linear_system_order": len(multiplier_indices),
        "support_digest": _index_digest(tuple(support_flat)),
        "active_rows_digest": _index_digest(tuple(active_rows_flat)),
        "multiplier_statistics": exact_coefficient_statistics(
            multiplier_coefficients, multiplier_degree
        ),
        "residual_statistics": exact_coefficient_statistics(
            residual_coefficients, target_degree
        ),
    }


def _polynomial_record(poly: sp.Poly) -> dict[str, object]:
    return {
        "degree": list(map(int, poly.degree_list())),
        "power_terms": len(poly.terms()),
        "sha256_power_expansion": _power_digest(poly),
    }


def _power_positive_record(poly: sp.Poly) -> dict[str, object]:
    smallest = min(_fraction(coefficient) for _, coefficient in poly.terms())
    return {
        **_polynomial_record(poly),
        "positive_power_coefficients": len(poly.terms()),
        "zero_stored_power_coefficients": 0,
        "negative_power_coefficients": 0,
        "smallest_positive": _fraction_record(smallest),
    }


def build_certificate() -> dict[str, object]:
    derived = derive_block_face_polynomials()
    cases = derived["cases"]

    pair23 = cases["two_below_top_pair"]
    quotient, remainder = sp.div(
        pair23["cube_target"], pair23["cube_order_factor"]
    )
    if quotient * pair23["cube_order_factor"] + remainder != pair23[
        "cube_target"
    ]:
        raise AssertionError("top-pair division identity failed")
    pair23_verified = {
        "quotient_statistics": exact_bernstein_statistics(quotient),
        "remainder_statistics": exact_bernstein_statistics(remainder),
        "quotient": quotient,
        "remainder": remainder,
    }

    pair32 = cases["three_below_top_pair"]
    pair32_verified = _sparse_exact_multiplier(
        pair32["cube_target"],
        pair32["cube_order_factor"],
        THREE_BELOW_TOP_PAIR_MULTIPLIER_DEGREE,
        THREE_BELOW_TOP_PAIR_SUPPORT_FLAT,
        THREE_BELOW_TOP_PAIR_ACTIVE_ROWS_FLAT,
    )
    if tuple(map(int, pair32["cube_target"].degree_list())) != tuple(
        THREE_BELOW_TOP_PAIR_TARGET_DEGREE
    ):
        raise AssertionError("three-below target witness degree changed")

    triple23 = cases["two_below_top_triple"]
    triple23_verified = _sparse_exact_multiplier(
        triple23["cube_target"],
        triple23["cube_order_factor"],
        TWO_BELOW_TOP_TRIPLE_MULTIPLIER_DEGREE,
        TWO_BELOW_TOP_TRIPLE_SUPPORT_FLAT,
        TWO_BELOW_TOP_TRIPLE_ACTIVE_ROWS_FLAT,
    )
    if tuple(map(int, triple23["cube_target"].degree_list())) != tuple(
        TWO_BELOW_TOP_TRIPLE_TARGET_DEGREE
    ):
        raise AssertionError("top-triple target witness degree changed")

    compactness = derived["compactness"]
    return {
        "schema_version": 1,
        "mathematical_status": (
            "ordinary face reduction with exact rational polynomial certificates"
        ),
        "scope": {
            "hypothesis": (
                'The section has five uniform coordinates and one of the three listed repeated-top knot shapes, an ordered centroid, and the reduced gap inequalities of its chart. The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative.'
            ),
            "conclusion": (
                "In each of the three nontrivial repeated-top configurations, "
                "p<=I_{gamma_0}(1,4)=1-(1-gamma_0)^4."
            ),
            "global_consequence": (
                'The listed local first-coordinate inequalities are certified. The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.'
            ),
        },
        "face_exhaustion": {
            "unique_largest_knot": (
                'The identities extend algebraically, but collision limits require preservation of the reduced centroid-gap inequalities; ordered centroid alone does not justify division by a vanished factor.'
            ),
            "three_below_top_pair": "(-1,-u,-u*s,t,t)",
            "two_below_top_pair": "(-1,-u,t*v,t,t)",
            "two_below_top_triple": "(-1,-u,t,t,t)",
            "one_below": (
                "Orderedness forces all four upper knots to coincide; the cap "
                "then equals I_{gamma_0}(1,4)."
            ),
            "one_above": "AM--GM gives p<=I_{M_4}(4,1).",
        },
        "three_below_top_pair": {
            "variables": ["u", "s", "tau"],
            "domain": "0<u,s,tau<1; t=(3/2)*tau",
            "hypothesis_used": "gamma_3=gamma_4>=gamma_2",
            "scale_bound": "G>=0 and exact positivity of -G(u,s,3/2+rho) imply t<3/2",
            "order_identity": pair32["order_identity"],
            "comparison_identity": (
                "I_{gamma_0}(1,4)-p=t*T/Delta_margin, with Delta_margin>0"
            ),
            "certificate_identity": "T_tilde=S*G_tilde+H",
            "polynomials": {
                "target_T_before_rescaling": _polynomial_record(pair32["target"]),
                "order_G_before_rescaling": _polynomial_record(
                    pair32["order_factor"]
                ),
                "shifted_minus_G": _polynomial_record(
                    compactness["three_below_top_pair"]["shifted_polynomial"]
                ),
            },
            "compactness_sign_audit": {
                "basis": "Bernstein in (u,s), power basis in rho",
                **compactness["three_below_top_pair"]["statistics"],
            },
            "exact_linear_witness": {
                "matrix_order": pair32_verified["linear_system_order"],
                "support_flat_index_sha256": pair32_verified["support_digest"],
                "active_target_rows_flat_index_sha256": pair32_verified[
                    "active_rows_digest"
                ],
            },
            "bernstein_basis": {
                "multiplier_S": pair32_verified["multiplier_statistics"],
                "residual_H": pair32_verified["residual_statistics"],
            },
        },
        "two_below_top_pair": {
            "variables": ["u", "v", "tau"],
            "domain": "0<u,v,tau<1; t=tau",
            "hypothesis_used": "gamma_3=gamma_4>=gamma_2",
            "scale_bound": (
                "R=t^2*u^3-u^4+Q with Q coefficientwise positive; "
                "R<=0 implies t^2<u<1"
            ),
            "order_identity": pair23["order_identity"],
            "comparison_identity": (
                "I_{gamma_0}(1,4)-p=t*T/Delta_margin, with Delta_margin>0"
            ),
            "certificate_identity": "T=S*G+H by exact lexicographic division",
            "polynomials": {
                "target_T": _polynomial_record(pair23["target"]),
                "order_G": _polynomial_record(pair23["order_factor"]),
                "confinement_Q": _power_positive_record(pair23["confinement"]),
                "quotient_S": _polynomial_record(pair23_verified["quotient"]),
                "remainder_H": _polynomial_record(pair23_verified["remainder"]),
            },
            "bernstein_basis": {
                "quotient_S": pair23_verified["quotient_statistics"],
                "remainder_H": pair23_verified["remainder_statistics"],
            },
        },
        "two_below_top_triple": {
            "variables": ["u", "tau"],
            "domain": "0<u,tau<1; t=(2/3)*tau",
            "hypothesis_used": "gamma_2=gamma_3=gamma_4>=gamma_1",
            "scale_bound": "G>=0 and exact positivity of -G(u,2/3+rho) imply t<2/3",
            "order_identity": triple23["order_identity"],
            "comparison_identity": (
                "I_{gamma_0}(1,4)-p=t*T/Delta_margin, with Delta_margin>0"
            ),
            "certificate_identity": "T_tilde=S*G_tilde+H",
            "polynomials": {
                "target_T_before_rescaling": _polynomial_record(
                    triple23["target"]
                ),
                "order_G_before_rescaling": _polynomial_record(
                    triple23["order_factor"]
                ),
                "shifted_minus_G": _polynomial_record(
                    compactness["two_below_top_triple"]["shifted_polynomial"]
                ),
            },
            "compactness_sign_audit": {
                "basis": "Bernstein in u, power basis in rho",
                **compactness["two_below_top_triple"]["statistics"],
            },
            "exact_linear_witness": {
                "matrix_order": triple23_verified["linear_system_order"],
                "support_flat_index_sha256": triple23_verified["support_digest"],
                "active_target_rows_flat_index_sha256": triple23_verified[
                    "active_rows_digest"
                ],
            },
            "bernstein_basis": {
                "multiplier_S": triple23_verified["multiplier_statistics"],
                "residual_H": triple23_verified["residual_statistics"],
            },
        },
        "arithmetic": {
            "derivation_and_division": "SymPy exact rational polynomial arithmetic",
            "sparse_linear_witnesses": "FLINT exact rational matrix solves",
            "bernstein_conversion_and_signs": "Python Fraction arithmetic",
            "floating_point_role": "none",
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    pair32 = certificate["three_below_top_pair"]["bernstein_basis"]
    pair23 = certificate["two_below_top_pair"]["bernstein_basis"]
    triple23 = certificate["two_below_top_triple"]["bernstein_basis"]
    print(
        "five-coordinate repeated-knot faces certified: "
        f"3--2 top pair {pair32['multiplier_S']['total_coefficients']}/"
        f"{pair32['residual_H']['total_coefficients']}; "
        f"2--3 top pair {pair23['quotient_S']['total_coefficients']}/"
        f"{pair23['remainder_H']['total_coefficients']}; "
        f"top triple {triple23['multiplier_S']['total_coefficients']}/"
        f"{triple23['residual_H']['total_coefficients']} exact coefficients"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
