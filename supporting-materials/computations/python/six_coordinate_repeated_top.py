"""Exact local certificates for six repeated-top knot strata.

The six-coordinate adjacent-block hierarchy has six non-elementary faces on
which the largest knot is repeated.  This module proves the first-step beta
cap bound under the reduced centroid-gap constraints for:

    (-1,-u,-v,-w,t,t),       (-1,-u,-v,t,t,t),
    (-1,-u,t,t,t,t),         (-1,-u,t*a,t,t,t),
    (-1,-u,t*a,t*b,t,t),   (-1,-u,-v,t*a,t,t).

The first three have one positive knot block.  Edge-intercept coordinates
give low-degree cap and centroid formulae, and sparse exact rational systems
construct nonnegative tensor-Bernstein multipliers.  The last two are derived
from confluent divided differences.  Exact multivariate division gives
``T = S*G + H``; a Bernstein homogenization proves ``S,H >= 0`` on the unit
cube.  The final ``32`` face uses exact projective blow-ups and a finite
dyadic Bernstein cover. On each relative-open stratum, orderedness gives
``G >= 0``, and the reduced constraints imply

    p <= I_{gamma_0}(1,5) = 1 - (1-gamma_0)^5.

The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative.
The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.

All derivations, matrix solves, divisions, and sign tests use exact rational
arithmetic.  Floating point is not used in the public certificate.
"""

from __future__ import annotations

import argparse
import gc
import json
from fractions import Fraction
from itertools import combinations_with_replacement
from math import prod
from pathlib import Path

import numpy as np
import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx

from five_coordinate_block_faces import _sparse_exact_multiplier
from five_coordinate_mixed_section import (
    _fraction,
    _fraction_record,
    _power_digest,
    exact_bernstein_coefficients,
    exact_coefficient_statistics,
)
from six_coordinate_repeated_top_witness import (
    TOP_BLOCK_2_ACTIVE_ROWS_FLAT,
    TOP_BLOCK_2_MULTIPLIER_DEGREE,
    TOP_BLOCK_2_ORDER_DEGREE,
    TOP_BLOCK_2_SUPPORT_FLAT,
    TOP_BLOCK_2_TARGET_DEGREE,
    TOP_BLOCK_3_ACTIVE_ROWS_FLAT,
    TOP_BLOCK_3_MULTIPLIER_DEGREE,
    TOP_BLOCK_3_ORDER_DEGREE,
    TOP_BLOCK_3_SUPPORT_FLAT,
    TOP_BLOCK_3_TARGET_DEGREE,
    TOP_BLOCK_4_ACTIVE_ROWS_FLAT,
    TOP_BLOCK_4_MULTIPLIER_DEGREE,
    TOP_BLOCK_4_ORDER_DEGREE,
    TOP_BLOCK_4_SUPPORT_FLAT,
    TOP_BLOCK_4_TARGET_DEGREE,
)
from six_coordinate_face32 import (
    derive_face32,
    face32_certificate,
    verify_face32,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent
    / "certificates"
    / "six-coordinate-repeated-top-certificate.json"
)


TOP_BLOCK_WITNESSES = {
    2: {
        "order_degree": TOP_BLOCK_2_ORDER_DEGREE,
        "target_degree": TOP_BLOCK_2_TARGET_DEGREE,
        "multiplier_degree": TOP_BLOCK_2_MULTIPLIER_DEGREE,
        "support": TOP_BLOCK_2_SUPPORT_FLAT,
        "rows": TOP_BLOCK_2_ACTIVE_ROWS_FLAT,
    },
    3: {
        "order_degree": TOP_BLOCK_3_ORDER_DEGREE,
        "target_degree": TOP_BLOCK_3_TARGET_DEGREE,
        "multiplier_degree": TOP_BLOCK_3_MULTIPLIER_DEGREE,
        "support": TOP_BLOCK_3_SUPPORT_FLAT,
        "rows": TOP_BLOCK_3_ACTIVE_ROWS_FLAT,
    },
    4: {
        "order_degree": TOP_BLOCK_4_ORDER_DEGREE,
        "target_degree": TOP_BLOCK_4_TARGET_DEGREE,
        "multiplier_degree": TOP_BLOCK_4_MULTIPLIER_DEGREE,
        "support": TOP_BLOCK_4_SUPPORT_FLAT,
        "rows": TOP_BLOCK_4_ACTIVE_ROWS_FLAT,
    },
}


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
    """Remove the greatest monomial dividing every term."""
    minima = tuple(
        min(powers[index] for powers, _ in poly.terms())
        for index in range(len(poly.gens))
    )
    monomial = _poly(
        sp.prod(poly.gens[index] ** minima[index] for index in range(len(minima))),
        poly.gens,
    )
    return sp.exquo(poly, monomial), minima


def _complete_homogeneous(variables, degree: int):
    if degree == 0:
        return sp.Integer(1)
    return sp.Add(
        *(
            sp.prod(variables[index] for index in indices)
            for indices in combinations_with_replacement(
                range(len(variables)), degree
            )
        )
    )


def derive_one_positive_block(multiplicity: int) -> dict[str, object]:
    """Derive the low-degree face with one positive top block of size ``m``."""
    if multiplicity not in (2, 3, 4):
        raise ValueError("top multiplicity must be 2, 3, or 4")
    lower_count = 6 - multiplicity
    r_all = sp.symbols("r0:4")
    x_all = sp.symbols("x0:4")
    r = r_all[:lower_count]
    x = x_all[:lower_count]

    # After the top knot is scaled to one, the negative knots are
    # -(1-r_i)/r_i.  Confluent interpolation simplifies to this complete
    # homogeneous-polynomial formula for the cap probability.
    cap = sp.expand(
        sp.prod(r)
        * _complete_homogeneous(
            (sp.Integer(1), *(1 - coordinate for coordinate in r)),
            multiplicity - 1,
        )
    )
    derivatives = tuple(sp.diff(cap, coordinate) for coordinate in r)
    normalizer = sp.expand(
        sum(coordinate * derivative for coordinate, derivative in zip(r, derivatives))
    )
    gamma0_numerator = sp.expand(r[0] ** 2 * derivatives[0])

    # The individual top-block centroid is
    # sum r_i(1-r_i)P_i/(m Q).  Thus the last displayed expression is
    # mQ times (gamma_top-gamma_last_lower).
    raw_order = sp.expand(
        sum(
            coordinate * (1 - coordinate) * derivative
            for coordinate, derivative in zip(r, derivatives)
        )
        - multiplicity * r[-1] ** 2 * derivatives[-1]
    )
    raw_target = sp.expand(
        normalizer**5
        - (normalizer - gamma0_numerator) ** 5
        - cap * normalizer**5
    )

    ordered_substitution = {
        r[index]: sp.prod(x[index:]) for index in range(lower_count)
    }
    order, order_monomial = _remove_positive_monomial(
        _poly(raw_order.subs(ordered_substitution), x)
    )
    target, target_monomial = _remove_positive_monomial(
        _poly(raw_target.subs(ordered_substitution), x)
    )

    witness = TOP_BLOCK_WITNESSES[multiplicity]
    if tuple(map(int, order.degree_list())) != witness["order_degree"]:
        raise AssertionError("one-positive-block order degree changed")
    if tuple(map(int, target.degree_list())) != witness["target_degree"]:
        raise AssertionError("one-positive-block target degree changed")
    expected_terms = {2: (15, 505), 3: (20, 660), 4: (15, 210)}
    if (len(order.terms()), len(target.terms())) != expected_terms[multiplicity]:
        raise AssertionError("one-positive-block sparse structure changed")

    return {
        "multiplicity": multiplicity,
        "lower_count": lower_count,
        "r_variables": r,
        "cube_variables": x,
        "cap": _poly(cap, r),
        "normalizer": _poly(normalizer, r),
        "gamma0_numerator": _poly(gamma0_numerator, r),
        "raw_order": _poly(raw_order, r),
        "raw_target": _poly(raw_target, r),
        "order": order,
        "target": target,
        "order_monomial": order_monomial,
        "target_monomial": target_monomial,
    }


def verify_one_positive_block(derived: dict[str, object]) -> dict[str, object]:
    multiplicity = int(derived["multiplicity"])
    witness = TOP_BLOCK_WITNESSES[multiplicity]
    result = _sparse_exact_multiplier(
        derived["target"],
        derived["order"],
        witness["multiplier_degree"],
        witness["support"],
        witness["rows"],
    )
    return {
        "linear_system_order": result["linear_system_order"],
        "support_flat_index_sha256": result["support_digest"],
        "active_target_rows_flat_index_sha256": result["active_rows_digest"],
        "multiplier": result["multiplier_statistics"],
        "residual": result["residual_statistics"],
    }


def _confluent_divided_difference(knots, multiplicities, positive, power):
    """Confluent divided difference of ``z**power`` over selected residues."""
    z = sp.symbols("z")
    answer = 0
    for index in positive:
        summand = z**power / sp.prod(
            (z - knots[j]) ** multiplicities[j]
            for j in range(len(knots))
            if j != index
        )
        answer += sp.diff(summand, z, multiplicities[index] - 1).subs(
            z, knots[index]
        ) / sp.factorial(multiplicities[index] - 1)
    return sp.cancel(answer)


def _derive_mixed_block_case(name: str) -> dict[str, object]:
    """Derive the 2--3 or 2--2 repeated-top face from interpolation."""
    u, a, b, t = sp.symbols("u a b t")
    if name == "23":
        specialized_knots = (-1, -u, t * a, t)
        multiplicities = (1, 1, 1, 3)
        variables = (u, a, t)
        positive = range(2, 4)
        ratio = a
        bound = sp.Rational(2, 3)
    elif name == "22":
        specialized_knots = (-1, -u, t * a, t * b, t)
        multiplicities = (1, 1, 1, 1, 2)
        variables = (u, a, b, t)
        positive = range(2, 5)
        ratio = b
        bound = sp.Integer(1)
    else:
        raise ValueError("mixed block case must be '23' or '22'")

    abstract_knots = sp.symbols(f"Y0:{len(specialized_knots)}")
    substitution = dict(zip(abstract_knots, specialized_knots, strict=True))
    cap_abstract = _confluent_divided_difference(
        abstract_knots, multiplicities, positive, 5
    )
    density_abstract = 5 * _confluent_divided_difference(
        abstract_knots, multiplicities, positive, 4
    )
    cap = sp.cancel(cap_abstract.subs(substitution))
    density = sp.cancel(density_abstract.subs(substitution))
    derivatives = tuple(
        sp.cancel(sp.diff(cap_abstract, knot).subs(substitution))
        for knot in abstract_knots
    )
    gamma0 = sp.cancel(derivatives[0] / density)
    top = len(abstract_knots) - 1
    previous = top - 1
    order_gap = sp.cancel(
        (
            derivatives[top] / multiplicities[top]
            - derivatives[previous] / multiplicities[previous]
        )
        / density
    )
    order_numerator, order_denominator = map(sp.factor, sp.fraction(order_gap))
    _coefficientwise_positive(_poly(order_denominator, variables))

    # The collision factor ratio-1 is nonpositive.  With the following sign,
    # the ordered adjacent-centroid inequality is exactly G>=0 after positive
    # monomial factors are removed.
    order_before_ordering = -sp.exquo(
        _poly(order_numerator, variables), _poly(ratio - 1, variables)
    )

    gamma_numerator, gamma_denominator = (
        _poly(expression, variables) for expression in sp.fraction(gamma0)
    )
    cap_numerator, cap_denominator = (
        _poly(expression, variables) for expression in sp.fraction(cap)
    )
    _coefficientwise_positive(gamma_denominator)
    _coefficientwise_positive(cap_denominator)

    raw_target = (
        (
            gamma_denominator**5
            - (gamma_denominator - gamma_numerator) ** 5
        )
        * cap_denominator
        - cap_numerator * gamma_denominator**5
    )
    positive_denominator = gamma_denominator**5 * cap_denominator
    common = sp.gcd(raw_target, positive_denominator)
    reduced_target = sp.exquo(raw_target, common)
    reduced_target, first_target_monomial = _remove_positive_monomial(
        reduced_target
    )
    content, reduced_target = reduced_target.primitive()
    if content <= 0:
        raise AssertionError("comparison target acquired a negative content")

    U, A, B, tau, rho = sp.symbols("U A B tau rho")
    if name == "23":
        cube_base = (U, A, t)
        ordered_substitution = {u: U, a: A}
    else:
        cube_base = (U, A, B, t)
        ordered_substitution = {u: U, a: A * B, b: B}
    order, order_monomial = _remove_positive_monomial(
        _poly(order_before_ordering.as_expr().subs(ordered_substitution), cube_base)
    )
    target, target_monomial = _remove_positive_monomial(
        _poly(reduced_target.as_expr().subs(ordered_substitution), cube_base)
    )

    expected = {
        "23": ((5, 2, 6), 32, (19, 11, 28), 2969),
        "22": ((5, 1, 3, 6), 50, (19, 11, 22, 33), 29801),
    }[name]
    if (
        tuple(map(int, order.degree_list())),
        len(order.terms()),
        tuple(map(int, target.degree_list())),
        len(target.terms()),
    ) != expected:
        raise AssertionError(f"mixed block {name} polynomials changed")

    shifted = _poly(
        (-order.as_expr()).subs(t, bound + rho), (*cube_base[:-1], rho)
    )
    compact_coefficients = exact_bernstein_coefficients(
        shifted, axes=tuple(range(len(cube_base) - 1))
    )
    compact_statistics = exact_coefficient_statistics(
        compact_coefficients, tuple(map(int, shifted.degree_list()))
    )
    if not any(
        compact_coefficients[(*index, 0)] > 0
        for index in np.ndindex(compact_coefficients.shape[:-1])
    ):
        raise AssertionError("compactness polynomial is not strict in the interior")

    cube_variables = (*cube_base[:-1], tau)
    cube_order = _poly(
        order.as_expr().subs(t, bound * tau), cube_variables
    )
    cube_target = _poly(
        target.as_expr().subs(t, bound * tau), cube_variables
    )

    return {
        "name": name,
        "variables": variables,
        "cube_variables": cube_variables,
        "cap": cap,
        "gamma0": gamma0,
        "order_gap": order_gap,
        "order_denominator": _poly(order_denominator, variables),
        "target": target,
        "order": order,
        "cube_target": cube_target,
        "cube_order": cube_order,
        "bound": bound,
        "shifted_minus_order": shifted,
        "compact_statistics": compact_statistics,
        "first_target_monomial": first_target_monomial,
        "order_monomial": order_monomial,
        "target_monomial": target_monomial,
    }


def derive_mixed_block_cases() -> dict[str, dict[str, object]]:
    return {name: _derive_mixed_block_case(name) for name in ("23", "22")}


def _to_flint_polynomial(poly: sp.Poly, context):
    return context.from_dict(
        {powers: fmpq(str(coefficient)) for powers, coefficient in poly.terms()}
    )


def _flint_fraction_record(value) -> dict[str, str]:
    return {"numerator": str(int(value.p)), "denominator": str(int(value.q))}


def _flint_polynomial_record(poly) -> dict[str, object]:
    return {
        "degree": list(map(int, poly.degrees())),
        "power_terms": len(poly),
    }


def _flint_bernstein_sign_statistics(poly, names) -> dict[str, object]:
    """Audit Bernstein signs through exact multihomogenization.

    If ``d`` is the coordinatewise degree, the coefficient of ``z**nu`` in

        prod_i (1+z_i)**d_i P(z_i/(1+z_i))

    is the tensor-Bernstein coefficient at ``nu`` multiplied by the positive
    number ``prod_i binomial(d_i,nu_i)``.  The two coefficient arrays therefore
    have exactly the same signs.
    """
    degrees = tuple(map(int, poly.degrees()))
    count = len(names)
    homogeneous_context = fmpq_mpoly_ctx.get(
        tuple(f"bz{index}" for index in range(count))
        + tuple(f"bh{index}" for index in range(count)),
        "lex",
    )
    homogeneous = homogeneous_context.from_dict(
        {
            tuple(powers)
            + tuple(degrees[index] - powers[index] for index in range(count)):
            coefficient
            for powers, coefficient in poly.terms()
        }
    )
    variables = poly.context().gens()
    transformed = homogeneous.compose(
        *(tuple(variables) + tuple(1 + variable for variable in variables))
    )
    del homogeneous

    positive = 0
    negative = 0
    smallest_positive = None
    largest = None
    for _, coefficient in transformed.terms():
        if coefficient > 0:
            positive += 1
            if smallest_positive is None or coefficient < smallest_positive:
                smallest_positive = coefficient
            if largest is None or coefficient > largest:
                largest = coefficient
        else:
            negative += 1
    total = prod(degree + 1 for degree in degrees)
    zero = total - len(transformed)
    del transformed
    gc.collect()
    if negative:
        raise AssertionError(f"found {negative} negative Bernstein coefficients")
    if smallest_positive is None or largest is None:
        raise AssertionError("Bernstein expansion has no positive coefficient")
    return {
        "degree": list(degrees),
        "total_coefficients": total,
        "positive_coefficients": positive,
        "zero_coefficients": zero,
        "negative_coefficients": negative,
        "smallest_positive_scaled": _flint_fraction_record(smallest_positive),
        "largest_scaled": _flint_fraction_record(largest),
    }


def verify_mixed_block_case(derived: dict[str, object]) -> dict[str, object]:
    names = tuple(map(str, derived["cube_variables"]))
    context = fmpq_mpoly_ctx.get(names, "lex")
    target = _to_flint_polynomial(derived["cube_target"], context)
    order = _to_flint_polynomial(derived["cube_order"], context)
    quotient, remainder = divmod(target, order)
    if quotient * order + remainder != target:
        raise AssertionError("exact multivariate division identity failed")
    quotient_record = _flint_polynomial_record(quotient)
    remainder_record = _flint_polynomial_record(remainder)
    quotient_signs = _flint_bernstein_sign_statistics(quotient, names)
    remainder_signs = _flint_bernstein_sign_statistics(remainder, names)
    return {
        "quotient": quotient_record,
        "remainder": remainder_record,
        "quotient_bernstein_signs": quotient_signs,
        "remainder_bernstein_signs": remainder_signs,
    }


def _top_block_certificate(derived, verified) -> dict[str, object]:
    multiplicity = int(derived["multiplicity"])
    lower_count = int(derived["lower_count"])
    face = "(" + ",".join(
        ["negative knot"] * lower_count + ["t"] * multiplicity
    ) + ")"
    return {
        "face_type": f"{lower_count}{multiplicity}",
        "normalized_shape": face,
        "variables": [str(variable) for variable in derived["cube_variables"]],
        "coordinate_change": (
            "top knot=1; lower knots=-(1-r_i)/r_i; "
            "r_i=product_{j=i}^{L-1} x_j"
        ),
        "cap_formula": (
            "P=(product r_i) h_{m-1}(1,1-r_0,...,1-r_{L-1})"
        ),
        "order_identity": (
            "m*Q*(gamma_top-gamma_last_lower) equals a positive monomial times G"
        ),
        "comparison_identity": (
            "Q^5*(I_{gamma_0}(1,5)-P) equals a positive monomial times T"
        ),
        "certificate_identity": "T=S*G+H in the tensor-Bernstein basis",
        "removed_positive_monomials": {
            "order": list(derived["order_monomial"]),
            "target": list(derived["target_monomial"]),
        },
        "polynomials": {
            "order_G": _polynomial_record(derived["order"]),
            "target_T": _polynomial_record(derived["target"]),
        },
        "exact_linear_witness": {
            "matrix_order": verified["linear_system_order"],
            "support_flat_index_sha256": verified[
                "support_flat_index_sha256"
            ],
            "active_target_rows_flat_index_sha256": verified[
                "active_target_rows_flat_index_sha256"
            ],
        },
        "bernstein_basis": {
            "multiplier_S": verified["multiplier"],
            "residual_H": verified["residual"],
        },
    }


def _mixed_certificate(derived, verified) -> dict[str, object]:
    name = str(derived["name"])
    if name == "23":
        face = "(-1,-u,t*a,t,t,t)"
        order_hypothesis = "gamma_top>=gamma_at_t*a"
    else:
        face = "(-1,-u,t*a,t*b,t,t)"
        order_hypothesis = "gamma_top>=gamma_at_t*b"
    return {
        "face_type": name,
        "normalized_shape": face,
        "variables": [str(variable) for variable in derived["cube_variables"]],
        "domain": (
            "unit cube after t=" + str(derived["bound"]) + "*tau"
        ),
        "hypothesis_used": order_hypothesis,
        "order_identity": (
            "the adjacent centroid-gap numerator is (ratio-1)*(-G) "
            "times a positive monomial, with positive denominator"
        ),
        "comparison_identity": (
            "I_{gamma_0}(1,5)-p is a positive rational factor times T"
        ),
        "scale_bound": (
            "G>=0 and the exact mixed Bernstein/power expansion of -G at "
            f"t={derived['bound']}+rho confine the ordered section to "
            f"t<{derived['bound']}"
        ),
        "certificate_identity": "T=S*G+H by exact lexicographic division",
        "removed_positive_monomials": {
            "initial_target": list(derived["first_target_monomial"]),
            "order_after_parameterization": list(derived["order_monomial"]),
            "target_after_parameterization": list(derived["target_monomial"]),
        },
        "polynomials": {
            "order_G_before_rescaling": _polynomial_record(derived["order"]),
            "target_T_before_rescaling": _polynomial_record(derived["target"]),
            "shifted_minus_G": _polynomial_record(
                derived["shifted_minus_order"]
            ),
            "quotient_S": verified["quotient"],
            "remainder_H": verified["remainder"],
        },
        "compactness_sign_audit": {
            "basis": "Bernstein in bounded variables and power basis in rho",
            **derived["compact_statistics"],
        },
        "bernstein_homogenization": {
            "coefficient_scaling": (
                "the transformed z^nu coefficient is the tensor-Bernstein "
                "coefficient times product_i binomial(d_i,nu_i)>0"
            ),
            "quotient_S": verified["quotient_bernstein_signs"],
            "remainder_H": verified["remainder_bernstein_signs"],
        },
    }


def build_certificate() -> dict[str, object]:
    top_cases = {
        multiplicity: derive_one_positive_block(multiplicity)
        for multiplicity in (2, 3, 4)
    }
    top_verified = {
        multiplicity: verify_one_positive_block(derived)
        for multiplicity, derived in top_cases.items()
    }
    mixed_cases = derive_mixed_block_cases()
    mixed_verified = {
        name: verify_mixed_block_case(derived)
        for name, derived in mixed_cases.items()
    }
    final_repeated_top = derive_face32()
    final_repeated_top_verified = verify_face32(final_repeated_top)

    return {
        "schema_version": 1,
        "mathematical_status": (
            "ordinary interpolation and face reductions with exact rational "
            "polynomial certificates"
        ),
        "scope": {
            "hypothesis": (
                'The section has six uniform coordinates, one of the six listed repeated-top knot shapes, an ordered centroid, and the reduced gap inequalities of its chart. The certified geometric implication applies on the listed relative-open knot strata with all reduced centroid-gap inequalities imposed, and to limits that preserve those inequalities. A vanished collision factor does not imply that the divided centroid-gap polynomial remains nonnegative.'
            ),
            "proved_faces": ["42", "33", "24", "23", "22", "32"],
            "open_repeated_top_faces": [],
            "conclusion_on_proved_faces": (
                "p<=I_{gamma_0}(1,5)=1-(1-gamma_0)^5"
            ),
            "global_limitation": (
                'The six repeated-top families and companion first-coordinate certificates address all nine local families, not the missing active-prefix implication. The first-coordinate beta bound is not an active-prefix bound on repeated-lowest faces. The global monotone cap theorem remains unresolved.'
            ),
        },
        "one_positive_top_block": {
            str(multiplicity): _top_block_certificate(
                top_cases[multiplicity], top_verified[multiplicity]
            )
            for multiplicity in (2, 3, 4)
        },
        "mixed_positive_blocks": {
            name: _mixed_certificate(mixed_cases[name], mixed_verified[name])
            for name in ("23", "22")
        },
        "final_repeated_top_face": face32_certificate(
            final_repeated_top, final_repeated_top_verified
        ),
        "arithmetic": {
            "symbolic_derivation": "SymPy exact rational arithmetic",
            "sparse_multiplier_systems": "FLINT exact rational matrix solves",
            "multivariate_division": "FLINT exact rational polynomial division",
            "bernstein_signs": (
                "Python Fraction arrays for sparse multipliers; FLINT exact "
                "Bernstein homogenization for division quotients, remainders, "
                "and dyadic projective-chart leaves"
            ),
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
    print(
        "all six-coordinate repeated-top faces certified exactly: "
        + ", ".join(certificate["scope"]["proved_faces"])
        + "; companion certificates prove all three unique-top families"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
