"""Exact certificate for the final repeated-top face in six coordinates.

The normalized knots are ``(-1,-u,-v,t*a,t,t)`` with
``1 >= u >= v > 0``, ``0 < a < 1``, and ``t > 0``.  Confluent divided
differences derive the cap and all four adjacent section-centroid gaps.
After positive factors are removed, orderedness is ``G_i >= 0`` and the
desired first beta-cap comparison is ``T >= 0``.

The proof uses two scale charts.  The chart ``t <= u`` closes by a small
dyadic tensor-Bernstein subdivision.  In the chart ``u <= t``, the final
centroid gap first gives the exact quadratic confinement ``t <= 3(u/t)^2``.
Two projective charts then compare the remaining positive knot ratio with
``u/t``.  Every leaf is closed by ``T``, by an infeasible sign ``G_i < 0``,
or by ``T-mu*G_i >= 0`` for an explicit nonnegative rational ``mu``.

Every symbolic derivation, affine chart, multiplier, and Bernstein sign is
checked with exact rational arithmetic.  Floating point is not used.
"""

from __future__ import annotations

import gc
import hashlib
from fractions import Fraction
from math import prod

import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx

from five_coordinate_mixed_section import (
    _fraction_record,
    _power_digest,
    exact_bernstein_coefficients,
    exact_coefficient_statistics,
)
from six_coordinate_face32_witness import (
    A_LE_Z_LEAVES,
    T_LE_U_LEAVES,
    Z_LE_A_LEAVES,
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


def _confluent_divided_difference(knots, multiplicities, positive, power):
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


def _bernstein_statistics(
    poly: sp.Poly,
    axes: tuple[int, ...] | None = None,
) -> dict[str, object]:
    coefficients = exact_bernstein_coefficients(poly, axes=axes)
    return exact_coefficient_statistics(
        coefficients, tuple(map(int, poly.degree_list()))
    )


def derive_face32() -> dict[str, object]:
    """Derive the cap target and four ordered-centroid gap polynomials."""
    u, v, a, t = sp.symbols("u v a t")
    variables = (u, v, a, t)
    specialized_knots = (-1, -u, -v, t * a, t)
    multiplicities = (1, 1, 1, 1, 2)
    abstract_knots = sp.symbols("Y0:5")
    specialization = dict(zip(abstract_knots, specialized_knots, strict=True))

    cap_abstract = _confluent_divided_difference(
        abstract_knots, multiplicities, range(3, 5), 5
    )
    density_abstract = 5 * _confluent_divided_difference(
        abstract_knots, multiplicities, range(3, 5), 4
    )
    cap = sp.cancel(cap_abstract.subs(specialization))
    density = sp.cancel(density_abstract.subs(specialization))
    derivatives = tuple(
        sp.cancel(sp.diff(cap_abstract, knot).subs(specialization))
        for knot in abstract_knots
    )
    centroids = tuple(
        sp.cancel(derivative / (multiplicity * density))
        for derivative, multiplicity in zip(
            derivatives, multiplicities, strict=True
        )
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
    reduced_target, initial_target_monomial = _remove_positive_monomial(
        reduced_target
    )
    target_content, reduced_target = reduced_target.primitive()
    if target_content <= 0:
        raise AssertionError("face 32 target acquired negative content")

    expected_target = ((13, 13, 17, 35), 40_732, (0, 0, 0, 1))
    if (
        tuple(map(int, reduced_target.degree_list())),
        len(reduced_target.terms()),
        initial_target_monomial,
    ) != expected_target:
        raise AssertionError("face 32 target structure changed")

    U, V, A, q = sp.symbols("U V A q")
    cube_variables = (U, V, A, q)
    ordered_substitution = {u: U, v: U * V, a: A, t: q}
    cube_target, cube_target_monomial = _remove_positive_monomial(
        _poly(reduced_target.as_expr().subs(ordered_substitution), cube_variables)
    )
    if cube_target_monomial != (0, 0, 0, 0):
        raise AssertionError("unexpected face 32 cube target monomial")
    if (
        tuple(map(int, cube_target.degree_list())),
        len(cube_target.terms()),
    ) != ((26, 13, 17, 35), 40_732):
        raise AssertionError("face 32 cube target structure changed")

    # Each numerator is an elementary positive collision factor times G_i.
    # Thus the four adjacent centroid inequalities are exactly G_i >= 0 in
    # the relative interior of the ordered cube.
    gap_cofactors = (
        _poly(q * (1 - U), cube_variables),
        _poly(U * q * (1 - V), cube_variables),
        _poly(1, cube_variables),
        _poly(1 - A, cube_variables),
    )
    expected_gaps = (
        ((5, 2, 4, 7), 102),
        ((6, 3, 4, 7), 102),
        ((7, 4, 5, 8), 169),
        ((8, 4, 4, 9), 213),
    )
    gaps = []
    gap_denominators = []
    cofactor_statistics = []
    for index, cofactor in enumerate(gap_cofactors):
        gap = sp.cancel(centroids[index + 1] - centroids[index])
        numerator, denominator = sp.fraction(gap)
        numerator_poly = _poly(numerator.subs(ordered_substitution), cube_variables)
        denominator_poly = _poly(
            denominator.subs(ordered_substitution), cube_variables
        )
        _coefficientwise_positive(denominator_poly)
        gap_poly = sp.exquo(numerator_poly, cofactor)
        gap_content, gap_poly = gap_poly.primitive()
        if gap_content <= 0:
            raise AssertionError("face 32 gap acquired negative content")
        if (
            tuple(map(int, gap_poly.degree_list())),
            len(gap_poly.terms()),
        ) != expected_gaps[index]:
            raise AssertionError(f"face 32 gap {index} structure changed")
        if numerator_poly != gap_content * cofactor * gap_poly:
            raise AssertionError("face 32 gap factorization changed")
        cofactor_statistics.append(_bernstein_statistics(cofactor))
        gaps.append(gap_poly)
        gap_denominators.append(denominator_poly)

    rho = sp.symbols("rho")
    shifted_minus_last_gap = _poly(
        (-gaps[3].as_expr()).subs(q, 1 + rho), (U, V, A, rho)
    )
    scale_statistics = _bernstein_statistics(
        shifted_minus_last_gap, axes=(0, 1, 2)
    )
    if (
        scale_statistics["positive_coefficients"],
        scale_statistics["zero_coefficients"],
    ) != (1_347, 903):
        raise AssertionError("face 32 scale-confinement signs changed")

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
        "shifted_minus_last_gap": shifted_minus_last_gap,
        "scale_statistics": scale_statistics,
    }


def _chart_polynomial(
    poly: sp.Poly,
    substitution: dict[sp.Symbol, sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Poly, tuple[int, ...]]:
    return _remove_positive_monomial(
        _poly(poly.as_expr().subs(substitution), variables)
    )


def derive_face32_charts(derived: dict[str, object]) -> dict[str, object]:
    """Construct the three compact projective charts used by the proof."""
    U, V, A, q = derived["cube_variables"]
    target = derived["target"]
    gaps = derived["gaps"]
    Z, R, B, rho = sp.symbols("Z R B rho")

    t_le_u_variables = (U, V, A, Z)
    t_le_u_substitution = {q: U * Z}
    t_le_u_target, t_le_u_target_monomial = _chart_polynomial(
        target, t_le_u_substitution, t_le_u_variables
    )
    t_le_u_gaps_with_monomials = tuple(
        _chart_polynomial(gap, t_le_u_substitution, t_le_u_variables)
        for gap in gaps
    )
    t_le_u_gaps = tuple(value[0] for value in t_le_u_gaps_with_monomials)

    u_le_t_variables = (q, V, A, Z)
    u_le_t_substitution = {U: q * Z}
    u_le_t_target, u_le_t_target_monomial = _chart_polynomial(
        target, u_le_t_substitution, u_le_t_variables
    )
    u_le_t_gaps_with_monomials = tuple(
        _chart_polynomial(gap, u_le_t_substitution, u_le_t_variables)
        for gap in gaps
    )
    u_le_t_gaps = tuple(value[0] for value in u_le_t_gaps_with_monomials)

    quadratic_confinement = _poly(
        (-u_le_t_gaps[3].as_expr()).subs(q, 3 * Z**2 + rho),
        (Z, V, A, rho),
    )
    quadratic_statistics = _bernstein_statistics(quadratic_confinement)
    if (
        tuple(map(int, quadratic_confinement.degree_list())),
        len(quadratic_confinement.terms()),
        quadratic_statistics["positive_coefficients"],
        quadratic_statistics["zero_coefficients"],
    ) != ((16, 4, 4, 4), 464, 1_630, 495):
        raise AssertionError("face 32 quadratic confinement changed")

    qz2_variables = (Z, V, A, R)
    qz2_substitution = {q: 3 * Z**2 * R}
    qz2_target, qz2_target_monomial = _chart_polynomial(
        u_le_t_target, qz2_substitution, qz2_variables
    )
    qz2_gaps_with_monomials = tuple(
        _chart_polynomial(gap, qz2_substitution, qz2_variables)
        for gap in u_le_t_gaps
    )
    qz2_gaps = tuple(value[0] for value in qz2_gaps_with_monomials)

    a_le_z_variables = (Z, B, V, R)
    a_le_z_target, a_le_z_target_monomial = _chart_polynomial(
        qz2_target, {A: Z * B}, a_le_z_variables
    )
    a_le_z_gaps_with_monomials = tuple(
        _chart_polynomial(gap, {A: Z * B}, a_le_z_variables)
        for gap in qz2_gaps
    )
    a_le_z_gaps = tuple(value[0] for value in a_le_z_gaps_with_monomials)

    z_le_a_variables = (A, B, V, R)
    z_le_a_target, z_le_a_target_monomial = _chart_polynomial(
        qz2_target, {Z: A * B}, z_le_a_variables
    )
    z_le_a_gaps_with_monomials = tuple(
        _chart_polynomial(gap, {Z: A * B}, z_le_a_variables)
        for gap in qz2_gaps
    )
    z_le_a_gaps = tuple(value[0] for value in z_le_a_gaps_with_monomials)

    expected = {
        "t_le_u": ((19, 13, 17, 35), (4, 4, 4, 9)),
        "a_le_z": ((65, 17, 13, 19), (12, 4, 4, 4)),
        "z_le_a": ((65, 62, 13, 19), (12, 14, 4, 4)),
    }
    actual = {
        "t_le_u": (
            tuple(map(int, t_le_u_target.degree_list())),
            tuple(map(int, t_le_u_gaps[3].degree_list())),
        ),
        "a_le_z": (
            tuple(map(int, a_le_z_target.degree_list())),
            tuple(map(int, a_le_z_gaps[3].degree_list())),
        ),
        "z_le_a": (
            tuple(map(int, z_le_a_target.degree_list())),
            tuple(map(int, z_le_a_gaps[3].degree_list())),
        ),
    }
    if actual != expected:
        raise AssertionError("face 32 projective chart degrees changed")

    return {
        "t_le_u": {
            "variables": t_le_u_variables,
            "target": t_le_u_target,
            "gaps": t_le_u_gaps,
            "target_monomial": t_le_u_target_monomial,
            "gap_monomials": tuple(
                value[1] for value in t_le_u_gaps_with_monomials
            ),
        },
        "quadratic_confinement": quadratic_confinement,
        "quadratic_statistics": quadratic_statistics,
        "u_le_t_monomials": {
            "target": u_le_t_target_monomial,
            "gaps": tuple(value[1] for value in u_le_t_gaps_with_monomials),
        },
        "qz2_monomials": {
            "target": qz2_target_monomial,
            "gaps": tuple(value[1] for value in qz2_gaps_with_monomials),
        },
        "a_le_z": {
            "variables": a_le_z_variables,
            "target": a_le_z_target,
            "gaps": a_le_z_gaps,
            "target_monomial": a_le_z_target_monomial,
            "gap_monomials": tuple(
                value[1] for value in a_le_z_gaps_with_monomials
            ),
        },
        "z_le_a": {
            "variables": z_le_a_variables,
            "target": z_le_a_target,
            "gaps": z_le_a_gaps,
            "target_monomial": z_le_a_target_monomial,
            "gap_monomials": tuple(
                value[1] for value in z_le_a_gaps_with_monomials
            ),
        },
    }


def _validate_leaf_partition(leaves) -> dict[str, int]:
    trie: dict[str, object] = {}
    maximum_depth = 0
    for path, closure, multiplier in leaves:
        if not path or any(bit not in "01" for bit in path):
            raise AssertionError("invalid binary subdivision path")
        if closure not in {"T", "G2", "G3", "M2", "M3"}:
            raise AssertionError("invalid face 32 leaf closure")
        if closure.startswith("M"):
            if multiplier is None or Fraction(*multiplier) < 0:
                raise AssertionError("invalid face 32 multiplier")
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
            raise AssertionError("subdivision tree does not cover the cube")
        left_nodes, left_leaves = check_complete(node["0"])
        right_nodes, right_leaves = check_complete(node["1"])
        return 1 + left_nodes + right_nodes, left_leaves + right_leaves

    nodes, leaf_count = check_complete(trie)
    if leaf_count != len(leaves):
        raise AssertionError("leaf count mismatch")
    return {
        "nodes": nodes,
        "leaves": leaf_count,
        "maximum_depth": maximum_depth,
    }


def _dyadic_box(path: str) -> tuple[tuple[Fraction, Fraction], ...]:
    intervals = [(Fraction(0), Fraction(1)) for _ in range(4)]
    for depth, bit in enumerate(path):
        axis = depth % 4
        lower, upper = intervals[axis]
        middle = (lower + upper) / 2
        intervals[axis] = (middle, upper) if bit == "1" else (lower, middle)
    return tuple(intervals)


def _to_flint_polynomial(poly: sp.Poly, context):
    return context.from_dict(
        {powers: fmpq(str(coefficient)) for powers, coefficient in poly.terms()}
    )


def _flint_fraction_record(value) -> dict[str, str]:
    return {"numerator": str(int(value.p)), "denominator": str(int(value.q))}


def _local_bernstein_sign_record(
    poly,
    variables,
    path: str,
    degree: tuple[int, ...],
) -> dict[str, object]:
    context = poly.context()
    generators = context.gens()
    substitutions = []
    for generator, (lower, upper) in zip(
        generators, _dyadic_box(path), strict=True
    ):
        lo = fmpq(lower.numerator, lower.denominator)
        width_fraction = upper - lower
        width = fmpq(width_fraction.numerator, width_fraction.denominator)
        substitutions.append(lo + width * generator)
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
    positive = len(coefficients)
    smallest_positive = min(coefficients) if coefficients else None
    largest = max(coefficients) if coefficients else None
    negative = (
        sum(coefficient < 0 for coefficient in coefficients)
        if smallest_positive is not None and smallest_positive < 0
        else 0
    )
    total = prod(value + 1 for value in degree)
    del coefficients, transformed
    gc.collect()
    if negative:
        raise AssertionError(
            f"face 32 leaf {path} has {negative} negative Bernstein coefficients"
        )
    if smallest_positive is None or largest is None:
        raise AssertionError("face 32 leaf certificate is not strict")
    return {
        "degree": list(degree),
        "total_coefficients": total,
        "positive_coefficients": positive,
        "zero_coefficients": total - positive,
        "negative_coefficients": negative,
        "smallest_positive_scaled": _flint_fraction_record(smallest_positive),
        "largest_scaled": _flint_fraction_record(largest),
    }


def verify_chart(chart: dict[str, object], leaves) -> dict[str, object]:
    """Replay one exact dyadic leaf partition and aggregate its signs."""
    tree = _validate_leaf_partition(leaves)
    variables = chart["variables"]
    target = chart["target"]
    gaps = chart["gaps"]
    context = fmpq_mpoly_ctx.get(tuple(map(str, variables)), "lex")
    target_flint = _to_flint_polynomial(target, context)
    gap_flint = tuple(_to_flint_polynomial(gap, context) for gap in gaps)
    target_degree = tuple(map(int, target.degree_list()))

    closure_counts: dict[str, int] = {}
    total_coefficients = 0
    positive_coefficients = 0
    zero_coefficients = 0
    smallest_positive = None
    largest = None
    leaf_digest = hashlib.sha256()
    multiplier_records = []
    for path, closure, multiplier_pair in leaves:
        closure_counts[closure] = closure_counts.get(closure, 0) + 1
        if closure == "T":
            source = target_flint
            degree = target_degree
        elif closure.startswith("G"):
            gap_index = int(closure[1:])
            source = -gap_flint[gap_index]
            degree = tuple(map(int, gaps[gap_index].degree_list()))
        else:
            gap_index = int(closure[1:])
            multiplier = Fraction(*multiplier_pair)
            source = target_flint - fmpq(
                multiplier.numerator, multiplier.denominator
            ) * gap_flint[gap_index]
            degree = target_degree
            multiplier_records.append(
                {
                    "path": path,
                    "gap": gap_index,
                    "multiplier": _fraction_record(multiplier),
                }
            )
        record = _local_bernstein_sign_record(
            source, variables, path, degree
        )
        leaf_digest.update(
            (
                f"{path}:{closure}:{multiplier_pair}:"
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
        "subdivision": {
            "rule": "bisect axis depth mod 4 at 1/2",
            **tree,
            "closure_counts": closure_counts,
        },
        "leaf_closures": (
            "T means T>=0; Gj means -G_j>=0 and excludes the ordered "
            "relative interior; Mj means T-mu*G_j>=0 with mu>=0"
        ),
        "exact_multipliers": multiplier_records,
        "bernstein_signs": {
            "basis": (
                "dyadic local tensor-Bernstein basis; coefficients are "
                "recorded after positive binomial scaling"
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


def verify_face32(derived: dict[str, object]) -> dict[str, object]:
    charts = derive_face32_charts(derived)
    return {
        "charts": charts,
        "t_le_u_audit": verify_chart(charts["t_le_u"], T_LE_U_LEAVES),
        "a_le_z_audit": verify_chart(charts["a_le_z"], A_LE_Z_LEAVES),
        "z_le_a_audit": verify_chart(charts["z_le_a"], Z_LE_A_LEAVES),
    }


def face32_certificate(
    derived: dict[str, object],
    verified: dict[str, object],
) -> dict[str, object]:
    charts = verified["charts"]
    return {
        "face_type": "32",
        "normalized_shape": "(-1,-u,-v,t*a,t,t)",
        "base_variables": [str(value) for value in derived["cube_variables"]],
        "base_parameterization": "u=U, v=U*V, a=A, physical scale t=q",
        "hypotheses_used": "all four adjacent centroid gaps G_0,...,G_3 are nonnegative",
        "comparison_identity": (
            "I_{gamma_0}(1,5)-p is a positive rational factor times T"
        ),
        "gap_identity": (
            "gamma_{i+1}-gamma_i is a positive collision factor times "
            "G_i over a coefficientwise-positive denominator"
        ),
        "scale_bound": (
            "G_3>=0 and the exact Bernstein/power expansion of -G_3 at "
            "q=1+rho imply q<1 on regular ordered sections"
        ),
        "quadratic_bound": (
            "in the chart U=q*Z, G_3>=0 and the exact Bernstein expansion "
            "of -G_3 at q=3*Z^2+rho imply q<=3*Z^2"
        ),
        "projective_cover": (
            "split q<=U versus U<=q; in the latter write q=3*Z^2*R "
            "and split A=Z*B versus Z=A*B"
        ),
        "removed_positive_monomials": {
            "initial_target": list(derived["initial_target_monomial"]),
            "t_le_u_target": list(charts["t_le_u"]["target_monomial"]),
            "t_le_u_gaps": [list(value) for value in charts["t_le_u"]["gap_monomials"]],
            "u_le_t": {
                "target": list(charts["u_le_t_monomials"]["target"]),
                "gaps": [list(value) for value in charts["u_le_t_monomials"]["gaps"]],
            },
            "q_equals_3_Z_squared_R": {
                "target": list(charts["qz2_monomials"]["target"]),
                "gaps": [list(value) for value in charts["qz2_monomials"]["gaps"]],
            },
            "A_le_Z": {
                "target": list(charts["a_le_z"]["target_monomial"]),
                "gaps": [list(value) for value in charts["a_le_z"]["gap_monomials"]],
            },
            "Z_le_A": {
                "target": list(charts["z_le_a"]["target_monomial"]),
                "gaps": [list(value) for value in charts["z_le_a"]["gap_monomials"]],
            },
        },
        "polynomials": {
            "target_T": _polynomial_record(derived["target"]),
            "centroid_gaps": [
                _polynomial_record(gap) for gap in derived["gaps"]
            ],
            "shifted_minus_G3": _polynomial_record(
                derived["shifted_minus_last_gap"]
            ),
            "quadratic_confinement_minus_G3": _polynomial_record(
                charts["quadratic_confinement"]
            ),
            "t_le_u_target": _polynomial_record(charts["t_le_u"]["target"]),
            "A_le_Z_target": _polynomial_record(charts["a_le_z"]["target"]),
            "Z_le_A_target": _polynomial_record(charts["z_le_a"]["target"]),
        },
        "positive_gap_factors": list(derived["cofactor_statistics"]),
        "scale_confinement_signs": derived["scale_statistics"],
        "quadratic_confinement_signs": charts["quadratic_statistics"],
        "dyadic_charts": {
            "t_le_u": verified["t_le_u_audit"],
            "A_le_Z": verified["a_le_z_audit"],
            "Z_le_A": verified["z_le_a_audit"],
        },
    }
