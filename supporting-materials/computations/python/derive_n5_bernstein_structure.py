"""Derive the exact polynomial structure used by the ``n=5`` proof.

The normalized five-observation comparison has five uniform-simplex cap
regions.  This script derives the cleared residual polynomials, verifies all
boundary-factor divisions exactly, records a fixed four-simplex triangulation,
and proves the structural-zero Bernstein patterns by polynomial ideal
membership over ``QQ(b,c,d,e,f)``.

A face condition ``sum(i_j for j in S) <= q`` means that the homogenized
residual vanishes to order at least ``q+1`` on the corresponding face of the
four-simplex.  The script verifies this by reducing the residual modulo the
``(q+1)``-st power of the affine face ideal.  Thus the zero indices are exact
identities, not values inferred from floating-point tolerances.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp


REGION_SIMPLICES = {
    "A": [["v0", "h4_01", "h4_02", "h4_03", "h4_04"]],
    "B": [
        ["h4_04", "h3_12", "h4_01", "h4_03", "h4_02"],
        ["h3_13", "h4_04", "h3_12", "h4_01", "h3_14"],
        ["h3_13", "h4_04", "h3_12", "h4_01", "h4_03"],
    ],
    "C": [
        ["h2_23", "h3_04", "h3_02", "h3_12", "h3_03"],
        ["h2_23", "h3_04", "h3_02", "h3_12", "h2_24"],
        ["h2_23", "h3_13", "h3_04", "h3_12", "h3_03"],
        ["h3_14", "h2_23", "h3_04", "h3_12", "h2_24"],
        ["h3_14", "h2_23", "h3_13", "h3_04", "h3_12"],
    ],
    "D": [
        ["h2_13", "h2_04", "h2_23", "h1_34", "h2_03"],
        ["h2_14", "h2_13", "h2_04", "h2_23", "h1_34"],
        ["h2_24", "h2_14", "h2_04", "h2_23", "h1_34"],
    ],
    "E": [["h1_04", "h1_14", "h1_24", "h1_34", "v4"]],
}

# Bernstein indices are ordered (lambda_1,...,lambda_4,lambda_0), where
# lambda_0 belongs to the first listed simplex vertex.  Each condition is a
# subset of these five index positions and the largest permitted sum.
FACE_CONDITIONS = {
    "A": [[([0, 1, 2, 3], 0)]],
    "B": [
        [([0, 1], 0), ([0, 1, 3], 1), ([0, 1, 2, 3], 2)],
        [([1, 2, 3, 4], 2)],
        [([1, 2, 4], 1), ([1, 2, 3, 4], 2)],
    ],
    "C": [
        [([1, 2, 4], 1), ([1, 2, 3, 4], 3)],
        [([1, 2, 3, 4], 3)],
        [([2, 4], 0), ([0, 2, 3, 4], 3), ([0, 2, 4], 1)],
        [([0, 2, 3], 1), ([0, 2, 3, 4], 3)],
        [([0, 3], 0), ([0, 1, 3, 4], 3), ([0, 1, 3], 1)],
    ],
    "D": [
        [([1, 2, 3, 4], 2)],
        [([0, 2, 3], 1), ([0, 2, 3, 4], 2)],
        [([2, 3], 0), ([2, 3, 4], 1), ([0, 2, 3, 4], 2)],
    ],
    "E": [[([0, 1, 2, 4], 0)]],
}

EXPECTED_ZERO_COUNTS = {
    "A": [1],
    "B": [66, 15, 39],
    "C": [59, 35, 92, 59, 92],
    "D": [15, 39, 66],
    "E": [1],
}


def _symbolic_regions():
    x, y, z, w, b, c, d, e, f = sp.symbols("x y z w b c d e f")
    coordinates = (x, y, z, w)
    weights = (b, c, d, e, f)
    s = b * x + c * y + d * z + e * w + f
    common = (1 - x) * (1 - y) * (1 - z) * (1 - w)
    factor_w = w * (w - x) * (w - y) * (w - z)
    factor_z = z * (z - x) * (z - y)
    factor_y = y * (y - x)

    q1 = e + f
    q2 = d + e + f
    q3 = c + d + e + f
    g = b + c + d + e + f

    polynomial_a = (1 - f) ** 5 * common - (1 - s) ** 5

    numerator_b = (
        (1 - q1) ** 4 * (1 + 4 * q1) * common * factor_w
        - (1 - s) ** 5 * factor_w
        + (w - s) ** 5 * (1 - x) * (1 - y) * (1 - z)
    )
    quotient_b, remainder_b = sp.div(
        sp.Poly(numerator_b, w, z, y, x),
        sp.Poly(1 - w, w, z, y, x),
    )
    if not remainder_b.is_zero:
        raise AssertionError("region-B factorization failed")

    numerator_c = (
        (1 - q2) ** 3 * (1 + 3 * q2 + 6 * q2 ** 2)
        * common * factor_w * factor_z
        - (1 - s) ** 5 * factor_w * factor_z
        + (w - s) ** 5 * (1 - x) * (1 - y) * (1 - z) * factor_z
        - (z - s) ** 5 * (1 - x) * (1 - y) * (1 - w)
        * w * (w - x) * (w - y)
    )
    quotient_c = sp.Poly(numerator_c, w, z, y, x)
    for factor in (1 - w, 1 - z, w - z):
        quotient_c, remainder = sp.div(
            quotient_c, sp.Poly(factor, w, z, y, x))
        if not remainder.is_zero:
            raise AssertionError("region-C factorization failed")

    numerator_d = (
        (1 - q3) ** 2 * (1 + 2 * q3 + 3 * q3 ** 2 + 4 * q3 ** 3)
        * common * factor_w * factor_z * factor_y
        - (1 - s) ** 5 * factor_w * factor_z * factor_y
        + (w - s) ** 5 * (1 - x) * (1 - y) * (1 - z)
        * factor_z * factor_y
        - (z - s) ** 5 * (1 - x) * (1 - y) * (1 - w)
        * w * (w - x) * (w - y) * factor_y
        + (y - s) ** 5 * (1 - x) * (1 - z) * (1 - w)
        * w * (w - x) * (w - z) * z * (z - x)
    )
    quotient_d = sp.Poly(numerator_d, w, z, y, x)
    for factor in (1 - w, 1 - z, 1 - y, w - z, w - y, z - y):
        quotient_d, remainder = sp.div(
            quotient_d, sp.Poly(factor, w, z, y, x))
        if not remainder.is_zero:
            raise AssertionError("region-D factorization failed")

    polynomial_e = s ** 5 - g ** 5 * x * y * z * w

    return coordinates, weights, s, {
        "A": {
            "degree": 5,
            "alpha_identity": "alpha=(1-f)^5",
            "extracted_nonnegative_factors": [],
            "polynomial": polynomial_a,
        },
        "B": {
            "degree": 8,
            "alpha_identity": "alpha=(1-e-f)^4*(1+4*(e+f))",
            "extracted_nonnegative_factors": ["1-w"],
            "polynomial": quotient_b.as_expr(),
        },
        "C": {
            "degree": 9,
            "alpha_identity": (
                "alpha=(1-d-e-f)^3*(1+3*(d+e+f)+6*(d+e+f)^2)"
            ),
            "extracted_nonnegative_factors": ["1-w", "1-z", "w-z"],
            "polynomial": quotient_c.as_expr(),
        },
        "D": {
            "degree": 8,
            "alpha_identity": (
                "alpha=(1-c-d-e-f)^2*(1+2*(c+d+e+f)+"
                "3*(c+d+e+f)^2+4*(c+d+e+f)^3)"
            ),
            "extracted_nonnegative_factors": [
                "1-w", "1-z", "1-y", "w-z", "w-y", "z-y"
            ],
            "polynomial": quotient_d.as_expr(),
        },
        "E": {
            "degree": 5,
            "alpha_identity": "1-alpha=(b+c+d+e+f)^5",
            "extracted_nonnegative_factors": [],
            "polynomial": polynomial_e,
        },
    }


def _ordered_simplex_vertices():
    return [
        tuple([sp.Integer(0)] * (4 - index)
              + [sp.Integer(1)] * index)
        for index in range(5)
    ]


def _vertices(coordinates, weights, s):
    b, c, d, e, f = weights
    del b, c, d, e
    base = _ordered_simplex_vertices()
    output = {"v0": base[0], "v4": base[4]}
    needed = {name for simplices in REGION_SIMPLICES.values()
              for simplex in simplices for name in simplex}
    for name in sorted(needed):
        if name in output:
            continue
        region = int(name[1])
        left, right = map(int, name.split("_")[1])
        boundary = s - coordinates[region - 1]
        left_value = sp.expand(boundary.subs(dict(zip(coordinates, base[left]))))
        right_value = sp.expand(boundary.subs(dict(zip(coordinates, base[right]))))
        parameter = sp.cancel(left_value / (left_value - right_value))
        output[name] = tuple(
            sp.cancel(base[left][j]
                      + parameter * (base[right][j] - base[left][j]))
            for j in range(4)
        )
    return output


def _bernstein_indices(degree):
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            for k in range(degree + 1 - i - j):
                for ell in range(degree + 1 - i - j - k):
                    yield (i, j, k, ell, degree - i - j - k - ell)


def _structural_zero_indices(degree, conditions):
    return sorted(
        index for index in _bernstein_indices(degree)
        if any(sum(index[position] for position in subset) <= maximum
               for subset, maximum in conditions)
    )


def _shared_face_boundaries(simplex, condition, vertices,
                            coordinates, s):
    subset, _ = condition
    # Convert (lambda_1,...,lambda_4,lambda_0) positions to the natural
    # listed-vertex order (lambda_0,...,lambda_4).
    position_map = {0: 1, 1: 2, 2: 3, 3: 4, 4: 0}
    excluded = {position_map[position] for position in subset}
    face_names = [name for index, name in enumerate(simplex)
                  if index not in excluded]
    x, y, z, w = coordinates
    boundaries = {
        "x=0": x,
        "x=y": y - x,
        "y=z": z - y,
        "z=w": w - z,
        "w=1": 1 - w,
        "s=x": s - x,
        "s=y": s - y,
        "s=z": s - z,
        "s=w": s - w,
    }
    shared = []
    for label, expression in boundaries.items():
        if all(sp.cancel(expression.subs(dict(zip(
                coordinates, vertices[name])))) == 0
               for name in face_names):
            shared.append((label, expression))
    return face_names, shared


def _verify_face_condition(region_name, region_poly, simplex, condition,
                           vertices, coordinates, weights, s, cache):
    face_names, shared = _shared_face_boundaries(
        simplex, condition, vertices, coordinates, s)
    subset, maximum = condition
    field = sp.QQ.frac_field(*weights)
    equations = [expression for _, expression in shared]
    # The common active boundaries must have the expected affine codimension.
    generic = {symbol: sp.Rational(index + 1, 17)
               for index, symbol in enumerate(weights)}
    coefficient_rows = []
    for expression in equations:
        specialized = sp.Poly(expression.subs(generic), *coordinates)
        coefficient_rows.append([
            specialized.coeff_monomial(variable) for variable in coordinates
        ])
    # Select an independent set of affine normals.  The other active
    # boundaries are redundant on the same face.  Using only a basis leaves
    # the face ideal unchanged and makes its higher powers far cheaper to
    # verify than a generating set containing every redundant equation.
    independent = []
    current_rows = []
    current_rank = 0
    for index, row in enumerate(coefficient_rows):
        candidate = current_rows + [row]
        rank = sp.Matrix(candidate).rank()
        if rank > current_rank:
            independent.append(index)
            current_rows = candidate
            current_rank = rank
    if current_rank != len(subset):
        raise AssertionError(
            f"face codimension mismatch: rank={current_rank}, "
            f"subset={subset}")
    equations = [equations[index] for index in independent]
    independent_labels = [shared[index][0] for index in independent]

    order = maximum + 1
    cache_key = (region_name, tuple(independent_labels), order)
    if cache_key not in cache:
        print(
            f"region {region_name}: exact face-ideal check "
            f"{independent_labels}, order {order}",
            flush=True,
        )
        generators = []
        for selection in itertools.combinations_with_replacement(
                range(len(equations)), order):
            generator = sp.Integer(1)
            for index in selection:
                generator *= equations[index]
            generators.append(generator)
        basis = sp.groebner(generators, *coordinates, domain=field)
        remainder = basis.reduce(region_poly.as_expr())[1]
        if remainder != 0:
            raise AssertionError(
                f"face-order check failed for {simplex}, {condition}")
        cache[cache_key] = len(generators)
    return {
        "bernstein_index_subset": subset,
        "maximum_subset_sum": maximum,
        "vanishing_order": order,
        "face_vertices": face_names,
        "active_boundaries": [label for label, _ in shared],
        "independent_ideal_generators": independent_labels,
        "ideal_power_generator_count": cache[cache_key],
    }


def derive(verify_faces=True):
    coordinates, weights, s, regions = _symbolic_regions()
    vertices = _vertices(coordinates, weights, s)
    output = {
        "schema_version": 1,
        "claim": (
            "Exact polynomial, four-simplex, and structural-zero structure "
            "for the n=5 Gaffke-domination certificate."
        ),
        "variables": [str(value) for value in weights],
        "coordinate_variables": [str(value) for value in coordinates],
        "vertex_formulas": {
            name: [str(sp.factor(value)) for value in point]
            for name, point in sorted(vertices.items())
        },
        "regions": {},
    }
    verification_cache = {}

    for region_name, region in regions.items():
        polynomial = sp.Poly(region["polynomial"], *coordinates)
        if polynomial.total_degree() != region["degree"]:
            raise AssertionError(f"unexpected degree in region {region_name}")
        polynomial_coefficients = [
            {
                "powers": list(powers),
                "expression": str(sp.factor(coefficient)),
            }
            for powers, coefficient in polynomial.terms()
        ]
        simplex_records = []
        for simplex_index, simplex in enumerate(
                REGION_SIMPLICES[region_name]):
            conditions = FACE_CONDITIONS[region_name][simplex_index]
            zero_indices = _structural_zero_indices(
                region["degree"], conditions)
            expected = EXPECTED_ZERO_COUNTS[region_name][simplex_index]
            if len(zero_indices) != expected:
                raise AssertionError(
                    f"region {region_name} simplex {simplex_index}: "
                    f"expected {expected} zeros, found {len(zero_indices)}")
            face_proofs = []
            if verify_faces:
                for condition in conditions:
                    face_proofs.append(_verify_face_condition(
                        region_name, polynomial, simplex, condition, vertices,
                        coordinates, weights, s, verification_cache))
            else:
                face_proofs = [
                    {
                        "bernstein_index_subset": subset,
                        "maximum_subset_sum": maximum,
                        "verification_skipped": True,
                    }
                    for subset, maximum in conditions
                ]
            simplex_records.append({
                "vertices": simplex,
                "degree": region["degree"],
                "coefficient_count": int(sp.binomial(
                    region["degree"] + 4, 4)),
                "structural_zero_count": len(zero_indices),
                "structural_zero_indices": [list(value)
                                            for value in zero_indices],
                "face_order_proofs": face_proofs,
            })
            print(
                f"region {region_name}, four-simplex {simplex_index}: "
                f"{len(zero_indices)} structural zeros"
            )
        output["regions"][region_name] = {
            "degree": region["degree"],
            "alpha_identity": region["alpha_identity"],
            "extracted_nonnegative_factors":
                region["extracted_nonnegative_factors"],
            "polynomial_coefficients": polynomial_coefficients,
            "simplices": simplex_records,
        }
    output["face_order_verification"] = (
        "exact_polynomial_ideal_membership" if verify_faces
        else "skipped_for_development"
    )
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--skip-face-proof", action="store_true",
        help="development only: derive formulas without ideal checks",
    )
    args = parser.parse_args(argv)
    args.out.write_text(json.dumps(
        derive(verify_faces=not args.skip_face_proof), indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
