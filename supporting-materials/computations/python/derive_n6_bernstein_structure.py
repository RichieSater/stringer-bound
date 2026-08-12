"""Derive the exact polynomial structure used by the ``n=6`` proof.

The normalized six-observation comparison has six uniform-simplex cap
regions.  This script derives their cleared residual polynomials, records a
fixed five-simplex triangulation, and proves every structural-zero Bernstein
pattern by exact polynomial ideal membership over ``QQ(b,c,d,e,f,g)``.

The ideal-power checks use the characteristic-zero differential criterion for
an affine linear face ideal ``I``: a polynomial belongs to ``I**q`` if and
only if every ordinary partial derivative of total order below ``q`` belongs
to ``I``.  Singular reduces those derivatives modulo the generic face ideal
over the rational-function field ``QQ(b,c,d,e,f,g)``.  This avoids the much
larger standard basis of ``I**q`` without weakening the exact check.  A
generic rational point is used only to select an independent subset of
redundant face equations; no weight specialization or floating-point
tolerance is used to establish a polynomial identity or zero.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import itertools
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import sympy as sp


REGION_NAMES = "ABCDEF"

REGION_SIMPLICES = {
    "A": [["v0", "h1_05", "h2_04", "h3_03", "h4_02", "h5_01"]],
    "B": [
        ["h1_05", "h5_01", "v1", "h4_02", "h3_03", "h2_04"],
        ["h1_05", "h4_12", "v1", "h4_02", "h3_03", "h2_04"],
        ["h1_05", "h4_12", "h3_13", "v1", "h3_03", "h2_04"],
        ["h2_14", "h1_05", "h4_12", "h3_13", "v1", "h2_04"],
        ["h2_14", "h1_05", "h4_12", "h3_13", "v1", "h1_15"],
    ],
    "C": [
        ["h2_04", "h2_14", "h1_05", "v2", "h3_23", "h2_24"],
        ["h2_04", "h3_03", "h1_05", "v2", "h4_02", "h4_12"],
        ["h2_04", "h3_03", "h3_13", "h1_05", "v2", "h3_23"],
        ["h2_04", "h2_14", "h3_13", "h1_05", "v2", "h3_23"],
        ["h2_04", "h3_03", "h3_13", "h1_05", "v2", "h4_12"],
        ["h2_04", "h2_14", "h3_13", "h1_05", "v2", "h4_12"],
        ["h1_15", "h2_14", "h3_13", "h1_05", "v2", "h3_23"],
        ["h1_15", "h1_05", "v2", "h3_23", "h2_24", "h1_25"],
        ["h1_15", "h2_14", "h1_05", "v2", "h3_23", "h2_24"],
        ["h1_15", "h2_14", "h3_13", "h1_05", "v2", "h4_12"],
    ],
    "D": [
        ["h2_04", "h1_05", "v3", "h2_14", "h2_24", "h3_23"],
        ["h2_04", "h3_13", "h1_05", "v3", "h3_03", "h3_23"],
        ["h2_04", "h3_13", "h1_05", "v3", "h2_14", "h3_23"],
        ["h2_34", "h2_04", "h1_05", "v3", "h2_14", "h2_24"],
        ["h1_15", "h1_25", "h2_34", "h1_05", "v3", "h1_35"],
        ["h1_15", "h3_13", "h1_05", "v3", "h2_14", "h3_23"],
        ["h1_15", "h2_34", "h1_05", "v3", "h2_14", "h2_24"],
        ["h1_15", "h1_25", "h2_34", "h1_05", "v3", "h2_24"],
        ["h1_15", "h1_05", "v3", "h2_14", "h2_24", "h3_23"],
        ["h1_15", "h1_25", "h1_05", "v3", "h2_24", "h3_23"],
    ],
    "E": [
        ["h1_35", "h1_05", "h1_15", "h1_25", "h1_45", "v4"],
        ["h1_35", "h1_05", "h1_15", "h1_25", "h2_34", "v4"],
        ["h2_24", "h2_14", "h1_05", "h2_34", "h2_04", "v4"],
        ["h2_24", "h1_05", "h1_15", "h1_25", "h2_34", "v4"],
        ["h2_24", "h2_14", "h1_05", "h1_15", "h2_34", "v4"],
    ],
    "F": [["v5", "h1_05", "h1_15", "h1_25", "h1_35", "h1_45"]],
}


def symbolic_regions():
    coordinates = sp.symbols("x y z w u")
    weights = sp.symbols("b c d e f g")
    s = sum(weight * coordinate for weight, coordinate
            in zip(weights[:5], coordinates)) + weights[5]
    common = sp.prod(1 - coordinate for coordinate in coordinates)
    node_factors = [
        sp.expand(coordinate * sp.prod(
            coordinate - coordinates[j] for j in range(i)))
        for i, coordinate in enumerate(coordinates)
    ]
    expected_degrees = (6, 10, 12, 12, 10)
    alpha_identities = (
        "alpha=(1-g)^6",
        "alpha=(1-f-g)^5*(1+5*(f+g))",
        "alpha=(1-e-f-g)^4*(1+4*(e+f+g)+10*(e+f+g)^2)",
        "alpha=(1-d-e-f-g)^3*(1+3*(d+e+f+g)+"
        "6*(d+e+f+g)^2+10*(d+e+f+g)^3)",
        "alpha=(1-c-d-e-f-g)^2*(1+2*(c+d+e+f+g)+"
        "3*(c+d+e+f+g)^2+4*(c+d+e+f+g)^3+"
        "5*(c+d+e+f+g)^4)",
    )

    regions = {}
    for k, name in enumerate(REGION_NAMES[:5]):
        upper = list(range(5 - k, 5))
        q = sum(weights[5 - k:])
        alpha_expression = sum(
            sp.binomial(6, j) * q ** j * (1 - q) ** (6 - j)
            for j in range(k + 1)
        )
        denominator = common * sp.prod(node_factors[i] for i in upper)
        numerator = (
            alpha_expression * denominator
            - (1 - s) ** 6 * sp.prod(node_factors[i] for i in upper)
        )
        for i in upper:
            ratio = sp.prod(
                1 - coordinates[j] for j in range(5) if j != i)
            for j in upper:
                if j == i:
                    continue
                factor = node_factors[j]
                if j > i:
                    factor = sp.cancel(
                        factor / (coordinates[j] - coordinates[i]))
                ratio *= factor
            numerator -= ((-1) ** (5 - i)) * (
                coordinates[i] - s) ** 6 * ratio

        # Expanding with the weights treated as an expression domain is very
        # slow at n=6.  Work instead in the exact eleven-variable ring over
        # QQ, divide out the coordinate-only factors there, and regroup the
        # resulting monomials by their five coordinate exponents.
        all_variables = coordinates + weights
        full_polynomial = sp.Poly(
            numerator, *all_variables, domain=sp.QQ)
        factor_expressions = (
            [1 - coordinates[i] for i in upper]
            + [coordinates[j] - coordinates[i]
               for offset, i in enumerate(upper)
               for j in upper[offset + 1:]]
        )
        for factor_expression in factor_expressions:
            full_polynomial, remainder = sp.div(
                full_polynomial,
                sp.Poly(factor_expression, *all_variables, domain=sp.QQ),
            )
            if not remainder.is_zero:
                raise AssertionError(f"region-{name} factorization failed")
        grouped = {}
        for powers, coefficient in full_polynomial.terms():
            coordinate_powers = powers[:5]
            weight_powers = powers[5:]
            term = coefficient * sp.prod(
                weight ** exponent
                for weight, exponent in zip(weights, weight_powers))
            grouped[coordinate_powers] = (
                grouped.get(coordinate_powers, sp.Integer(0)) + term)
        polynomial = sp.Poly.from_dict(
            grouped, coordinates, domain=sp.EX)
        factors = [str(value) for value in factor_expressions]
        if polynomial.total_degree() != expected_degrees[k]:
            raise AssertionError(f"unexpected degree in region {name}")
        regions[name] = {
            "degree": expected_degrees[k],
            "alpha_identity": alpha_identities[k],
            "extracted_nonnegative_factors": factors,
            "polynomial": polynomial.as_expr(),
        }
        print(
            f"derived region {name}: degree {polynomial.total_degree()}, "
            f"{len(polynomial.terms())} coordinate terms",
            flush=True,
        )

    total = sum(weights)
    final_expression = s ** 6 - total ** 6 * sp.prod(coordinates)
    full_final = sp.Poly(
        final_expression, *(coordinates + weights), domain=sp.QQ)
    grouped_final = {}
    for powers, coefficient in full_final.terms():
        coordinate_powers = powers[:5]
        weight_powers = powers[5:]
        term = coefficient * sp.prod(
            weight ** exponent
            for weight, exponent in zip(weights, weight_powers))
        grouped_final[coordinate_powers] = (
            grouped_final.get(coordinate_powers, sp.Integer(0)) + term)
    final_polynomial = sp.Poly.from_dict(
        grouped_final, coordinates, domain=sp.EX)
    regions["F"] = {
        "degree": 6,
        "alpha_identity": "1-alpha=(b+c+d+e+f+g)^6",
        "extracted_nonnegative_factors": [],
        "polynomial": final_polynomial.as_expr(),
    }
    print("derived region F: degree 6", flush=True)
    return coordinates, weights, s, regions


def ordered_simplex_vertices():
    return [tuple([sp.Integer(0)] * (5 - i) + [sp.Integer(1)] * i)
            for i in range(6)]


def vertices(coordinates, weights, s):
    base = ordered_simplex_vertices()
    output = {f"v{i}": point for i, point in enumerate(base)}
    needed = {name for simplices in REGION_SIMPLICES.values()
              for simplex in simplices for name in simplex}
    for name in sorted(needed):
        if name in output:
            continue
        region = int(name[1])
        left, right = map(int, name.split("_")[1])
        boundary = s - coordinates[region - 1]
        left_value = sp.expand(
            boundary.subs(dict(zip(coordinates, base[left]))))
        right_value = sp.expand(
            boundary.subs(dict(zip(coordinates, base[right]))))
        parameter = sp.cancel(left_value / (left_value - right_value))
        output[name] = tuple(
            sp.cancel(base[left][j] + parameter * (
                base[right][j] - base[left][j]))
            for j in range(5)
        )
    return output


def vertex_boundary_identities(verts, coordinates, s):
    """Record every ordered-domain or threshold boundary met identically."""

    x, y, z, w, u = coordinates
    boundaries = {
        "x=0": x,
        "x=y": y - x,
        "y=z": z - y,
        "z=w": w - z,
        "w=u": u - w,
        "u=1": 1 - u,
        "s=x": s - x,
        "s=y": s - y,
        "s=z": s - z,
        "s=w": s - w,
        "s=u": s - u,
    }
    return {
        name: [
            label for label, expression in boundaries.items()
            if sp.cancel(expression.subs(dict(zip(coordinates, point)))) == 0
        ]
        for name, point in sorted(verts.items())
    }


def bernstein_indices(degree):
    for prefix in itertools.product(range(degree + 1), repeat=5):
        if sum(prefix) <= degree:
            yield prefix + (degree - sum(prefix),)


def structural_zero_indices(degree, conditions):
    return sorted(index for index in bernstein_indices(degree)
                  if any(sum(index[position] for position in subset) <= maximum
                         for subset, maximum in conditions))


def shared_face_boundaries(simplex, condition, verts, coordinates, s):
    subset, _ = condition
    # Bernstein indices are (lambda_1,...,lambda_5,lambda_0).
    position_map = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 0}
    excluded = {position_map[position] for position in subset}
    face_names = [name for index, name in enumerate(simplex)
                  if index not in excluded]
    x, y, z, w, u = coordinates
    boundaries = {
        "x=0": x,
        "x=y": y - x,
        "y=z": z - y,
        "z=w": w - z,
        "w=u": u - w,
        "u=1": 1 - u,
        "s=x": s - x,
        "s=y": s - y,
        "s=z": s - z,
        "s=w": s - w,
        "s=u": s - u,
    }
    shared = []
    for label, expression in boundaries.items():
        if all(sp.cancel(expression.subs(dict(zip(
                coordinates, verts[name])))) == 0 for name in face_names):
            shared.append((label, expression))
    return face_names, shared


def independent_face_equations(simplex, condition, verts, coordinates,
                               weights, s):
    face_names, shared = shared_face_boundaries(
        simplex, condition, verts, coordinates, s)
    subset, _ = condition
    # This point selects a basis from redundant linear normals only.  A
    # nonzero minor here is a nonzero rational function, and the subsequent
    # membership calculation is performed generically over QQ(weights).
    generic = {symbol: sp.Rational(i + 1, 23)
               for i, symbol in enumerate(weights)}
    independent = []
    rows = []
    rank = 0
    for index, (_, expression) in enumerate(shared):
        specialized = sp.Poly(expression.subs(generic), *coordinates)
        row = [specialized.coeff_monomial(variable)
               for variable in coordinates]
        new_rank = sp.Matrix(rows + [row]).rank()
        if new_rank > rank:
            independent.append(index)
            rows.append(row)
            rank = new_rank
    if rank != len(subset):
        return None
    return (face_names,
            [shared[i][0] for i in independent],
            [shared[i][1] for i in independent],
            [label for label, _ in shared])


def verify_face_condition(region_name, region_poly, simplex, condition,
                          verts, coordinates, weights, s, cache=None):
    cache = {} if cache is None else cache
    result = independent_face_equations(
        simplex, condition, verts, coordinates, weights, s)
    if result is None:
        raise AssertionError(f"face codimension mismatch: {condition}")
    face_names, labels, equations, all_labels = result
    order = condition[1] + 1
    key = (region_name, tuple(labels), order)
    if key not in cache:
        generators = []
        for selection in itertools.combinations_with_replacement(
                range(len(equations)), order):
            generator = sp.Integer(1)
            for index in selection:
                generator *= equations[index]
            generators.append(generator)
        field = sp.QQ.frac_field(*weights)
        basis = sp.groebner(generators, *coordinates, domain=field)
        remainder = basis.reduce(region_poly.as_expr())[1]
        cache[key] = (remainder == 0, len(generators))
    valid, generator_count = cache[key]
    return valid, {
        "bernstein_index_subset": list(condition[0]),
        "maximum_subset_sum": condition[1],
        "vanishing_order": order,
        "face_vertices": face_names,
        "active_boundaries": all_labels,
        "independent_ideal_generators": labels,
        "ideal_power_generator_count": generator_count,
    }
FACE_CONDITIONS = {'A': [[((0, 1, 2, 3, 4), 0)]],
 'B': [[((0, 1), 0), ((0, 1, 2), 1), ((0, 1, 2, 3), 2), ((0, 1, 2, 3, 4), 3), ((0, 2, 3, 4, 5), 0)],
       [((0, 1), 0), ((0, 1, 2), 1), ((0, 1, 2, 3), 2), ((0, 1, 2, 3, 4), 3), ((0, 2, 3, 4, 5), 0)],
       [((0, 1, 2), 1), ((0, 1, 2, 3), 2), ((0, 1, 2, 3, 4), 3), ((0, 1, 3, 4, 5), 0)],
       [((1, 2, 3, 5), 2), ((1, 2, 3, 4, 5), 3), ((0, 1, 2, 4, 5), 0)],
       [((1, 2, 3, 4, 5), 3), ((0, 1, 2, 4, 5), 0)]],
 'C': [[((0, 2, 3, 4), 3), ((2, 3, 4), 1), ((0, 2, 3, 4, 5), 5), ((0, 1, 3, 4, 5), 0)],
       [((0, 2, 3, 4), 3), ((2, 3, 4), 1), ((0, 2, 3, 4, 5), 5), ((0, 1, 3, 4, 5), 0)],
       [((3, 4), 0), ((0, 1, 3, 4), 3), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 2, 4, 5), 0)],
       [((3, 4), 0), ((0, 1, 3, 4), 3), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 2, 4, 5), 0)],
       [((3, 4), 0), ((0, 1, 3, 4), 3), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 2, 4, 5), 0)],
       [((3, 4), 0), ((0, 1, 3, 4), 3), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 2, 4, 5), 0)],
       [((3, 4), 0), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 3, 4), 2), ((0, 1, 2, 4, 5), 0)],
       [((1, 2, 3, 4, 5), 5), ((1, 2, 3, 4), 2), ((0, 2, 3, 4, 5), 0)],
       [((2, 3, 4), 1), ((0, 2, 3, 4, 5), 5), ((0, 2, 3, 4), 2), ((0, 1, 3, 4, 5), 0)],
       [((3, 4), 0),
        ((1, 3, 4), 1),
        ((0, 1, 3, 4, 5), 5),
        ((0, 1, 3, 4), 2),
        ((0, 1, 2, 4, 5), 0)]],
 'D': [[((1, 4), 0), ((1, 3, 4), 1), ((1, 2, 3, 4, 5), 5), ((1, 2, 3, 4), 2), ((0, 2, 3, 4, 5), 0)],
       [((0, 2, 3, 4, 5), 5), ((0, 2, 3, 4), 2), ((0, 1, 3, 4, 5), 0)],
       [((0, 2, 4), 1), ((0, 2, 3, 4, 5), 5), ((0, 2, 3, 4), 2), ((0, 1, 3, 4, 5), 0)],
       [((2, 5), 0), ((2, 4, 5), 1), ((0, 2, 3, 4, 5), 5), ((2, 3, 4, 5), 2), ((0, 1, 3, 4, 5), 0)],
       [((0, 1, 3, 4), 3), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 2, 4, 5), 0)],
       [((0, 2, 3, 4), 3), ((0, 2, 4), 1), ((0, 2, 3, 4, 5), 5), ((0, 1, 3, 4, 5), 0)],
       [((0, 2), 0), ((0, 2, 3, 4), 3), ((0, 2, 4), 1), ((0, 2, 3, 4, 5), 5), ((0, 1, 3, 4, 5), 0)],
       [((1, 3), 0), ((0, 1, 3, 4), 3), ((1, 3, 4), 1), ((0, 1, 3, 4, 5), 5), ((0, 1, 2, 4, 5), 0)],
       [((1, 4), 0), ((1, 2, 3, 4), 3), ((1, 3, 4), 1), ((1, 2, 3, 4, 5), 5), ((0, 2, 3, 4, 5), 0)],
       [((2, 4), 0),
        ((0, 2, 3, 4), 3),
        ((2, 3, 4), 1),
        ((0, 2, 3, 4, 5), 5),
        ((0, 1, 3, 4, 5), 0)]],
 'E': [[((3, 4), 0), ((3, 4, 5), 1), ((2, 3, 4, 5), 2), ((1, 2, 3, 4, 5), 3), ((0, 1, 2, 3, 5), 0)],
       [((3, 4), 0), ((3, 4, 5), 1), ((2, 3, 4, 5), 2), ((1, 2, 3, 4, 5), 3), ((0, 1, 2, 3, 5), 0)],
       [((0, 2, 3, 4, 5), 3), ((0, 1, 2, 3, 5), 0)],
       [((3, 4, 5), 1), ((2, 3, 4, 5), 2), ((1, 2, 3, 4, 5), 3), ((0, 1, 2, 3, 5), 0)],
       [((0, 3, 4, 5), 2), ((0, 2, 3, 4, 5), 3), ((0, 1, 2, 3, 5), 0)]],
 'F': [[((0, 1, 2, 3, 4), 0)]]}
EXPECTED_ZERO_COUNTS = {'A': [1],
 'B': [456, 456, 302, 162, 57],
 'C': [666, 666, 911, 911, 911, 911, 813, 358, 554, 813],
 'D': [813, 358, 554, 813, 666, 666, 911, 911, 911, 911],
 'E': [456, 456, 57, 302, 162],
 'F': [1]}


def _singular_polynomial(polynomial, coordinates):
    names = [str(value) for value in coordinates]
    terms = []
    for powers, coefficient in polynomial.terms():
        coefficient_text = str(coefficient).replace("**", "^")
        monomial = "*".join(
            f"{name}^{exponent}" if exponent != 1 else name
            for name, exponent in zip(names, powers) if exponent
        )
        terms.append(
            f"({coefficient_text})" + (f"*{monomial}" if monomial else ""))
    return "+".join(terms)


def _singular_boundary(label):
    linear = "b*x+c*y+d*z+e*w+f*u+g"
    return {
        "x=0": "x", "x=y": "y-x", "y=z": "z-y",
        "z=w": "w-z", "w=u": "u-w", "u=1": "1-u",
        "s=x": f"({linear})-x", "s=y": f"({linear})-y",
        "s=z": f"({linear})-z", "s=w": f"({linear})-w",
        "s=u": f"({linear})-u",
    }[label]


def _run_singular_check(task):
    script_path, result_path, key = task
    completed = subprocess.run(
        ["Singular", "-q", str(script_path)],
        check=False, capture_output=True, text=True,
    )
    rendered = completed.stdout + completed.stderr
    result_path.write_text(rendered)
    if completed.returncode != 0 or rendered.strip() != "PASS":
        raise AssertionError(
            f"Singular face-ideal check failed for {key}: {rendered}")
    return key


def _singular_differential_checks(order):
    """Render an exact check of every derivative below ``order``.

    Nondecreasing variable-index tuples visit each derivative multi-index
    once.  Each derivative is constructed directly from ``P``, reduced, and
    killed before the next one.  Recomputing the short differentiation chain
    is substantially faster here than retaining a tree of expanded generic
    rational-function polynomials.
    """

    coordinates = ("x", "y", "z", "w", "u")
    lines = ["int failures=0;"]
    count = 0
    for degree in range(order):
        for path in itertools.combinations_with_replacement(
                range(len(coordinates)), degree):
            child = "D" + ("".join(str(value) for value in path)
                           if path else "root")
            expression = "P"
            for index in path:
                expression = f"diff({expression},{coordinates[index]})"
            lines.append(f"poly {child}={expression};")
            lines.append(
                f'if (reduce({child},G)!=0) {{ "FAIL {path}"; '
                "failures=failures+1; }")
            lines.append(f"kill {child};")
            count += 1
    lines.append('if (failures==0) { "PASS"; }')
    return "\n".join(lines), count


def _verify_faces_with_singular(regions, tasks):
    if shutil.which("Singular") is None:
        raise RuntimeError(
            "Singular is required for the n=6 exact face-ideal checks")
    unique = {}
    for task in tasks:
        unique.setdefault(task[0], task)
    verified = set()
    with tempfile.TemporaryDirectory(prefix="stringer-n6-") as directory:
        directory = Path(directory)
        run_tasks = []
        polynomial_text = {
            name: _singular_polynomial(
                sp.Poly(record["polynomial"], *record["coordinates"]),
                record["coordinates"],
            )
            for name, record in regions.items()
        }
        for index, (key, task) in enumerate(sorted(unique.items())):
            region_name, labels, order = key
            script_path = directory / f"{index:02d}.sing"
            result_path = directory / f"{index:02d}.result"
            exact_checks, condition_count = (
                _singular_differential_checks(order))
            expected_count = math.comb(5 + order - 1, 5)
            if condition_count != expected_count:
                raise AssertionError(
                    f"wrong derivative count for {key}: "
                    f"{condition_count} != {expected_count}")
            boundaries = [_singular_boundary(label) for label in labels]
            script_path.write_text(
                "option(redSB);\n"
                "ring r=(0,b,c,d,e,f,g),(x,y,z,w,u),lp;\n"
                f"poly P={polynomial_text[region_name]};\n"
                "ideal I=" + ",".join(boundaries) + ";\n"
                "ideal G=std(I);\n"
                f"{exact_checks}\n"
                "quit;\n"
            )
            run_tasks.append((script_path, result_path, key))
        # Generic rational-function arithmetic is CPU intensive.  A modest
        # worker cap keeps memory use predictable on both laptops and CI.
        workers = min(len(run_tasks), os.cpu_count() or 1, 6)
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=workers) as executor:
            for key in executor.map(_run_singular_check, run_tasks):
                verified.add(key)
                print(f"exact generic Singular face check passed: {key}",
                      flush=True)
    expected = set(unique)
    if verified != expected:
        raise AssertionError("not every n=6 face condition was verified")
    return verified


def derive(verify_faces=True):
    coordinates, weights, s, symbolic = symbolic_regions()
    verts = vertices(coordinates, weights, s)
    # Retain the coordinates in a private field for the Singular renderer.
    singular_regions = {
        name: {**record, "coordinates": coordinates}
        for name, record in symbolic.items()
    }
    face_tasks = []
    condition_metadata = {}
    for region_name in REGION_NAMES:
        polynomial = sp.Poly(symbolic[region_name]["polynomial"], *coordinates)
        for simplex_index, simplex in enumerate(
                REGION_SIMPLICES[region_name]):
            for condition in FACE_CONDITIONS[region_name][simplex_index]:
                result = independent_face_equations(
                    simplex, condition, verts, coordinates, weights, s)
                if result is None:
                    raise AssertionError(
                        f"face codimension mismatch: {region_name}, "
                        f"{simplex_index}, {condition}")
                face_names, independent_labels, equations, active_labels = result
                key = (region_name, tuple(independent_labels), condition[1] + 1)
                face_tasks.append((key, region_name, simplex_index, condition))
                condition_metadata[(region_name, simplex_index, condition)] = {
                    "bernstein_index_subset": list(condition[0]),
                    "maximum_subset_sum": condition[1],
                    "vanishing_order": condition[1] + 1,
                    "face_vertices": face_names,
                    "active_boundaries": active_labels,
                    "independent_ideal_generators": independent_labels,
                    "ideal_power_generator_count": int(sp.binomial(
                        len(equations) + condition[1], condition[1] + 1)),
                }
    verified = (_verify_faces_with_singular(singular_regions, face_tasks)
                if verify_faces else set())

    output = {
        "schema_version": 1,
        "claim": (
            "Exact polynomial, five-simplex, and structural-zero structure "
            "for the n=6 Gaffke-domination certificate."
        ),
        "variables": [str(value) for value in weights],
        "coordinate_variables": [str(value) for value in coordinates],
        "vertex_formulas": {
            name: [str(sp.factor(value)) for value in point]
            for name, point in sorted(verts.items())
        },
        "vertex_boundary_identities": vertex_boundary_identities(
            verts, coordinates, s),
        "regions": {},
        "face_order_verification": (
            "exact_generic_differential_ideal_membership"
            if verify_faces else "skipped_for_development"
        ),
    }
    for region_name in REGION_NAMES:
        record = symbolic[region_name]
        polynomial = sp.Poly(record["polynomial"], *coordinates)
        simplex_records = []
        for simplex_index, simplex in enumerate(
                REGION_SIMPLICES[region_name]):
            conditions = FACE_CONDITIONS[region_name][simplex_index]
            zero_indices = structural_zero_indices(
                record["degree"], conditions)
            expected = EXPECTED_ZERO_COUNTS[region_name][simplex_index]
            if len(zero_indices) != expected:
                raise AssertionError(
                    f"region {region_name}, simplex {simplex_index}: "
                    f"expected {expected} zeros, found {len(zero_indices)}")
            proofs = []
            for condition in conditions:
                metadata = dict(condition_metadata[
                    (region_name, simplex_index, condition)])
                if verify_faces:
                    key = (region_name,
                           tuple(metadata["independent_ideal_generators"]),
                           metadata["vanishing_order"])
                    if key not in verified:
                        raise AssertionError(f"unverified face condition {key}")
                    metadata["verification"] = (
                        "exact_generic_Singular_differential_ideal_"
                        "membership_over_QQ(b,c,d,e,f,g)")
                    metadata["derivative_reductions_checked"] = int(
                        sp.binomial(5 + metadata["vanishing_order"] - 1, 5))
                else:
                    metadata["verification_skipped"] = True
                proofs.append(metadata)
            simplex_records.append({
                "vertices": simplex,
                "degree": record["degree"],
                "coefficient_count": math.comb(record["degree"] + 5, 5),
                "structural_zero_count": len(zero_indices),
                "structural_zero_indices": [list(value)
                                            for value in zero_indices],
                "face_order_proofs": proofs,
            })
        output["regions"][region_name] = {
            "degree": record["degree"],
            "alpha_identity": record["alpha_identity"],
            "extracted_nonnegative_factors":
                record["extracted_nonnegative_factors"],
            "polynomial_coefficients": [
                {"powers": list(powers), "expression": str(coefficient)}
                for powers, coefficient in polynomial.terms()
            ],
            "simplices": simplex_records,
        }
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--skip-face-proof", action="store_true",
        help="development only: derive formulas without exact ideal checks",
    )
    args = parser.parse_args(argv)
    args.out.write_text(json.dumps(
        derive(verify_faces=not args.skip_face_proof), indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
