"""Derive and verify the exact polynomial structure for ``n=7``.

The seven uniform-simplex cap regions have 64 six-simplices.  Regions E,
F, and G are obtained from C, B, and A by the exact involution that reverses
the eight full Stringer increments and sends
``(x,y,z,w,u,v)`` to ``(1-v,1-u,1-w,1-z,1-y,1-x)``.  Only the four source
polynomials A--D therefore need to be stored.

Structural Bernstein zeros are proved over ``QQ(b,c,d,e,f,g,h)``.  For an
affine face ideal ``I`` and order ``q``, a generic affine change of variables
makes the generators of ``I`` normal coordinates.  Sparse multivariate
Horner evaluation takes place in the exact quotient by all normal monomials
of total degree ``q``.  A zero remainder is exactly the assertion
``P in I**q``.  The high-order checks are split into independently runnable
groups for continuous integration.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
import re
import shutil
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import sympy as sp

from n7_gaffke_structure_data import (
    DIMENSION,
    FACE_CONDITIONS,
    REFLECTED_BOUNDARY_LABEL,
    REFLECTED_REGIONS,
    REGION_DEGREES,
    REGION_NAMES,
    REGION_SIMPLICES,
    SOURCE_REGION,
    canonical_boundary_labels,
    structural_zero_indices,
)


HERE = Path(__file__).resolve().parent
DEFAULT_STRUCTURE_PATH = (
    HERE.parent / "certificates" / "n7-gaffke-bernstein-structure.json.gz"
)
COORDINATE_NAMES = ("x", "y", "z", "w", "u", "v")
WEIGHT_NAMES = ("b", "c", "d", "e", "f", "g", "h")
SOURCE_REGIONS = "ABCD"


def _product(expressions):
    return "*".join(f"({value})" for value in expressions) \
        if expressions else "1"


def _node(index):
    coordinate = COORDINATE_NAMES[index]
    return _product(
        [coordinate]
        + [f"{coordinate}-{COORDINATE_NAMES[j]}" for j in range(index)]
    )


def _singular_residual_source(region):
    """Render one factored cleared residual in Singular syntax."""

    k = ord(region) - ord("A")
    if not 0 <= k <= 3:
        raise ValueError(f"region {region} is not a stored source")
    n = 7
    d = DIMENSION
    s = "+".join(
        f"{weight}*{coordinate}"
        for weight, coordinate in zip(WEIGHT_NAMES[:d], COORDINATE_NAMES)
    ) + f"+{WEIGHT_NAMES[d]}"
    common = _product(f"1-{coordinate}" for coordinate in COORDINATE_NAMES)
    upper = list(range(d-k, d))
    q = "+".join(WEIGHT_NAMES[d-k:])
    alpha = "+".join(
        f"{math.comb(n, j)}*({q})^{j}*(1-({q}))^{n-j}"
        for j in range(k+1)
    )
    nodes = _product(_node(i) for i in upper)
    numerator = (
        f"({alpha})*({common})*({nodes})"
        f"-(1-({s}))^{n}*({nodes})"
    )
    for i in upper:
        ratio = [f"1-{COORDINATE_NAMES[j]}"
                 for j in range(d) if j != i]
        for j in upper:
            if j == i:
                continue
            pieces = [COORDINATE_NAMES[j]] + [
                f"{COORDINATE_NAMES[j]}-{COORDINATE_NAMES[t]}"
                for t in range(j) if t != i
            ]
            ratio.append(_product(pieces))
        sign = (-1) ** (d-i)
        term = (
            f"({COORDINATE_NAMES[i]}-({s}))^{n}"
            f"*({_product(ratio)})"
        )
        numerator += ("-" if sign == 1 else "+") + term
    factors = (
        [f"1-{COORDINATE_NAMES[i]}" for i in upper]
        + [f"{COORDINATE_NAMES[j]}-{COORDINATE_NAMES[i]}"
           for offset, i in enumerate(upper)
           for j in upper[offset+1:]]
    )
    lines = [f"poly P={numerator};"]
    lines.extend(f"P=P/({factor});" for factor in factors)
    lines.extend([
        f"poly check=({numerator})-P*({_product(factors)});",
        'if(check!=0){"DIVISION_FAILED";quit;}',
    ])
    return "\n".join(lines)


def _split_singular_terms(text):
    output = []
    start = 0
    depth = 0
    for index, character in enumerate(text):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif index > start and depth == 0 and character in "+-":
            output.append(text[start:index])
            start = index
    output.append(text[start:])
    return output


def _parse_singular_term(term):
    powers = [0] * DIMENSION
    depth = 0
    cut = len(term)
    for index, character in enumerate(term):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif (depth == 0 and character in COORDINATE_NAMES
              and (index == 0 or term[index-1] in "+-*")):
            cut = index
            break
    prefix = term[:cut]
    monomial = term[cut:]
    if prefix in ("", "+"):
        coefficient = "1"
    elif prefix == "-":
        coefficient = "-1"
    else:
        coefficient = prefix
        if coefficient.endswith("*"):
            coefficient = coefficient[:-1]
        if coefficient.startswith("+"):
            coefficient = coefficient[1:]
        if coefficient.startswith("(") and coefficient.endswith(")"):
            coefficient = coefficient[1:-1]
    for name, exponent in re.findall(
            r"(x|y|z|w|u|v)(?:\^(\d+))?", monomial):
        powers[COORDINATE_NAMES.index(name)] = int(exponent or 1)
    stripped = re.sub(
        r"\*?(?:x|y|z|w|u|v)(?:\^\d+)?", "", monomial)
    if stripped:
        raise ValueError((term, prefix, monomial, stripped))
    return {
        "powers": powers,
        "expression": coefficient.replace("^", "**"),
    }


def derive_source_region(region):
    if shutil.which("Singular") is None:
        raise RuntimeError("Singular is required for the n=7 derivation")
    script = (
        "option(redSB);\n"
        f"ring r=(0,{','.join(WEIGHT_NAMES)}),"
        f"({','.join(COORDINATE_NAMES)}),dp;\n"
        "short=0;\n"
        f"{_singular_residual_source(region)}\n"
        '"META|"+string(size(P))+"|"+string(deg(P));\n'
        "print(P);\nquit;\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".sing", delete=False) as f:
        f.write(script)
        path = Path(f.name)
    try:
        started = time.time()
        completed = subprocess.run(
            ["Singular", "-q", str(path)],
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        path.unlink(missing_ok=True)
    if completed.stderr.strip():
        raise RuntimeError(completed.stderr)
    lines = completed.stdout.splitlines()
    if len(lines) != 2 or not lines[0].startswith("META|"):
        raise RuntimeError(completed.stdout[:1000])
    _, term_count, degree = lines[0].split("|")
    records = [
        _parse_singular_term(term)
        for term in _split_singular_terms(lines[1])
    ]
    if len(records) != int(term_count):
        raise AssertionError((region, len(records), term_count))
    if int(degree) != REGION_DEGREES[region]:
        raise AssertionError((region, degree))
    print(
        f"derived source region {region}: degree {degree}, "
        f"{term_count} terms, {time.time()-started:.2f}s",
        flush=True,
    )
    return {
        "degree": int(degree),
        "term_count": int(term_count),
        "polynomial_coefficients": records,
    }


def _load_cached_source(region, directory):
    payload = json.loads((directory / f"n7{region}.json").read_text())
    if (payload.get("region") != region
            or payload.get("degree") != REGION_DEGREES[region]
            or payload.get("term_count") != len(
                payload.get("polynomial_coefficients", []))):
        raise ValueError(f"invalid cached source region {region}")
    return {
        "degree": payload["degree"],
        "term_count": payload["term_count"],
        "polynomial_coefficients": payload["polynomial_coefficients"],
    }


def symbolic_setup():
    coordinates = sp.symbols(" ".join(COORDINATE_NAMES))
    weights = sp.symbols(" ".join(WEIGHT_NAMES))
    s = sum(weight*coordinate for weight, coordinate
            in zip(weights[:DIMENSION], coordinates)) + weights[DIMENSION]
    return coordinates, weights, s


def ordered_simplex_vertices():
    return [
        tuple([sp.Integer(0)] * (DIMENSION-i) + [sp.Integer(1)] * i)
        for i in range(DIMENSION+1)
    ]


def symbolic_vertices(coordinates, weights, s):
    base = ordered_simplex_vertices()
    output = {f"v{i}": point for i, point in enumerate(base)}
    needed = {
        name for simplices in REGION_SIMPLICES.values()
        for simplex in simplices for name in simplex
    }
    for name in sorted(needed):
        if name in output:
            continue
        region = int(name[1])
        left, right = map(int, name.split("_")[1])
        boundary = s-coordinates[region-1]
        left_value = sp.expand(boundary.subs(dict(zip(coordinates, base[left]))))
        right_value = sp.expand(boundary.subs(dict(zip(coordinates, base[right]))))
        parameter = sp.cancel(left_value/(left_value-right_value))
        output[name] = tuple(
            sp.cancel(base[left][j] + parameter*(base[right][j]-base[left][j]))
            for j in range(DIMENSION)
        )
    return output


def boundary_expressions(coordinates, s):
    x, y, z, w, u, v = coordinates
    return {
        "x=0": x,
        "x=y": y-x,
        "y=z": z-y,
        "z=w": w-z,
        "w=u": u-w,
        "u=v": v-u,
        "v=1": 1-v,
        "s=x": s-x,
        "s=y": s-y,
        "s=z": s-z,
        "s=w": s-w,
        "s=u": s-u,
        "s=v": s-v,
    }


def vertex_boundary_identities(vertices, coordinates, s):
    boundaries = boundary_expressions(coordinates, s)
    return {
        name: [
            label for label, expression in boundaries.items()
            if sp.cancel(expression.subs(dict(zip(coordinates, point)))) == 0
        ]
        for name, point in sorted(vertices.items())
    }


def _matrix_rank(rows):
    return sp.Matrix(rows).rank() if rows else 0


def independent_face_equations(simplex, condition, vertices, coordinates,
                               weights, s):
    subset, _ = condition
    position_map = {i: i+1 for i in range(DIMENSION)}
    position_map[DIMENSION] = 0
    excluded = {position_map[position] for position in subset}
    face_names = [name for index, name in enumerate(simplex)
                  if index not in excluded]
    boundaries = boundary_expressions(coordinates, s)
    shared = []
    for label, expression in boundaries.items():
        if all(sp.cancel(expression.subs(dict(zip(
                coordinates, vertices[name])))) == 0 for name in face_names):
            shared.append((label, expression))
    generic = {symbol: sp.Rational(index+2, 37)
               for index, symbol in enumerate(weights)}
    independent = []
    rows = []
    rank = 0
    for label, expression in shared:
        specialized = sp.Poly(expression.subs(generic), *coordinates)
        row = [specialized.coeff_monomial(variable)
               for variable in coordinates]
        new_rank = _matrix_rank(rows+[row])
        if new_rank > rank:
            independent.append(label)
            rows.append(row)
            rank = new_rank
    if rank != len(subset):
        raise AssertionError(
            f"face codimension mismatch: {simplex}, {condition}, {rank}")
    return face_names, [label for label, _ in shared], independent


def _canonical_independent_labels(labels, coordinates, weights, s):
    expressions = boundary_expressions(coordinates, s)
    generic = {symbol: sp.Rational(index+2, 37)
               for index, symbol in enumerate(weights)}
    selected = []
    rows = []
    rank = 0
    for label in canonical_boundary_labels(labels):
        polynomial = sp.Poly(expressions[label].subs(generic), *coordinates)
        row = [polynomial.coeff_monomial(variable)
               for variable in coordinates]
        new_rank = _matrix_rank(rows+[row])
        if new_rank > rank:
            selected.append(label)
            rows.append(row)
            rank = new_rank
    return canonical_boundary_labels(selected)


def _proof_source(region, active_labels, coordinates, weights, s):
    labels = active_labels
    if region in REFLECTED_REGIONS:
        labels = [REFLECTED_BOUNDARY_LABEL[label] for label in labels]
    return (SOURCE_REGION[region],
            _canonical_independent_labels(
                labels, coordinates, weights, s))


def alpha_identity(region):
    index = ord(region)-ord("A")
    if index == DIMENSION:
        return "1-alpha=(b+c+d+e+f+g+h)^7"
    q = "+".join(WEIGHT_NAMES[DIMENSION-index:])
    terms = [
        f"{math.comb(7, j)}*({q})^{j}*(1-({q}))^{7-j}"
        for j in range(index+1)
    ]
    return "alpha=" + "+".join(terms)


def extracted_factors(region):
    index = ord(region)-ord("A")
    if index in (0, DIMENSION):
        return []
    upper = list(range(DIMENSION-index, DIMENSION))
    return (
        [f"1-{COORDINATE_NAMES[i]}" for i in upper]
        + [f"{COORDINATE_NAMES[j]}-{COORDINATE_NAMES[i]}"
           for offset, i in enumerate(upper)
           for j in upper[offset+1:]]
    )


def derive(source_cache=None):
    coordinates, weights, s = symbolic_setup()
    vertices = symbolic_vertices(coordinates, weights, s)
    sources = {
        region: (_load_cached_source(region, source_cache)
                 if source_cache else derive_source_region(region))
        for region in SOURCE_REGIONS
    }
    output = {
        "schema_version": 1,
        "claim": (
            "Exact polynomial, six-simplex, symmetry, and structural-zero "
            "structure for the n=7 Gaffke-domination certificate."
        ),
        "variables": list(WEIGHT_NAMES),
        "coordinate_variables": list(COORDINATE_NAMES),
        "polynomial_sources": sources,
        "reflection_identity": {
            "weight_map": (
                "(a,b,c,d,e,f,g,h)->(h,g,f,e,d,c,b,a), "
                "a=1-b-c-d-e-f-g-h"
            ),
            "coordinate_map": "(x,y,z,w,u,v)->(1-v,1-u,1-w,1-z,1-y,1-x)",
            "residual_map": "P_G=-R(P_A), P_F=-R(P_B), P_E=-R(P_C)",
            "justification": (
                "Reversing Dirichlet spacings sends the cap probability V "
                "to 1-V and the binomial tail identity alpha to 1-alpha."
            ),
        },
        "vertex_formulas": {
            name: [str(sp.factor(value)) for value in point]
            for name, point in sorted(vertices.items())
        },
        "vertex_boundary_identities": vertex_boundary_identities(
            vertices, coordinates, s),
        "regions": {},
        "face_order_verification": (
            "exact_generic_Singular_quotient_Horner_ideal_power_membership"
        ),
    }
    proof_source_pivots = {}
    for region in REGION_NAMES:
        simplex_records = []
        total_zeros = 0
        for simplex_index, simplex in enumerate(REGION_SIMPLICES[region]):
            conditions = FACE_CONDITIONS[region][simplex_index]
            zeros = structural_zero_indices(REGION_DEGREES[region], conditions)
            total_zeros += len(zeros)
            proofs = []
            for condition in conditions:
                face_names, active, independent = independent_face_equations(
                    simplex, condition, vertices, coordinates, weights, s)
                source_region, source_labels = _proof_source(
                    region, active, coordinates, weights, s)
                source_labels = tuple(source_labels)
                if source_labels not in proof_source_pivots:
                    proof_source_pivots[source_labels] = (
                        _inverse_face_map(source_labels)[1])
                source_pivot = proof_source_pivots[source_labels]
                proofs.append({
                    "bernstein_index_subset": list(condition[0]),
                    "maximum_subset_sum": condition[1],
                    "vanishing_order": condition[1]+1,
                    "face_vertices": face_names,
                    "active_boundaries": active,
                    "independent_ideal_generators": independent,
                    "proof_source_region": source_region,
                    "proof_source_ideal_generators": list(source_labels),
                    "proof_source_inverse_pivot_columns": [
                        COORDINATE_NAMES[index] for index in source_pivot
                    ],
                    "verification": (
                        "exact_generic_Singular_quotient_Horner_"
                        "ideal_power_membership_over_QQ(b,c,d,e,f,g,h)"
                    ),
                })
            simplex_records.append({
                "vertices": simplex,
                "degree": REGION_DEGREES[region],
                "coefficient_count": math.comb(
                    REGION_DEGREES[region]+DIMENSION, DIMENSION),
                "structural_zero_count": len(zeros),
                "face_conditions": [
                    [list(subset), maximum] for subset, maximum in conditions
                ],
                "face_order_proofs": proofs,
            })
        output["regions"][region] = {
            "degree": REGION_DEGREES[region],
            "source_region": SOURCE_REGION[region],
            "reflection_applied": region in REFLECTED_REGIONS,
            "alpha_identity": alpha_identity(region),
            "extracted_nonnegative_factors": extracted_factors(region),
            "simplex_count": len(simplex_records),
            "total_structural_zero_count": total_zeros,
            "simplices": simplex_records,
        }
    return output


def _json_bytes(payload):
    return (json.dumps(payload, sort_keys=True, separators=(",", ":"))
            + "\n").encode()


def write_gzip_json(path, payload):
    raw = _json_bytes(payload)
    with path.open("wb") as handle:
        with gzip.GzipFile(
                filename="", mode="wb", fileobj=handle,
                compresslevel=9, mtime=0) as compressed:
            compressed.write(raw)


def read_gzip_json(path):
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def structure_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _boundary_expression_strings():
    linear = "b*x+c*y+d*z+e*w+f*u+g*v+h"
    return {
        "x=0": "x", "x=y": "y-x", "y=z": "z-y",
        "z=w": "w-z", "w=u": "u-w", "u=v": "v-u",
        "v=1": "1-v",
        **{f"s={name}": f"({linear})-{name}"
           for name in COORDINATE_NAMES},
    }


def _inverse_face_map(labels):
    coordinates, weights, s = symbolic_setup()
    expressions = boundary_expressions(coordinates, s)
    normal_count = len(labels)
    normals = sp.symbols(" ".join(f"r{i}" for i in range(normal_count)))
    if normal_count == 1:
        normals = (normals,)
    tangent_count = DIMENSION-normal_count
    tangents = (sp.symbols(" ".join(f"t{i}" for i in range(tangent_count)))
                if tangent_count else ())
    if tangent_count == 1:
        tangents = (tangents,)
    if all(label.startswith("s=") for label in labels):
        coordinate_by_name = dict(zip(COORDINATE_NAMES, coordinates))
        del coordinate_by_name  # names, rather than symbols, index this map
        weight_by_name = dict(zip(COORDINATE_NAMES, weights[:DIMENSION]))
        labeled_names = [label.split("=", 1)[1] for label in labels]
        free_names = [name for name in COORDINATE_NAMES
                      if name not in labeled_names]
        tangent_by_name = dict(zip(free_names, tangents))
        denominator = 1-sum(weight_by_name[name]
                            for name in labeled_names)
        threshold = weights[DIMENSION]
        threshold += sum(weight_by_name[name]*tangent_by_name[name]
                         for name in free_names)
        threshold -= sum(weight_by_name[name]*normal
                         for name, normal in zip(labeled_names, normals))
        threshold = sp.factor(threshold/denominator)
        normal_by_name = dict(zip(labeled_names, normals))
        images = [
            threshold-normal_by_name[name]
            if name in normal_by_name else tangent_by_name[name]
            for name in COORDINATE_NAMES
        ]
        pivot = tuple(COORDINATE_NAMES.index(name)
                      for name in labeled_names)
    else:
        matrix = sp.Matrix([
            [sp.diff(expressions[label], coordinate)
             for coordinate in coordinates]
            for label in labels
        ])
        constant = sp.Matrix([
            expressions[label].subs(dict.fromkeys(coordinates, 0))
            for label in labels
        ])
        pivot = None
        for columns in itertools.combinations(range(DIMENSION), normal_count):
            if sp.factor(matrix[:, columns].det()) != 0:
                pivot = columns
                break
        if pivot is None:
            raise AssertionError(labels)
        free = tuple(index for index in range(DIMENSION)
                     if index not in pivot)
        right = sp.Matrix(normals)-constant
        if free:
            right -= matrix[:, free]*sp.Matrix(tangents)
        solution = matrix[:, pivot].inv()*right
        images = [None]*DIMENSION
        for column, value in zip(pivot, solution):
            images[column] = sp.factor(value)
        for column, value in zip(free, tangents):
            images[column] = value
    return (
        tuple(str(sp.factor(value)).replace("**", "^") for value in images),
        tuple(pivot),
    )


class _HornerEmitter:
    def __init__(self, polynomial):
        self.polynomial = polynomial
        self.lines = []
        self.counter = 0

    def _new(self):
        self.counter += 1
        return f"p{self.counter}"

    def emit(self, data=None, axis=0):
        data = self.polynomial if data is None else data
        if axis == DIMENSION:
            name = self._new()
            expression = "+".join(f"({value})" for value in data.values())
            self.lines.append(f"poly {name}={expression or '0'};")
            return name
        groups = defaultdict(dict)
        for powers, coefficient in data.items():
            exponent = powers[axis]
            child = list(powers)
            child[axis] = 0
            groups[exponent][tuple(child)] = coefficient
        result = self.emit(groups[max(groups)], axis+1)
        for exponent in range(max(groups)-1, -1, -1):
            self.lines.append(f"{result}=reduce({result}*X{axis},Z);")
            if exponent in groups:
                child = self.emit(groups[exponent], axis+1)
                self.lines.append(f"{result}=reduce({result}+{child},Z);")
                self.lines.append(f"kill {child};")
        return result


def render_face_proof_script(structure, key):
    region, labels, order = key
    source = structure["polynomial_sources"][region]
    polynomial = {
        tuple(item["powers"]): item["expression"].replace("**", "^")
        for item in source["polynomial_coefficients"]
    }
    images, pivot = _inverse_face_map(labels)
    stored_pivots = {
        tuple(proof["proof_source_inverse_pivot_columns"])
        for region_record in structure["regions"].values()
        for simplex in region_record["simplices"]
        for proof in simplex["face_order_proofs"]
        if (
            proof["proof_source_region"],
            tuple(proof["proof_source_ideal_generators"]),
            proof["vanishing_order"],
        ) == key
    }
    expected_pivot = tuple(COORDINATE_NAMES[index] for index in pivot)
    if stored_pivots != {expected_pivot}:
        raise AssertionError(
            f"stored inverse pivot differs for {key}: {stored_pivots}")
    emitter = _HornerEmitter(polynomial)
    root = emitter.emit()
    variables = [f"r{i}" for i in range(len(labels))]
    variables.extend(
        f"t{i}" for i in range(DIMENSION-len(labels)))
    normal_ideal = ",".join(variables[:len(labels)])
    image_lines = ";".join(
        f"poly X{i}={image}" for i, image in enumerate(images))
    script = f'''option(redSB);
ring raw=(0,{','.join(WEIGHT_NAMES)}),({','.join(variables)}),dp;
ideal J={normal_ideal};J=J^{order};
qring target=std(J);
ideal Z=std(0);
{image_lines};
{chr(10).join(emitter.lines)}
if({root}==0){{"PASS";}}else{{"FAIL|"+string(size({root}));}}
quit;
'''
    return script, pivot, emitter.counter


def face_proof_keys(structure):
    keys = set()
    for region in structure["regions"].values():
        for simplex in region["simplices"]:
            for proof in simplex["face_order_proofs"]:
                keys.add((
                    proof["proof_source_region"],
                    tuple(proof["proof_source_ideal_generators"]),
                    proof["vanishing_order"],
                ))
    return sorted(keys, key=lambda key: (key[0], key[2], key[1]))


FACE_GROUPS = {
    "a": lambda key: key[0] == "A",
    "b-basic": lambda key: key[0] == "B" and key[2] < 5,
    "b5": lambda key: key[0] == "B" and key[2] == 5,
    "c-basic": lambda key: key[0] == "C" and key[2] < 6,
    "c6": lambda key: key[0] == "C" and key[2] == 6,
    "c8": lambda key: key[0] == "C" and key[2] == 8,
    "d-basic": lambda key: key[0] == "D" and key[2] < 6,
    "d6": lambda key: key[0] == "D" and key[2] == 6,
    "d9": lambda key: key[0] == "D" and key[2] == 9,
}


def verify_face_group(structure_path, group):
    if group not in FACE_GROUPS:
        raise ValueError(f"unknown n=7 face group: {group}")
    if shutil.which("Singular") is None:
        raise RuntimeError("Singular is required for n=7 face proofs")
    structure = read_gzip_json(structure_path)
    keys = [key for key in face_proof_keys(structure)
            if FACE_GROUPS[group](key)]
    if not keys:
        raise AssertionError(f"empty n=7 face group: {group}")
    with tempfile.TemporaryDirectory(prefix="stringer-n7-face-") as directory:
        directory = Path(directory)
        for index, key in enumerate(keys):
            script, pivot, temporary_count = render_face_proof_script(
                structure, key)
            path = directory / f"{index:02d}.sing"
            path.write_text(script)
            started = time.time()
            completed = subprocess.run(
                ["Singular", "-q", str(path)],
                capture_output=True,
                text=True,
            )
            rendered = completed.stdout.strip()
            if completed.returncode != 0 or rendered != "PASS":
                raise AssertionError(
                    f"n=7 face proof failed for {key}: "
                    f"{rendered}\n{completed.stderr}")
            print(
                f"n=7 exact generic face proof passed: {key}; "
                f"pivot={pivot}; temporaries={temporary_count}; "
                f"seconds={time.time()-started:.2f}",
                flush=True,
            )
    print(f"n=7 face-proof group {group}: PASS ({len(keys)} keys)")


def _proof_independent_projection(structure):
    # All stored structure fields are proof inputs or claims; execution of
    # the split face jobs produces no machine-dependent data in the artifact.
    return structure


def check_structure_data(structure_path):
    expected = read_gzip_json(structure_path)
    observed = derive()
    if (_proof_independent_projection(observed)
            != _proof_independent_projection(expected)):
        raise AssertionError("regenerated n=7 structure differs")
    print("n=7 proof-independent structure regeneration: PASS")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--source-cache-dir", type=Path)
    parser.add_argument("--check-data-against", type=Path)
    parser.add_argument("--face-proof-group", choices=sorted(FACE_GROUPS))
    parser.add_argument("--structure", type=Path,
                        default=DEFAULT_STRUCTURE_PATH)
    args = parser.parse_args(argv)
    selected = sum(bool(value) for value in (
        args.out, args.check_data_against, args.face_proof_group))
    if selected != 1:
        parser.error("select exactly one action")
    if args.source_cache_dir and not args.out:
        parser.error("--source-cache-dir requires --out")
    if args.out:
        structure = derive(args.source_cache_dir)
        write_gzip_json(args.out, structure)
        print(
            f"wrote {args.out}; compressed_sha256="
            f"{structure_sha256(args.out)}",
            flush=True,
        )
    elif args.check_data_against:
        check_structure_data(args.check_data_against)
    else:
        verify_face_group(args.structure, args.face_proof_group)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
