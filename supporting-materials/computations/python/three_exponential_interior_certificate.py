"""Rigorous interior certificate for the three-exponential inequality.

The theory note ``THREE-EXPONENTIAL-QUANTILE.md`` reduces convexity at the
sharp candidate level ``alpha = 4*exp(-3)`` to

    exp(-h) * (1 + A*h + Z*h**2) >= 4*exp(-3),

where ``h`` is the largest generalized eigenvalue of the boundary and bulk
Gram matrices ``(L,K)``.  This module proves the inequality on the closed
interior square ``1 <= z,w <= 4``.

For each parameter box it proposes a rational witness ``x`` numerically and
then proves, using outward-rounded Arb balls, both

    x*K - L positive definite,
    exp(-x) * (1 + A*x + Z*x**2) > 4*exp(-3).

The first inequality gives ``h < x``; monotonicity of the tail gives the
target.  Numerical eigensolvers and root finding only propose ``x`` and make
no sign decision.

Moment enclosures use fourth-order Taylor expansions about each box center.
Every coefficient is an elementary tilted-simplex or boundary moment.  The
remainder follows from ``|delta_z*u + delta_w*v| <= H`` and the exponential
Taylor bound ``exp(H)*H**5/5!``.  Centering the Gram polynomials at a fixed
rational approximation to the tilted mean avoids subtracting broad raw
moment intervals; the centered identities are exact for any chosen center.

The computation is serial and memory bounded.  On the development machine it
takes about one minute and peaks below 85 MB resident memory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import flint
import numpy as np
from flint import arb, ctx, fmpq
from scipy.linalg import eigvalsh
from scipy.optimize import brentq


PRECISION_BITS = 256
SERIES_TERMS = 120
TAYLOR_DEGREE = 4
WITNESS_BITS = 30
DOMAIN_LOWER = Fraction(1)
DOMAIN_UPPER = Fraction(4)
COARSE_STEP = Fraction(1, 50)
FINE_STEP = COARSE_STEP / 2
GRID_SIZE = int((DOMAIN_UPPER - DOMAIN_LOWER) / COARSE_STEP)

ctx.prec = PRECISION_BITS


def _const(value) -> arb:
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def _radius_ball(radius) -> arb:
    radius = Fraction(radius)
    return arb(0, fmpq(radius.numerator, radius.denominator))


def _arb_power(value: arb, exponent: int) -> arb:
    """Nonnegative integer power without a zero-midpoint ``arb**k`` edge."""
    result = arb(1)
    for _ in range(exponent):
        result *= value
    return result


def _moment_f_series(order: int, argument: arb,
                     absolute_argument: Fraction) -> arb:
    """Enclose integral_0^1 t^order exp(-argument*t) dt by a series."""
    value = arb(0)
    for k in range(SERIES_TERMS, -1, -1):
        coefficient = Fraction(
            -1 if k % 2 else 1,
            factorial(k) * (order + k + 1),
        )
        value = value * argument + _const(coefficient)
    radius = _const(absolute_argument)
    tail = (radius.exp() * _arb_power(radius, SERIES_TERMS + 1)
            / factorial(SERIES_TERMS + 1))
    return value + arb(0, tail.upper())


def _moment_f(order: int, argument: arb,
              absolute_argument: Fraction) -> arb:
    """Enclose ``F_order(r)=integral_0^1 t^order exp(-r*t)dt``."""
    if absolute_argument <= Fraction(1, 2):
        return _moment_f_series(order, argument, absolute_argument)
    exponential = (-argument).exp()
    partial = arb(0)
    term = arb(1)
    for k in range(order + 1):
        if k:
            term = term * argument / k
        partial += term
    return (factorial(order) * (1 - exponential * partial)
            / _arb_power(argument, order + 1))


def _center_moments(z_value: Fraction, w_value: Fraction, max_degree: int):
    """Return exact-ball bulk and hypotenuse monomial moments at a point."""
    z = _const(z_value)
    w = _const(w_value)
    difference = z - w
    difference_size = abs(z_value - w_value)
    exp_minus_w = (-w).exp()
    difference_moments = [
        _moment_f(order, difference, difference_size)
        for order in range(max_degree + 1)
    ]
    z_moments = [
        _moment_f(order, z, abs(z_value))
        for order in range(max_degree + 1)
    ]
    bulk = {}
    boundary = {}
    for i in range(max_degree + 1):
        for j in range(max_degree + 1 - i):
            inner = arb(0)
            for r in range(j + 1):
                polynomial = arb(0)
                for q in range(r + 1):
                    polynomial += (
                        (-1) ** q * comb(r, q)
                        * difference_moments[i + q]
                    )
                inner += _arb_power(w, r) * polynomial / factorial(r)
            bulk[i, j] = (
                factorial(j)
                / _arb_power(w, j + 1)
                * (z_moments[i] - exp_minus_w * inner)
            )
            polynomial = arb(0)
            for q in range(j + 1):
                polynomial += (
                    (-1) ** q * comb(j, q)
                    * difference_moments[i + q]
                )
            boundary[i, j] = exp_minus_w * polynomial
    return bulk, boundary


def _linear_center(coefficients, moments, a=0, b=0) -> arb:
    value = arb(0)
    for (i, j), coefficient in coefficients.items():
        value += _const(coefficient) * moments[i + a, j + b]
    return value


def _linear_enclosure(coefficients, moments, half_z: Fraction,
                      half_w: Fraction, supremum=1) -> arb:
    """Taylor-enclose a polynomial-weighted moment on one parameter box."""
    delta_z = _radius_ball(half_z)
    delta_w = _radius_ball(half_w)
    value = arb(0)
    for a in range(TAYLOR_DEGREE + 1):
        for b in range(TAYLOR_DEGREE + 1 - a):
            value += (
                (-1) ** (a + b)
                * _arb_power(delta_z, a)
                * _arb_power(delta_w, b)
                / factorial(a) / factorial(b)
                * _linear_center(coefficients, moments, a, b)
            )
    total_radius = half_z + half_w
    radius = _const(total_radius)
    tail = (
        _const(supremum)
        * moments[0, 0].upper()
        * radius.exp()
        * _arb_power(radius, TAYLOR_DEGREE + 1)
        / factorial(TAYLOR_DEGREE + 1)
    )
    return value + arb(0, tail.upper())


def _gram_enclosures(box):
    """Enclose ``Z,A,K,L`` throughout a rational parameter box."""
    z_lower, z_upper, w_lower, w_upper = map(Fraction, box)
    z_center = (z_lower + z_upper) / 2
    w_center = (w_lower + w_upper) / 2
    half_z = (z_upper - z_lower) / 2
    half_w = (w_upper - w_lower) / 2
    bulk, boundary = _center_moments(
        z_center, w_center, 2 + TAYLOR_DEGREE)

    # Any fixed center gives the exact covariance identities below.  A
    # numerical mean only improves enclosure width and is rounded to a
    # dyadic rational before entering the proof computation.
    center_denominator = 1 << 50
    mean_u = Fraction(
        round(float((bulk[1, 0] / bulk[0, 0]).mid()) * center_denominator),
        center_denominator,
    )
    mean_v = Fraction(
        round(float((bulk[0, 1] / bulk[0, 0]).mid()) * center_denominator),
        center_denominator,
    )
    polynomials = {
        "one": {(0, 0): Fraction(1)},
        "q_u": {(1, 0): Fraction(1), (0, 0): -mean_u},
        "q_v": {(0, 1): Fraction(1), (0, 0): -mean_v},
        "q_uu": {
            (2, 0): Fraction(1),
            (1, 0): -2 * mean_u,
            (0, 0): mean_u * mean_u,
        },
        "q_uv": {
            (1, 1): Fraction(1),
            (1, 0): -mean_v,
            (0, 1): -mean_u,
            (0, 0): mean_u * mean_v,
        },
        "q_vv": {
            (0, 2): Fraction(1),
            (0, 1): -2 * mean_v,
            (0, 0): mean_v * mean_v,
        },
    }

    z_partition = _linear_enclosure(
        polynomials["one"], bulk, half_z, half_w)
    q_u = _linear_enclosure(polynomials["q_u"], bulk, half_z, half_w)
    q_v = _linear_enclosure(polynomials["q_v"], bulk, half_z, half_w)
    q_uu = _linear_enclosure(polynomials["q_uu"], bulk, half_z, half_w)
    q_uv = _linear_enclosure(polynomials["q_uv"], bulk, half_z, half_w)
    q_vv = _linear_enclosure(polynomials["q_vv"], bulk, half_z, half_w)
    k_uu = q_uu - q_u * q_u / z_partition
    k_uv = q_uv - q_u * q_v / z_partition
    k_vv = q_vv - q_v * q_v / z_partition

    boundary_mass = _linear_enclosure(
        polynomials["one"], boundary, half_z, half_w)
    boundary_u = _linear_enclosure(
        polynomials["q_u"], boundary, half_z, half_w)
    boundary_v = _linear_enclosure(
        polynomials["q_v"], boundary, half_z, half_w)
    boundary_uu = _linear_enclosure(
        polynomials["q_uu"], boundary, half_z, half_w)
    boundary_uv = _linear_enclosure(
        polynomials["q_uv"], boundary, half_z, half_w)
    boundary_vv = _linear_enclosure(
        polynomials["q_vv"], boundary, half_z, half_w)
    mean_delta_u = q_u / z_partition
    mean_delta_v = q_v / z_partition
    l_uu = (
        boundary_uu - 2 * mean_delta_u * boundary_u
        + mean_delta_u * mean_delta_u * boundary_mass
    )
    l_uv = (
        boundary_uv
        - mean_delta_u * boundary_v
        - mean_delta_v * boundary_u
        + mean_delta_u * mean_delta_v * boundary_mass
    )
    l_vv = (
        boundary_vv - 2 * mean_delta_v * boundary_v
        + mean_delta_v * mean_delta_v * boundary_mass
    )

    # A=B(z)+zZ.  Use a one-variable centered Taylor enclosure for B.
    z_midpoint = _const(z_center)
    delta_z = _radius_ball(half_z)
    z_ball = z_midpoint + delta_z
    b_value = arb(0)
    for a in range(TAYLOR_DEGREE + 1):
        b_value += (
            (-1) ** a * _arb_power(delta_z, a) / factorial(a)
            * _moment_f(a, z_midpoint, abs(z_center))
        )
    radius = _const(half_z)
    b_tail = (
        _moment_f(0, z_midpoint, abs(z_center)).upper()
        * radius.exp()
        * _arb_power(radius, TAYLOR_DEGREE + 1)
        / factorial(TAYLOR_DEGREE + 1)
    )
    b_value += arb(0, b_tail.upper())
    a_value = b_value + z_ball * z_partition
    return (
        z_partition,
        a_value,
        (k_uu, k_uv, k_vv),
        (l_uu, l_uv, l_vv),
    )


def _midpoint_float(value: arb) -> float:
    return float(value.mid())


def certify_cell(box):
    """Certify one box and return its rigorous margins and rational witness."""
    z_partition, a_value, k_matrix, l_matrix = _gram_enclosures(box)
    k_midpoint = np.array([
        [_midpoint_float(k_matrix[0]), _midpoint_float(k_matrix[1])],
        [_midpoint_float(k_matrix[1]), _midpoint_float(k_matrix[2])],
    ])
    l_midpoint = np.array([
        [_midpoint_float(l_matrix[0]), _midpoint_float(l_matrix[1])],
        [_midpoint_float(l_matrix[1]), _midpoint_float(l_matrix[2])],
    ])
    proposed_h = eigvalsh(l_midpoint, k_midpoint)[-1]
    proposed_a = _midpoint_float(a_value)
    proposed_z = _midpoint_float(z_partition)
    target = 4 * np.exp(-3)
    proposed_tail_root = brentq(
        lambda x: np.exp(-x) * (
            1 + proposed_a * x + proposed_z * x * x) - target,
        max(proposed_h, 0),
        10,
    )
    proposed_witness = (proposed_h + proposed_tail_root) / 2
    denominator = 1 << WITNESS_BITS
    witness = Fraction(round(proposed_witness * denominator), denominator)
    x = _const(witness)

    m_uu = x * k_matrix[0] - l_matrix[0]
    m_uv = x * k_matrix[1] - l_matrix[1]
    m_vv = x * k_matrix[2] - l_matrix[2]
    determinant = m_uu * m_vv - m_uv * m_uv
    tail_margin = (
        (-x).exp() * (1 + a_value * x + z_partition * x * x)
        - 4 * arb(-3).exp()
    )
    margins = (m_uu, determinant, tail_margin)
    certified = all(margin > 0 for margin in margins)
    return certified, margins, witness


def _lower_text(value: arb) -> str:
    return str(value.lower())


def build_certificate():
    coarse_leaves = 0
    fine_leaves = 0
    minima = [(float("inf"), None) for _ in range(3)]
    trace = hashlib.sha256()
    for i in range(GRID_SIZE):
        for j in range(GRID_SIZE):
            z_lower = DOMAIN_LOWER + i * COARSE_STEP
            w_lower = DOMAIN_LOWER + j * COARSE_STEP
            box = (
                z_lower, z_lower + COARSE_STEP,
                w_lower, w_lower + COARSE_STEP,
            )
            certified, margins, witness = certify_cell(box)
            leaves = []
            if certified:
                coarse_leaves += 1
                leaves.append((0, 0, margins, witness))
            else:
                for z_half in range(2):
                    for w_half in range(2):
                        subbox = (
                            z_lower + z_half * FINE_STEP,
                            z_lower + (z_half + 1) * FINE_STEP,
                            w_lower + w_half * FINE_STEP,
                            w_lower + (w_half + 1) * FINE_STEP,
                        )
                        subcertified, submargins, subwitness = certify_cell(
                            subbox)
                        if not subcertified:
                            raise AssertionError(
                                "three-exponential cell is not certified: "
                                f"coarse=({i},{j}), half=({z_half},{w_half})")
                        fine_leaves += 1
                        leaves.append((
                            z_half, w_half, submargins, subwitness))

            for z_half, w_half, leaf_margins, leaf_witness in leaves:
                lower_texts = []
                for margin_index, margin in enumerate(leaf_margins):
                    lower = margin.lower()
                    lower_text = str(lower)
                    lower_texts.append(lower_text)
                    lower_float = float(lower)
                    if lower_float < minima[margin_index][0]:
                        minima[margin_index] = (
                            lower_float,
                            {
                                "coarse_indices": [i, j],
                                "half_indices": [z_half, w_half],
                                "lower_enclosure": lower_text,
                                "witness": (
                                    f"{leaf_witness.numerator}/"
                                    f"{leaf_witness.denominator}"
                                ),
                            },
                        )
                trace.update((
                    f"{i},{j},{z_half},{w_half},"
                    f"{leaf_witness.numerator}/{leaf_witness.denominator}|"
                    + "|".join(lower_texts) + "\n"
                ).encode())

    labels = (
        "leading_principal_minor",
        "matrix_determinant",
        "tail_margin",
    )
    minimum_records = {
        label: minima[index][1] for index, label in enumerate(labels)
    }
    return {
        "schema_version": 1,
        "claim": (
            "The explicit three-exponential inequality (16) holds strictly "
            "throughout the closed gap-coordinate square 1<=z,w<=4."
        ),
        "mathematical_status": "rigorous computer-assisted theorem",
        "implication": (
            "At alpha<=4*exp(-3), the three-exponential quantile Hessian is "
            "positive semidefinite at every coefficient vector whose "
            "normalized gap coordinates lie in this square."
        ),
        "domain": {
            "z": ["1", "4"],
            "w": ["1", "4"],
        },
        "proof_interface": (
            "Each leaf has a rational x with x*K-L positive definite and "
            "S_(z,w)(x)>4*exp(-3). Therefore h(z,w)<x and monotonicity gives "
            "S_(z,w)(h)>4*exp(-3)."
        ),
        "arithmetic": {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
            "arb_precision_bits": PRECISION_BITS,
            "moment_series_terms": SERIES_TERMS,
            "parameter_taylor_degree": TAYLOR_DEGREE,
            "witness_dyadic_bits": WITNESS_BITS,
            "floating_point_role": (
                "Only proposes each rational witness and selects which "
                "already-positive lower margin to summarize."
            ),
        },
        "partition": {
            "coarse_step": "1/50",
            "fine_step_on_failed_coarse_cells": "1/100",
            "coarse_grid_shape": [GRID_SIZE, GRID_SIZE],
            "coarse_leaf_count": coarse_leaves,
            "fine_leaf_count": fine_leaves,
            "total_leaf_count": coarse_leaves + fine_leaves,
            "coverage": "closed square, with shared box boundaries",
        },
        "minimum_rigorous_lower_enclosures": minimum_records,
        "leaf_trace_sha256": trace.hexdigest(),
        "conclusion": (
            "All leafwise Arb positivity checks passed; no unresolved box "
            "remains inside 1<=z,w<=4."
        ),
        "limitation": (
            "This closes a compact interior square only. It does not prove "
            "three-exponential convexity on the remaining off-symmetry "
            "domain or imply a general-n Stringer theorem."
        ),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        help="write the machine-readable certificate")
    args = parser.parse_args(argv)
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    partition = certificate["partition"]
    print(
        "three-exponential interior: certified "
        f"{partition['total_leaf_count']} leaves on 1<=z,w<=4 "
        f"({partition['coarse_leaf_count']} coarse, "
        f"{partition['fine_leaf_count']} refined)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
