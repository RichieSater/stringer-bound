"""Exact certificates for the five-coordinate distinct-knot cap barrier.

Let ``D`` be uniform on the five-coordinate simplex.  The two nontrivial
mixed crossing types normalize to

    (-1, -u, t*v, t*w, t),
    (-1, -u, -v, t*w, t).

where ``0 < u,v,w < 1``, ``v < w``, and ``t > 0``.  Write ``p`` for the
upper-cap volume and ``gamma`` for the section centroid.  This module checks
the exact implication

    gamma[4] >= gamma[3]  ==>  p <= 1 - (1-gamma[0])**4.

In both types the final-centroid order condition gives a polynomial
``G >= 0`` and confines the scale to a compact interval.  For the first type,
exact multivariate division gives

    T = S*G + H,

For the second, a sparse support and active-row witness determines ``S`` by
an exact rational linear solve, followed by the same identity.  In each case
``t*T`` is the numerator of the target first-step margin and all
tensor-Bernstein coefficients of ``S`` and ``H`` are nonnegative on the unit
cube.  Together with the elementary endpoint crossings, this checks every
regular distinct-knot section.  Repeated-knot block faces remain outside the
claim.  Floating point is not used by the certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import comb, prod
from pathlib import Path

import numpy as np
import sympy as sp
from flint import fmpq, fmpq_mat

from five_coordinate_mixed_witness import (
    THREE_TWO_ACTIVE_ROWS_FLAT,
    THREE_TWO_MULTIPLIER_DEGREE,
    THREE_TWO_SUPPORT_FLAT,
    THREE_TWO_TARGET_DEGREE,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    HERE.parent / "certificates"
    / "five-coordinate-mixed-section-certificate.json"
)

VARIABLE_NAMES = ("u", "v", "w", "t")
COMPARISON_DEGREE = (11, 9, 9, 22)
ORDER_DEGREE = (4, 1, 2, 5)
QUOTIENT_DEGREE = (7, 15, 22, 54)
REMAINDER_DEGREE = (3, 16, 24, 59)
THREE_TWO_ORDER_DEGREE = (3, 3, 4, 7)
THREE_TWO_COMPARISON_TERMS = 4015


def _poly(expression, variables) -> sp.Poly:
    return sp.Poly(expression, *variables, domain=sp.QQ)


def _fraction(value: sp.Rational) -> Fraction:
    return Fraction(int(value.p), int(value.q))


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _power_digest(poly: sp.Poly) -> str:
    """Hash a sparse exact power-basis expansion in canonical term order."""
    digest = hashlib.sha256()
    digest.update(("degree=" + ",".join(map(str, poly.degree_list())) + "\n").encode())
    for powers, coefficient in poly.terms():
        value = _fraction(coefficient)
        digest.update(
            (
                ",".join(map(str, powers))
                + f":{value.numerator}/{value.denominator}\n"
            ).encode()
        )
    return digest.hexdigest()


def derive_polynomials():
    """Derive the two-below/three-above certificate polynomials exactly."""
    u, v, w, t = sp.symbols("u v w t")
    variables = (u, v, w, t)
    z = sp.symbols("z0:5")

    # Divided-difference formula for the cap above zero.  Translation of all
    # knots differentiates in the opposite direction to the cutting level,
    # so sum_i partial_i(cap) is the section density.  Its normalized gradient
    # is therefore the section centroid.
    cap = sum(
        z[i] ** 4 / sp.prod(z[i] - z[j] for j in range(5) if j != i)
        for i in (2, 3, 4)
    )
    substitution = {
        z[0]: -1,
        z[1]: -u,
        z[2]: t * v,
        z[3]: t * w,
        z[4]: t,
    }
    density = sum(sp.diff(cap, coordinate) for coordinate in z)
    p = sp.factor(cap.subs(substitution))
    gamma0 = sp.factor(sp.diff(cap, z[0]).subs(substitution)
                       / density.subs(substitution))
    gamma3 = sp.factor(sp.diff(cap, z[3]).subs(substitution)
                       / density.subs(substitution))
    gamma4 = sp.factor(sp.diff(cap, z[4]).subs(substitution)
                       / density.subs(substitution))

    p_numerator, p_denominator = map(sp.factor, sp.fraction(p))
    g0_numerator, g0_denominator = map(sp.factor, sp.fraction(gamma0))
    linear_positive = (t + 1) * (t * v + 1) * (t * w + 1)
    other_positive = (t + u) * (t * v + u) * (t * w + u)
    if sp.factor(p_denominator - linear_positive * other_positive) != 0:
        raise AssertionError("cap denominator changed")
    if sp.rem(_poly(p_numerator, variables), _poly(t**2, variables)) != 0:
        raise AssertionError("cap numerator lost its t^2 factor")
    if sp.rem(_poly(g0_numerator, variables), _poly(t, variables)) != 0:
        raise AssertionError("first centroid numerator lost its t factor")

    # B=4*linear_positive*F is the positive denominator of gamma0.  Keeping
    # the powers as Poly objects avoids an expensive generic cancellation of
    # the 6,895-term comparison numerator.
    if sp.rem(
        _poly(g0_denominator, variables),
        _poly(4 * linear_positive, variables),
    ) != 0:
        raise AssertionError("first centroid denominator changed")
    F = sp.exquo(
        _poly(g0_denominator, variables),
        _poly(4 * linear_positive, variables),
    )
    if any(coefficient <= 0 for _, coefficient in F.terms()):
        raise AssertionError("positive density factor lost coefficient positivity")

    A = sp.exquo(_poly(g0_numerator, variables), _poly(t, variables))
    B = _poly(g0_denominator, variables)
    P = sp.exquo(_poly(p_numerator, variables), _poly(t**2, variables))
    U = _poly(other_positive, variables)
    t_poly = _poly(t, variables)

    # If gamma0=t*A/B and p=t^2*P/(linear_positive*U), cancellation of the
    # common factor linear_positive gives
    #
    #   1-(1-gamma0)^4-p = t*T/(B^4*U).
    first_part = (
        4 * B**3 * A
        - 6 * B**2 * t_poly * A**2
        + 4 * B * t_poly**2 * A**3
        - t_poly**3 * A**4
    ) * U
    b4_over_linear = sp.exquo(B**4, _poly(linear_positive, variables))
    T = first_part - t_poly * P * b4_over_linear
    if tuple(T.degree_list()) != COMPARISON_DEGREE:
        raise AssertionError("comparison numerator degree changed")

    # Cross multiplication independently checks the displayed margin identity.
    raw_margin_numerator = (
        (B**4 - (B - _poly(g0_numerator, variables))**4)
        * _poly(p_denominator, variables)
        - _poly(p_numerator, variables) * B**4
    )
    if raw_margin_numerator != t_poly * T * _poly(linear_positive, variables):
        raise AssertionError("first-boundary margin identity changed")

    order_gap = sp.factor(gamma4 - gamma3)
    order_numerator, order_denominator = map(
        sp.factor, sp.fraction(order_gap)
    )
    R = sp.exquo(
        _poly(order_numerator, variables),
        _poly(w - 1, variables),
    )
    G = -R
    if tuple(G.degree_list()) != ORDER_DEGREE or len(G.terms()) != 29:
        raise AssertionError("last-centroid order factor changed")
    expected_order_denominator = (
        4 * (t + 1) * (t + u) * (t * w + 1) * (t * w + u)
        * F.as_expr()
    )
    if sp.factor(order_denominator - expected_order_denominator) != 0:
        raise AssertionError("last-centroid order denominator changed")

    # R=(t^2*u^3-u^4)+Q with Q coefficientwise nonnegative.  Thus R<=0
    # forces t^2<=u<1 and puts the remaining scale parameter in the unit cube.
    confinement = R - _poly(t**2 * u**3 - u**4, variables)
    if len(confinement.terms()) != 27 or any(
        coefficient <= 0 for _, coefficient in confinement.terms()
    ):
        raise AssertionError("unit-cube confinement identity changed")

    quotient, remainder = sp.div(T, G)
    if quotient * G + remainder != T:
        raise AssertionError("multivariate division identity failed")
    if tuple(quotient.degree_list()) != QUOTIENT_DEGREE:
        raise AssertionError("quotient degree changed")
    if tuple(remainder.degree_list()) != REMAINDER_DEGREE:
        raise AssertionError("remainder degree changed")

    return {
        "variables": variables,
        "cap": p,
        "gamma0": gamma0,
        "order_denominator": order_denominator,
        "margin_denominator": (B**4 * U).as_expr(),
        "comparison": T,
        "order_factor": G,
        "confinement": confinement,
        "quotient": quotient,
        "remainder": remainder,
    }


def derive_three_below_polynomials():
    """Derive the three-below/two-above certificate polynomials exactly."""
    u, v, w, t = sp.symbols("u v w t")
    variables = (u, v, w, t)
    z = sp.symbols("z0:5")
    cap = sum(
        z[i] ** 4 / sp.prod(z[i] - z[j] for j in range(5) if j != i)
        for i in (3, 4)
    )
    substitution = {
        z[0]: -1,
        z[1]: -u,
        z[2]: -v,
        z[3]: t * w,
        z[4]: t,
    }
    density = sum(sp.diff(cap, coordinate) for coordinate in z)
    p = sp.factor(cap.subs(substitution))
    gamma0 = sp.factor(sp.diff(cap, z[0]).subs(substitution)
                       / density.subs(substitution))
    gamma3 = sp.factor(sp.diff(cap, z[3]).subs(substitution)
                       / density.subs(substitution))
    gamma4 = sp.factor(sp.diff(cap, z[4]).subs(substitution)
                       / density.subs(substitution))

    p_numerator, p_denominator = map(sp.factor, sp.fraction(p))
    g0_numerator, g0_denominator = map(sp.factor, sp.fraction(gamma0))
    linear_positive = (t + 1) * (t * w + 1)
    other_positive = (
        (t + u) * (t + v) * (t * w + u) * (t * w + v)
    )
    if sp.factor(p_denominator - linear_positive * other_positive) != 0:
        raise AssertionError("three-below cap denominator changed")
    if sp.rem(_poly(p_numerator, variables), _poly(t**3, variables)) != 0:
        raise AssertionError("three-below cap numerator lost its t^3 factor")
    if sp.rem(_poly(g0_numerator, variables), _poly(t, variables)) != 0:
        raise AssertionError("three-below centroid numerator lost its t factor")
    if sp.rem(
        _poly(g0_denominator, variables),
        _poly(4 * linear_positive, variables),
    ) != 0:
        raise AssertionError("three-below centroid denominator changed")
    F = sp.exquo(
        _poly(g0_denominator, variables),
        _poly(4 * linear_positive, variables),
    )
    if any(coefficient <= 0 for _, coefficient in F.terms()):
        raise AssertionError("three-below density factor lost positivity")

    A = sp.exquo(_poly(g0_numerator, variables), _poly(t, variables))
    B = _poly(g0_denominator, variables)
    U = _poly(other_positive, variables)
    t_poly = _poly(t, variables)
    first_part = (
        4 * B**3 * A
        - 6 * B**2 * t_poly * A**2
        + 4 * B * t_poly**2 * A**3
        - t_poly**3 * A**4
    ) * U
    b4_over_linear = sp.exquo(B**4, _poly(linear_positive, variables))
    T = first_part - sp.exquo(
        _poly(p_numerator, variables), t_poly
    ) * b4_over_linear
    if (tuple(T.degree_list()) != THREE_TWO_TARGET_DEGREE
            or len(T.terms()) != THREE_TWO_COMPARISON_TERMS):
        raise AssertionError("three-below comparison numerator changed")

    raw_margin_numerator = (
        (B**4 - (B - _poly(g0_numerator, variables))**4)
        * _poly(p_denominator, variables)
        - _poly(p_numerator, variables) * B**4
    )
    if raw_margin_numerator != t_poly * T * _poly(linear_positive, variables):
        raise AssertionError("three-below margin identity changed")

    order_gap = sp.factor(gamma4 - gamma3)
    order_numerator, order_denominator = map(
        sp.factor, sp.fraction(order_gap)
    )
    R = sp.exquo(
        _poly(order_numerator, variables),
        _poly(w - 1, variables),
    )
    G = -R
    if tuple(G.degree_list()) != THREE_TWO_ORDER_DEGREE or len(G.terms()) != 87:
        raise AssertionError("three-below order factor changed")
    expected_order_denominator = (
        4 * linear_positive * other_positive * F.as_expr()
    )
    if sp.factor(order_denominator - expected_order_denominator) != 0:
        raise AssertionError("three-below order denominator changed")

    # Orderedness gives R<=0.  To prove t<2, substitute t=2+rho.  The
    # coefficients are Bernstein-nonnegative in (u,v,w) and power-nonnegative
    # in rho; at least one rho^0 coefficient is positive on the open cube.
    rho = sp.symbols("rho")
    shifted_R = _poly(R.as_expr().subs(t, 2 + rho), (u, v, w, rho))

    tau = sp.symbols("tau")
    cube_variables = (u, v, w, tau)
    cube_T = _poly(T.as_expr().subs(t, 2 * tau), cube_variables)
    cube_G = _poly(G.as_expr().subs(t, 2 * tau), cube_variables)
    if tuple(cube_T.degree_list()) != THREE_TWO_TARGET_DEGREE:
        raise AssertionError("rescaled three-below target degree changed")
    if tuple(
        cube_T.degree_list()[index] - cube_G.degree_list()[index]
        for index in range(4)
    ) != THREE_TWO_MULTIPLIER_DEGREE:
        raise AssertionError("three-below multiplier degree changed")

    return {
        "variables": variables,
        "cube_variables": cube_variables,
        "cap": p,
        "gamma0": gamma0,
        "order_denominator": order_denominator,
        "margin_denominator": (B**4 * U).as_expr(),
        "comparison": T,
        "order_R": R,
        "order_factor": G,
        "shifted_R": shifted_R,
        "cube_comparison": cube_T,
        "cube_order_factor": cube_G,
    }


def exact_bernstein_coefficients(
    poly: sp.Poly,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray:
    """Return an exact Bernstein array, optionally transforming selected axes."""
    degree = tuple(int(value) for value in poly.degree_list())
    shape = tuple(value + 1 for value in degree)
    coefficients = np.empty(shape, dtype=object)
    coefficients.fill(Fraction(0))
    for powers, coefficient in poly.terms():
        coefficients[powers] = _fraction(coefficient)

    # In one variable, a_k*x^k has degree-d Bernstein coefficients
    # a_k*C(i,k)/C(d,k), i>=k.  Apply this exact transform on each chosen axis.
    selected_axes = tuple(range(len(degree))) if axes is None else axes
    for axis in reversed(selected_axes):
        axis_degree = degree[axis]
        moved = np.moveaxis(coefficients, axis, -1)
        rows = moved.reshape((-1, axis_degree + 1))
        transformed = np.empty_like(rows)
        ratios = tuple(
            tuple(
                Fraction(comb(index, power), comb(axis_degree, power))
                for power in range(index + 1)
            )
            for index in range(axis_degree + 1)
        )
        for row_index in range(rows.shape[0]):
            row = rows[row_index]
            for index in range(axis_degree + 1):
                transformed[row_index, index] = sum(
                    (row[power] * ratios[index][power]
                     for power in range(index + 1)),
                    Fraction(0),
                )
        coefficients = np.moveaxis(
            transformed.reshape(moved.shape), -1, axis
        )
    return coefficients


def exact_coefficient_statistics(
    coefficients: np.ndarray,
    degree: tuple[int, ...],
) -> dict[str, object]:
    """Audit and hash a dense exact coefficient array."""
    shape = tuple(value + 1 for value in degree)
    if coefficients.shape != shape:
        raise AssertionError("coefficient array has the wrong shape")

    positive = 0
    zero = 0
    negative = 0
    smallest_positive: Fraction | None = None
    largest: Fraction | None = None
    digest = hashlib.sha256()
    digest.update(("degree=" + ",".join(map(str, degree)) + "\n").encode())
    for index in np.ndindex(shape):
        value = coefficients[index]
        digest.update(
            (
                ",".join(map(str, index))
                + f":{value.numerator}/{value.denominator}\n"
            ).encode()
        )
        if value > 0:
            positive += 1
            if smallest_positive is None or value < smallest_positive:
                smallest_positive = value
            if largest is None or value > largest:
                largest = value
        elif value == 0:
            zero += 1
        else:
            negative += 1

    if negative:
        raise AssertionError(f"found {negative} negative Bernstein coefficients")
    if smallest_positive is None or largest is None:
        raise AssertionError("Bernstein expansion has no positive coefficient")
    if positive + zero != prod(shape):
        raise AssertionError("Bernstein coefficient count mismatch")
    return {
        "degree": list(degree),
        "total_coefficients": prod(shape),
        "positive_coefficients": positive,
        "zero_coefficients": zero,
        "negative_coefficients": negative,
        "smallest_positive": _fraction_record(smallest_positive),
        "largest": _fraction_record(largest),
        "sha256_all_coefficients": digest.hexdigest(),
    }


def exact_bernstein_statistics(poly: sp.Poly) -> dict[str, object]:
    """Convert to the exact tensor-Bernstein basis and audit every sign."""
    degree = tuple(int(value) for value in poly.degree_list())
    coefficients = exact_bernstein_coefficients(poly)
    return exact_coefficient_statistics(coefficients, degree)


def _product_entry(
    order_coefficients: np.ndarray,
    target_degree: tuple[int, ...],
    multiplier_degree: tuple[int, ...],
    target_index: tuple[int, ...],
    multiplier_index: tuple[int, ...],
) -> Fraction:
    """Coefficient of one multiplier basis term in one product basis term."""
    order_degree = tuple(value - multiplier_degree[index]
                         for index, value in enumerate(target_degree))
    order_index = tuple(
        target_index[index] - multiplier_index[index]
        for index in range(len(target_degree))
    )
    if any(
        value < 0 or value > order_degree[index]
        for index, value in enumerate(order_index)
    ):
        return Fraction(0)
    order_value = order_coefficients[order_index]
    if order_value == 0:
        return Fraction(0)
    factor = Fraction(1)
    for index in range(len(target_degree)):
        factor *= Fraction(
            comb(multiplier_degree[index], multiplier_index[index])
            * comb(order_degree[index], order_index[index]),
            comb(target_degree[index], target_index[index]),
        )
    return order_value * factor


def _to_fmpq(value: Fraction) -> fmpq:
    return fmpq(value.numerator, value.denominator)


def _from_fmpq(value: fmpq) -> Fraction:
    return Fraction(str(value))


def _index_digest(values: tuple[int, ...]) -> str:
    digest = hashlib.sha256()
    digest.update((",".join(map(str, values)) + "\n").encode())
    return digest.hexdigest()


def verify_three_below_certificate(derived) -> dict[str, object]:
    """Solve the sparse exact witness and check the complete residual."""
    shifted_R = derived["shifted_R"]
    shifted_degree = tuple(map(int, shifted_R.degree_list()))
    shifted_coefficients = exact_bernstein_coefficients(
        shifted_R, axes=(0, 1, 2)
    )
    shifted_statistics = exact_coefficient_statistics(
        shifted_coefficients, shifted_degree
    )
    if not any(
        shifted_coefficients[index] > 0
        for index in np.ndindex(shifted_coefficients.shape[:-1])
        for index in [(*index, 0)]
    ):
        raise AssertionError("t=2 confinement boundary lost strict positivity")

    target = derived["cube_comparison"]
    order_factor = derived["cube_order_factor"]
    target_degree = tuple(map(int, target.degree_list()))
    order_degree = tuple(map(int, order_factor.degree_list()))
    multiplier_degree = tuple(THREE_TWO_MULTIPLIER_DEGREE)
    if tuple(
        order_degree[index] + multiplier_degree[index]
        for index in range(4)
    ) != target_degree:
        raise AssertionError("three-below product degree mismatch")

    target_coefficients = exact_bernstein_coefficients(target)
    order_coefficients = exact_bernstein_coefficients(order_factor)
    target_shape = tuple(value + 1 for value in target_degree)
    multiplier_shape = tuple(value + 1 for value in multiplier_degree)
    if len(THREE_TWO_SUPPORT_FLAT) != len(THREE_TWO_ACTIVE_ROWS_FLAT):
        raise AssertionError("three-below witness is not square")
    if len(set(THREE_TWO_SUPPORT_FLAT)) != len(THREE_TWO_SUPPORT_FLAT):
        raise AssertionError("three-below multiplier support is not unique")
    if len(set(THREE_TWO_ACTIVE_ROWS_FLAT)) != len(THREE_TWO_ACTIVE_ROWS_FLAT):
        raise AssertionError("three-below active rows are not unique")
    if not all(0 <= value < prod(multiplier_shape)
               for value in THREE_TWO_SUPPORT_FLAT):
        raise AssertionError("three-below multiplier index is out of range")
    if not all(0 <= value < prod(target_shape)
               for value in THREE_TWO_ACTIVE_ROWS_FLAT):
        raise AssertionError("three-below row index is out of range")

    multiplier_indices = tuple(
        tuple(map(int, np.unravel_index(value, multiplier_shape)))
        for value in THREE_TWO_SUPPORT_FLAT
    )
    target_indices = tuple(
        tuple(map(int, np.unravel_index(value, target_shape)))
        for value in THREE_TWO_ACTIVE_ROWS_FLAT
    )
    matrix = [
        [
            _product_entry(
                order_coefficients,
                target_degree,
                multiplier_degree,
                target_index,
                multiplier_index,
            )
            for multiplier_index in multiplier_indices
        ]
        for target_index in target_indices
    ]
    right_hand_side = [target_coefficients[index] for index in target_indices]
    exact_matrix = fmpq_mat(
        [[_to_fmpq(value) for value in row] for row in matrix]
    )
    exact_rhs = fmpq_mat([[_to_fmpq(value)] for value in right_hand_side])
    solution = exact_matrix.solve(exact_rhs)
    multiplier_values = tuple(
        _from_fmpq(solution[index, 0])
        for index in range(len(multiplier_indices))
    )
    if any(value <= 0 for value in multiplier_values):
        raise AssertionError("three-below exact multiplier lost positivity")

    multiplier_coefficients = np.empty(multiplier_shape, dtype=object)
    multiplier_coefficients.fill(Fraction(0))
    for index, value in zip(multiplier_indices, multiplier_values, strict=True):
        multiplier_coefficients[index] = value

    product_coefficients = np.empty(target_shape, dtype=object)
    product_coefficients.fill(Fraction(0))
    nonzero_order = tuple(
        (index, order_coefficients[index])
        for index in np.ndindex(order_coefficients.shape)
        if order_coefficients[index] != 0
    )
    for multiplier_index, multiplier_value in zip(
        multiplier_indices, multiplier_values, strict=True
    ):
        for order_index, order_value in nonzero_order:
            target_index = tuple(
                multiplier_index[axis] + order_index[axis]
                for axis in range(4)
            )
            factor = Fraction(1)
            for axis in range(4):
                factor *= Fraction(
                    comb(multiplier_degree[axis], multiplier_index[axis])
                    * comb(order_degree[axis], order_index[axis]),
                    comb(target_degree[axis], target_index[axis]),
                )
            product_coefficients[target_index] += (
                multiplier_value * order_value * factor
            )

    residual_coefficients = np.empty(target_shape, dtype=object)
    for index in np.ndindex(target_shape):
        residual_coefficients[index] = (
            target_coefficients[index] - product_coefficients[index]
        )
    if any(
        residual_coefficients.flat[index] != 0
        for index in THREE_TWO_ACTIVE_ROWS_FLAT
    ):
        raise AssertionError("three-below active subsystem was not solved exactly")

    return {
        "shifted_order_statistics": shifted_statistics,
        "multiplier_statistics": exact_coefficient_statistics(
            multiplier_coefficients, multiplier_degree
        ),
        "residual_statistics": exact_coefficient_statistics(
            residual_coefficients, target_degree
        ),
        "linear_system_order": len(multiplier_indices),
        "support_digest": _index_digest(tuple(THREE_TWO_SUPPORT_FLAT)),
        "active_rows_digest": _index_digest(tuple(THREE_TWO_ACTIVE_ROWS_FLAT)),
    }


def _polynomial_record(poly: sp.Poly) -> dict[str, object]:
    return {
        "degree": list(map(int, poly.degree_list())),
        "power_terms": len(poly.terms()),
        "sha256_power_expansion": _power_digest(poly),
    }


def build_certificate() -> dict[str, object]:
    two_below = derive_polynomials()
    quotient = two_below["quotient"]
    remainder = two_below["remainder"]
    quotient_bernstein = exact_bernstein_statistics(quotient)
    remainder_bernstein = exact_bernstein_statistics(remainder)
    three_below = derive_three_below_polynomials()
    three_below_verified = verify_three_below_certificate(three_below)
    return {
        "schema_version": 2,
        "mathematical_status": (
            "ordinary reduction with an exact rational polynomial certificate"
        ),
        "scope": {
            "hypothesis": (
                "Every regular section has five distinct ordered knots and "
                "an ordered section centroid."
            ),
            "conclusion": (
                "The section satisfies the local maximum barrier. In both "
                "mixed crossing types, the stronger inequality "
                "p<=I_{gamma_0}(1,4)=1-(1-gamma_0)^4 holds."
            ),
            "nonclaim": (
                "Repeated-knot block faces are not covered, so this is not "
                "yet a complete five-coordinate vertex-maximization theorem."
            ),
        },
        "endpoint_crossings": {
            "four_below_one_above": "AM--GM gives p<=I_{M_4}(4,1).",
            "one_below_four_above": (
                "Knot order and centroid order force the four upper knots "
                "to coincide, so this crossing has no distinct-knot instance."
            ),
        },
        "two_below_three_above": {
            "normalized_knots": ["-1", "-u", "t*v", "t*w", "t"],
            "domain": "0<u<1, 0<v<w<1, t>0",
            "hypothesis_used": "gamma_4>=gamma_3",
            "identities": {
                "order_gap": (
                    "gamma_4-gamma_3=(w-1)*R/Delta_order, with "
                    "Delta_order>0 and G=-R"
                ),
                "unit_cube": (
                    "R=t^2*u^3-u^4+Q, where every power coefficient of Q "
                    "is positive; hence R<=0 implies t^2<=u<1"
                ),
                "comparison": (
                    "I_{gamma_0}(1,4)-p=t*T/Delta_margin, with "
                    "Delta_margin>0"
                ),
                "division": "T=S*G+H",
                "logic": "G,S,H>=0 on the unit cube imply T>=0",
            },
            "polynomials": {
                "comparison_T": _polynomial_record(two_below["comparison"]),
                "order_G": _polynomial_record(two_below["order_factor"]),
                "confinement_Q": _polynomial_record(two_below["confinement"]),
                "quotient_S": _polynomial_record(quotient),
                "remainder_H": _polynomial_record(remainder),
            },
            "bernstein_basis": {
                "variables": list(VARIABLE_NAMES),
                "domain": "[0,1]^4",
                "quotient_S": quotient_bernstein,
                "remainder_H": remainder_bernstein,
            },
        },
        "three_below_two_above": {
            "normalized_knots": ["-1", "-u", "-v", "t*w", "t"],
            "domain": "0<v<u<1, 0<w<1, t>0",
            "hypothesis_used": "gamma_4>=gamma_3",
            "identities": {
                "order_gap": (
                    "gamma_4-gamma_3=(w-1)*R/Delta_order, with "
                    "Delta_order>0 and G=-R"
                ),
                "compactness": (
                    "The mixed Bernstein/power expansion of R(u,v,w,2+rho) "
                    "is nonnegative and strictly positive in the open "
                    "parameter domain; hence R<=0 implies t<2."
                ),
                "comparison": (
                    "I_{gamma_0}(1,4)-p=t*T/Delta_margin, with "
                    "Delta_margin>0"
                ),
                "rescaling": "t=2*tau maps the ordered domain into [0,1]^4",
                "certificate": "T_tilde=S*G_tilde+H",
                "logic": "G_tilde,S,H>=0 imply T_tilde>=0",
            },
            "polynomials": {
                "comparison_T_before_rescaling": _polynomial_record(
                    three_below["comparison"]
                ),
                "order_G_before_rescaling": _polynomial_record(
                    three_below["order_factor"]
                ),
                "shifted_R_for_t_at_least_2": _polynomial_record(
                    three_below["shifted_R"]
                ),
            },
            "compactness_sign_audit": {
                "basis": "Bernstein in (u,v,w), power basis in rho",
                **three_below_verified["shifted_order_statistics"],
            },
            "exact_linear_witness": {
                "matrix_order": three_below_verified["linear_system_order"],
                "multiplier_support_flat_index_sha256": (
                    three_below_verified["support_digest"]
                ),
                "active_target_rows_flat_index_sha256": (
                    three_below_verified["active_rows_digest"]
                ),
                "solution_method": (
                    "exact nonsingular rational solve; the resulting "
                    "multiplier is then checked against every target row"
                ),
            },
            "bernstein_basis": {
                "variables": ["u", "v", "w", "tau"],
                "domain": "[0,1]^4",
                "multiplier_S": three_below_verified["multiplier_statistics"],
                "residual_H": three_below_verified["residual_statistics"],
            },
        },
        "arithmetic": {
            "derivation_and_division": "SymPy exact rational polynomial arithmetic",
            "sparse_linear_witness": "FLINT exact rational matrix solve",
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
    first = certificate["two_below_three_above"]["bernstein_basis"]
    second = certificate["three_below_two_above"]["bernstein_basis"]
    print(
        "five-coordinate distinct-knot sections certified: "
        f"2--3 S/H {first['quotient_S']['total_coefficients']}/"
        f"{first['remainder_H']['total_coefficients']}; "
        f"3--2 S/H {second['multiplier_S']['total_coefficients']}/"
        f"{second['residual_H']['total_coefficients']} exact coefficients"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
