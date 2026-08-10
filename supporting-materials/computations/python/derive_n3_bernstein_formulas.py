"""Derive the symbolic Bernstein formulas used by the n=3 certificate.

This is intentionally separate from :mod:`n3_gaffke_certificate`: the latter
performs quick exact sign checks, while this script starts from the
middle-region simplex-cap volume, differentiates it symbolically, carries out
the two affine triangle maps, and regenerates every coefficient formula.

The dependency versions are pinned, so comparing this output with the
committed formula file checks that no unverified transcription has entered
the certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp


def _bernstein_coefficients(polynomial, U, W, degree):
    power = sp.Poly(sp.cancel(polynomial), U, W)
    output = []
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            k = degree - i - j
            value = 0
            for (p, q), coefficient in power.terms():
                if i >= p and j >= q:
                    value += coefficient * sp.Rational(
                        math.comb(i, p) * math.comb(j, q),
                        math.comb(degree, p + q) * math.comb(p + q, p),
                    )
            output.append({
                "index": [i, j, k],
                "expression": str(sp.factor(value)),
            })
    return output


def derive():
    x, y, b, c, d, rho, tau = sp.symbols(
        "x y b c d rho tau")
    A = 1 - b - c
    a = 1 - b - c - d

    threshold = b * x + c * y + d
    volume = (
        (1 - threshold) ** 3 / ((1 - x) * (1 - y))
        - (y - threshold) ** 3 / (y * (y - x) * (1 - y))
    )
    numerator, denominator = sp.together(sp.diff(volume, y)).as_numer_denom()
    # denominator = y^2(1-x)(y-x)^2(1-y)^2.  Removing (1-y)^2
    # leaves Q with the positive denominator y^2(1-x)(y-x)^2.
    Q = sp.cancel(numerator / (1 - y) ** 2)

    x_rt = (d - (1 - c) * rho + c * tau) / A
    y_rt = x_rt + rho + tau
    # A is positive.  Multiplication by A^3 clears the remaining
    # denominators and leaves the sign unchanged.
    Q_rt = sp.Poly(
        sp.cancel(Q.subs({x: x_rt, y: y_rt}) * A ** 3),
        rho, tau,
    ).as_expr()

    U, W = sp.symbols("U W")
    rho_1 = d / (1 - c)
    tau_3 = a / (1 - b)
    rho_2 = c + d
    tau_2 = a + b
    triangle_maps = [
        (rho_1 * U + rho_2 * W, tau_2 * W),
        (rho_2 * U, tau_2 * U + tau_3 * W),
    ]
    triangles = [
        _bernstein_coefficients(
            sp.cancel(Q_rt.subs({rho: rho_map, tau: tau_map})),
            U, W, 5,
        )
        for rho_map, tau_map in triangle_maps
    ]
    return {
        "schema_version": 1,
        "quantity": (
            "Degree-five triangular Bernstein coefficients of (1-b-c)^3 "
            "Q, where Q=y^2(1-x)(y-x)^2 times the y-derivative of the "
            "middle-region simplex-cap volume."
        ),
        "triangle_maps": [
            "(rho,tau)=(d/(1-c) U + (c+d) W, (1-c-d) W)",
            "(rho,tau)=((c+d) U, (1-c-d) U + "
            "(1-b-c-d)/(1-b) W)",
        ],
        "degree": 5,
        "triangles": triangles,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    args.out.write_text(json.dumps(derive(), indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
