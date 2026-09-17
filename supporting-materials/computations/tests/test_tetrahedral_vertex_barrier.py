"""Tests for the exact four-coordinate section-centroid barrier."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent / "certificates"
    / "tetrahedral-vertex-barrier-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from tetrahedral_vertex_barrier import (  # noqa: E402
    build_certificate,
    derive_polynomials,
)


class TetrahedralVertexBarrierTests(unittest.TestCase):
    def test_checked_in_certificate_is_exactly_regenerated(self):
        expected = json.loads(CERTIFICATE.read_text())
        self.assertEqual(build_certificate(), expected)

    def test_certificate_closes_both_bernstein_sign_problems(self):
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["multiplier"]["positive_coefficients"], 79)
        self.assertEqual(data["multiplier"]["negative_coefficients"], 0)
        self.assertEqual(data["residual"]["positive_coefficients"], 631)
        self.assertEqual(data["residual"]["zero_coefficients"], 569)
        boundary = data["repeated_upper_boundary"]
        self.assertEqual(boundary["residual_positive_coefficients"], 19)
        self.assertEqual(boundary["residual_zero_coefficients"], 11)
        self.assertEqual(boundary["residual_negative_coefficients"], 0)
        self.assertEqual(
            data["n3_stringer_application"]["alpha_star"],
            "((19+sqrt(21))/34)^3",
        )

    def test_exact_rational_section_fixture(self):
        variables, p, gamma, R, _, _ = derive_polynomials()
        u, v, t = variables
        point = {u: Fraction(1, 2), v: Fraction(2, 3), t: Fraction(1, 5)}
        cap = Fraction(p.subs(point))
        centroid = [Fraction(value.subs(point)) for value in gamma]
        order_factor = Fraction(R.as_expr().subs(point))
        middle_mass = centroid[0] + centroid[1]
        facet_cap = 3 * middle_mass**2 - 2 * middle_mass**3

        self.assertEqual(sum(centroid), 1)
        self.assertLess(order_factor, 0)
        self.assertLess(cap, facet_cap)
        self.assertLessEqual(centroid[2], centroid[3])


if __name__ == "__main__":
    unittest.main()
