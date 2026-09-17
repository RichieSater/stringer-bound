"""Tests for the exact five-coordinate mixed-section certificate."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent / "certificates"
    / "five-coordinate-mixed-section-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from five_coordinate_mixed_section import (  # noqa: E402
    build_certificate,
    derive_polynomials,
    derive_three_below_polynomials,
)


class FiveCoordinateMixedSectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected = json.loads(CERTIFICATE.read_text())
        cls.regenerated = build_certificate()

    def test_checked_in_certificate_is_exactly_regenerated(self):
        self.assertEqual(self.regenerated, self.expected)

    def test_every_bernstein_coefficient_has_certified_sign(self):
        basis = self.expected["two_below_three_above"]["bernstein_basis"]
        quotient = basis["quotient_S"]
        remainder = basis["remainder_H"]
        self.assertEqual(quotient["positive_coefficients"], 151230)
        self.assertEqual(quotient["zero_coefficients"], 10690)
        self.assertEqual(quotient["negative_coefficients"], 0)
        self.assertEqual(remainder["positive_coefficients"], 80067)
        self.assertEqual(remainder["zero_coefficients"], 21933)
        self.assertEqual(remainder["negative_coefficients"], 0)
        self.assertEqual(
            quotient["sha256_all_coefficients"],
            "c469dca6c885b078c3229de75c0a27902a181c3fd10a3090c7ed09416c03f0fd",
        )
        self.assertEqual(
            remainder["sha256_all_coefficients"],
            "fb5412edbad29673a576c5d52e6876fa32e336cffe0f69aa03b267a078e60454",
        )

        three_below = self.expected["three_below_two_above"]
        compactness = three_below["compactness_sign_audit"]
        multiplier = three_below["bernstein_basis"]["multiplier_S"]
        residual = three_below["bernstein_basis"]["residual_H"]
        self.assertEqual(compactness["positive_coefficients"], 366)
        self.assertEqual(compactness["negative_coefficients"], 0)
        self.assertEqual(multiplier["positive_coefficients"], 363)
        self.assertEqual(multiplier["zero_coefficients"], 1925)
        self.assertEqual(multiplier["negative_coefficients"], 0)
        self.assertEqual(residual["positive_coefficients"], 7526)
        self.assertEqual(residual["zero_coefficients"], 7174)
        self.assertEqual(residual["negative_coefficients"], 0)
        self.assertEqual(
            multiplier["sha256_all_coefficients"],
            "a35201138de220041db5e1e01a57c4a95f75ffe0310ee07ae5403e27cfe75494",
        )
        self.assertEqual(
            residual["sha256_all_coefficients"],
            "7c45cb9822cb7f4b2403803af6747994fd56056d346fedeafad397448238b34c",
        )

    def test_exact_rational_section_fixture(self):
        derived = derive_polynomials()
        u, v, w, t = derived["variables"]
        point = {
            u: Fraction(1, 2),
            v: Fraction(1, 5),
            w: Fraction(3, 5),
            t: Fraction(1, 5),
        }
        cap = Fraction(derived["cap"].subs(point))
        gamma0 = Fraction(derived["gamma0"].subs(point))
        G = Fraction(derived["order_factor"].as_expr().subs(point))
        order_denominator = Fraction(derived["order_denominator"].subs(point))
        last_gap = (1 - point[w]) * G / order_denominator
        first_step_cap = 1 - (1 - gamma0) ** 4

        self.assertEqual(cap, Fraction(8923, 87048))
        self.assertGreater(G, 0)
        self.assertGreater(last_gap, 0)
        self.assertLess(cap, first_step_cap)

    def test_exact_rational_three_below_fixture(self):
        derived = derive_three_below_polynomials()
        u, v, w, t = derived["variables"]
        point = {
            u: Fraction(1, 2),
            v: Fraction(1, 5),
            w: Fraction(3, 5),
            t: Fraction(1, 5),
        }
        cap = Fraction(derived["cap"].subs(point))
        gamma0 = Fraction(derived["gamma0"].subs(point))
        G = Fraction(derived["order_factor"].as_expr().subs(point))
        order_denominator = Fraction(derived["order_denominator"].subs(point))
        last_gap = (1 - point[w]) * G / order_denominator

        self.assertEqual(cap, Fraction(997, 20832))
        self.assertEqual(gamma0, Fraction(2743, 56112))
        self.assertGreater(G, 0)
        self.assertGreater(last_gap, 0)
        self.assertLess(cap, 1 - (1 - gamma0) ** 4)


if __name__ == "__main__":
    unittest.main()
