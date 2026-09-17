"""Tests for the sharp n=4 Stringer-weight monotonicity certificate."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent
    / "certificates"
    / "n4-weight-monotonicity-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from n4_weight_monotonicity import (  # noqa: E402
    build_certificate,
    derive_certificate_data,
)


class N4WeightMonotonicityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.expected = json.loads(CERTIFICATE.read_text())
        cls.regenerated = build_certificate()
        cls.derived = derive_certificate_data()

    def test_checked_in_certificate_is_exactly_regenerated(self):
        self.assertEqual(self.regenerated, self.expected)

    def test_terminal_endpoint_bracket(self):
        endpoint = self.expected["critical_endpoint"]
        self.assertEqual(
            endpoint["alpha_display_interval"],
            "[0.3214292869970077, 0.3214292869970078]",
        )
        lower = self.derived["r_interval"][0]
        upper = self.derived["r_interval"][1]
        cubic = self.derived["nontrivial"]
        self.assertLess(cubic.eval(lower), 0)
        self.assertGreater(cubic.eval(upper), 0)
        self.assertLess(upper**4 - lower**4, Fraction(1, 10**16))
        alpha_lower = Fraction(
            int(endpoint["alpha_lower"]["numerator"]),
            int(endpoint["alpha_lower"]["denominator"]),
        )
        alpha_upper = Fraction(
            int(endpoint["alpha_upper"]["numerator"]),
            int(endpoint["alpha_upper"]["denominator"]),
        )
        self.assertLess(Fraction(3214292869970077, 10**16), alpha_lower)
        self.assertLess(alpha_upper, Fraction(3214292869970078, 10**16))
        self.assertTrue(endpoint["alpha_display_interval_verified_exactly"])

    def test_resultant_factor_and_exact_root_counts(self):
        result = self.expected["middle_spacing_resultant"]
        self.assertEqual(result["degree"], 8)
        self.assertEqual(result["roots_in_open_unit_interval"], 2)
        self.assertEqual(
            result["sha256"],
            "61803877c3ca1d8d91fd9e3251a733c5f6194874646cb0ed57b33301aeca2983",
        )
        self.assertEqual(
            [entry["root_count"] for entry in result["isolating_intervals"]],
            [1, 1],
        )

    def test_domain_exclusion_identities(self):
        self.assertEqual(
            str(self.derived["q_one_identity"]),
            "(x - 1)**3*(19*x - 11)/16",
        )
        self.assertEqual(
            self.derived["resultant"],
            144
            * self.derived["relation"].gens[0] ** 2
            * (self.derived["relation"].gens[0] - 1) ** 6
            * self.derived["resultant_factor"].as_expr(),
        )


if __name__ == "__main__":
    unittest.main()
