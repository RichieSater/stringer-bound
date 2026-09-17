"""Tests for the sharp n=5 Stringer-weight monotonicity certificate."""

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
    / "n5-weight-monotonicity-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from n5_weight_monotonicity import (  # noqa: E402
    build_certificate,
    derive_certificate_data,
)


class N5WeightMonotonicityTests(unittest.TestCase):
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
            "[0.3141689898050253, 0.3141689898050254]",
        )
        lower, upper = self.derived["r_interval"]
        quartic = self.derived["nontrivial"]
        self.assertLess(quartic.eval(lower), 0)
        self.assertGreater(quartic.eval(upper), 0)
        self.assertLess(upper**5 - lower**5, Fraction(1, 10**16))
        alpha_lower = Fraction(
            int(endpoint["alpha_lower"]["numerator"]),
            int(endpoint["alpha_lower"]["denominator"]),
        )
        alpha_upper = Fraction(
            int(endpoint["alpha_upper"]["numerator"]),
            int(endpoint["alpha_upper"]["denominator"]),
        )
        self.assertLess(Fraction(3141689898050253, 10**16), alpha_lower)
        self.assertLess(alpha_upper, Fraction(3141689898050254, 10**16))
        self.assertTrue(endpoint["alpha_display_interval_verified_exactly"])

    def test_d3_resultant_factor_and_root_count(self):
        result = self.expected["d3_resultant"]
        self.assertEqual(result["degree"], 11)
        self.assertEqual(result["roots_in_open_unit_interval"], 1)
        self.assertEqual(
            result["sha256"],
            "88cd52e2df765e4e144e8e3a3af663e9401472c4d0d63a712ad50a91dde64000",
        )
        self.assertEqual(result["isolating_interval"]["root_count"], 1)

    def test_d2_resultant_factors_and_root_count(self):
        result = self.expected["d2_resultant"]
        self.assertEqual(result["half_factor"]["degree"], 5)
        self.assertEqual(result["extraneous_factor"]["degree"], 8)
        self.assertEqual(
            result["half_factor"]["sha256"],
            "a6b1e112f5e9615b992f8f862812d45d2cc6017f35e7502cb8d22880378a400b",
        )
        self.assertEqual(
            result["extraneous_factor"]["sha256"],
            "73fa013ec096ef0580a30b227d5239d6207d95d73a42f1c1024c22796f66ac06",
        )
        self.assertEqual(
            result["extraneous_factor"]["roots_in_open_unit_interval"], 1
        )
        self.assertEqual(
            result["extraneous_factor"]["isolating_interval"]["root_count"],
            1,
        )

    def test_exact_elimination_and_boundary_identities(self):
        x = self.derived["d3_relation"].gens[0]
        self.assertEqual(
            self.derived["d3_resultant"],
            400000
            * x**2
            * (x - 1) ** 12
            * self.derived["d3_factor"].as_expr(),
        )
        self.assertEqual(
            self.derived["d2_resultant"],
            1600000
            * x**6
            * (x - 1) ** 6
            * self.derived["half_factor"].as_expr()
            * self.derived["extraneous_factor"].as_expr(),
        )
        self.assertEqual(
            str(self.derived["q_one_identity"]),
            "(x - 1)**3*(67*x**2 - 39*x - 8)/16",
        )

    def test_rational_sign_witnesses(self):
        self.assertEqual(
            tuple(value > 0 for value in self.derived["low_checks"].values()),
            (True, False, False, True, False),
        )
        self.assertEqual(
            tuple(value > 0 for value in self.derived["high_checks"].values()),
            (False, True),
        )


if __name__ == "__main__":
    unittest.main()
