"""Regression tests for the one-upper-knot certificate."""

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from one_cap_certificate import certify_alpha  # noqa: E402


class OneCapCertificateTests(unittest.TestCase):
    def test_small_exact_range(self):
        result = certify_alpha("0.05", n_max=4, factor_bits=32)
        self.assertEqual(
            result["strict_nonterminal_inequalities_checked"], 6
        )
        self.assertEqual(result["analytic_terminal_equalities"], 4)
        worst = result["worst_nonterminal_case"]
        upper = Fraction(
            int(worst["certified_upper_bound"]["numerator"]),
            int(worst["certified_upper_bound"]["denominator"]),
        )
        self.assertLess(upper, Fraction(1, 20))
        width = Fraction(
            int(result["maximum_factor_bracket_width"]["numerator"]),
            int(result["maximum_factor_bracket_width"]["denominator"]),
        )
        self.assertLessEqual(width, Fraction(1, 1 << 32))


if __name__ == "__main__":
    unittest.main()
