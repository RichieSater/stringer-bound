"""Regression tests for the analytic all-n one-cap proof checks."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from one_cap_all_n_check import (  # noqa: E402
    exact_constant_checks,
    numerical_regression,
    symbolic_identity_checks,
)


class OneCapAllNCheckTests(unittest.TestCase):
    def test_exact_and_symbolic_checks(self):
        self.assertIn("e_upper", exact_constant_checks())
        self.assertIn("boundary_tail", symbolic_identity_checks())

    def test_small_numerical_regression(self):
        rows = numerical_regression(n_max=40)
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(
            row["minimum_tail_divided_by_alpha"] >= 1 for row in rows
        ))


if __name__ == "__main__":
    unittest.main()
