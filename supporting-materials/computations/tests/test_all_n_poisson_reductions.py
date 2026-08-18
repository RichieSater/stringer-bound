"""Tests for the exact analytic reductions in the all-n Poisson program."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from all_n_poisson_reductions import (  # noqa: E402
    check_equal_weight_hessian_thresholds,
    check_three_exponential_equal_smaller_line,
    check_three_exponential_reduction,
    check_three_exponential_repeated_max_boundary,
    check_two_exponential_global_convexity,
    check_two_exponential_obstruction,
)


class AllNPoissonReductionTests(unittest.TestCase):
    def test_two_exponential_local_threshold_identity(self):
        check_two_exponential_obstruction()

    def test_two_exponential_global_convexity_factorization(self):
        check_two_exponential_global_convexity()

    def test_equal_weight_hessian_thresholds(self):
        check_equal_weight_hessian_thresholds()

    def test_three_exponential_reduction(self):
        check_three_exponential_reduction()

    def test_three_exponential_repeated_max_boundary(self):
        check_three_exponential_repeated_max_boundary()

    def test_three_exponential_equal_smaller_line(self):
        check_three_exponential_equal_smaller_line()


if __name__ == "__main__":
    unittest.main()
