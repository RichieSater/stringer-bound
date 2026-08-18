"""Tests for the exact analytic reductions in the all-n Poisson program."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from all_n_poisson_reductions import (  # noqa: E402
    check_equal_weight_hessian_thresholds,
    check_three_exponential_axis_transversality,
    check_three_exponential_diagonal_transversality,
    check_three_exponential_equal_smaller_line,
    check_three_exponential_infinite_gap_boundary,
    check_three_exponential_reduction,
    check_three_exponential_repeated_max_boundary,
    check_three_exponential_sharp_corner_expansion,
    check_three_exponential_small_gap_region,
    check_three_exponential_trace_identity,
    check_three_exponential_two_large_gap_region,
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

    def test_three_exponential_trace_identity(self):
        check_three_exponential_trace_identity()

    def test_three_exponential_repeated_max_boundary(self):
        check_three_exponential_repeated_max_boundary()

    def test_three_exponential_axis_transversality(self):
        check_three_exponential_axis_transversality()

    def test_three_exponential_equal_smaller_line(self):
        check_three_exponential_equal_smaller_line()

    def test_three_exponential_diagonal_transversality(self):
        check_three_exponential_diagonal_transversality()

    def test_three_exponential_small_gap_region(self):
        check_three_exponential_small_gap_region()

    def test_three_exponential_two_large_gap_region(self):
        check_three_exponential_two_large_gap_region()

    def test_three_exponential_infinite_gap_boundary(self):
        check_three_exponential_infinite_gap_boundary()

    def test_three_exponential_sharp_corner_expansion(self):
        check_three_exponential_sharp_corner_expansion()


if __name__ == "__main__":
    unittest.main()
