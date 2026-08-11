"""Tests for the open Dirichlet--Poissonization research certificate."""

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from dirichlet_poissonization import (  # noqa: E402
    build_certificate,
    equal_block_check,
    saffine_localization_obstruction,
    verify_saffine_symbolic_identities,
)


class DirichletPoissonizationTests(unittest.TestCase):
    def test_saffine_closed_forms_are_symbolically_derived(self):
        verify_saffine_symbolic_identities()

    def test_saffine_obstruction_preserves_mean_and_has_opposite_tail_signs(self):
        result = saffine_localization_obstruction()
        mean_record = result["law"]["mean"]
        mean = Fraction(
            int(mean_record["numerator"]),
            int(mean_record["denominator"]),
        )
        tail_record = result["tail_probability_P_Y_gt_1"]
        tail = Fraction(
            int(tail_record["numerator"]),
            int(tail_record["denominator"]),
        )
        self.assertLess(mean, Fraction(10, 11))
        self.assertGreater(tail, Fraction(7, 15))
        self.assertEqual(
            result["exact_comparisons"]["P_S_Y_gt_10"],
            "< 7/15",
        )

    def test_equal_block_boundary_and_strict_cases(self):
        boundary = equal_block_check(10, 4, Fraction(4))
        strict = equal_block_check(10, 4, Fraction(7))
        self.assertEqual(boundary["lambda"]["numerator"], "4")
        self.assertEqual(strict["lambda"]["numerator"], "7")
        self.assertFalse(
            boundary["certified_lower_margin"]["decimal"].startswith("-")
        )
        self.assertFalse(
            strict["certified_lower_margin"]["decimal"].startswith("-")
        )

    def test_certificate_is_explicitly_non_theorem_research_support(self):
        certificate = build_certificate()
        self.assertIn("not an all-sample-size coverage certificate", certificate["status"])
        self.assertEqual(
            len(certificate["equal_block_profiles"]["regression_checks"]),
            75,
        )


if __name__ == "__main__":
    unittest.main()
