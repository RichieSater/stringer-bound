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
    two_level_profile_regression,
    verify_n2_three_knot_theorem,
    verify_n3_four_knot_theorem,
    verify_radial_symbolic_identities,
    verify_saffine_symbolic_identities,
    verify_three_positive_convex_core,
    verify_three_positive_far_cap,
    verify_three_positive_middle_knot_region,
    verify_two_positive_knot_theorem,
    verify_two_level_symbolic_identities,
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

    def test_radial_zero_knot_algebra(self):
        verify_radial_symbolic_identities()

    def test_two_level_algebra_and_strict_sum_reduction(self):
        verify_two_level_symbolic_identities()
        result = two_level_profile_regression(
            5,
            3,
            Fraction(1, 3),
            Fraction(5, 4),
        )
        self.assertGreater(
            float(result["coefficient_sum_slack"]["decimal"]),
            0.0,
        )
        self.assertGreater(
            float(result["radial_path"]["endpoint_sum_slack"]["decimal"]),
            0.0,
        )

    def test_complete_n2_three_knot_algebra(self):
        verify_n2_three_knot_theorem()

    def test_complete_n3_four_knot_algebra(self):
        verify_n3_four_knot_theorem()

    def test_all_n_three_positive_convex_core_algebra(self):
        verify_three_positive_convex_core()

    def test_all_n_three_positive_far_cap_algebra(self):
        verify_three_positive_far_cap()

    def test_all_n_three_positive_middle_knot_algebra(self):
        verify_three_positive_middle_knot_region()

    def test_all_n_two_positive_knot_algebra(self):
        verify_two_positive_knot_theorem()

    def test_certificate_is_explicitly_non_theorem_research_support(self):
        certificate = build_certificate()
        self.assertEqual(certificate["schema_version"], 10)
        self.assertIn("not an all-sample-size coverage certificate", certificate["status"])
        self.assertEqual(
            len(certificate["equal_block_profiles"]["regression_checks"]),
            75,
        )
        self.assertEqual(
            certificate["radial_zero_knot_reduction"][
                "symbolic_identity_check"
            ],
            "passed",
        )
        self.assertEqual(
            len(
                certificate["two_level_profiles"][
                    "exact_rational_regression_checks"
                ]
            ),
            10,
        )
        self.assertEqual(
            certificate["n2_three_knot_theorem"][
                "symbolic_identity_and_constant_check"
            ],
            "passed",
        )
        self.assertEqual(
            certificate["n3_four_knot_theorem"][
                "symbolic_identity_and_constant_check"
            ],
            "passed",
        )
        self.assertEqual(
            certificate["three_positive_convex_core_all_n"][
                "symbolic_identity_check"
            ],
            "passed",
        )
        self.assertEqual(
            certificate["three_positive_far_cap_all_n"][
                "symbolic_identity_and_constant_check"
            ],
            "passed",
        )
        middle_knot = certificate["three_positive_middle_knot_all_n"]
        self.assertEqual(
            middle_knot["symbolic_identity_and_constant_check"],
            "passed",
        )
        self.assertIn("b<=n/3", middle_knot["scope"])
        self.assertEqual(
            middle_knot["exact_rational_bounds"][
                "refined_lambda_3_bernstein_minimum"
            ]["numerator"],
            "264755763361",
        )
        self.assertEqual(
            certificate["two_positive_knots_all_n"][
                "symbolic_identity_and_constant_check"
            ],
            "passed",
        )


if __name__ == "__main__":
    unittest.main()
