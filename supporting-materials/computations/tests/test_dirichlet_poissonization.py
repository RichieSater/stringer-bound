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
    verify_four_positive_far_cap,
    verify_n2_three_knot_theorem,
    verify_n3_four_knot_theorem,
    verify_n4_five_knot_theorem,
    verify_n5_four_positive_face,
    verify_radial_symbolic_identities,
    verify_saffine_symbolic_identities,
    verify_sparse_convex_core,
    verify_three_positive_convex_core,
    verify_three_positive_far_cap,
    verify_three_positive_full_face,
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

    def test_complete_n4_five_knot_certificate(self):
        record = verify_n4_five_knot_theorem()
        self.assertEqual(record["arb_precision_bits"], 160)
        self.assertEqual(
            record["main_four_positive_face"]["total_branch_calls"],
            77401,
        )
        self.assertEqual(
            record["secondary_constrained_convexity"][
                "lipschitz_terminal_boxes"
            ],
            12187,
        )
        self.assertEqual(
            record["secondary_constrained_convexity"][
                "terminal_transcript_sha256"
            ],
            "baf76e5da205718ac2f7e7037bde03a4d02e5cacfa9351c14254c24f1ca31dfe",
        )
        self.assertEqual(
            record["main_four_positive_face"][
                "terminal_transcript_sha256"
            ],
            "66dfba78573895b27c4ba6dd7616b68ed96a377f0c767107d30a5f08de8fbf94",
        )
        self.assertEqual(
            record["scalar_pruning_monotonicity"][
                "f_4_4_minus_f_4_3_lower_margin"
            ]["numerator"],
            "4159877",
        )

    def test_complete_n5_four_positive_face_certificate(self):
        record = verify_n5_four_positive_face()
        self.assertEqual(record["arb_precision_bits"], 160)
        self.assertEqual(
            record["four_positive_face"]["total_branch_calls"],
            81703,
        )
        self.assertEqual(
            record["four_positive_face"]["terminal_transcript_sha256"],
            "f4bae2dcb30244286c97891596caffa5c8b1580f41aa80c600fba63c5af7d5a3",
        )
        self.assertEqual(
            record["scalar_terminal_transcript_sha256"],
            "b5cc3d8cbda175722021249a795d275b0f29adbb6b76c658dada7a3903b2db2e",
        )

    def test_all_n_three_positive_convex_core_algebra(self):
        verify_three_positive_convex_core()

    def test_all_n_sparse_convex_core_algebra(self):
        record = verify_sparse_convex_core()
        self.assertEqual(
            record["symbolic_n_derivative_orders_checked"],
            [1, 2, 3, 4, 5, 6],
        )
        self.assertEqual(record["integer_parameter_pair_count"], 45)

    def test_all_n_three_positive_far_cap_algebra(self):
        verify_three_positive_far_cap()

    def test_all_n_four_positive_far_cap_algebra(self):
        record = verify_four_positive_far_cap()
        self.assertEqual(
            record["n4_left_endpoint_margin"]["numerator"],
            "29899236229",
        )
        self.assertEqual(
            record["n4_right_endpoint_margin"]["numerator"],
            "6461863",
        )

    def test_all_n_three_positive_middle_knot_algebra(self):
        verify_three_positive_middle_knot_region()

    def test_all_n_complete_three_positive_face(self):
        record = verify_three_positive_full_face()
        self.assertEqual(record["arb_precision_bits"], 128)
        self.assertEqual(
            len(record["finite_fprime_upper_bound"]),
            18,
        )

    def test_all_n_two_positive_knot_algebra(self):
        verify_two_positive_knot_theorem()

    def test_certificate_is_explicitly_non_theorem_research_support(self):
        certificate = build_certificate()
        self.assertEqual(certificate["schema_version"], 15)
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
            certificate["n4_five_knot_theorem"][
                "symbolic_exact_and_interval_checks"
            ],
            "passed",
        )
        self.assertEqual(
            certificate["n5_four_positive_face"][
                "symbolic_exact_and_interval_checks"
            ],
            "passed",
        )
        self.assertEqual(
            certificate["three_positive_convex_core_all_n"][
                "symbolic_identity_check"
            ],
            "passed",
        )
        sparse_core = certificate["sparse_convex_core_all_n"]
        self.assertEqual(
            sparse_core["symbolic_identity_and_threshold_check"],
            "passed",
        )
        self.assertIn("four-positive", sparse_core["scope"])
        self.assertEqual(
            certificate["three_positive_far_cap_all_n"][
                "symbolic_identity_and_constant_check"
            ],
            "passed",
        )
        self.assertEqual(
            certificate["four_positive_far_cap_all_n"][
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
        full_face = certificate["three_positive_full_face_all_n"]
        self.assertEqual(
            full_face["symbolic_identity_and_constant_check"],
            "passed",
        )
        self.assertIn(
            "Every n>=4",
            full_face["scope"],
        )
        self.assertEqual(
            certificate["two_positive_knots_all_n"][
                "symbolic_identity_and_constant_check"
            ],
            "passed",
        )


if __name__ == "__main__":
    unittest.main()
