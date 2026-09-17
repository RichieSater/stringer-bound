"""Tests for the six-coordinate unique-top ``L=4`` certificate."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent
    / "certificates"
    / "six-coordinate-unique-top-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from six_coordinate_unique_top import _validate_leaf_partition  # noqa: E402
from six_coordinate_unique_top_witness import (  # noqa: E402
    B_GE_V_FINAL_LEAVES,
    LOWER_LOW_Z_LEAVES,
    LOWER_MID_CORNER_LEAVES,
    MIDDLE_LOW_Z_LEAVES,
    Z_LE_A_FINAL_LEAVES,
)


class SixCoordinateUniqueTopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(CERTIFICATE.read_text())

    def test_certificate_scope_is_exact(self):
        scope = self.certificate["scope"]
        self.assertEqual(scope["proved_unique_top_families"], ["L=4"])
        self.assertEqual(
            scope["companion_unique_top_families"], ["L=2", "L=3"]
        )
        self.assertIn("all nine", scope["global_consequence"])
        self.assertIn("unresolved", scope["global_consequence"])
        self.assertIn("active-prefix", scope["global_consequence"])
        self.assertIn("monotone cap theorem", scope["global_consequence"])

    def test_polynomial_structures(self):
        expected = {
            "target_T": ([21, 14, 7, 23, 30], 87_047),
            "q_le_u_target": ([14, 14, 7, 23, 30], 87_047),
            "z_le_a_final_target": ([36, 61, 14, 7, 14], 87_047),
            "b_ge_v_final_target": ([14, 7, 36, 23, 14], 87_047),
            "middle_low_target": ([14, 23, 36, 23, 7], 87_047),
            "lower_low_target": ([14, 23, 23, 36, 23], 87_047),
        }
        for name, (degree, terms) in expected.items():
            with self.subTest(polynomial=name):
                record = self.certificate["polynomials"][name]
                self.assertEqual(record["degree"], degree)
                self.assertEqual(record["power_terms"], terms)

        gaps = self.certificate["polynomials"]["centroid_gaps"]
        self.assertEqual(
            [(record["degree"], record["power_terms"]) for record in gaps],
            [
                ([4, 2, 1, 5, 6], 74),
                ([5, 3, 1, 5, 6], 74),
                ([5, 4, 2, 5, 6], 74),
                ([7, 5, 3, 7, 8], 218),
                ([9, 6, 3, 6, 10], 572),
            ],
        )

    def test_exact_partition_shapes_and_closures(self):
        expected = {
            "z_le_a": (
                Z_LE_A_FINAL_LEAVES,
                53,
                27,
                {"T": 8, "M4": 11, "G4": 8},
            ),
            "b_ge_v": (
                B_GE_V_FINAL_LEAVES,
                17,
                9,
                {"T": 5, "G4": 4},
            ),
            "middle_low": (
                MIDDLE_LOW_Z_LEAVES,
                19,
                10,
                {"T": 6, "G4": 4},
            ),
            "lower_low": (
                LOWER_LOW_Z_LEAVES,
                53,
                27,
                {"T": 13, "G4": 14},
            ),
            "lower_mid_corner": (
                LOWER_MID_CORNER_LEAVES,
                5,
                3,
                {"T": 2, "G4": 1},
            ),
        }
        covers = self.certificate["bernstein_covers"]
        for name, (leaves, nodes, count, closures) in expected.items():
            with self.subTest(chart=name):
                tree = _validate_leaf_partition(
                    leaves, {"T", "G4", "M4"}
                )
                self.assertEqual(tree["nodes"], nodes)
                self.assertEqual(tree["leaves"], count)
                self.assertEqual(covers[name]["closure_counts"], closures)
                self.assertEqual(
                    covers[name]["bernstein_signs"]["negative_coefficients"],
                    0,
                )

    def test_complete_cover_has_83_boxes(self):
        covers = self.certificate["bernstein_covers"]
        self.assertEqual(
            {name: record["box_count"] for name, record in covers.items()},
            {
                "q_le_u": 1,
                "z_le_a": 27,
                "b_ge_v": 9,
                "middle_low": 10,
                "middle_high": 2,
                "lower_low": 27,
                "lower_large": 4,
                "lower_mid_corner": 3,
            },
        )
        self.assertEqual(
            sum(record["box_count"] for record in covers.values()), 83
        )
        self.assertEqual(
            covers["q_le_u"]["bernstein_signs"]["positive_coefficients"],
            776_265,
        )
        self.assertEqual(
            covers["q_le_u"]["bernstein_signs"]["zero_coefficients"],
            562_935,
        )

    def test_scale_coefficient_counts(self):
        expected = {
            "global": (10_279, 11_281),
            "z_le_a_q": (9_380, 1_820),
            "z_le_a_r": (12_852, 252),
            "z_le_a_s": (14_868, 252),
            "z_le_a_c": (16_884, 252),
            "b_ge_v": (6_356, 1_484),
            "middle_low": (5_816, 2_024),
            "lower_low": (9_460, 4_260),
        }
        for name, (positive, zero) in expected.items():
            with self.subTest(scale_bound=name):
                signs = self.certificate["scale_bounds"][name][
                    "bernstein_signs"
                ]
                self.assertEqual(signs["positive_coefficients"], positive)
                self.assertEqual(signs["zero_coefficients"], zero)
                self.assertEqual(signs["negative_coefficients"], 0)

    def test_all_scale_and_box_signs_are_exactly_nonnegative(self):
        for name, record in self.certificate["scale_bounds"].items():
            with self.subTest(scale_bound=name):
                self.assertEqual(
                    record["bernstein_signs"]["negative_coefficients"], 0
                )
        for name, record in self.certificate["bernstein_covers"].items():
            with self.subTest(cover=name):
                self.assertEqual(
                    record["bernstein_signs"]["negative_coefficients"], 0
                )

    def test_multipliers_are_nonnegative_integers(self):
        records = self.certificate["bernstein_covers"]["z_le_a"][
            "exact_multipliers"
        ]
        self.assertEqual(len(records), 11)
        for record in records:
            self.assertEqual(record["multiplier"]["denominator"], "1")
            self.assertGreaterEqual(int(record["multiplier"]["numerator"]), 0)


if __name__ == "__main__":
    unittest.main()
