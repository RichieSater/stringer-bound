"""Tests for the unique-top ``L=2,3`` certificate."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import sympy as sp


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent
    / "certificates"
    / "six-coordinate-unique-top-l2-l3-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from certified_binary64_bernstein import (  # noqa: E402
    _witness_trie,
    verify_binary64_bernstein_tree,
)
from six_coordinate_unique_top import _validate_leaf_partition  # noqa: E402
from six_coordinate_unique_top_l2_l3_witness import (  # noqa: E402
    L2_C_LE_Z,
    L2_Q_LE_U,
    L2_Z_LE_C_LOWER,
    L2_Z_LE_C_MIDDLE,
    L2_Z_LE_C_UPPER,
    L3_B_LE_Z_LOWER,
    L3_B_LE_Z_MIDDLE,
    L3_B_LE_Z_UPPER,
    L3_Q_LE_U,
    L3_Z_LE_B_LOWER,
    L3_Z_LE_B_MIDDLE,
    L3_Z_LE_B_UPPER,
)


class SixCoordinateUniqueTopL2L3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(CERTIFICATE.read_text())

    def test_scope_closes_the_two_remaining_rows(self):
        scope = self.certificate["scope"]
        self.assertEqual(scope["proved_unique_top_families"], ["L=2", "L=3"])
        self.assertEqual(scope["companion_exact_family"], "L=4")
        self.assertIn("all nine", scope["global_consequence"])
        self.assertIn("unresolved", scope["global_consequence"])
        self.assertIn("active-prefix", scope["global_consequence"])

    def test_sharp_rational_minorant(self):
        comparison = self.certificate["comparison"]
        self.assertEqual(
            comparison["stronger_conclusion"],
            "p<=5*gamma_0/(1+(10974/4651)*gamma_0)",
        )
        minorant = comparison["minorant"]
        self.assertEqual(minorant["interval"], ["0", "1/6"])
        self.assertEqual(
            minorant["sharp_constant"],
            {"numerator": "10974", "denominator": "4651"},
        )
        self.assertEqual(minorant["endpoint_equality"], "g=1/6")

    def test_base_polynomial_structures(self):
        families = self.certificate["families"]
        self.assertEqual(families["L=2"]["target"]["degree"], [7, 4, 8, 12, 14])
        self.assertEqual(families["L=2"]["target"]["power_terms"], 2_779)
        self.assertEqual(families["L=3"]["target"]["degree"], [10, 5, 6, 12, 15])
        self.assertEqual(families["L=3"]["target"]["power_terms"], 4_687)
        self.assertEqual(
            [
                (record["degree"], record["power_terms"])
                for record in families["L=2"]["centroid_gaps"]
            ],
            [
                ([6, 3, 6, 9, 10], 572),
                ([7, 3, 5, 7, 8], 218),
                ([5, 2, 4, 5, 6], 74),
                ([5, 1, 3, 5, 6], 74),
                ([5, 1, 2, 4, 6], 74),
            ],
        )

    def test_all_exact_scale_and_leaf_signs_are_nonnegative(self):
        for family in self.certificate["families"].values():
            for name, record in family["scale_bounds"].items():
                with self.subTest(scale=name):
                    self.assertEqual(
                        record["bernstein_signs"]["negative_coefficients"], 0
                    )
            for name, record in family["exact_rational_covers"].items():
                with self.subTest(cover=name):
                    self.assertEqual(
                        record["bernstein_signs"]["negative_coefficients"], 0
                    )

    def test_cover_sizes_and_closures(self):
        l2 = self.certificate["families"]["L=2"]
        self.assertEqual(
            {
                name: (record["box_count"], record["closure_counts"])
                for name, record in l2["exact_rational_covers"].items()
            },
            {
                "q_le_u": (37, {"G4": 11, "T": 26}),
                "c_le_z": (78, {"G4": 26, "T": 52}),
                "z_le_c_middle": (54, {"G4": 21, "T": 33}),
                "z_le_c_upper": (30, {"G4": 12, "T": 18}),
            },
        )
        enclosed = l2["certified_binary64_cover"]
        self.assertEqual(enclosed["subdivision"]["leaves"], 11)
        self.assertEqual(
            enclosed["closure_counts"], {"G4": 1, "M4": 8, "T": 2}
        )
        margin = Fraction(
            int(enclosed["minimum_certified_margin"]["numerator"]),
            int(enclosed["minimum_certified_margin"]["denominator"]),
        )
        self.assertGreater(margin, 0)

        l3 = self.certificate["families"]["L=3"]["exact_rational_covers"]
        self.assertEqual(sum(record["box_count"] for record in l3.values()), 56)
        self.assertEqual(
            {name: record["box_count"] for name, record in l3.items()},
            {
                "q_le_u": 8,
                "z_le_b_upper": 27,
                "z_le_b_middle": 8,
                "z_le_b_lower": 2,
                "b_le_z_upper": 7,
                "b_le_z_middle": 2,
                "b_le_z_lower": 2,
            },
        )

    def test_witness_partitions_are_prefix_free_and_complete(self):
        witnesses = (
            L2_Q_LE_U,
            L2_C_LE_Z,
            L2_Z_LE_C_MIDDLE,
            L2_Z_LE_C_UPPER,
            L3_Q_LE_U,
            L3_Z_LE_B_UPPER,
            L3_Z_LE_B_MIDDLE,
            L3_Z_LE_B_LOWER,
            L3_B_LE_Z_UPPER,
            L3_B_LE_Z_MIDDLE,
            L3_B_LE_Z_LOWER,
        )
        for witness in witnesses:
            with self.subTest(leaves=len(witness)):
                record = _validate_leaf_partition(
                    witness, {"T", "G4", "M4"}
                )
                self.assertEqual(record["leaves"], len(witness))
        _, record = _witness_trie(L2_Z_LE_C_LOWER)
        self.assertEqual(record, {"nodes": 21, "leaves": 11, "maximum_depth": 7})

    def test_binary64_enclosure_rejects_mutated_multiplier(self):
        x = sp.symbols("x")
        target = sp.Poly(2, x, domain=sp.QQ)
        gap = sp.Poly(1, x, domain=sp.QQ)
        good = verify_binary64_bernstein_tree(
            target, gap, (("", "M4", (1, 1)),), (0,)
        )
        self.assertEqual(good["subdivision"]["leaves"], 1)
        with self.assertRaises(AssertionError):
            verify_binary64_bernstein_tree(
                target, gap, (("", "M4", (3, 1)),), (0,)
            )

    def test_incomplete_binary64_tree_is_rejected(self):
        with self.assertRaises(AssertionError):
            _witness_trie((("0", "T", None),))

    def test_arithmetic_boundary_is_explicit(self):
        arithmetic = self.certificate["arithmetic"]
        self.assertIn("exact rational", arithmetic["symbolic_derivation"])
        self.assertIn("FLINT exact rational", arithmetic["exact_bernstein_signs"])
        self.assertIn("forward-error", arithmetic["large_l2_lower_chart"])


if __name__ == "__main__":
    unittest.main()
