"""Tests for the six-coordinate repeated-top local theorem."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import sympy as sp


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
CERTIFICATE = (
    PYTHON_DIR.parent
    / "certificates"
    / "six-coordinate-repeated-top-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from six_coordinate_repeated_top import (  # noqa: E402
    _derive_mixed_block_case,
    derive_one_positive_block,
)
from six_coordinate_face32 import _validate_leaf_partition  # noqa: E402
from six_coordinate_face32_witness import (  # noqa: E402
    A_LE_Z_LEAVES,
    T_LE_U_LEAVES,
    Z_LE_A_LEAVES,
)


class SixCoordinateRepeatedTopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(CERTIFICATE.read_text())

    def test_certificate_has_exactly_the_proved_and_open_scope(self):
        scope = self.certificate["scope"]
        self.assertEqual(
            scope["proved_faces"], ["42", "33", "24", "23", "22", "32"]
        )
        self.assertEqual(scope["open_repeated_top_faces"], [])
        self.assertIn("all nine", scope["global_limitation"])
        self.assertIn("unresolved", scope["global_limitation"])
        self.assertIn("active-prefix", scope["global_limitation"])
        self.assertIn("companion", scope["global_limitation"])

    def test_face32_exact_partition_scope(self):
        expected = {
            "t_le_u": (T_LE_U_LEAVES, 17, 9, {"T": 7, "G3": 2}),
            "A_le_Z": (
                A_LE_Z_LEAVES,
                127,
                64,
                {"T": 35, "M3": 7, "G3": 18, "G2": 4},
            ),
            "Z_le_A": (
                Z_LE_A_LEAVES,
                149,
                75,
                {"M3": 14, "G3": 19, "T": 22, "M2": 7, "G2": 13},
            ),
        }
        recorded = self.certificate["final_repeated_top_face"]["dyadic_charts"]
        for name, (leaves, nodes, count, closures) in expected.items():
            with self.subTest(chart=name):
                tree = _validate_leaf_partition(leaves)
                self.assertEqual(tree["nodes"], nodes)
                self.assertEqual(tree["leaves"], count)
                self.assertEqual(
                    recorded[name]["subdivision"]["closure_counts"], closures
                )
                self.assertEqual(
                    recorded[name]["bernstein_signs"]["negative_coefficients"],
                    0,
                )

        confinement = self.certificate["final_repeated_top_face"]
        self.assertEqual(
            confinement["scale_confinement_signs"]["negative_coefficients"],
            0,
        )
        self.assertEqual(
            confinement["quadratic_confinement_signs"]["negative_coefficients"],
            0,
        )

    def test_one_positive_block_derivations_and_witness_counts(self):
        expected = {
            2: ((2, 2, 2, 2), (9, 9, 9, 9), 415, 9505, 495),
            3: ((3, 3, 3), (14, 14, 14), 549, 2588, 787),
            4: ((4, 4), (19, 19), 126, 154, 246),
        }
        for multiplicity, values in expected.items():
            with self.subTest(multiplicity=multiplicity):
                derived = derive_one_positive_block(multiplicity)
                order_degree, target_degree, multiplier_positive, residual_positive, residual_zero = values
                self.assertEqual(tuple(derived["order"].degree_list()), order_degree)
                self.assertEqual(tuple(derived["target"].degree_list()), target_degree)
                recorded = self.certificate["one_positive_top_block"][
                    str(multiplicity)
                ]["bernstein_basis"]
                self.assertEqual(
                    recorded["multiplier_S"]["positive_coefficients"],
                    multiplier_positive,
                )
                self.assertEqual(
                    recorded["multiplier_S"]["negative_coefficients"], 0
                )
                self.assertEqual(
                    recorded["residual_H"]["positive_coefficients"],
                    residual_positive,
                )
                self.assertEqual(
                    recorded["residual_H"]["zero_coefficients"], residual_zero
                )
                self.assertEqual(recorded["residual_H"]["negative_coefficients"], 0)

    def test_exact_rational_one_positive_block_fixtures(self):
        points = {
            2: (sp.Rational(1, 2),) * 4,
            3: (sp.Rational(1, 2),) * 3,
            4: (sp.Rational(1, 10),) * 2,
        }
        for multiplicity, values in points.items():
            with self.subTest(multiplicity=multiplicity):
                derived = derive_one_positive_block(multiplicity)
                cube_point = dict(zip(derived["cube_variables"], values, strict=True))
                r_values = tuple(
                    sp.prod(derived["cube_variables"][index:]).subs(cube_point)
                    for index in range(derived["lower_count"])
                )
                r_point = dict(zip(derived["r_variables"], r_values, strict=True))
                cap = derived["cap"].as_expr().subs(r_point)
                normalizer = derived["normalizer"].as_expr().subs(r_point)
                gamma0 = (
                    derived["gamma0_numerator"].as_expr().subs(r_point)
                    / normalizer
                )
                self.assertGreater(derived["order"].as_expr().subs(cube_point), 0)
                self.assertLess(cap, 1 - (1 - gamma0) ** 5)

    def test_23_interpolation_derivation_and_certificate_signs(self):
        # The 22 face is intentionally left to the full make target because
        # deriving and homogenizing its 3.4 million exact coefficients is the
        # expensive part of the standalone certificate replay.
        derived = _derive_mixed_block_case("23")
        self.assertEqual(tuple(derived["order"].degree_list()), (5, 2, 6))
        self.assertEqual(tuple(derived["target"].degree_list()), (19, 11, 28))
        recorded = self.certificate["mixed_positive_blocks"]["23"]
        self.assertEqual(recorded["compactness_sign_audit"]["negative_coefficients"], 0)
        self.assertEqual(
            recorded["bernstein_homogenization"]["quotient_S"][
                "negative_coefficients"
            ],
            0,
        )
        self.assertEqual(
            recorded["bernstein_homogenization"]["remainder_H"][
                "negative_coefficients"
            ],
            0,
        )


if __name__ == "__main__":
    unittest.main()
