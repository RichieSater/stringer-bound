from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = (
    ROOT
    / "computations"
    / "certificates"
    / "systematic-minimax-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from systematic_minimax import (  # noqa: E402
    ComplexityLimitError,
    _minimal_probability_coalitions,
    build_certificate,
    exact_minimax_frame_bound,
    solve_payload,
    systematic_frame_atoms,
)


class SystematicMinimaxTests(unittest.TestCase):
    def test_minimal_coalitions_match_brute_force(self) -> None:
        probabilities = (Fraction(1, 2), Fraction(1, 3), Fraction(1, 6))
        for alpha in (Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)):
            actual, _ = _minimal_probability_coalitions(
                probabilities,
                alpha,
                max_coalitions=100,
                max_search_nodes=1000,
            )
            expected = []
            for size in range(1, len(probabilities) + 1):
                for coalition in combinations(range(len(probabilities)), size):
                    total = sum(
                        (probabilities[index] for index in coalition),
                        Fraction(0),
                    )
                    if total <= alpha:
                        continue
                    if all(
                        total - probabilities[index] <= alpha
                        for index in coalition
                    ):
                        expected.append(coalition)
            self.assertEqual(actual, tuple(sorted(expected)))

    def test_four_item_nonbinary_result_is_exact_and_dominates_hybrid(self) -> None:
        weights = (1, 1, 1, 1)
        atoms, _, raw_items = systematic_frame_atoms(weights, 1)
        self.assertEqual(len(atoms), 4)
        observed = {raw_items[0][0]: Fraction(1, 2)}
        result = exact_minimax_frame_bound(
            weights, 1, 0, observed, Fraction(1, 4)
        )
        self.assertEqual(result.bound, Fraction(3, 4))
        self.assertEqual(result.markov_bound, Fraction(7, 8))
        self.assertEqual(result.completion_bound, Fraction(7, 8))
        self.assertEqual(result.witness_lower_tail_probability, Fraction(1, 2))
        self.assertLess(result.bound, result.prior_hybrid_bound)

    def test_irregular_frame_accepts_integral_exact_simplex_coordinates(self) -> None:
        weights = (3, 3, 4, 2)
        result = exact_minimax_frame_bound(
            weights,
            2,
            0,
            {0: Fraction(0), 2: Fraction(1, 2)},
            Fraction(1, 5),
        )
        self.assertEqual(result.bound, Fraction(7, 12))
        self.assertEqual(result.witness_lower_tail_probability, Fraction(1, 2))

    def test_tiny_frame_has_exact_design_coverage_on_a_rational_grid(self) -> None:
        weights = (1, 1, 1, 1)
        atoms, raw_to_atom, raw_items = systematic_frame_atoms(weights, 1)
        self.assertEqual(tuple(raw_to_atom), (0, 1, 2, 3))
        alpha = Fraction(1, 4)
        values = (Fraction(0), Fraction(1, 2), Fraction(1))
        cache = {}
        for phase_index, selected in enumerate(raw_items):
            item_index = selected[0]
            for value in values:
                cache[(phase_index, value)] = exact_minimax_frame_bound(
                    weights,
                    1,
                    phase_index,
                    {item_index: value},
                    alpha,
                ).bound

        for t0 in values:
            for t1 in values:
                for t2 in values:
                    for t3 in values:
                        taints = (t0, t1, t2, t3)
                        target = sum(taints, Fraction(0)) / 4
                        failure = sum(
                            (
                                atom.probability
                                for phase_index, atom in enumerate(atoms)
                                if cache[(phase_index, taints[phase_index])]
                                < target
                            ),
                            Fraction(0),
                        )
                        self.assertLessEqual(failure, alpha)

    def test_practical_frame_aggregates_to_twenty_exact_atoms(self) -> None:
        sample_size = 100
        weights = (1,) * 2000
        atoms, _, raw_items = systematic_frame_atoms(weights, sample_size)
        self.assertEqual(len(atoms), 20)
        self.assertTrue(all(atom.probability == Fraction(1, 20) for atom in atoms))
        observed = {index: 0 for index in set(raw_items[0])}
        result = exact_minimax_frame_bound(
            weights, sample_size, 0, observed, Fraction(1, 20)
        )
        self.assertEqual(result.bound, Fraction(9, 10))
        self.assertEqual(result.prior_hybrid_bound, Fraction(19, 20))
        self.assertEqual(result.minimal_coalition_count, 190)
        self.assertEqual(result.maximizing_coalition, (0, 1))
        self.assertEqual(result.maximizing_coalition_probability, Fraction(1, 10))
        self.assertEqual(result.witness_lower_tail_probability, Fraction(1, 10))

    def test_search_limits_fail_closed(self) -> None:
        weights = (1, 1, 1, 1)
        _, _, raw_items = systematic_frame_atoms(weights, 1)
        observed = {raw_items[0][0]: 0}
        with self.assertRaises(ComplexityLimitError):
            exact_minimax_frame_bound(
                weights,
                1,
                0,
                observed,
                Fraction(1, 4),
                max_coalitions=1,
            )

    def test_generic_json_interface_preserves_exact_rationals(self) -> None:
        payload = solve_payload(
            {
                "weights": ["1", "1", "1", "1"],
                "sample_size": 1,
                "observed_phase_index": 0,
                "observed_taints": {"0": "1/2"},
                "alpha": "1/4",
            }
        )
        result = payload["result"]
        self.assertEqual(result["exact_minimax_bound"]["numerator"], "3")
        self.assertEqual(result["exact_minimax_bound"]["denominator"], "4")
        self.assertEqual(result["witness_taints"], ["1/2", "1/2", "1/1", "1/1"])

    def test_certificate_matches_generator(self) -> None:
        self.assertEqual(json.loads(CERTIFICATE.read_text()), build_certificate())


if __name__ == "__main__":
    unittest.main()
