from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = (
    ROOT
    / "computations"
    / "certificates"
    / "systematic-disjoint-dp-certificate.json"
)
sys.path.insert(0, str(PYTHON_DIR))

from systematic_disjoint_dp import (  # noqa: E402
    DPComplexityLimitError,
    build_certificate,
    equal_disjoint_phase_formula,
    equal_multiple_start_all_zero_formula,
    exact_disjoint_phase_bound,
    exact_disjoint_systematic_frame_bound,
    exact_multiple_start_all_zero_bound,
    minimum_strict_tail_coalition,
    solve_payload,
)
from systematic_minimax import (  # noqa: E402
    exact_minimax_frame_bound,
    systematic_frame_atoms,
)


def _brute_minimum_mass(
    probabilities: tuple[Fraction, ...],
    required: tuple[int, ...],
    alpha: Fraction,
    starts: int,
) -> Fraction:
    required_set = set(required)
    best = Fraction(2)
    for mask in range(1 << len(probabilities)):
        coalition = {index for index in range(len(probabilities)) if mask >> index & 1}
        if not required_set <= coalition:
            continue
        mass = sum((probabilities[index] for index in coalition), Fraction(0))
        if mass**starts > alpha:
            best = min(best, mass)
    if best > 1:
        raise AssertionError("all phases should always be a feasible coalition")
    return best


class SystematicDisjointDPTests(unittest.TestCase):
    def test_dp_matches_exhaustive_coalitions(self) -> None:
        mass_vectors = (
            (1, 1, 1, 1),
            (1, 2, 3, 4),
            (2, 3, 5, 7),
            (1, 1, 2, 3, 5),
        )
        alphas = (Fraction(1, 20), Fraction(1, 4), Fraction(2, 5))
        for masses in mass_vectors:
            total = sum(masses)
            probabilities = tuple(Fraction(mass, total) for mass in masses)
            for starts in (1, 2, 3):
                for alpha in alphas:
                    for required in ((0,), (len(masses) - 1,), (0, 1)):
                        result = minimum_strict_tail_coalition(
                            probabilities,
                            alpha,
                            required,
                            independent_starts=starts,
                        )
                        expected = _brute_minimum_mass(
                            probabilities, required, alpha, starts
                        )
                        self.assertEqual(result.coalition_probability, expected)
                        self.assertGreater(
                            result.product_lower_tail_probability, alpha
                        )

    def test_equal_formula_matches_general_rational_lp_solver(self) -> None:
        for phase_count in range(2, 7):
            weights = (1,) * phase_count
            for alpha in (Fraction(1, 20), Fraction(1, 4), Fraction(2, 5)):
                for observed_mean in (Fraction(0), Fraction(1, 3), Fraction(1)):
                    closed = equal_disjoint_phase_formula(
                        phase_count, observed_mean, alpha
                    )
                    general = exact_minimax_frame_bound(
                        weights,
                        1,
                        0,
                        {0: observed_mean},
                        alpha,
                    )
                    self.assertEqual(closed.bound, general.bound)

    def test_unequal_formula_matches_general_rational_lp_solver(self) -> None:
        fixtures = (
            ((1, 2, 3, 4), 0, Fraction(1, 3), Fraction(7, 20)),
            ((2, 7, 5, 3), 2, Fraction(1, 2), Fraction(1, 5)),
            ((11, 2, 6, 8, 3), 4, Fraction(2, 5), Fraction(9, 20)),
        )
        for weights, observed, mean, alpha in fixtures:
            total = sum(weights)
            probabilities = tuple(Fraction(weight, total) for weight in weights)
            specialized = exact_disjoint_phase_bound(
                probabilities, observed, mean, alpha
            )
            general = exact_minimax_frame_bound(
                weights, 1, observed, {observed: mean}, alpha
            )
            self.assertEqual(specialized.bound, general.bound)
            self.assertEqual(
                specialized.coalition.coalition_probability,
                general.maximizing_coalition_probability,
            )

    def test_frame_wrapper_verifies_partition_and_matches_general_solver(self) -> None:
        weights = (1,) * 12
        sample_size = 3
        _, _, raw_items = systematic_frame_atoms(weights, sample_size)
        observed = {index: Fraction(1, 3) for index in set(raw_items[0])}
        specialized = exact_disjoint_systematic_frame_bound(
            weights,
            sample_size,
            0,
            observed,
            Fraction(1, 4),
        )
        general = exact_minimax_frame_bound(
            weights,
            sample_size,
            0,
            observed,
            Fraction(1, 4),
        )
        self.assertEqual(specialized.phase_atom_count, 4)
        self.assertEqual(specialized.observed_mean, Fraction(1, 3))
        self.assertEqual(specialized.bound, general.bound)

    def test_frame_wrapper_rejects_overlapping_atoms(self) -> None:
        with self.assertRaisesRegex(ValueError, "pairwise disjoint"):
            exact_disjoint_systematic_frame_bound(
                (3, 3, 4, 2),
                2,
                0,
                {0: 0, 2: 0},
                Fraction(1, 5),
            )

    def test_subset_sum_reduction_distinguishes_yes_and_no_instances(self) -> None:
        def reduced(instance: tuple[int, ...], target: int):
            weights = (1,) + tuple(2 * value for value in instance)
            total = sum(weights)
            result = exact_disjoint_phase_bound(
                tuple(Fraction(weight, total) for weight in weights),
                0,
                0,
                Fraction(2 * target, total),
            )
            decoded = (
                result.coalition.coalition_scaled_mass - 1
            ) // 2
            return result, decoded

        yes, yes_sum = reduced((3, 5, 6), 8)
        no, no_sum = reduced((3, 5, 9), 7)
        self.assertEqual(yes_sum, 8)
        self.assertEqual(yes.bound, Fraction(12, 29))
        self.assertEqual(no_sum, 8)
        self.assertNotEqual(no_sum, 7)
        self.assertEqual(no.bound, Fraction(18, 35))

    def test_multiple_start_equal_formula_matches_exact_dp(self) -> None:
        for phase_count in range(2, 9):
            probabilities = (Fraction(1, phase_count),) * phase_count
            for starts in (1, 2, 3):
                for distinct in range(1, min(phase_count, starts) + 1):
                    observed = tuple(range(distinct)) + (0,) * (starts - distinct)
                    for alpha in (
                        Fraction(1, 20),
                        Fraction(1, 4),
                        Fraction(2, 5),
                    ):
                        formula = equal_multiple_start_all_zero_formula(
                            phase_count, starts, distinct, alpha
                        )
                        dp = exact_multiple_start_all_zero_bound(
                            probabilities, observed, starts, alpha
                        )
                        self.assertEqual(formula.bound, dp.bound)
                        self.assertEqual(
                            formula.required_zero_phase_count,
                            dp.coalition.coalition_scaled_mass,
                        )

    def test_strict_threshold_handles_exact_discrete_boundary(self) -> None:
        result = minimum_strict_tail_coalition(
            (Fraction(1, 4),) * 4,
            Fraction(1, 4),
            (0,),
        )
        self.assertEqual(result.strict_scaled_threshold, 2)
        self.assertEqual(result.coalition_probability, Fraction(1, 2))

    def test_scaled_state_limit_fails_closed(self) -> None:
        with self.assertRaises(DPComplexityLimitError):
            minimum_strict_tail_coalition(
                (Fraction(1, 101), Fraction(100, 101)),
                Fraction(1, 20),
                (0,),
                max_scaled_total=100,
            )

    def test_generic_json_interfaces_preserve_exact_values(self) -> None:
        one_start = solve_payload(
            {
                "probabilities": ["1/10", "1/5", "3/10", "2/5"],
                "observed_phase_index": 0,
                "observed_mean": "1/3",
                "alpha": "7/20",
            }
        )
        self.assertEqual(one_start["bound"]["numerator"], "11")
        self.assertEqual(one_start["bound"]["denominator"], "15")

        multiple = solve_payload(
            {
                "probabilities": ["1/20"] * 20,
                "observed_phase_indices": [0, 1],
                "independent_starts": 2,
                "alpha": "1/20",
            }
        )
        self.assertEqual(multiple["bound"]["numerator"], "3")
        self.assertEqual(multiple["bound"]["denominator"], "4")

        frame = solve_payload(
            {
                "weights": ["1"] * 4,
                "sample_size": 2,
                "observed_phase_index": 0,
                "observed_taints": {"0": "0", "2": "0"},
                "alpha": "1/4",
            }
        )
        self.assertEqual(
            frame["mode"], "verified phase-partitioned systematic-PPS frame"
        )
        self.assertEqual(frame["bound"]["numerator"], "1")
        self.assertEqual(frame["bound"]["denominator"], "2")

    def test_certificate_matches_generator(self) -> None:
        self.assertEqual(json.loads(CERTIFICATE.read_text()), build_certificate())


if __name__ == "__main__":
    unittest.main()
