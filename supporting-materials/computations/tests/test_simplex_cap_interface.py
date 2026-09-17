"""Tests for the solution-free simplex-cap challenge interface."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
ROOT = PYTHON_DIR.parents[2]
INDEX_PATH = (
    PYTHON_DIR.parent / "certificates" / "simplex-cap-solution-index.json")
sys.path.insert(0, str(PYTHON_DIR))

from simplex_cap_interface import build_solution_index  # noqa: E402


class SimplexCapInterfaceTests(unittest.TestCase):
    def test_solution_index_regenerates_byte_for_byte(self):
        generated = build_solution_index()
        committed = json.loads(INDEX_PATH.read_text())
        self.assertEqual(generated, committed)
        self.assertEqual(len(generated["instances"]), 6)
        unresolved = {"binomial-n4-confidence-range", "binomial-n5-confidence-range"}
        for item in generated["instances"]:
            expected = "unresolved" if item["id"] in unresolved else "proved"
            self.assertEqual(item["mathematical_status"], expected)
            if expected == "unresolved":
                self.assertIn("active-prefix", item["proof_obligation"])
                self.assertIn("repeated-lowest", item["proof_obligation"])

    def test_every_bound_source_exists_and_matches_digest(self):
        index = build_solution_index()
        import hashlib
        for instance in index["instances"]:
            self.assertTrue(instance["verification_commands"])
            for source in instance["sources"]:
                path = ROOT / source["path"]
                self.assertTrue(path.is_file())
                self.assertEqual(
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    source["sha256"],
                )


if __name__ == "__main__":
    unittest.main()
