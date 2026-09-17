import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MODULE_PATH = ROOT / "supporting-materials/computations/python/terminal_edge_bifurcation.py"
SPEC = importlib.util.spec_from_file_location("terminal_edge_bifurcation", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalEdgeBifurcationTests(unittest.TestCase):
    def test_exact_special_cases_and_limit_bracket(self):
        self.assertEqual(MODULE.critical_polynomial(2, MODULE.Fraction(3, 5)), 0)
        certificate = MODULE.build_certificate()
        finite = certificate["finite_root_brackets_n2_through_n30"]
        self.assertEqual(len(finite), 29)
        self.assertEqual(finite[1]["n"], 3)
        self.assertIn("0.333685211867271", finite[1]["alpha_display_interval"])
        self.assertIn("0.284668137040838", certificate["limiting_root"]["display_interval"])

    def test_committed_certificate_is_reproducible(self):
        committed = MODULE.DEFAULT_OUTPUT.read_text()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "certificate.json"
            path.write_text(json.dumps(MODULE.build_certificate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(path.read_text(), committed)


if __name__ == "__main__":
    unittest.main()
