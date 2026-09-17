#!/usr/bin/env python3
"""Run the canonical fail-closed public-text check from the CI entry point."""
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[2]
runpy.run_path(
    str(ROOT / "supporting-materials/computations/python/check_public_corpus.py"),
    run_name="__main__",
)
