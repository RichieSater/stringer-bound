"""Build the manuscript counterexample table from certificate JSON files.

Only version-2 records marked ``confirmed_counterexample`` are included.
The minimum is selected by comparing the printed exact rational coverages,
not their floating-point convenience fields.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _decimal(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits + 10
        rendered = Decimal(value.numerator) / Decimal(value.denominator)
    return format(rendered, f".{digits}f")


def load_records(paths):
    for path_text in paths:
        path = Path(path_text)
        with path.open() as handle:
            payload = json.load(handle)
        records = payload if isinstance(payload, list) else [payload]
        for record in records:
            if record.get("certificate_version") != 2:
                raise ValueError(f"{path}: expected certificate_version 2")
            if record.get("confirmed_counterexample"):
                yield path, record


def summarize(paths):
    groups = defaultdict(list)
    for path, record in load_records(paths):
        candidate = record["input"]
        alpha = Fraction(str(candidate["alpha"]))
        key = (Fraction(1) - alpha, int(candidate["n"]))
        groups[key].append((path, record))

    rows = []
    for (confidence, n), members in sorted(groups.items()):
        minimum_path, minimum_record = min(
            members, key=lambda item: Fraction(item[1]["exact_coverage"])
        )
        minimum = Fraction(minimum_record["exact_coverage"])
        percent = confidence * 100
        rows.append({
            "nominal_confidence": _fraction_text(confidence),
            "nominal_percent": (
                str(percent.numerator) if percent.denominator == 1
                else _decimal(percent, 6).rstrip("0").rstrip(".")),
            "n": n,
            "smallest_exact_coverage": _fraction_text(minimum),
            "smallest_coverage_decimal": _decimal(minimum),
            "table_display": _decimal(minimum, 5),
            "certified_examples": len(members),
            "minimum_source": minimum_path.name,
            "sources": sorted({path.name for path, _ in members}),
        })
    return rows


def default_paths():
    certificate_dir = Path(__file__).resolve().parents[1] / "certificates"
    return sorted(glob.glob(str(certificate_dir / "certified-alpha*.json")))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("certificates", nargs="*",
                        help="version-2 certificate JSON files")
    parser.add_argument("--out", help="write the machine-readable summary")
    args = parser.parse_args(argv)

    paths = args.certificates or default_paths()
    rows = summarize(paths)
    payload = {
        "generated_by": "summarize_certificates.py",
        "certificate_version": 2,
        "rows": rows,
    }
    if args.out:
        with open(args.out, "w") as handle:
            json.dump(payload, handle, indent=1)

    for row in rows:
        print(f"${row['nominal_percent']}\\%$ & ${row['n']}$ & "
              f"${row['table_display']}$ & "
              f"${row['certified_examples']}$\\\\")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
