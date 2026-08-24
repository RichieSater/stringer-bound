"""Prove that one hard primary ``n=7`` root exceeds the fixed call cutoff."""

from __future__ import annotations

import argparse
import json
import os
import pathlib

from n7_dirichlet_poissonization import (
    FACE7_EXTRA_REFINEMENT_CALL_CUTOFF,
    FACE7_WORKER_CALL_LIMIT,
    FaceCallLimitExceeded,
    _face7_primary_tasks,
    face7_worker,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=int, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    if FACE7_WORKER_CALL_LIMIT != FACE7_EXTRA_REFINEMENT_CALL_CUTOFF:
        raise RuntimeError(
            "set N7_FACE7_CALL_LIMIT to the fixed extra-refinement cutoff"
        )
    tasks = {task[0]: task for task in _face7_primary_tasks()}
    if args.root not in tasks:
        raise ValueError("unknown primary root")

    try:
        record = face7_worker(tasks[args.root])
    except FaceCallLimitExceeded as exc:
        if exc.root_index != args.root or exc.calls != FACE7_WORKER_CALL_LIMIT + 1:
            raise
        result = {
            "primary_root": args.root,
            "fixed_call_cutoff": FACE7_EXTRA_REFINEMENT_CALL_CUTOFF,
            "observed_calls_before_termination": exc.calls,
            "strictly_exceeds_cutoff": True,
            "source_commit": os.environ.get("GITHUB_SHA", "local"),
        }
    else:
        raise AssertionError(
            f"root {args.root} terminated in {record['calls']} calls; "
            "it does not satisfy the fixed strict-cutoff policy"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
