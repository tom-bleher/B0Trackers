#!/usr/bin/env python3
"""Compare two B0 report JSON files and fail on configured regressions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from .b0_metrics import evaluate_regression_policy, nested_get
except ImportError:
    from b0_metrics import evaluate_regression_policy, nested_get


DEFAULT_CHECKS = [
    {
        "path": "stages.filtered_over_eligible.value",
        "mode": "max_drop",
        "tolerance": 0.01,
    },
    {
        "path": "stages.truth_matched_over_eligible.value",
        "mode": "max_drop",
        "tolerance": 0.01,
    },
    {
        "path": "truth_matched.relative_delta_p.std",
        "mode": "max_fractional_increase",
        "tolerance": 0.10,
    },
    {
        "path": "sensor_mapping.fallback_fraction_of_resolved",
        "mode": "max_increase",
        "tolerance": 0.001,
    },
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument(
        "--policy",
        type=Path,
        help="JSON file containing a `checks` array; defaults are conservative examples",
    )
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text())
    candidate = json.loads(args.candidate.read_text())
    checks = DEFAULT_CHECKS
    if args.policy:
        checks = json.loads(args.policy.read_text())["checks"]

    print(f"{'metric':56} {'baseline':>12} {'candidate':>12}")
    print("-" * 84)
    for check in checks:
        path = check["path"]
        try:
            base = nested_get(baseline, path)
            cand = nested_get(candidate, path)
        except KeyError:
            print(f"{path:56} {'MISSING':>12} {'MISSING':>12}")
            continue
        print(f"{path:56} {str(base):>12} {str(cand):>12}")

    try:
        failures = evaluate_regression_policy(baseline, candidate, checks)
    except KeyError as exc:
        print(f"policy references missing metric: {exc}")
        return 2

    if failures:
        print("\nREGRESSION CHECK FAILED")
        for failure in failures:
            print(json.dumps(failure, sort_keys=True))
        return 1

    print("\nAll configured B0 regression checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
