#!/usr/bin/env python3
"""Compare two B0 report JSON files and fail on configured regressions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from .b0_metrics import evaluate_regression_policy, nested_get, provenance_mismatches
except ImportError:
    from b0_metrics import evaluate_regression_policy, nested_get, provenance_mismatches


DEFAULT_CHECKS = [
    {
        "path": "event_stage_presence.truth_matched_over_eligible.value",
        "mode": "max_drop",
        "tolerance": 0.01,
    },
    {
        "path": "truth_matched.relative_delta_p.std",
        "mode": "max_fractional_increase",
        "tolerance": 0.10,
    },
    {
        "path": "sensor_mapping.stub_ckf.fallback_fraction_of_resolved",
        "mode": "max_increase",
        "tolerance": 0.001,
    },
]

DEFAULT_COMPATIBILITY_PATHS = [
    "provenance.tree",
    "provenance.schema_versions",
    "provenance.geometry_names",
    "selection.min_stations",
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
    parser.add_argument(
        "--allow-incompatible",
        action="store_true",
        help="compare even if schema/geometry/selection provenance differs",
    )
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text())
    candidate = json.loads(args.candidate.read_text())
    checks = DEFAULT_CHECKS
    compatibility_paths = list(DEFAULT_COMPATIBILITY_PATHS)
    if args.policy:
        policy = json.loads(args.policy.read_text())
        checks = policy["checks"]
        compatibility_paths.extend(policy.get("compatibility_paths", []))

    compatibility_paths = list(dict.fromkeys(compatibility_paths))
    mismatches = provenance_mismatches(baseline, candidate, compatibility_paths)
    if mismatches and not args.allow_incompatible:
        print("INCOMPATIBLE B0 REPORTS")
        print("Regression thresholds are not meaningful until these provenance fields agree:")
        for mismatch in mismatches:
            print(json.dumps(mismatch, sort_keys=True))
        print("Use --allow-incompatible only for an intentional cross-configuration comparison.")
        return 3

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
