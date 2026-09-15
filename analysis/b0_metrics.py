"""Pure metric helpers for B0 tracking validation.

This module deliberately has no ROOT/uproot dependency so the numerical
regression logic can be unit-tested in a minimal Python environment.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from statistics import fmean, median


@dataclass(frozen=True)
class Efficiency:
    numerator: int
    denominator: int
    value: float
    low: float
    high: float

    def as_dict(self) -> dict[str, float | int]:
        return {
            "numerator": self.numerator,
            "denominator": self.denominator,
            "value": self.value,
            "low": self.low,
            "high": self.high,
        }


def wilson_efficiency(numerator: int, denominator: int, z: float = 1.0) -> Efficiency:
    """Return a binomial efficiency with a Wilson interval."""
    if denominator < 0 or numerator < 0 or numerator > denominator:
        raise ValueError("require 0 <= numerator <= denominator")
    if denominator == 0:
        nan = math.nan
        return Efficiency(numerator, denominator, nan, nan, nan)

    p = numerator / denominator
    z2 = z * z
    denom = 1.0 + z2 / denominator
    center = (p + z2 / (2.0 * denominator)) / denom
    half = (
        z
        * math.sqrt(p * (1.0 - p) / denominator + z2 / (4.0 * denominator**2))
        / denom
    )
    return Efficiency(numerator, denominator, p, max(0.0, center - half), min(1.0, center + half))


def finite(values: Iterable[float]) -> list[float]:
    return [float(v) for v in values if math.isfinite(float(v))]


def quantile(values: Sequence[float], q: float) -> float:
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must be in [0, 1]")
    vals = sorted(finite(values))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]
    pos = q * (len(vals) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac


def distribution_summary(values: Iterable[float]) -> dict[str, float | int]:
    vals = finite(values)
    if not vals:
        return {
            "n": 0,
            "mean": math.nan,
            "std": math.nan,
            "median": math.nan,
            "q16": math.nan,
            "q84": math.nan,
            "rms": math.nan,
        }

    mean = fmean(vals)
    variance = fmean((x - mean) ** 2 for x in vals)
    return {
        "n": len(vals),
        "mean": mean,
        "std": math.sqrt(variance),
        "median": median(vals),
        "q16": quantile(vals, 0.16),
        "q84": quantile(vals, 0.84),
        "rms": math.sqrt(fmean(x * x for x in vals)),
    }


def pull_summary(values: Iterable[float]) -> dict[str, float | int]:
    vals = finite(values)
    out = distribution_summary(vals)
    if vals:
        out["frac_abs_lt_1"] = sum(abs(x) < 1.0 for x in vals) / len(vals)
        out["frac_abs_lt_2"] = sum(abs(x) < 2.0 for x in vals) / len(vals)
    else:
        out["frac_abs_lt_1"] = math.nan
        out["frac_abs_lt_2"] = math.nan
    return out


def stage_efficiencies(
    eligible: Sequence[bool],
    seeded: Sequence[bool],
    unfiltered: Sequence[bool],
    filtered: Sequence[bool],
    truth_matched: Sequence[bool],
) -> dict[str, dict[str, float | int]]:
    """Summarize event-level B0 reconstruction flow for a selected truth particle."""
    arrays = [eligible, seeded, unfiltered, filtered, truth_matched]
    n = len(eligible)
    if any(len(a) != n for a in arrays):
        raise ValueError("stage arrays must have identical length")

    elig = [bool(x) for x in eligible]
    stages = {
        "seeded": [bool(s) and e for s, e in zip(seeded, elig)],
        "unfiltered": [bool(s) and e for s, e in zip(unfiltered, elig)],
        "filtered": [bool(s) and e for s, e in zip(filtered, elig)],
        "truth_matched": [bool(s) and e for s, e in zip(truth_matched, elig)],
    }
    denom = sum(elig)

    result: dict[str, dict[str, float | int]] = {
        "eligible": {"events": denom, "total_events": n}
    }
    for name, mask in stages.items():
        result[f"{name}_over_eligible"] = wilson_efficiency(sum(mask), denom).as_dict()

    sequential = [
        ("unfiltered_over_seeded", stages["unfiltered"], stages["seeded"]),
        ("filtered_over_unfiltered", stages["filtered"], stages["unfiltered"]),
        ("truth_matched_over_filtered", stages["truth_matched"], stages["filtered"]),
    ]
    for name, num_mask, den_mask in sequential:
        den = sum(den_mask)
        num = sum(nv and dv for nv, dv in zip(num_mask, den_mask))
        result[name] = wilson_efficiency(num, den).as_dict()
    return result


def binned_efficiency(
    values: Sequence[float],
    passed: Sequence[bool],
    eligible: Sequence[bool],
    edges: Sequence[float],
):
    """Return Wilson efficiencies in explicit half-open bins [lo, hi), last bin inclusive."""
    if not (len(values) == len(passed) == len(eligible)):
        raise ValueError("values/passed/eligible must have identical length")
    if len(edges) < 2 or any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("edges must be strictly increasing")

    rows = []
    last = len(edges) - 2
    for i, (lo, hi) in enumerate(zip(edges, edges[1:])):
        den = num = 0
        for value, ok, use in zip(values, passed, eligible):
            if not use or not math.isfinite(float(value)):
                continue
            inside = lo <= value <= hi if i == last else lo <= value < hi
            if inside:
                den += 1
                num += bool(ok)
        eff = wilson_efficiency(num, den).as_dict()
        rows.append(
            {
                "bin_low": float(lo),
                "bin_high": float(hi),
                "numerator": eff["numerator"],
                "denominator": eff["denominator"],
                "value": eff["value"],
                "interval_low": eff["low"],
                "interval_high": eff["high"],
            }
        )
    return rows


def nested_get(mapping: dict, path: str):
    value = mapping
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(path)
        value = value[part]
    return value


def provenance_mismatches(baseline: dict, candidate: dict, paths: Sequence[str]) -> list[dict]:
    """Return exact-value provenance mismatches for paths that must agree."""
    mismatches = []
    for path in paths:
        base_missing = cand_missing = False
        try:
            base = nested_get(baseline, path)
        except KeyError:
            base_missing = True
            base = "<missing>"
        try:
            cand = nested_get(candidate, path)
        except KeyError:
            cand_missing = True
            cand = "<missing>"
        if base_missing or cand_missing or base != cand:
            mismatches.append({"path": path, "baseline": base, "candidate": cand})
    return mismatches


def evaluate_regression_policy(
    baseline: dict,
    candidate: dict,
    checks: Sequence[dict],
) -> list[dict]:
    """Evaluate generic regression checks and return detailed failures."""
    failures: list[dict] = []
    for check in checks:
        path = check["path"]
        mode = check["mode"]
        tolerance = float(check["tolerance"])
        base = float(nested_get(baseline, path))
        cand = float(nested_get(candidate, path))
        if not (math.isfinite(base) and math.isfinite(cand)):
            failures.append({"path": path, "reason": "non-finite", "baseline": base, "candidate": cand})
            continue

        failed = False
        if mode == "max_drop":
            failed = cand < base - tolerance
        elif mode == "max_increase":
            failed = cand > base + tolerance
        elif mode == "max_fractional_increase":
            failed = cand > tolerance if base == 0.0 else cand > base * (1.0 + tolerance)
        elif mode == "max_abs_change":
            failed = abs(cand - base) > tolerance
        else:
            raise ValueError(f"unknown regression mode: {mode}")

        if failed:
            failures.append(
                {
                    "path": path,
                    "mode": mode,
                    "tolerance": tolerance,
                    "baseline": base,
                    "candidate": cand,
                }
            )
    return failures
