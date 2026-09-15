#!/usr/bin/env python3
"""Create an interpretable B0 tracking validation report from B0Trackers ROOT output."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

try:
    from .b0_metrics import distribution_summary, pull_summary, stage_efficiencies, wilson_efficiency
except ImportError:
    from b0_metrics import distribution_summary, pull_summary, stage_efficiencies, wilson_efficiency


def _sanitize(value):
    if isinstance(value, dict):
        return {key: _sanitize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize(item) for item in value]
    if isinstance(value, tuple):
        return [_sanitize(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _require_dependencies():
    try:
        import awkward as ak
        import numpy as np
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "b0_report.py requires numpy, awkward, and uproot. "
            "Install them in the analysis environment before running this tool."
        ) from exc
    return ak, np, uproot


def _scalar(tree, name, np):
    return np.asarray(tree[name].array(library="np"))


def _flat(tree, name, ak, np):
    array = tree[name].array(library="ak")
    return np.asarray(ak.to_numpy(ak.flatten(array, axis=1)))


def _optional_scalar(tree, name, np, default):
    if name not in tree.keys():
        return np.full(tree.num_entries, default)
    return _scalar(tree, name, np)


def _plot_stage_flow(outdir, stages, np):
    import matplotlib.pyplot as plt

    names = ["Seed", "CKF candidate", "After ambiguity", "Truth matched"]
    keys = [
        "seeded_over_eligible",
        "unfiltered_over_eligible",
        "filtered_over_eligible",
        "truth_matched_over_eligible",
    ]
    values = np.array([stages[key]["value"] for key in keys], dtype=float)
    lows = np.array([stages[key]["low"] for key in keys], dtype=float)
    highs = np.array([stages[key]["high"] for key in keys], dtype=float)
    yerr = np.vstack([values - lows, highs - values])

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.errorbar(range(len(names)), values, yerr=yerr, fmt="o", capsize=4)
    ax.set_xticks(range(len(names)), names, rotation=15, ha="right")
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Fraction of truth-eligible events")
    ax.set_title("B0 reconstruction stage flow")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "stage_flow.png", dpi=160)
    plt.close(fig)


def _plot_hist(outdir, values, filename, xlabel, title, np, bins=60, hist_range=None):
    import matplotlib.pyplot as plt

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.hist(values, bins=bins, range=hist_range, histtype="step", linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    ax.set_title(title)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=160)
    plt.close(fig)


def _plot_residual_rms(outdir, residuals, np):
    import matplotlib.pyplot as plt

    stations = sorted(int(station) for station in residuals)
    if not stations:
        return
    loc0 = [residuals[str(s)]["loc0"]["rms"] for s in stations]
    loc1 = [residuals[str(s)]["loc1"]["rms"] for s in stations]
    x = np.arange(len(stations), dtype=float)
    width = 0.36

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.bar(x - width / 2, loc0, width, label="loc0")
    ax.bar(x + width / 2, loc1, width, label="loc1")
    ax.set_xticks(x, [str(s) for s in stations])
    ax.set_xlabel("B0 station")
    ax.set_ylabel("Residual RMS [mm]")
    ax.set_title("Measurement-state residual RMS")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "measurement_residual_rms_by_station.png", dpi=160)
    plt.close(fig)


def build_report(path: Path, tree_name: str, min_stations: int, make_plots: bool, outdir: Path):
    ak, np, uproot = _require_dependencies()

    with uproot.open(path) as root_file:
        if tree_name not in root_file:
            raise SystemExit(f"tree {tree_name!r} not found in {path}")
        tree = root_file[tree_name]
        keys = set(tree.keys())

        required = {
            "schema_version",
            "n_stations_primary",
            "n_stub_seeds",
            "n_ckf_unfiltered",
            "n_ckf_filtered",
            "ckf_truth_matched_trk_index",
            "sel_primary_p",
            "ckf_truth_matched_trk_delta_p",
            "ckf_truth_matched_trk_pull_qOverP",
            "ckf_truth_matched_trk_pull_theta",
            "ckf_truth_matched_trk_pull_phi",
        }
        availability_required = {
            "has_stub_seeds",
            "has_ckf_tracks_unfiltered",
            "has_ckf_tracks",
            "has_ckf_assocs",
        }
        missing = sorted((required | availability_required) - keys)
        if missing:
            raise SystemExit("missing required B0Trackers branches: " + ", ".join(missing))

        availability = {
            name: bool(np.all(_scalar(tree, name, np).astype(bool)))
            for name in sorted(availability_required)
        }
        unavailable = [name for name, present in availability.items() if not present]
        if unavailable:
            raise SystemExit(
                "cannot interpret stage efficiencies because these factories/collections were unavailable: "
                + ", ".join(unavailable)
                + ". Missing is configuration absence, not reconstruction inefficiency."
            )

        schema = _scalar(tree, "schema_version", np)
        stations = _scalar(tree, "n_stations_primary", np)
        n_seeds = _scalar(tree, "n_stub_seeds", np)
        n_unfiltered = _scalar(tree, "n_ckf_unfiltered", np)
        n_filtered = _scalar(tree, "n_ckf_filtered", np)
        matched_index = _scalar(tree, "ckf_truth_matched_trk_index", np)

        eligible = stations >= min_stations
        seeded = n_seeds > 0
        unfiltered = n_unfiltered > 0
        filtered = n_filtered > 0
        matched = matched_index >= 0
        stages = stage_efficiencies(
            eligible.tolist(),
            seeded.tolist(),
            unfiltered.tolist(),
            filtered.tolist(),
            matched.tolist(),
        )

        truth_p = _scalar(tree, "sel_primary_p", np).astype(float)
        delta_p = _scalar(tree, "ckf_truth_matched_trk_delta_p", np).astype(float)
        matched_mask = eligible & matched & np.isfinite(truth_p) & (truth_p > 0.0) & np.isfinite(delta_p)
        rel_delta_p = np.full(tree.num_entries, np.nan, dtype=float)
        rel_delta_p[matched_mask] = delta_p[matched_mask] / truth_p[matched_mask]

        truth_metrics = {
            "events": int(np.count_nonzero(matched_mask)),
            "delta_p_GeV": distribution_summary(delta_p[matched_mask].tolist()),
            "relative_delta_p": distribution_summary(rel_delta_p[matched_mask].tolist()),
            "pull_qOverP": pull_summary(
                _scalar(tree, "ckf_truth_matched_trk_pull_qOverP", np)[matched_mask].tolist()
            ),
            "pull_theta": pull_summary(
                _scalar(tree, "ckf_truth_matched_trk_pull_theta", np)[matched_mask].tolist()
            ),
            "pull_phi": pull_summary(
                _scalar(tree, "ckf_truth_matched_trk_pull_phi", np)[matched_mask].tolist()
            ),
        }

        seed_metrics = {}
        if {"seed_made_unfiltered_track", "seed_survived_ambiguity"} <= keys:
            made_unfiltered = _flat(tree, "seed_made_unfiltered_track", ak, np).astype(int)
            survived = _flat(tree, "seed_survived_ambiguity", ak, np).astype(int)
            known_ckf = made_unfiltered >= 0
            seed_metrics["made_unfiltered"] = wilson_efficiency(
                int(np.count_nonzero(made_unfiltered[known_ckf] == 1)),
                int(np.count_nonzero(known_ckf)),
            ).as_dict()
            ambiguity_den = made_unfiltered == 1
            seed_metrics["survived_ambiguity_given_ckf"] = wilson_efficiency(
                int(np.count_nonzero(survived[ambiguity_den] == 1)),
                int(np.count_nonzero(ambiguity_den)),
            ).as_dict()

        mapping = {}
        if {"n_sensor_map_exact", "n_sensor_map_fallback", "n_sensor_map_failed"} <= keys:
            exact = int(np.sum(_scalar(tree, "n_sensor_map_exact", np)))
            fallback = int(np.sum(_scalar(tree, "n_sensor_map_fallback", np)))
            failed = int(np.sum(_scalar(tree, "n_sensor_map_failed", np)))
            resolved = exact + fallback
            all_states = resolved + failed
            mapping = {
                "exact": exact,
                "fallback": fallback,
                "failed_raw": failed,
                "exact_fraction_of_resolved": exact / resolved if resolved else math.nan,
                "fallback_fraction_of_resolved": fallback / resolved if resolved else math.nan,
                "failed_fraction_all_states": failed / all_states if all_states else math.nan,
                "note": (
                    "schema-2 n_sensor_map_failed also counts non-physics states that are not expected to map "
                    "to a B0 sensor; use fallback_fraction_of_resolved for regression gating."
                ),
            }

        residuals: dict[str, dict] = {}
        state_required = {
            "ckf_trk_state_resid_loc0",
            "ckf_trk_state_resid_loc1",
            "ckf_trk_state_type",
            "ckf_trk_aclgad_station",
        }
        if state_required <= keys:
            resid0 = _flat(tree, "ckf_trk_state_resid_loc0", ak, np).astype(float)
            resid1 = _flat(tree, "ckf_trk_state_resid_loc1", ak, np).astype(float)
            state_type = _flat(tree, "ckf_trk_state_type", ak, np).astype(int)
            state_station = _flat(tree, "ckf_trk_aclgad_station", ak, np).astype(int)
            measurement_state = (state_type & 1) != 0
            for station in sorted(set(state_station[measurement_state].tolist())):
                if station <= 0:
                    continue
                mask = measurement_state & (state_station == station)
                residuals[str(int(station))] = {
                    "loc0": distribution_summary(resid0[mask].tolist()),
                    "loc1": distribution_summary(resid1[mask].tolist()),
                }

        diagnostics = {
            name: int(np.sum(_optional_scalar(tree, name, np, 0)))
            for name in [
                "n_simhits_unresolved_cellid",
                "n_missing_mc_relation",
                "n_pixel_snap_failed",
            ]
        }

        report = {
            "input": str(path),
            "tree": tree_name,
            "events": int(tree.num_entries),
            "schema_versions": sorted(set(int(x) for x in schema.tolist())),
            "availability": availability,
            "selection": {
                "definition": f"selected primary crosses at least {min_stations} truth B0 stations",
                "min_stations": min_stations,
                "note": (
                    "This is a truth-reachability denominator, not a measurement-level reconstructability definition. "
                    "Use it primarily for controlled single-primary samples."
                ),
            },
            "stages": stages,
            "truth_matched": truth_metrics,
            "seed_survival": seed_metrics,
            "sensor_mapping": mapping,
            "measurement_residuals_by_station": residuals,
            "residual_interpretation": (
                "Current schema-2 state residuals use the best available state in the order "
                "smoothed, filtered, predicted. They are fit diagnostics, not intrinsic detector-resolution measurements."
            ),
            "diagnostic_failures": diagnostics,
        }

        if make_plots:
            outdir.mkdir(parents=True, exist_ok=True)
            _plot_stage_flow(outdir, stages, np)
            _plot_hist(
                outdir,
                rel_delta_p[matched_mask],
                "relative_momentum_residual.png",
                r"$(p_{reco}-p_{truth})/p_{truth}$",
                "Truth-matched B0 momentum residual",
                np,
            )
            for branch, label, filename in [
                ("ckf_truth_matched_trk_pull_qOverP", "q/p pull", "pull_qoverp.png"),
                ("ckf_truth_matched_trk_pull_theta", "theta pull", "pull_theta.png"),
                ("ckf_truth_matched_trk_pull_phi", "phi pull", "pull_phi.png"),
            ]:
                _plot_hist(
                    outdir,
                    _scalar(tree, branch, np)[matched_mask],
                    filename,
                    label,
                    f"Truth-matched B0 {label}",
                    np,
                    hist_range=(-5.0, 5.0),
                )
            _plot_residual_rms(outdir, residuals, np)

    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="B0Trackers ROOT histogram file")
    parser.add_argument("--tree", default="B0Trackers/hits")
    parser.add_argument("--output-dir", type=Path, default=Path("b0-report"))
    parser.add_argument("--min-stations", type=int, default=3)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = build_report(
        args.input,
        args.tree,
        args.min_stations,
        not args.no_plots,
        args.output_dir,
    )
    output = args.output_dir / "summary.json"
    output.write_text(json.dumps(_sanitize(report), indent=2, sort_keys=True) + "\n")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
