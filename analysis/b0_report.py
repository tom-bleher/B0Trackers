#!/usr/bin/env python3
"""Create an interpretable B0 tracking validation report from B0Trackers ROOT output."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

try:
    from .b0_metrics import (
        binned_efficiency,
        distribution_summary,
        pull_summary,
        stage_efficiencies,
        wilson_efficiency,
    )
except ImportError:
    from b0_metrics import (
        binned_efficiency,
        distribution_summary,
        pull_summary,
        stage_efficiencies,
        wilson_efficiency,
    )


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


def _jagged(tree, name, ak):
    return tree[name].array(library="ak")


def _flat(tree, name, ak, np):
    array = _jagged(tree, name, ak)
    return np.asarray(ak.to_numpy(ak.flatten(array, axis=1)))


def _optional_scalar(tree, name, np, default):
    if name not in tree.keys():
        return np.full(tree.num_entries, default)
    return _scalar(tree, name, np)


def _unique_strings(tree, name, ak):
    if name not in tree.keys():
        return []
    values = _jagged(tree, name, ak).to_list()
    return sorted({str(value) for value in values})


def _parse_edges(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if len(values) < 2 or any(b <= a for a, b in zip(values, values[1:])):
        raise argparse.ArgumentTypeError("bin edges must be a strictly increasing comma-separated list")
    return values


def _mapping_summary(methods, types, np):
    methods = np.asarray(methods, dtype=int)
    types = np.asarray(types, dtype=int)
    physics = (types & (1 | 4 | 8)) != 0  # measurement | outlier | hole
    exact = int(np.count_nonzero(physics & (methods == 1)))
    fallback = int(np.count_nonzero(physics & (methods == 2)))
    unresolved = int(np.count_nonzero(physics & (methods == 0)))
    resolved = exact + fallback
    total = resolved + unresolved
    return {
        "exact": exact,
        "fallback": fallback,
        "unresolved_physics_states": unresolved,
        "exact_fraction_of_resolved": exact / resolved if resolved else math.nan,
        "fallback_fraction_of_resolved": fallback / resolved if resolved else math.nan,
        "unresolved_fraction_of_physics_states": unresolved / total if total else math.nan,
    }


def _station_residual_summary(station, state_type, resid0, resid1, np, extra_mask=None):
    station = np.asarray(station, dtype=int)
    state_type = np.asarray(state_type, dtype=int)
    resid0 = np.asarray(resid0, dtype=float)
    resid1 = np.asarray(resid1, dtype=float)
    mask = (state_type & 1) != 0
    if extra_mask is not None:
        mask &= np.asarray(extra_mask, dtype=bool)

    out: dict[str, dict] = {}
    for value in sorted(set(station[mask].tolist())):
        if value <= 0:
            continue
        use = mask & (station == value)
        out[str(int(value))] = {
            "loc0": distribution_summary(resid0[use].tolist()),
            "loc1": distribution_summary(resid1[use].tolist()),
        }
    return out


def _plot_stage_flow(outdir, stages, np, *, truth_specific: bool):
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
    ax.set_ylabel("Fraction of truth-reachable events")
    if truth_specific:
        ax.set_title("Selected-primary B0 reconstruction flow")
        filename = "selected_primary_stage_flow.png"
    else:
        ax.set_title("Event-level B0 object presence (not primary-specific)")
        filename = "event_stage_presence.png"
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=160)
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


def _plot_residual_rms(outdir, residuals, np, filename, title):
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
    ax.set_title(title)
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=160)
    plt.close(fig)


def build_report(
    path: Path,
    tree_name: str,
    min_stations: int,
    make_plots: bool,
    outdir: Path,
    *,
    momentum_bins: list[float],
    angle_bins_mrad: list[float],
    manifest: dict | None = None,
):
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
            "sel_primary_px",
            "sel_primary_py",
            "sel_primary_thscat_mrad",
            "sel_primary_mcIndex",
            "sel_primary_mcCollectionID",
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
                "cannot interpret B0 reconstruction metrics because these factories/collections were unavailable: "
                + ", ".join(unavailable)
                + ". Missing is configuration absence, not reconstruction inefficiency."
            )

        schema = _scalar(tree, "schema_version", np)
        schema_versions = sorted(set(int(x) for x in schema.tolist()))
        is_schema3 = bool(schema_versions and min(schema_versions) >= 3)
        if is_schema3 and "has_ckf_assocs_unfiltered" in keys:
            availability["has_ckf_assocs_unfiltered"] = bool(
                np.all(_scalar(tree, "has_ckf_assocs_unfiltered", np).astype(bool))
            )
            if not availability["has_ckf_assocs_unfiltered"]:
                raise SystemExit(
                    "schema-3 selected-primary stage efficiencies require "
                    "B0TrackerCKFTrackUnfilteredAssociations; the collection was unavailable"
                )
        elif is_schema3:
            raise SystemExit("schema 3 is missing has_ckf_assocs_unfiltered")
        stations = _scalar(tree, "n_stations_primary", np)
        n_seeds = _scalar(tree, "n_stub_seeds", np)
        n_unfiltered = _scalar(tree, "n_ckf_unfiltered", np)
        n_filtered = _scalar(tree, "n_ckf_filtered", np)
        matched_index = _scalar(tree, "ckf_truth_matched_trk_index", np).astype(int)

        eligible = stations >= min_stations
        matched = matched_index >= 0

        # Schema 2 only knows whether *any* seed/candidate exists in an event.
        # Keep these useful event-presence diagnostics, but do not call the
        # first three stages selected-primary efficiencies.
        event_stage_presence = stage_efficiencies(
            eligible.tolist(),
            (n_seeds > 0).tolist(),
            (n_unfiltered > 0).tolist(),
            (n_filtered > 0).tolist(),
            matched.tolist(),
        )

        selected_primary_stages = None
        truth_stage_branches = {
            "sel_primary_has_seed",
            "sel_primary_has_unfiltered_track",
            "sel_primary_has_filtered_track",
            "sel_primary_has_truth_matched_track",
        }
        if truth_stage_branches <= keys:
            selected_primary_stages = stage_efficiencies(
                eligible.tolist(),
                _scalar(tree, "sel_primary_has_seed", np).astype(bool).tolist(),
                _scalar(tree, "sel_primary_has_unfiltered_track", np).astype(bool).tolist(),
                _scalar(tree, "sel_primary_has_filtered_track", np).astype(bool).tolist(),
                _scalar(tree, "sel_primary_has_truth_matched_track", np).astype(bool).tolist(),
            )

        truth_p = _scalar(tree, "sel_primary_p", np).astype(float)
        truth_theta_mrad = _scalar(tree, "sel_primary_thscat_mrad", np).astype(float)
        truth_px = _scalar(tree, "sel_primary_px", np).astype(float)
        truth_py = _scalar(tree, "sel_primary_py", np).astype(float)
        truth_phi = np.arctan2(truth_py, truth_px)
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
            "efficiency_vs_truth_p_GeV": binned_efficiency(
                truth_p.tolist(), matched.tolist(), eligible.tolist(), momentum_bins
            ),
            "efficiency_vs_scattering_angle_mrad": binned_efficiency(
                truth_theta_mrad.tolist(), matched.tolist(), eligible.tolist(), angle_bins_mrad
            ),
            "truth_phi": distribution_summary(truth_phi[eligible].tolist()),
        }

        # Resolution is meaningful only for truth-matched tracks. Report robust
        # width and bias versus truth momentum/angle in the same explicit bins.
        truth_metrics["relative_delta_p_vs_truth_p_GeV"] = []
        for lo, hi in zip(momentum_bins, momentum_bins[1:]):
            use = matched_mask & (truth_p >= lo) & (truth_p < hi)
            truth_metrics["relative_delta_p_vs_truth_p_GeV"].append(
                {"low": lo, "high": hi, **distribution_summary(rel_delta_p[use].tolist())}
            )
        truth_metrics["relative_delta_p_vs_scattering_angle_mrad"] = []
        for lo, hi in zip(angle_bins_mrad, angle_bins_mrad[1:]):
            use = matched_mask & (truth_theta_mrad >= lo) & (truth_theta_mrad < hi)
            truth_metrics["relative_delta_p_vs_scattering_angle_mrad"].append(
                {"low": lo, "high": hi, **distribution_summary(rel_delta_p[use].tolist())}
            )

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

        mapping: dict[str, dict] = {}
        for label, prefix in (("stub_ckf", "ckf_trk_"), ("truth_seeded", "trk_")):
            needed = {prefix + "state_mapping_method", prefix + "state_type"}
            if needed <= keys:
                mapping[label] = _mapping_summary(
                    _flat(tree, prefix + "state_mapping_method", ak, np),
                    _flat(tree, prefix + "state_type", ak, np),
                    np,
                )

        residuals: dict[str, dict] = {}
        state_required = {
            "ckf_trk_state_resid_loc0",
            "ckf_trk_state_resid_loc1",
            "ckf_trk_state_type",
            "ckf_trk_aclgad_station",
            "ckf_trk_state_track_index",
        }
        if state_required <= keys:
            resid0_j = _jagged(tree, "ckf_trk_state_resid_loc0", ak)
            resid1_j = _jagged(tree, "ckf_trk_state_resid_loc1", ak)
            type_j = _jagged(tree, "ckf_trk_state_type", ak)
            station_j = _jagged(tree, "ckf_trk_aclgad_station", ak)
            track_index_j = _jagged(tree, "ckf_trk_state_track_index", ak)

            resid0 = np.asarray(ak.to_numpy(ak.flatten(resid0_j, axis=1)), dtype=float)
            resid1 = np.asarray(ak.to_numpy(ak.flatten(resid1_j, axis=1)), dtype=float)
            state_type = np.asarray(ak.to_numpy(ak.flatten(type_j, axis=1)), dtype=int)
            state_station = np.asarray(ak.to_numpy(ak.flatten(station_j, axis=1)), dtype=int)
            state_track_index = np.asarray(ak.to_numpy(ak.flatten(track_index_j, axis=1)), dtype=int)

            residuals["all_stub_ckf_tracks"] = _station_residual_summary(
                state_station, state_type, resid0, resid1, np
            )

            stable_state_identity = {
                "ckf_trk_state_parent_track_index",
                "ckf_trk_state_parent_track_collectionID",
                "ckf_trk_state_parent_identity_valid",
                "ckf_truth_matched_trk_object_index",
                "ckf_truth_matched_trk_object_collectionID",
                "ckf_truth_matched_trk_identity_valid",
            }
            selected_state_mask = None
            if stable_state_identity <= keys:
                parent_index_j = _jagged(tree, "ckf_trk_state_parent_track_index", ak)
                parent_collection_j = _jagged(tree, "ckf_trk_state_parent_track_collectionID", ak)
                parent_valid_j = _jagged(tree, "ckf_trk_state_parent_identity_valid", ak)
                matched_object_index = _scalar(tree, "ckf_truth_matched_trk_object_index", np).astype(int)
                matched_object_collection = _scalar(
                    tree, "ckf_truth_matched_trk_object_collectionID", np
                ).astype(int)
                matched_identity_valid = _scalar(
                    tree, "ckf_truth_matched_trk_identity_valid", np
                ).astype(bool)

                parent_index = np.asarray(
                    ak.to_numpy(ak.flatten(parent_index_j, axis=1)), dtype=int
                )
                parent_collection = np.asarray(
                    ak.to_numpy(ak.flatten(parent_collection_j, axis=1)), dtype=int
                )
                parent_valid = np.asarray(
                    ak.to_numpy(ak.flatten(parent_valid_j, axis=1)), dtype=bool
                )
                repeated_match_index = np.asarray(
                    ak.to_numpy(
                        ak.flatten(ak.broadcast_arrays(parent_index_j, matched_object_index)[1], axis=1)
                    ),
                    dtype=int,
                )
                repeated_match_collection = np.asarray(
                    ak.to_numpy(
                        ak.flatten(
                            ak.broadcast_arrays(parent_collection_j, matched_object_collection)[1], axis=1
                        )
                    ),
                    dtype=int,
                )
                repeated_match_valid = np.asarray(
                    ak.to_numpy(
                        ak.flatten(ak.broadcast_arrays(parent_valid_j, matched_identity_valid)[1], axis=1)
                    ),
                    dtype=bool,
                )
                selected_state_mask = (
                    parent_valid
                    & repeated_match_valid
                    & (parent_index == repeated_match_index)
                    & (parent_collection == repeated_match_collection)
                )
                residuals["truth_matched_track_stable_identity"] = _station_residual_summary(
                    state_station, state_type, resid0, resid1, np, selected_state_mask
                )
            else:
                # Schema 2 fallback: retain the explicitly labelled positional view.
                repeated_match = np.asarray(
                    ak.to_numpy(ak.flatten(ak.broadcast_arrays(track_index_j, matched_index)[1], axis=1)),
                    dtype=int,
                )
                selected_state_mask = (repeated_match >= 0) & (state_track_index == repeated_match)
                residuals["nominal_truth_matched_track_schema2_positional"] = _station_residual_summary(
                    state_station, state_type, resid0, resid1, np, selected_state_mask
                )

            innovation_required = {
                "ckf_trk_state_innov_chi2",
                "ckf_trk_state_innov_pull0",
                "ckf_trk_state_innov_pull1",
            }
            if innovation_required <= keys:
                innov_chi2 = _flat(tree, "ckf_trk_state_innov_chi2", ak, np).astype(float)
                innov_pull0 = _flat(tree, "ckf_trk_state_innov_pull0", ak, np).astype(float)
                innov_pull1 = _flat(tree, "ckf_trk_state_innov_pull1", ak, np).astype(float)
                innovation = {
                    "all_stub_ckf_measurements": {
                        "chi2": distribution_summary(innov_chi2.tolist()),
                        "pull0": pull_summary(innov_pull0.tolist()),
                        "pull1": pull_summary(innov_pull1.tolist()),
                    },
                    "truth_matched_track": {
                        "chi2": distribution_summary(innov_chi2[selected_state_mask].tolist()),
                        "pull0": pull_summary(innov_pull0[selected_state_mask].tolist()),
                        "pull1": pull_summary(innov_pull1[selected_state_mask].tolist()),
                    },
                    "by_station_truth_matched": {},
                }
                for station in sorted(set(state_station[selected_state_mask].tolist())):
                    if station <= 0:
                        continue
                    use = selected_state_mask & (state_station == station)
                    innovation["by_station_truth_matched"][str(int(station))] = {
                        "chi2": distribution_summary(innov_chi2[use].tolist()),
                        "pull0": pull_summary(innov_pull0[use].tolist()),
                        "pull1": pull_summary(innov_pull1[use].tolist()),
                    }
            else:
                innovation = {}

        association_quality = {}
        assoc_needed = {
            "ckf_trk_assoc_mcIndex",
            "ckf_trk_assoc_mcCollectionID",
            "ckf_trk_assoc_weight",
        }
        if assoc_needed <= keys:
            assoc_idx = _jagged(tree, "ckf_trk_assoc_mcIndex", ak)
            assoc_col = _jagged(tree, "ckf_trk_assoc_mcCollectionID", ak)
            assoc_w = _jagged(tree, "ckf_trk_assoc_weight", ak)
            sel_idx = _scalar(tree, "sel_primary_mcIndex", np).astype(int)
            sel_col = _scalar(tree, "sel_primary_mcCollectionID", np).astype(int)

            n_tracks = np.asarray(ak.to_numpy(ak.num(assoc_idx, axis=1)), dtype=int)
            sel_idx_b = ak.broadcast_arrays(assoc_idx, sel_idx)[1]
            sel_col_b = ak.broadcast_arrays(assoc_col, sel_col)[1]
            primary_assoc = (assoc_idx == sel_idx_b) & (assoc_col == sel_col_b)
            n_primary_tracks = np.asarray(ak.to_numpy(ak.sum(primary_assoc, axis=1)), dtype=int)
            duplicate_events = eligible & (n_primary_tracks > 1)
            association_quality = {
                "eligible_events_with_duplicate_primary_tracks": wilson_efficiency(
                    int(np.count_nonzero(duplicate_events)), int(np.count_nonzero(eligible))
                ).as_dict(),
                "primary_associated_tracks_per_eligible_event": distribution_summary(
                    n_primary_tracks[eligible].tolist()
                ),
                "all_ckf_tracks_per_eligible_event": distribution_summary(n_tracks[eligible].tolist()),
            }
            flat_weight = np.asarray(ak.to_numpy(ak.flatten(assoc_w, axis=1)), dtype=float)
            flat_idx = np.asarray(ak.to_numpy(ak.flatten(assoc_idx, axis=1)), dtype=int)
            low_purity = (flat_idx < 0) | ~np.isfinite(flat_weight) | (flat_weight < 0.5)
            association_quality["low_purity_or_unmatched_track_fraction_proxy"] = (
                float(np.mean(low_purity)) if low_purity.size else math.nan
            )
            association_quality["note"] = (
                "This is an association-quality proxy, not a formal fake rate. A formal fake definition should be "
                "fixed for the analysis sample and validated after stable track identity is exported."
            )

        diagnostics = {
            name: int(np.sum(_optional_scalar(tree, name, np, 0)))
            for name in [
                "n_simhits_unresolved_cellid",
                "n_missing_mc_relation",
                "n_pixel_snap_failed",
            ]
        }

        provenance = {
            "tree": tree_name,
            "schema_versions": schema_versions,
            "geometry_names": _unique_strings(tree, "geometry_name", ak),
            "manifest": manifest or {},
        }

        report = {
            "input": str(path),
            "events": int(tree.num_entries),
            "provenance": provenance,
            "availability": availability,
            "selection": {
                "definition": f"selected primary crosses at least {min_stations} truth B0 stations",
                "min_stations": min_stations,
                "note": (
                    "This is a truth-reachability denominator, not a measurement-level reconstructability definition. "
                    "Use it primarily for controlled single-primary samples."
                ),
            },
            "event_stage_presence": event_stage_presence,
            "event_stage_interpretation": (
                "In schema 2, seed/unfiltered/filtered mean at least one such object exists anywhere in the event. "
                "Only truth_matched is selected-primary-specific. Do not quote the first three as primary efficiency."
            ),
            "selected_primary_stages": selected_primary_stages,
            "truth_matched": truth_metrics,
            "seed_survival": seed_metrics,
            "sensor_mapping": mapping,
            "measurement_residuals_by_station": residuals,
            "predicted_innovation": innovation,
            "residual_interpretation": (
                "Legacy state residuals are fit diagnostics, not intrinsic detector-resolution measurements. "
                "Schema 3 selects truth-matched states by stable parent-track identity and separately exports "
                "the pre-update predicted innovation and normalized innovation chi2 for CKF compatibility."
            ),
            "track_association_quality": association_quality,
            "diagnostic_failures": diagnostics,
        }

        if make_plots:
            outdir.mkdir(parents=True, exist_ok=True)
            _plot_stage_flow(outdir, event_stage_presence, np, truth_specific=False)
            if selected_primary_stages is not None:
                _plot_stage_flow(outdir, selected_primary_stages, np, truth_specific=True)
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
            if residuals.get("all_stub_ckf_tracks"):
                _plot_residual_rms(
                    outdir,
                    residuals["all_stub_ckf_tracks"],
                    np,
                    "measurement_residual_rms_all_tracks.png",
                    "All stub-CKF measurement-state residual RMS",
                )
            if residuals.get("nominal_truth_matched_track_schema2_positional"):
                _plot_residual_rms(
                    outdir,
                    residuals["nominal_truth_matched_track_schema2_positional"],
                    np,
                    "measurement_residual_rms_truth_matched_nominal.png",
                    "Nominal truth-matched track residual RMS (schema-2 positional)",
                )

    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="B0Trackers ROOT histogram file")
    parser.add_argument("--tree", default="B0Trackers/hits")
    parser.add_argument("--output-dir", type=Path, default=Path("b0-report"))
    parser.add_argument("--min-stations", type=int, default=3)
    parser.add_argument(
        "--momentum-bins",
        type=_parse_edges,
        default=_parse_edges("0,5,10,15,20,30,40,60,100"),
        help="truth-momentum bin edges in GeV",
    )
    parser.add_argument(
        "--angle-bins-mrad",
        type=_parse_edges,
        default=_parse_edges("0,2,4,6,8,10,14,18,24,35"),
        help="truth scattering-angle bin edges in mrad",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="optional JSON run manifest (dataset/commits/material-map hash/config) embedded in provenance",
    )
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text()) if args.manifest else None
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = build_report(
        args.input,
        args.tree,
        args.min_stations,
        not args.no_plots,
        args.output_dir,
        momentum_bins=args.momentum_bins,
        angle_bins_mrad=args.angle_bins_mrad,
        manifest=manifest,
    )
    output = args.output_dir / "summary.json"
    output.write_text(json.dumps(_sanitize(report), indent=2, sort_keys=True) + "\n")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
