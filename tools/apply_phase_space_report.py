#!/usr/bin/env python3
from pathlib import Path


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {count}")
    return text.replace(old, new, 1)


# Fix binned-efficiency key collision: bin edges and Wilson interval bounds
# previously both used `low`/`high`, so the interval silently overwrote the bin.
metric = Path("analysis/b0_metrics.py")
m = metric.read_text()
m = replace_once(
    m,
    '        rows.append({"low": float(lo), "high": float(hi), **wilson_efficiency(num, den).as_dict()})\n',
    '''        eff = wilson_efficiency(num, den).as_dict()
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
''',
    "binned efficiency fields",
)
metric.write_text(m)

path = Path("analysis/b0_report.py")
text = path.read_text()

text = replace_once(
    text,
    'def _plot_residual_rms(outdir, residuals, np, filename, title):\n',
    '''def _plot_binned_efficiency(outdir, rows, filename, xlabel, title, np):
    import matplotlib.pyplot as plt

    valid = [
        row
        for row in rows
        if row["denominator"] > 0 and math.isfinite(float(row["value"]))
    ]
    if not valid:
        return
    centers = np.array(
        [(row["bin_low"] + row["bin_high"]) / 2.0 for row in valid], dtype=float
    )
    half_width = np.array(
        [(row["bin_high"] - row["bin_low"]) / 2.0 for row in valid], dtype=float
    )
    values = np.array([row["value"] for row in valid], dtype=float)
    interval_low = np.array([row["interval_low"] for row in valid], dtype=float)
    interval_high = np.array([row["interval_high"] for row in valid], dtype=float)

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.errorbar(
        centers,
        values,
        xerr=half_width,
        yerr=np.vstack([values - interval_low, interval_high - values]),
        fmt="o",
        capsize=3,
    )
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Truth-matched efficiency")
    ax.set_title(title)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=160)
    plt.close(fig)


def _selected_primary_patterns(tree, eligible, matched, ak):
    needed = {"stationP", "sideP", "isSelPrimaryP"}
    if not needed <= set(tree.keys()):
        return {}
    stations = _jagged(tree, "stationP", ak).to_list()
    sides = _jagged(tree, "sideP", ak).to_list()
    selected = _jagged(tree, "isSelPrimaryP", ak).to_list()
    counts: dict[str, dict[str, int]] = {}
    for event, (event_stations, event_sides, event_selected) in enumerate(
        zip(stations, sides, selected)
    ):
        if not eligible[event]:
            continue
        pairs = sorted(
            {
                (int(station), int(side))
                for station, side, flag in zip(event_stations, event_sides, event_selected)
                if int(flag) == 1 and int(station) > 0
            }
        )
        key = (
            "none"
            if not pairs
            else ",".join(
                f"S{station}{'F' if side == 1 else 'B' if side == 0 else 'U'}"
                for station, side in pairs
            )
        )
        row = counts.setdefault(key, {"events": 0, "truth_matched": 0})
        row["events"] += 1
        row["truth_matched"] += int(bool(matched[event]))

    out = {}
    for key, row in sorted(counts.items()):
        out[key] = {
            **row,
            "efficiency": wilson_efficiency(
                row["truth_matched"], row["events"]
            ).as_dict(),
        }
    return out


def _plot_residual_rms(outdir, residuals, np, filename, title):
''',
    "phase-space plotting helpers",
)

text = replace_once(
    text,
    '''    momentum_bins: list[float],
    angle_bins_mrad: list[float],
    manifest: dict | None = None,
):
''',
    '''    momentum_bins: list[float],
    angle_bins_mrad: list[float],
    phi_bins_rad: list[float] | None = None,
    min_truth_weight: float = 0.5,
    manifest: dict | None = None,
):
''',
    "report signature",
)
text = replace_once(
    text,
    '    ak, np, uproot = _require_dependencies()\n\n    with uproot.open(path) as root_file:\n',
    '    ak, np, uproot = _require_dependencies()\n'
    '    if phi_bins_rad is None:\n'
    '        phi_bins_rad = [-math.pi, -3 * math.pi / 4, -math.pi / 2, -math.pi / 4, 0.0, '\
    'math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi]\n\n'
    '    with uproot.open(path) as root_file:\n',
    "default phi bins",
)

text = replace_once(
    text,
    '''            "efficiency_vs_scattering_angle_mrad": binned_efficiency(
                truth_theta_mrad.tolist(), matched.tolist(), eligible.tolist(), angle_bins_mrad
            ),
            "truth_phi": distribution_summary(truth_phi[eligible].tolist()),
''',
    '''            "efficiency_vs_scattering_angle_mrad": binned_efficiency(
                truth_theta_mrad.tolist(), matched.tolist(), eligible.tolist(), angle_bins_mrad
            ),
            "efficiency_vs_truth_phi_rad": binned_efficiency(
                truth_phi.tolist(), matched.tolist(), eligible.tolist(), phi_bins_rad
            ),
            "truth_phi": distribution_summary(truth_phi[eligible].tolist()),
''',
    "phi efficiency",
)

text = replace_once(
    text,
    '        seed_metrics = {}\n',
    '        station_side_patterns = _selected_primary_patterns(tree, eligible, matched, ak)\n\n'
    '        seed_metrics = {}\n',
    "station pattern calculation",
)

text = replace_once(
    text,
    '''            primary_assoc = (assoc_idx == sel_idx_b) & (assoc_col == sel_col_b)
            n_primary_tracks = np.asarray(ak.to_numpy(ak.sum(primary_assoc, axis=1)), dtype=int)
''',
    '''            primary_assoc = (
                (assoc_idx == sel_idx_b)
                & (assoc_col == sel_col_b)
                & np.isfinite(assoc_w)
                & (assoc_w >= min_truth_weight)
            )
            n_primary_tracks = np.asarray(ak.to_numpy(ak.sum(primary_assoc, axis=1)), dtype=int)
''',
    "weighted primary association",
)

text = replace_once(
    text,
    '''            flat_weight = np.asarray(ak.to_numpy(ak.flatten(assoc_w, axis=1)), dtype=float)
            flat_idx = np.asarray(ak.to_numpy(ak.flatten(assoc_idx, axis=1)), dtype=int)
            low_purity = (flat_idx < 0) | ~np.isfinite(flat_weight) | (flat_weight < 0.5)
            association_quality["low_purity_or_unmatched_track_fraction_proxy"] = (
                float(np.mean(low_purity)) if low_purity.size else math.nan
            )
            association_quality["note"] = (
                "This is an association-quality proxy, not a formal fake rate. A formal fake definition should be "
                "fixed for the analysis sample and validated after stable track identity is exported."
            )
''',
    '''            eligible_assoc_w = assoc_w[eligible]
            eligible_assoc_idx = assoc_idx[eligible]
            flat_weight = np.asarray(
                ak.to_numpy(ak.flatten(eligible_assoc_w, axis=1)), dtype=float
            )
            flat_idx = np.asarray(
                ak.to_numpy(ak.flatten(eligible_assoc_idx, axis=1)), dtype=int
            )
            association_fake = (
                (flat_idx < 0)
                | ~np.isfinite(flat_weight)
                | (flat_weight < min_truth_weight)
            )
            association_quality["association_defined_fake_fraction"] = (
                float(np.mean(association_fake)) if association_fake.size else math.nan
            )
            association_quality["truth_weight_threshold"] = min_truth_weight
            association_quality["definition"] = (
                "Within truth-reachable events, a reconstructed track is association-matched when its dominant "
                f"MC association weight is at least {min_truth_weight:.3g}; otherwise it is an association-defined "
                "fake. A duplicate-primary event contains more than one matched track dominantly associated to "
                "the selected primary."
            )
''',
    "association fake definition",
)

text = replace_once(
    text,
    '            "truth_matched": truth_metrics,\n',
    '            "truth_matched": truth_metrics,\n'
    '            "selected_primary_station_side_patterns": station_side_patterns,\n',
    "pattern report field",
)

text = replace_once(
    text,
    '''            if selected_primary_stages is not None:
                _plot_stage_flow(outdir, selected_primary_stages, np, truth_specific=True)
            _plot_hist(
''',
    '''            if selected_primary_stages is not None:
                _plot_stage_flow(outdir, selected_primary_stages, np, truth_specific=True)
            _plot_binned_efficiency(
                outdir,
                truth_metrics["efficiency_vs_truth_p_GeV"],
                "efficiency_vs_truth_p.png",
                "truth p [GeV]",
                "B0 truth-matched efficiency vs momentum",
                np,
            )
            _plot_binned_efficiency(
                outdir,
                truth_metrics["efficiency_vs_scattering_angle_mrad"],
                "efficiency_vs_scattering_angle.png",
                "truth scattering angle [mrad]",
                "B0 truth-matched efficiency vs scattering angle",
                np,
            )
            _plot_binned_efficiency(
                outdir,
                truth_metrics["efficiency_vs_truth_phi_rad"],
                "efficiency_vs_truth_phi.png",
                "truth phi [rad]",
                "B0 truth-matched efficiency vs phi",
                np,
            )
            _plot_hist(
''',
    "phase-space plots",
)

text = replace_once(
    text,
    '''    parser.add_argument(
        "--manifest",
''',
    '''    parser.add_argument(
        "--phi-bins-rad",
        type=_parse_edges,
        default=_parse_edges(
            "-3.1415926536,-2.3561944902,-1.5707963268,-0.7853981634,0,"
            "0.7853981634,1.5707963268,2.3561944902,3.1415926536"
        ),
        help="truth phi bin edges in radians",
    )
    parser.add_argument(
        "--min-truth-weight",
        type=float,
        default=0.5,
        help="dominant track-to-MC association threshold used for fake/duplicate definitions",
    )
    parser.add_argument(
        "--manifest",
''',
    "phi/fake CLI",
)

text = replace_once(
    text,
    '''        momentum_bins=args.momentum_bins,
        angle_bins_mrad=args.angle_bins_mrad,
        manifest=manifest,
''',
    '''        momentum_bins=args.momentum_bins,
        angle_bins_mrad=args.angle_bins_mrad,
        phi_bins_rad=args.phi_bins_rad,
        min_truth_weight=args.min_truth_weight,
        manifest=manifest,
''',
    "build call args",
)
path.write_text(text)

# Update unit tests for distinct bin edge and interval fields.
test = Path("tests/test_b0_metrics.py")
t = test.read_text()
t = replace_once(
    t,
    '        self.assertEqual(result[1]["numerator"], 1)\n',
    '        self.assertEqual(result[1]["numerator"], 1)\n'
    '        self.assertEqual(result[0]["bin_low"], 0.0)\n'
    '        self.assertEqual(result[0]["bin_high"], 5.0)\n'
    '        self.assertLessEqual(result[0]["interval_low"], result[0]["value"])\n'
    '        self.assertGreaterEqual(result[0]["interval_high"], result[0]["value"])\n',
    "binned field tests",
)
test.write_text(t)
