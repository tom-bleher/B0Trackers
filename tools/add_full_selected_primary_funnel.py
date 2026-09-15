from pathlib import Path

path = Path("analysis/b0_report.py")
text = path.read_text()

anchor = '''def _plot_hist(outdir, values, filename, xlabel, title, np, bins=60, hist_range=None):
'''
plot_func = '''def _plot_full_selected_primary_funnel(outdir, funnel, np):
    import matplotlib.pyplot as plt

    stages = funnel.get("stages", [])
    if not stages:
        return
    names = [row["label"] for row in stages]
    values = np.array([row["over_generated"]["value"] for row in stages], dtype=float)
    lows = np.array([row["over_generated"]["low"] for row in stages], dtype=float)
    highs = np.array([row["over_generated"]["high"] for row in stages], dtype=float)
    yerr = np.vstack([values - lows, highs - values])

    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.errorbar(range(len(names)), values, yerr=yerr, fmt="o", capsize=4)
    ax.set_xticks(range(len(names)), names, rotation=20, ha="right")
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Fraction of generated selected primaries")
    ax.set_title("Selected-primary B0 reconstruction funnel")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "selected_primary_full_funnel.png", dpi=160)
    plt.close(fig)


'''
if anchor not in text:
    raise SystemExit("plot insertion anchor not found")
text = text.replace(anchor, plot_func + anchor, 1)

anchor = '''        truth_p = _scalar(tree, "sel_primary_p", np).astype(float)
'''
funnel = '''        selected_primary_full_funnel = None
        full_funnel_branches = truth_stage_branches | {
            "sel_primary_measurement_reconstructable",
            "sel_primary_mcIndex",
        }
        if full_funnel_branches <= keys:
            measurement_flag = _scalar(
                tree, "sel_primary_measurement_reconstructable", np
            ).astype(int)
            generated = _scalar(tree, "sel_primary_mcIndex", np).astype(int) >= 0
            measurement_known = measurement_flag >= 0
            base = generated & measurement_known
            reachable_stage = base & eligible
            reconstructable_stage = reachable_stage & (measurement_flag == 1)
            seed_stage = reconstructable_stage & _scalar(
                tree, "sel_primary_has_seed", np
            ).astype(bool)
            unfiltered_stage = seed_stage & _scalar(
                tree, "sel_primary_has_unfiltered_track", np
            ).astype(bool)
            filtered_stage = unfiltered_stage & _scalar(
                tree, "sel_primary_has_filtered_track", np
            ).astype(bool)
            matched_stage = filtered_stage & _scalar(
                tree, "sel_primary_has_truth_matched_track", np
            ).astype(bool)

            ordered_stages = [
                ("truth_reachable", "Truth reachable", reachable_stage),
                ("measurement_reconstructable", "Measurement reconstructable", reconstructable_stage),
                ("seeded", "Seeded", seed_stage),
                ("unfiltered", "CKF candidate", unfiltered_stage),
                ("filtered", "After ambiguity", filtered_stage),
                ("truth_matched", "Truth matched", matched_stage),
            ]
            base_count = int(np.count_nonzero(base))
            stage_rows = []
            for name, label, mask in ordered_stages:
                stage_rows.append(
                    {
                        "name": name,
                        "label": label,
                        "events": int(np.count_nonzero(mask)),
                        "over_generated": wilson_efficiency(
                            int(np.count_nonzero(mask)), base_count
                        ).as_dict(),
                    }
                )

            sequential = {}
            previous_name = "generated"
            previous_mask = base
            for name, _, mask in ordered_stages:
                denominator = int(np.count_nonzero(previous_mask))
                numerator = int(np.count_nonzero(mask & previous_mask))
                sequential[f"{name}_over_{previous_name}"] = wilson_efficiency(
                    numerator, denominator
                ).as_dict()
                previous_name = name
                previous_mask = mask

            selected_primary_full_funnel = {
                "definition": (
                    "Selected-primary funnel using dominant-truth measurement attribution before seeding"
                ),
                "generated_events_with_known_measurement_truth": base_count,
                "unknown_measurement_truth_events": int(
                    np.count_nonzero(generated & ~measurement_known)
                ),
                "measurement_station_requirement": (
                    int(_scalar(tree, "min_measurement_stations_required", np)[0])
                    if "min_measurement_stations_required" in keys and tree.num_entries > 0
                    else None
                ),
                "stages": stage_rows,
                "sequential": sequential,
            }

'''
if anchor not in text:
    raise SystemExit("funnel insertion anchor not found")
text = text.replace(anchor, funnel + anchor, 1)

old_report = '''            "selected_primary_stages": selected_primary_stages,
            "truth_matched": truth_metrics,
'''
new_report = '''            "selected_primary_stages": selected_primary_stages,
            "selected_primary_full_funnel": selected_primary_full_funnel,
            "truth_matched": truth_metrics,
'''
if old_report not in text:
    raise SystemExit("report insertion anchor not found")
text = text.replace(old_report, new_report, 1)

old_plot = '''            if selected_primary_stages is not None:
                _plot_stage_flow(outdir, selected_primary_stages, np, truth_specific=True)
            _plot_binned_efficiency(
'''
new_plot = '''            if selected_primary_stages is not None:
                _plot_stage_flow(outdir, selected_primary_stages, np, truth_specific=True)
            if selected_primary_full_funnel is not None:
                _plot_full_selected_primary_funnel(outdir, selected_primary_full_funnel, np)
            _plot_binned_efficiency(
'''
if old_plot not in text:
    raise SystemExit("plot call anchor not found")
text = text.replace(old_plot, new_plot, 1)

path.write_text(text)
