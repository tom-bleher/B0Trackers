from pathlib import Path

path = Path("analysis/b0_report.py")
text = path.read_text()
old = '''            # Schema 2 uses positional trajectory/ACTS-track indices. Preserve a
            # separately labelled view for the nominal truth-matched track, but
            # make the identity assumption explicit until schema 3 exports a
            # stable parent-track ObjectID relation.
            repeated_match = np.asarray(
                ak.to_numpy(ak.flatten(ak.broadcast_arrays(track_index_j, matched_index)[1], axis=1)),
                dtype=int,
            )
            nominal_selected = (repeated_match >= 0) & (state_track_index == repeated_match)
            residuals["nominal_truth_matched_track_schema2_positional"] = _station_residual_summary(
                state_station, state_type, resid0, resid1, np, nominal_selected
            )
'''
new = '''            stable_identity_branches = {
                "ckf_trk_state_parent_track_index",
                "ckf_trk_state_parent_track_collectionID",
                "ckf_trk_state_parent_identity_valid",
                "ckf_truth_matched_trk_object_index",
                "ckf_truth_matched_trk_object_collectionID",
                "ckf_truth_matched_trk_identity_valid",
            }
            if stable_identity_branches <= keys:
                parent_idx_j = _jagged(tree, "ckf_trk_state_parent_track_index", ak)
                parent_col_j = _jagged(tree, "ckf_trk_state_parent_track_collectionID", ak)
                parent_valid_j = _jagged(tree, "ckf_trk_state_parent_identity_valid", ak)
                match_obj_idx = _scalar(tree, "ckf_truth_matched_trk_object_index", np).astype(int)
                match_obj_col = _scalar(tree, "ckf_truth_matched_trk_object_collectionID", np).astype(int)
                match_obj_valid = _scalar(tree, "ckf_truth_matched_trk_identity_valid", np).astype(bool)

                parent_idx = np.asarray(ak.to_numpy(ak.flatten(parent_idx_j, axis=1)), dtype=int)
                parent_col = np.asarray(ak.to_numpy(ak.flatten(parent_col_j, axis=1)), dtype=int)
                parent_valid = np.asarray(
                    ak.to_numpy(ak.flatten(parent_valid_j, axis=1)), dtype=bool
                )
                repeated_idx = np.asarray(
                    ak.to_numpy(ak.flatten(ak.broadcast_arrays(parent_idx_j, match_obj_idx)[1], axis=1)),
                    dtype=int,
                )
                repeated_col = np.asarray(
                    ak.to_numpy(ak.flatten(ak.broadcast_arrays(parent_col_j, match_obj_col)[1], axis=1)),
                    dtype=int,
                )
                repeated_valid = np.asarray(
                    ak.to_numpy(
                        ak.flatten(ak.broadcast_arrays(parent_valid_j, match_obj_valid)[1], axis=1)
                    ),
                    dtype=bool,
                )
                stable_selected = (
                    parent_valid
                    & repeated_valid
                    & (parent_idx == repeated_idx)
                    & (parent_col == repeated_col)
                )
                residuals["truth_matched_track_stable_identity"] = _station_residual_summary(
                    state_station, state_type, resid0, resid1, np, stable_selected
                )
            else:
                # Schema 2 uses positional trajectory/ACTS-track indices. Preserve a
                # separately labelled view for the nominal truth-matched track, but
                # make the identity assumption explicit until schema 3 is available.
                repeated_match = np.asarray(
                    ak.to_numpy(
                        ak.flatten(ak.broadcast_arrays(track_index_j, matched_index)[1], axis=1)
                    ),
                    dtype=int,
                )
                nominal_selected = (repeated_match >= 0) & (state_track_index == repeated_match)
                residuals["nominal_truth_matched_track_schema2_positional"] = (
                    _station_residual_summary(
                        state_station, state_type, resid0, resid1, np, nominal_selected
                    )
                )
'''
if old not in text:
    raise SystemExit("schema2 residual selection block not found")
text = text.replace(old, new, 1)
old_note = '''            "residual_interpretation": (
                "Schema-2 state residuals use the best available state in the order smoothed, filtered, predicted. "
                "They are fit diagnostics, not intrinsic detector-resolution measurements. The nominal selected-track "
                "view additionally relies on schema-2 positional ACTS/trajectory alignment."
            ),
'''
new_note = '''            "residual_interpretation": (
                "State residuals are fit diagnostics, not intrinsic detector-resolution measurements. "
                "When schema-3 parent-track ObjectIDs are available, the selected truth-matched view uses "
                "stable PODIO identity; otherwise the explicitly labelled schema-2 nominal view relies on "
                "positional ACTS/trajectory alignment."
            ),
'''
if old_note not in text:
    raise SystemExit("residual interpretation block not found")
text = text.replace(old_note, new_note, 1)
old_plot = '''            if residuals.get("nominal_truth_matched_track_schema2_positional"):
                _plot_residual_rms(
                    outdir,
                    residuals["nominal_truth_matched_track_schema2_positional"],
                    np,
                    "measurement_residual_rms_truth_matched_nominal.png",
                    "Nominal truth-matched track residual RMS (schema-2 positional)",
                )
'''
new_plot = '''            if residuals.get("truth_matched_track_stable_identity"):
                _plot_residual_rms(
                    outdir,
                    residuals["truth_matched_track_stable_identity"],
                    np,
                    "measurement_residual_rms_truth_matched_stable.png",
                    "Truth-matched track residual RMS (stable ObjectID)",
                )
            elif residuals.get("nominal_truth_matched_track_schema2_positional"):
                _plot_residual_rms(
                    outdir,
                    residuals["nominal_truth_matched_track_schema2_positional"],
                    np,
                    "measurement_residual_rms_truth_matched_nominal.png",
                    "Nominal truth-matched track residual RMS (schema-2 positional)",
                )
'''
if old_plot not in text:
    raise SystemExit("residual plot block not found")
text = text.replace(old_plot, new_plot, 1)
path.write_text(text)
