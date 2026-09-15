#!/usr/bin/env python3
from pathlib import Path


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {count}")
    return text.replace(old, new, 1)


path = Path("analysis/b0_report.py")
text = path.read_text()

text = replace_once(
    text,
    '        schema = _scalar(tree, "schema_version", np)\n        schema_versions = sorted(set(int(x) for x in schema.tolist()))\n',
    '        schema = _scalar(tree, "schema_version", np)\n'
    '        schema_versions = sorted(set(int(x) for x in schema.tolist()))\n'
    '        is_schema3 = bool(schema_versions and min(schema_versions) >= 3)\n'
    '        if is_schema3 and "has_ckf_assocs_unfiltered" in keys:\n'
    '            availability["has_ckf_assocs_unfiltered"] = bool(\n'
    '                np.all(_scalar(tree, "has_ckf_assocs_unfiltered", np).astype(bool))\n'
    '            )\n'
    '            if not availability["has_ckf_assocs_unfiltered"]:\n'
    '                raise SystemExit(\n'
    '                    "schema-3 selected-primary stage efficiencies require "\n'
    '                    "B0TrackerCKFTrackUnfilteredAssociations; the collection was unavailable"\n'
    '                )\n'
    '        elif is_schema3:\n'
    '            raise SystemExit("schema 3 is missing has_ckf_assocs_unfiltered")\n',
    "schema3 availability",
)

old_residual = '''            # Schema 2 uses positional trajectory/ACTS-track indices. Preserve a
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
new_residual = '''            stable_state_identity = {
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
'''
text = replace_once(text, old_residual, new_residual, "stable residual/innovation")

text = replace_once(
    text,
    '            "measurement_residuals_by_station": residuals,\n',
    '            "measurement_residuals_by_station": residuals,\n'
    '            "predicted_innovation": innovation,\n',
    "innovation report field",
)

text = replace_once(
    text,
    '                "Schema-2 state residuals use the best available state in the order smoothed, filtered, predicted. "\n                "They are fit diagnostics, not intrinsic detector-resolution measurements. The nominal selected-track "\n                "view additionally relies on schema-2 positional ACTS/trajectory alignment."\n',
    '                "Legacy state residuals are fit diagnostics, not intrinsic detector-resolution measurements. "\n'
    '                "Schema 3 selects truth-matched states by stable parent-track identity and separately exports "\n'
    '                "the pre-update predicted innovation and normalized innovation chi2 for CKF compatibility."\n',
    "residual interpretation",
)

path.write_text(text)
