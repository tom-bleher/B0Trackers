import math
import tempfile
import unittest
from pathlib import Path

try:
    import awkward as ak
    import numpy as np
    import uproot
except ImportError:  # pragma: no cover - dependency-light local builds may omit ROOT Python IO
    ak = np = uproot = None

from analysis.b0_report import build_report


@unittest.skipIf(uproot is None, "requires numpy, awkward, and uproot")
class B0ReportIntegrationTest(unittest.TestCase):
    def test_synthetic_schema2_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root_path = tmp / "toy.root"
            with uproot.recreate(root_path) as root:
                root["B0Trackers/hits"] = {
                    "schema_version": np.array([2, 2, 2], dtype=np.int32),
                    "geometry_name": ak.Array(["toy-geom", "toy-geom", "toy-geom"]),
                    "n_stations_primary": np.array([4, 4, 2], dtype=np.int32),
                    "n_stub_seeds": np.array([1, 1, 0], dtype=np.int32),
                    "n_ckf_unfiltered": np.array([1, 0, 0], dtype=np.int32),
                    "n_ckf_filtered": np.array([1, 0, 0], dtype=np.int32),
                    "ckf_truth_matched_trk_index": np.array([0, -1, -1], dtype=np.int32),
                    "sel_primary_p": np.array([10.0, 20.0, 30.0]),
                    "sel_primary_px": np.array([0.1, 0.2, 0.3]),
                    "sel_primary_py": np.array([0.0, 0.1, 0.0]),
                    "sel_primary_thscat_mrad": np.array([5.0, 12.0, 25.0]),
                    "sel_primary_mcIndex": np.array([7, 8, 9], dtype=np.int32),
                    "sel_primary_mcCollectionID": np.array([1, 1, 1], dtype=np.uint32),
                    "ckf_truth_matched_trk_delta_p": np.array([0.2, np.nan, np.nan]),
                    "ckf_truth_matched_trk_pull_qOverP": np.array([0.5, np.nan, np.nan]),
                    "ckf_truth_matched_trk_pull_theta": np.array([-0.2, np.nan, np.nan]),
                    "ckf_truth_matched_trk_pull_phi": np.array([0.1, np.nan, np.nan]),
                    "has_stub_seeds": np.array([True, True, True]),
                    "has_ckf_tracks_unfiltered": np.array([True, True, True]),
                    "has_ckf_tracks": np.array([True, True, True]),
                    "has_ckf_assocs": np.array([True, True, True]),
                    "seed_made_unfiltered_track": ak.Array([[1], [0], []]),
                    "seed_survived_ambiguity": ak.Array([[1], [-1], []]),
                    "ckf_trk_state_resid_loc0": ak.Array([[0.01, -0.02], [0.03], []]),
                    "ckf_trk_state_resid_loc1": ak.Array([[0.02, 0.01], [-0.01], []]),
                    "ckf_trk_state_type": ak.Array([[1, 1], [1], []]),
                    "ckf_trk_aclgad_station": ak.Array([[1, 2], [1], []]),
                    "ckf_trk_state_track_index": ak.Array([[0, 0], [0], []]),
                    "ckf_trk_state_mapping_method": ak.Array([[1, 1], [2], []]),
                    "ckf_trk_assoc_mcIndex": ak.Array([[7], [-1], []]),
                    "ckf_trk_assoc_mcCollectionID": ak.Array([[1], [0], []]),
                    "ckf_trk_assoc_weight": ak.Array([[1.0], [0.0], []]),
                    "n_simhits_unresolved_cellid": np.array([0, 0, 0], dtype=np.int32),
                    "n_missing_mc_relation": np.array([0, 1, 0], dtype=np.int32),
                    "n_pixel_snap_failed": np.array([0, 0, 0], dtype=np.int32),
                }

            report = build_report(
                root_path,
                "B0Trackers/hits",
                3,
                False,
                tmp / "report",
                momentum_bins=[0.0, 15.0, 40.0],
                angle_bins_mrad=[0.0, 10.0, 30.0],
                manifest={"input_dataset": "toy"},
            )

            self.assertEqual(report["event_stage_presence"]["eligible"]["events"], 2)
            self.assertEqual(
                report["event_stage_presence"]["truth_matched_over_eligible"]["numerator"], 1
            )
            self.assertIsNone(report["selected_primary_stages"])
            self.assertEqual(report["sensor_mapping"]["stub_ckf"]["exact"], 2)
            self.assertEqual(report["sensor_mapping"]["stub_ckf"]["fallback"], 1)
            self.assertEqual(report["provenance"]["geometry_names"], ["toy-geom"])
            self.assertEqual(report["provenance"]["manifest"]["input_dataset"], "toy")
            self.assertEqual(report["truth_matched"]["efficiency_vs_truth_p_GeV"][0]["numerator"], 1)
            self.assertEqual(
                report["track_association_quality"][
                    "eligible_events_with_duplicate_primary_tracks"
                ]["numerator"],
                0,
            )
            self.assertIn("1", report["measurement_residuals_by_station"]["all_stub_ckf_tracks"])
            self.assertIn(
                "nominal_truth_matched_track_schema2_positional",
                report["measurement_residuals_by_station"],
            )

    def test_schema3_residuals_use_stable_parent_object_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root_path = tmp / "toy-schema3.root"
            with uproot.recreate(root_path) as root:
                root["B0Trackers/hits"] = {
                    "schema_version": np.array([3], dtype=np.int32),
                    "geometry_name": ak.Array(["toy-geom"]),
                    "n_stations_primary": np.array([4], dtype=np.int32),
                    "n_stub_seeds": np.array([1], dtype=np.int32),
                    "n_ckf_unfiltered": np.array([2], dtype=np.int32),
                    "n_ckf_filtered": np.array([2], dtype=np.int32),
                    "ckf_truth_matched_trk_index": np.array([0], dtype=np.int32),
                    "ckf_truth_matched_trk_object_index": np.array([42], dtype=np.int32),
                    "ckf_truth_matched_trk_object_collectionID": np.array([7], dtype=np.uint32),
                    "ckf_truth_matched_trk_identity_valid": np.array([1], dtype=np.int32),
                    "sel_primary_p": np.array([20.0]),
                    "sel_primary_px": np.array([0.2]),
                    "sel_primary_py": np.array([0.1]),
                    "sel_primary_thscat_mrad": np.array([12.0]),
                    "sel_primary_mcIndex": np.array([8], dtype=np.int32),
                    "sel_primary_mcCollectionID": np.array([1], dtype=np.uint32),
                    "sel_primary_has_seed": np.array([1], dtype=np.int32),
                    "sel_primary_has_unfiltered_track": np.array([1], dtype=np.int32),
                    "sel_primary_has_filtered_track": np.array([1], dtype=np.int32),
                    "sel_primary_has_truth_matched_track": np.array([1], dtype=np.int32),
                    "ckf_truth_matched_trk_delta_p": np.array([0.1]),
                    "ckf_truth_matched_trk_pull_qOverP": np.array([0.1]),
                    "ckf_truth_matched_trk_pull_theta": np.array([0.2]),
                    "ckf_truth_matched_trk_pull_phi": np.array([-0.1]),
                    "has_stub_seeds": np.array([True]),
                    "has_ckf_tracks_unfiltered": np.array([True]),
                    "has_ckf_tracks": np.array([True]),
                    "has_ckf_assocs": np.array([True]),
                    "has_ckf_assocs_unfiltered": np.array([True]),
                    "ckf_trk_state_resid_loc0": ak.Array([[0.01, 0.50, 0.02]]),
                    "ckf_trk_state_resid_loc1": ak.Array([[0.02, 0.40, -0.01]]),
                    "ckf_trk_state_type": ak.Array([[1, 1, 1]]),
                    "ckf_trk_aclgad_station": ak.Array([[1, 1, 2]]),
                    "ckf_trk_state_track_index": ak.Array([[0, 1, 0]]),
                    "ckf_trk_state_mapping_method": ak.Array([[1, 1, 1]]),
                    "ckf_trk_state_parent_track_index": ak.Array([[42, 99, 42]]),
                    "ckf_trk_state_parent_track_collectionID": ak.Array([[7, 7, 7]]),
                    "ckf_trk_state_parent_identity_valid": ak.Array([[1, 1, 1]]),
                    "n_simhits_unresolved_cellid": np.array([0], dtype=np.int32),
                    "n_missing_mc_relation": np.array([0], dtype=np.int32),
                    "n_pixel_snap_failed": np.array([0], dtype=np.int32),
                }

            report = build_report(
                root_path,
                "B0Trackers/hits",
                3,
                False,
                tmp / "report",
                momentum_bins=[0.0, 40.0],
                angle_bins_mrad=[0.0, 30.0],
                manifest=None,
            )

            residuals = report["measurement_residuals_by_station"]
            self.assertIn("truth_matched_track_stable_identity", residuals)
            self.assertNotIn("nominal_truth_matched_track_schema2_positional", residuals)
            # The second station-1 state belongs to track ObjectID (7,99), so it
            # must not contaminate the selected truth-matched track (7,42).
            self.assertEqual(
                residuals["truth_matched_track_stable_identity"]["1"]["loc0"]["n"], 1
            )
            self.assertEqual(
                residuals["all_stub_ckf_tracks"]["1"]["loc0"]["n"], 2
            )

    def test_synthetic_schema3_stable_identity_and_innovation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root_path = tmp / "toy-schema3.root"
            with uproot.recreate(root_path) as root:
                root["B0Trackers/hits"] = {
                    "schema_version": np.array([3, 3], dtype=np.int32),
                    "geometry_name": ak.Array(["toy-schema3", "toy-schema3"]),
                    "n_stations_primary": np.array([4, 4], dtype=np.int32),
                    "n_stub_seeds": np.array([1, 1], dtype=np.int32),
                    "n_ckf_unfiltered": np.array([1, 1], dtype=np.int32),
                    "n_ckf_filtered": np.array([1, 0], dtype=np.int32),
                    "ckf_truth_matched_trk_index": np.array([0, -1], dtype=np.int32),
                    "ckf_truth_matched_trk_object_index": np.array([5, -1], dtype=np.int32),
                    "ckf_truth_matched_trk_object_collectionID": np.array([10, 0], dtype=np.uint32),
                    "ckf_truth_matched_trk_identity_valid": np.array([1, 0], dtype=np.int32),
                    "sel_primary_p": np.array([20.0, 25.0]),
                    "sel_primary_px": np.array([0.2, 0.25]),
                    "sel_primary_py": np.array([0.0, 0.0]),
                    "sel_primary_thscat_mrad": np.array([8.0, 10.0]),
                    "sel_primary_mcIndex": np.array([7, 8], dtype=np.int32),
                    "sel_primary_mcCollectionID": np.array([1, 1], dtype=np.uint32),
                    "sel_primary_has_seed": np.array([1, 1], dtype=np.int32),
                    "sel_primary_has_unfiltered_track": np.array([1, 1], dtype=np.int32),
                    "sel_primary_has_filtered_track": np.array([1, 0], dtype=np.int32),
                    "sel_primary_has_truth_matched_track": np.array([1, 0], dtype=np.int32),
                    "ckf_truth_matched_trk_delta_p": np.array([0.1, np.nan]),
                    "ckf_truth_matched_trk_pull_qOverP": np.array([0.2, np.nan]),
                    "ckf_truth_matched_trk_pull_theta": np.array([-0.1, np.nan]),
                    "ckf_truth_matched_trk_pull_phi": np.array([0.05, np.nan]),
                    "has_stub_seeds": np.array([True, True]),
                    "has_ckf_tracks_unfiltered": np.array([True, True]),
                    "has_ckf_tracks": np.array([True, True]),
                    "has_ckf_assocs": np.array([True, True]),
                    "has_ckf_assocs_unfiltered": np.array([True, True]),
                    "seed_made_unfiltered_track": ak.Array([[1], [1]]),
                    "seed_survived_ambiguity": ak.Array([[1], [0]]),
                    "ckf_trk_state_resid_loc0": ak.Array([[0.01, -0.01], [0.03]]),
                    "ckf_trk_state_resid_loc1": ak.Array([[0.02, 0.00], [0.01]]),
                    "ckf_trk_state_type": ak.Array([[1, 1], [1]]),
                    "ckf_trk_aclgad_station": ak.Array([[1, 2], [1]]),
                    "ckf_trk_state_track_index": ak.Array([[0, 0], [0]]),
                    "ckf_trk_state_mapping_method": ak.Array([[1, 1], [1]]),
                    "ckf_trk_state_parent_track_index": ak.Array([[5, 5], [6]]),
                    "ckf_trk_state_parent_track_collectionID": ak.Array([[10, 10], [10]]),
                    "ckf_trk_state_parent_identity_valid": ak.Array([[1, 1], [1]]),
                    "ckf_trk_state_innov_chi2": ak.Array([[1.2, 0.8], [2.0]]),
                    "ckf_trk_state_innov_pull0": ak.Array([[0.5, -0.3], [1.0]]),
                    "ckf_trk_state_innov_pull1": ak.Array([[0.2, 0.1], [0.4]]),
                    "ckf_trk_assoc_mcIndex": ak.Array([[7], [8]]),
                    "ckf_trk_assoc_mcCollectionID": ak.Array([[1], [1]]),
                    "ckf_trk_assoc_weight": ak.Array([[1.0], [1.0]]),
                    "n_simhits_unresolved_cellid": np.array([0, 0], dtype=np.int32),
                    "n_missing_mc_relation": np.array([0, 0], dtype=np.int32),
                    "n_pixel_snap_failed": np.array([0, 0], dtype=np.int32),
                }

            report = build_report(
                root_path,
                "B0Trackers/hits",
                3,
                False,
                tmp / "report3",
                momentum_bins=[0.0, 22.0, 30.0],
                angle_bins_mrad=[0.0, 9.0, 15.0],
            )

            self.assertEqual(
                report["selected_primary_stages"]["truth_matched_over_eligible"]["numerator"], 1
            )
            stable = report["measurement_residuals_by_station"][
                "truth_matched_track_stable_identity"
            ]
            self.assertIn("1", stable)
            self.assertIn("2", stable)
            innovation = report["predicted_innovation"]
            self.assertEqual(innovation["truth_matched_track"]["chi2"]["n"], 2)
            self.assertIn("1", innovation["by_station_truth_matched"])
            self.assertIn("2", innovation["by_station_truth_matched"])

    def _schema11_branches(self):
        # Beam proton along the ion axis, rotated by the 25 mrad crossing angle
        # in the x-z plane. Every selected primary scatters by 10 mrad towards
        # +y around that axis: beam-frame phi = +pi/2, lab phi close to pi.
        angle = -0.025
        bx, bz = 100.0 * math.sin(angle), 100.0 * math.cos(angle)
        p = 100.0
        px = p * math.cos(0.010) * math.sin(angle)
        py = p * math.sin(0.010)
        pz = p * math.cos(0.010) * math.cos(angle)
        return {
            "schema_version": np.array([11, 11, 11], dtype=np.int32),
            "geometry_name": ak.Array(["toy-face", "toy-face", "toy-face"]),
            "detector_path": ak.Array(["/toy/epic", "/toy/epic", "/toy/epic"]),
            "station": ak.Array([[1, 1, 2, 2], [1, 2], []]),
            "plane": ak.Array([[1, 2, 3, 4], [1, 3], []]),
            "n_stations_primary": np.array([4, 4, 4], dtype=np.int32),
            "n_stub_seeds": np.array([2, 1, 0], dtype=np.int32),
            "n_ckf_unfiltered": np.array([1, 1, 0], dtype=np.int32),
            "n_ckf_filtered": np.array([1, 1, 0], dtype=np.int32),
            "ckf_truth_matched_trk_index": np.array([0, 0, -1], dtype=np.int32),
            # Event 1's "truth match" is only 40% primary: not efficient at 0.5.
            "ckf_truth_matched_trk_assoc_weight": np.array([1.0, 0.4, np.nan]),
            "sel_primary_p": np.array([p, p, 150.0]),
            "sel_primary_px": np.array([px, px, px]),
            "sel_primary_py": np.array([py, py, py]),
            "sel_primary_pz": np.array([pz, pz, pz]),
            "sel_primary_thscat_mrad": np.array([10.0, 10.0, np.nan]),
            "sel_primary_mcIndex": np.array([7, 8, 9], dtype=np.int32),
            "sel_primary_mcCollectionID": np.array([1, 1, 1], dtype=np.uint32),
            "beam_proton_px": ak.Array([[bx], [bx], [bx]]),
            "beam_proton_py": ak.Array([[0.0], [0.0], [0.0]]),
            "beam_proton_pz": ak.Array([[bz], [bz], [bz]]),
            "ckf_truth_matched_trk_delta_p": np.array([0.2, 0.3, np.nan]),
            "ckf_truth_matched_trk_pull_qOverP": np.array([0.5, 0.1, np.nan]),
            "ckf_truth_matched_trk_pull_theta": np.array([-0.2, 0.1, np.nan]),
            "ckf_truth_matched_trk_pull_phi": np.array([0.1, 0.1, np.nan]),
            "has_stub_seeds": np.array([True, True, True]),
            "has_ckf_tracks_unfiltered": np.array([True, True, True]),
            "has_ckf_tracks": np.array([True, True, True]),
            "has_ckf_assocs": np.array([True, True, True]),
            "has_ckf_assocs_unfiltered": np.array([True, True, True]),
            "seed_made_unfiltered_track": ak.Array([[1, 0], [1], []]),
            "seed_survived_ambiguity": ak.Array([[1, -1], [0], []]),
            "seed_assoc_mcIndex": ak.Array([[7, 3], [8], []]),
            "seed_assoc_mcCollectionID": ak.Array([[1, 1], [1], []]),
            "ckf_trk_state_resid_loc0": ak.Array([[0.01, -0.02], [0.03], []]),
            "ckf_trk_state_resid_loc1": ak.Array([[0.02, 0.01], [-0.01], []]),
            "ckf_trk_state_type": ak.Array([[1, 1], [1], []]),
            "ckf_trk_state_estimate_kind": ak.Array([[3, 3], [2], []]),
            "ckf_trk_aclgad_station": ak.Array([[1, 1], [2], []]),
            "ckf_trk_aclgad_side": ak.Array([[0, 1], [0], []]),
            "ckf_trk_state_track_index": ak.Array([[0, 0], [0], []]),
            "ckf_trk_state_mapping_method": ak.Array([[1, 1], [1], []]),
            "ckf_trk_assoc_mcIndex": ak.Array([[7], [8], []]),
            "ckf_trk_assoc_mcCollectionID": ak.Array([[1], [1], []]),
            "ckf_trk_assoc_weight": ak.Array([[1.0], [0.4], []]),
            "n_simhits_unresolved_cellid": np.array([0, 0, 0], dtype=np.int32),
            "n_missing_mc_relation": np.array([0, 0, 0], dtype=np.int32),
            "n_pixel_snap_failed": np.array([0, 0, 0], dtype=np.int32),
        }

    def _report(self, tmp, branches):
        root_path = tmp / "toy.root"
        with uproot.recreate(root_path) as root:
            root["B0Trackers/hits"] = branches
        return build_report(
            root_path,
            "B0Trackers/hits",
            3,
            False,
            tmp / "report",
            momentum_bins=[0.0, 120.0],
            angle_bins_mrad=[0.0, 30.0],
            phi_bins_rad=[-math.pi, 0.0, math.pi],
            min_truth_weight=0.5,
        )

    def test_schema11_truth_match_threshold_is_shared(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self._report(Path(tmp), self._schema11_branches())
            # The 0.4-weight track is a fake, so it must not also be efficient.
            efficiency = report["event_stage_presence"]["truth_matched_over_eligible"]
            self.assertEqual((efficiency["numerator"], efficiency["denominator"]), (1, 3))
            self.assertEqual(report["selection"]["truth_match_min_weight"], 0.5)
            self.assertIn("Wilson", report["selection"]["efficiency_interval"])

    def test_schema11_phi_is_beam_relative_and_overflow_is_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self._report(Path(tmp), self._schema11_branches())
            truth = report["truth_matched"]
            self.assertTrue(truth["phi_frame"].startswith("beam"))
            self.assertAlmostEqual(truth["phi"]["median"], math.pi / 2, places=3)
            self.assertEqual(truth["out_of_range"]["truth_p_GeV"]["overflow"], 1)
            self.assertEqual(truth["out_of_range"]["scattering_angle_mrad"]["nonfinite"], 1)
            self.assertEqual(report["provenance"]["detector_paths"], ["/toy/epic"])

    def test_schema11_residuals_by_face_and_primary_seeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = self._report(Path(tmp), self._schema11_branches())
            residuals = report["measurement_residuals_by_station"]
            self.assertEqual(set(residuals["all_stub_ckf_tracks_by_face"]), {"S1B", "S1F", "S2B"})
            self.assertEqual(
                residuals["estimate_kind_of_measurement_states"], {"filtered": 1, "smoothed": 2}
            )
            primary = report["seed_survival"]["selected_primary_seeds"]
            self.assertEqual(primary["seeds"], 2)
            self.assertEqual(primary["made_unfiltered"]["numerator"], 2)
            self.assertEqual(primary["survived_ambiguity_given_ckf"]["numerator"], 1)

    def test_refuses_pre_schema11_per_face_station_numbering(self):
        branches = self._schema11_branches()
        branches["schema_version"] = np.array([10, 10, 10], dtype=np.int32)
        # Before schema 11 the realistic geometry wrote station = cellID layer.
        branches["station"] = ak.Array([[1, 2, 7, 8], [1, 3], []])
        branches["plane"] = ak.Array([[1, 2, 7, 8], [1, 3], []])
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit) as raised:
                self._report(Path(tmp), branches)
            self.assertIn("schema >= 11", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
