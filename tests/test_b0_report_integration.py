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
                    "trk_state_mapping_method": ak.Array([[], [], []]),
                    "trk_state_type": ak.Array([[], [], []]),
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


if __name__ == "__main__":
    unittest.main()
