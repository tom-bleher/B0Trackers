# B0Trackers schema 5: per-state CKF measurement compatibility

Schema 5 adds one `ckfdiag_unf_*` row per physics state on a mapped B0
station of every unfiltered CKF candidate (accepted and rejected), so the
measurement-selection step can be diagnosed directly: why CKF sees a
measurement but rejects it as an outlier.

## Upstream requirement

Same as schema 4: the UNFILTERED ACTS containers
(`B0TrackerCKFActsTrackStatesUnfiltered` / `B0TrackerCKFActsTracksUnfiltered`),
written by `CKFTracking` only with `numB0StationsMin > 0` AND the opt-in
`B0KeepCKFDiagnostics` flag, read with `write_track_states=1`. Without them
the new branches are empty. Seeds with zero candidates (find-failed /
find-empty markers) are stateless, so they contribute no rows; the number of
source links tested per surface and the silent max-branches pruning remain
unobservable without ACTS actor instrumentation.

## Branches (flat across the event, all parallel)

- `ckfdiag_unf_seedpos`: position in the stub-seed input order (joins to the
  per-seed `ckfdiag_*` vectors, which are parallel to that order).
- `ckfdiag_unf_status`: candidate row outcome (`ckf_status`: 0 = accepted,
  1-5 = rejected with reason; see schema 4).
- `ckfdiag_unf_station`: B0 station (exact ACTS surface map only).
- `ckfdiag_unf_kind`: 0 = accepted measurement, 1 = outlier, 2 = hole.
- `ckfdiag_unf_chi2`: innovation chi2, r^T S^-1 r with S = V + H C_pred H^T
  (NaN for holes, which carry no calibrated measurement).
- `ckfdiag_unf_pull0` / `ckfdiag_unf_pull1`: normalized residuals
  r_i / sqrt(S_ii).
- `ckfdiag_unf_innov0` / `ckfdiag_unf_innov1`: raw residuals m - H x_pred
  in mm.

The compatibility definition (`b0trk::innovCompat`) is shared with the
filtered-chain `ckf_trk_state_innov_*` dump, so both see the same decision
variable the CKF measurement selector used. The reco `Chi2CutOff` (15 in the
campaign) can be overlaid directly on the chi2 distributions.
