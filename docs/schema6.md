# B0Trackers schema 6: prediction-vs-truth and outlier ownership

Schema 6 adds per-state truth residuals to the schema-5 `ckfdiag_unf_*` rows
(calibrated states only; holes keep NaNs and zero counts). For each state,
truth SimTrackerHits at the same station are projected into the state's
reference surface (`globalToLocal` along the predicted direction); the
selected-primary hit nearest the measurement anchors the residuals.

## Branches (parallel to the schema-5 rows)

- `ckfdiag_unf_pred0` / `ckfdiag_unf_pred1`: predicted local position [mm].
- `ckfdiag_unf_meas0` / `ckfdiag_unf_meas1`: calibrated measurement [mm].
- `ckfdiag_unf_dpred_t0` / `ckfdiag_unf_dpred_t1`: prediction minus truth [mm].
- `ckfdiag_unf_dmeas_t0` / `ckfdiag_unf_dmeas_t1`: measurement minus truth [mm].
- `ckfdiag_unf_tdist_prim`: local-plane distance from the measurement to the
  nearest selected-primary truth hit at the same station [mm] (NaN if none).
- `ckfdiag_unf_tdist_other`: same for the nearest non-primary truth hit [mm].
- `ckfdiag_unf_nprim_st` / `ckfdiag_unf_noth_st`: truth-hit counts from
  selected-primary / other particles at the same station.

## Diagnostic reading

- `dmeas_t ~= 0`, `dpred_t` large: prediction/propagation is wrong.
- `dpred_t ~= 0`, `dmeas_t` large: measurement/source-link association is wrong.
- Both small but the pull huge: surface/local-coordinate or covariance issue.
- Ownership split for outliers: `tdist_prim < 1 mm` (primary-owned) vs
  `tdist_other < 1 mm` (other-owned) vs neither (unmatched/ambiguous). If the
  pull shelf survives on primary-owned outliers, the predicted trajectory is
  bad rather than the hit assignment.

Limitations: truth positions are SimTrackerHit cell centers (pixel-scale
approximation, fine against multi-mm discrepancies); hole states carry no
residuals; the number of source links tested per surface remains unobservable
without ACTS actor instrumentation.
