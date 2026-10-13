# B0Trackers schema 8: reference surface and predicted global position

Schema 8 adds per-`ckfdiag_unf_*`-row surface identity and the predicted
global position, to separate wrong-sensor navigation from bad propagation:

- `ckfdiag_unf_geoid`: ACTS `geometryId().value()` of the state's reference
  surface (join against the opt-in `dump_surface_map` JSON for sensor layout:
  volumeID, center, normal, bounds, mapped station).
- `ckfdiag_unf_pgx` / `ckfdiag_unf_pgy` / `ckfdiag_unf_pgz`: predicted
  position in global coordinates [mm] (`localToGlobal` of the predicted
  locals along the predicted direction).

Diagnostic reading for a hole row with target truth hit `x_truth`:

- `|x_pred_global - x_truth| small` but no measurement: the prediction is
  fine and the loss is on the measurement side (digitization, clustering,
  or source-link/surface mapping).
- `|x_pred_global - x_truth| large`: the prediction itself missed, and the
  `geoid` tells whether CKF navigated to the wrong sensor instance
  (prediction far from every same-station truth hit) or propagated badly
  onto the right one.
