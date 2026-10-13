# B0Trackers schema 9: single-propagation arrivals

Schema 9 adds `tsprop_*` rows: every selected-primary truth seed propagated
once from its origin perigee to every mapped B0 surface, with the same
propagator construction CKFTracking uses (`EigenStepper` + `Navigator` +
`MaterialInteractor`, same field and tracking geometry). Opt-in via
`B0Trackers:write_propagation=true` (default off); no CKF branching or
measurement selection is involved.

## Branches (flat across the event, all parallel)

- `tsprop_seedpos`: truth-seed position; `tsprop_mcindex` / `tsprop_mccol`:
  the seed's MC particle (kinematic match).
- `tsprop_station`, `tsprop_geoid`: mapped station and surface identity.
- `tsprop_reached`: 1 when the `SurfaceReached` propagation succeeded with
  finite arrival locals.
- `tsprop_loc0` / `tsprop_loc1`: arrival locals [mm].
- `tsprop_gx` / `tsprop_gy` / `tsprop_gz`: arrival global position [mm].

## Diagnostic reading

Join arrivals to the particle's SimTrackerHits (same event/station/MC):
global `|arrival - truth|` separates the hypotheses with no CKF involved —
small means the propagator delivers and the loss is CKF navigation or
measurement mapping; large from an exact start means propagator, field, or
tracking-geometry description.
