# B0Trackers schema 11: realistic-geometry stations, seed perigees, ckfdiag availability

Schema 11 targets `eic/EICrecon@b0_trackers_triplet_seeding` and
`tom-bleher/epic@b0-realistic-geometry` in `eic_xl:nightly`. Files written
by schema 10 or earlier on the realistic geometry carry wrong station and
side values (below) and should not be mixed with schema 11 files.

## Station and side (behavior change)

The realistic geometry names its faces `B0Tracker_station{N}_{back|front}_P`
(cellID layer `1..8`). Up to schema 10 only `B0Tracker_layer{N}` names were
recognized, so every face became its own station (`station` = layer, 1-8)
and every module had `side = 0`. Schema 11 takes station and side from the
face name and never from layer parity:

- `station`, `stationP`, `stationEE`, `raw_station`, `rec_station`,
  `*_aclgad_station`, `ckfdiag_*station*`, `tsprop_station`,
  `tsloc_*station` are physical stations (1-4).
- `side` (and its variants) is 1 = front, 0 = back.
- `n_stations_primary`, `n_measurement_stations_selected_primary`,
  `sel_primary_measurement_reconstructable` and `ckfdiag_max_station` count
  physical stations.
- The schema-10 next-station test (`tsloc_test == 1`) now targets station
  `n+1`, not the other face of station `n`.

The official 4-layer geometry (`B0Tracker_layer{N}_P`) is unchanged. Init
logs the `layer -> station/side` table and throws if two layers map to the
same station and side.

## Seed perigees

- `seed_perigee_x/y/z`, `truth_seed_perigee_x/y/z` [mm]: the point the
  seed's perigee line passes through. Truth seeds use the origin.
  `B0TripletSeeding` puts it `anchorDistance` (10 mm) upstream of the
  seed's first hit, so `seed_loc0/loc1/theta/phi/qOverP` are B0-local
  quantities, not IP quantities. Fitted tracks (`trk_*`, `ckf_trk_*`) are
  still extrapolated to the origin perigee by CKFTracking.

## ckfdiag availability

- `ckfdiag_status_available` (per event): 1 when the unfiltered ACTS
  containers carry `ckf_status`, 0 when they do not, -1 when the
  containers are absent.

The `ckf_status`/`ckf_find_err` columns, the zero-state find-failure markers
and the rejected candidates come only from a CKFTracking with B0
diagnostics. Upstream CKFTracking stores just the surviving tracks and the
`seed` column. When `ckfdiag_status_available == 0`:

- `ckfdiag_n_candidates` and `ckfdiag_best_status` are -1 (unknown);
  `ckfdiag_n_accepted` counts the surviving tracks.
- `ckfdiag_stage` is 0 for a seed with a surviving track and -1 otherwise:
  a failed, empty or fully rejected search cannot be told apart.
- `ckfdiag_unf_status` is -1 on every row.

## Other changes

- `B0Trackers:write_nonphysics_states` (default true;
  `config_write_nonphysics_states`): with false, `*_state_*` rows are
  written only for measurement, outlier and hole states. Material-only
  states on passive surfaces are about 80% of the rows on the realistic
  geometry; they are still counted in `*_n_sensor_map_unmapped_nonphysics`.
- The schema-10 local propagation now starts at the SimHit time in ACTS
  units (it was passed in ns).
