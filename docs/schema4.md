# B0Trackers schema 4: per-seed CKF failure diagnostics

Schema 4 adds one `ckfdiag_*` entry per stub seed (including seeds that
produced zero output tracks), so `P(CKF failure reason | theta_true)` can be
built directly. Join on `ckfdiag_seed_index` (PODIO ObjectID); the vectors are
parallel to the stub-seed input order.

## Upstream requirement

Candidate and marker rows live in the UNFILTERED ACTS containers
(`B0TrackerCKFActsTrackStatesUnfiltered` / `B0TrackerCKFActsTracksUnfiltered`),
written by `CKFTracking` only when it runs with `numB0StationsMin > 0` AND the
opt-in `B0KeepCKFDiagnostics` factory flag (`-P...:B0KeepCKFDiagnostics=true`;
both B0 chains qualify; central tracking is untouched, and production B0 jobs
leave the flag off). They are read only when
`write_track_states=1` (the default). Without them, every finding quantity is
`-1` (unknown) while seed- and truth-level quantities are still filled; see
`has_ckf_acts_states_unfiltered` / `has_ckf_acts_tracks_unfiltered`.

Every processed seed owns at least one row there: accepted and CKF-rejected
candidates keep their states, and seeds with zero candidates leave a stateless
marker with the findTracks outcome. `AmbiguitySolver` and `ActsToTracks` skip
non-accepted rows, so filtered collections and all EDM products are unchanged.

## Branches

- `ckfdiag_seed_index` + `ckfdiag_seed_collectionID`: full PODIO ObjectID
  of the seed (an index alone is not an identity).
- `ckfdiag_seed_n_stations`: distinct B0 stations over the seed's TrackerHits.
- `ckfdiag_seed_assoc_weight`: dominant MC energy fraction of the seed
  (same definition as `seed_assoc_weight`).
- `ckfdiag_truth_p / _theta / _phi`: associated MC particle kinematics (NaN if
  the seed has no truth association). `_theta` is the GLOBAL polar angle;
  use `_thscat_mrad` below for physics.
- `ckfdiag_truth_thscat_mrad`: beam-relative scattering angle of the
  associated MC particle, exactly the `sel_primary_thscat_mrad` definition
  (angle to the status-4 proton beam, mrad). This is the theta for
  `P(CKF failure reason | theta)`.
- `ckfdiag_n_candidates`: non-marker CKF candidates for this seed
  (`-1` without the unfiltered ACTS containers).
- `ckfdiag_n_accepted`: candidates passing all CKF cuts (pre-ambiguity).
- `ckfdiag_stage`: `-1` unknown, `0` accepted, `1` all CKF-rejected,
  `2` findTracks failed, `3` findTracks ok but zero candidates.
- `ckfdiag_find_err_class / _value`: `0/0` none; class `1` = CKF actor error
  (update, measurement selection, max steps), `2` = propagation/navigation
  error; value is the raw `std::error_code` value. `-1` when unknown.
- `ckfdiag_best_status`: rejection reason of the best candidate
  (`0` accepted, `1` no valid measurement, `2` too few hits,
  `3` too few B0 stations, `4` smoothing failed, `5` extrapolation failed).
  Best = accepted first, then most measurements, fewest holes, lowest status.
- `ckfdiag_best_n_states`: states on the best candidate; proxy for N surfaces
  visited (passive surfaces are excluded by the CKF navigator configuration).
- `ckfdiag_best_n_meas / _n_holes / _n_outliers`: accepted measurements /
  holes / outliers on the best candidate, with OUTLIER PRECEDENCE: in
  ACTS 44 an outlier state carries both MeasurementFlag and OutlierFlag, so
  a state counts as an outlier first, then as an accepted measurement
  (`MeasurementFlag && !OutlierFlag`), then as a hole. Holes mark surfaces
  reached with no compatible measurement.
- `ckfdiag_best_last_station`: max station over outlier, measurement and hole
  states (where CKF got, including incompatible surfaces).
- `ckfdiag_best_first_hole_station`: min hole station, approximating the first
  station with no compatible measurement (`-1` if none).
- `ckfdiag_best_station_mask`: bitmask of measurement stations (bit `s`).
- `ckfdiag_cand_per_station[seed][station]`: number of candidates with a
  measurement at each station (station = the `station` branch numbering;
  inner size = `ckfdiag_max_station + 1`).
- `ckfdiag_max_station`: geometry-defined maximum station from the sensor
  map (sizes the inner vector), not stations observed in this event.
- Status-code sync note: the `ckfdiag` code mirrors in `B0TrackersHelpers.h`
  duplicate EICrecon's `b0counters::ckfdiag` values by hand (the standalone
  helper build has no EICrecon headers, so a cross-repo `static_assert` is
  impossible there); keep them in sync manually and check the helper test.

## Stated limits (require ACTS actor instrumentation, not implemented)

N source links tested per surface, the MeasurementSelector chi2 shortlist per
surface (N compatible), and silent max-branches pruning inside the CKF actor
are not observable post-hoc. Use `ckfdiag_best_n_outliers` and
`ckfdiag_best_status` as proxies, plus the run-level `B0TrackingCounters`
`ckfFindFailed / ckfFindEmpty / ckfFindErrCkf / ckfFindErrPropagation` counts.
