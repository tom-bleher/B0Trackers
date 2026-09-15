# B0Trackers

JANA2 analysis processor. Writes a `B0Trackers/hits` TTree into the
histogram file (`-Phistsfile`, default `eicrecon.root`), not the PODIO
output.

Intended geometry: the realistic 8-plane AC-LGAD B0 in `~/eic/dev/epic`
(`epic_ip6_extended`). Official 4-layer `CartesianGridXZ` B0 is supported
for hit decoding and station numbering; Init throws if the sensor map
comes out empty.

Portability caveat: hit decoding and station numbering are geometry-agnostic,
but the nearest-sensor *fallback* assumes rectangular sensor bounds
(`B0TrackerSensorWidth`/`Length`). Those bounds are correct for this fork's
16x16 mm tiles and wrong for the official trapezoidal B0 layout, so on
upstream geometry only exact ACTS surface mapping should be trusted -- check
`n_sensor_map_exact` against `n_sensor_map_fallback`, or set
`B0Trackers:fail_on_incomplete_surface_map=1` to require exact coverage.

## Build

```bash
source ~/eic/eic-shell   # from this directory, or enter the container first
cmake -S . -B build
cmake --build build -j$(nproc)
cmake --install build    # -> $EICrecon_MY/plugins/B0Trackers.so
```

Helpers (no EICrecon). The checks report failures explicitly instead of via
`assert`, so they stay live in a `-DNDEBUG` build. They can now be configured
through CMake without loading the EICrecon dependency stack:

```bash
cmake -S . -B build-helpers -DB0Trackers_BUILD_PLUGIN=OFF
cmake --build build-helpers
ctest --test-dir build-helpers --output-on-failure
```

The direct standalone build remains available:

```bash
c++ -std=c++17 -I. tests/test_helpers.cc -o build/B0TrackersHelpers_test
./build/B0TrackersHelpers_test
```

## Run

```bash
eicrecon -Pplugins=B0Trackers \
  -Phistsfile=b0trackers.root \
  sim.edm4hep.root
```

Do not blindly override `acts:MaterialMap` with a generic
`calibrations/materials-map.cbor`. The B0 CKF is sensitive to the material
model, so use the material map selected for the exact detector geometry, or
pass an explicit override only when that map was generated and validated for
the geometry being run.

For production validation runs, prefer exact ACTS surface coverage:

```bash
eicrecon -Pplugins=B0Trackers \
  -PB0Trackers:fail_on_incomplete_surface_map=1 \
  -Phistsfile=b0trackers.root \
  sim.edm4hep.root
```

With multiple JANA worker threads, event analysis may complete out of order.
Per-event analysis runs in thread-local buffers; only the final swap into the
ROOT branch-backed buffers and `TTree::Fill()` are serialized. TTree row order
is therefore not a stable event identity. Use the `eventNumber` branch when
comparing outputs across thread counts or runtime implementations.

Inputs are split into required and optional. `MCParticles`, `B0TrackerHits`,
`B0TrackerRawHits`, `B0TrackerRecHits` and `B0TrackerMeasurements` are
required: a missing factory or an upstream exception aborts the job rather
than writing a zero count. Everything else (seeds, both CKF chains, the Acts
containers) is optional, and its availability is reported per event in the
`has_*` branches -- `false` there means the factory is not registered at all,
which is a different fact from a registered factory producing nothing.

`ckf_*` is empty on stock EICrecon orthogonal seeding (`zMax=1700 mm`).
Use a B0 stub-seeder build for those branches.

## Schema 2 — what to use for performance

| Quantity | Branch | Notes |
|---|---|---|
| Truth-matched track | `truth_matched_trk_*` / `ckf_truth_matched_trk_*` | Assoc → selected MC, ranked by weight |
| Reco-only “best” | `reco_best_trk_*` | nMeasurements, then holes, then χ² |
| Oracle (cheat) | `oracle_best_trk_*` and `best_trk_*` | min \|p−p_sel\|; do not quote as resolution |
| Selected primary | `sel_primary_*`, `isSelPrimary` | Highest-p `primary_pdg`/`primary_status` |
| Selector flag | `matchesPrimarySelector` (`isPrimary`) | Every matching MC, not just selected |
| Stations on selected | `n_stations_primary` | Unique stations of `isSelPrimary` truth hits |
| RecHits / RawHits | `rec_*`, `raw_*` | B0 digitization uses a 10 keV threshold and 30 ps time resolution; RecHit *x* is the cell center |
| Cell truth purity | `{rec,raw}_nContribSim`, `_nContribMc`, `_dominantFrac`, `_totalEdep`, `_mixedCell` | `mcIndex` is the largest *summed* contributor, not the largest single step |
| Cell fired | `cell_fired` | The cell produced a RawHit -- **not** that this SimHit contributed (the digitizer links subthreshold SimHits to a fired cell) |
| Seed survival | `seed_became_track`, `seed_made_unfiltered_track`, `seed_survived_ambiguity`, `seed_n_{unfiltered,filtered}_tracks` | Separates a CKF failure from an ambiguity-solver rejection; `-1` = unknowable (no unfiltered collection) |
| Truth q/p pull | `sel_primary_charge` | Signed; the pull uses `q_truth/p_truth`, so a non-proton `primary_pdg` is handled |
| SimHit cell center | `xR/yR/zR` | Same converter as RecHit position |
| First/last SimHit | `xFirstHit` / `xLastHit` (`xEntry`/`xExit`) | Min/max time, not silicon faces |
| Track *p*, θ, loc | `trk_*` | **IP perigee**, not a B0-plane fit |
| On-plane state | `trk_state_*`, `trk_x_on_plane` | `state_index` 0 = innermost |
| Sensor map | `trk_state_mapping_method` | 0 unresolved, 1 exact, 2 gated fallback |
| PDG | `trk_pdg` | From `edm4eic::Track`; the current B0 CKF configuration uses the proton hypothesis (2212) |
| q/p = 0 | `p` is NaN, `momentum_resolved=0` | Not *p* = 0 |

Stage counters: `n_simhits`, `n_rawhits`, `n_rechits`, `n_measurements`,
`n_truth_seeds`, `n_stub_seeds`, `n_ts_unfiltered`, `n_ts_filtered`,
`n_ckf_unfiltered`, `n_ckf_filtered`, plus skip diagnostics
(`n_simhits_unresolved_cellid`, `n_missing_mc_relation`,
`n_sensor_map_*`, `n_pixel_snap_failed`).

Input availability: `has_raw_assocs`, `has_{stub,truth}_seeds`, and
`has_{ts,ckf}_{track_params,trajectories,trajectories_unfiltered,tracks,assocs,acts_states,acts_tracks,tracks_unfiltered}`.

`schema_version`, `geometry_name` (`$DETECTOR_CONFIG`), and
`detector_path` are written on every entry.
