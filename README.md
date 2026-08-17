# B0Trackers

JANA2 analysis processor. Writes a `B0Trackers/hits` TTree into the
histogram file (`-Phistsfile`, default `eicrecon.root`), not the PODIO
output.

Intended geometry: the realistic 8-plane AC-LGAD B0 in `~/eic/dev/epic`
(`epic_ip6_extended`). Official 4-layer `CartesianGridXZ` B0 is supported
for hit decoding and station numbering; Init throws if the sensor map
comes out empty.

## Build

```bash
source ~/eic/eic-shell   # from this directory, or enter the container first
cmake -S . -B build
cmake --build build -j$(nproc)
cmake --install build    # -> $EICrecon_MY/plugins/B0Trackers.so
```

Helpers (no EICrecon):

```bash
c++ -std=c++17 -I. tests/test_helpers.cc -o build/B0TrackersHelpers_test
./build/B0TrackersHelpers_test
```

## Run

```bash
eicrecon -Pplugins=B0Trackers \
  -Phistsfile=b0trackers.root \
  -Pacts:MaterialMap=$DETECTOR_PATH/calibrations/materials-map.cbor \
  sim.edm4hep.root
```

The processor Gets collections optionally: a missing stub-seeded or Acts
chain does not drop SimHits / RecHits / the other CKF.

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
| RecHits / RawHits | `rec_*`, `raw_*` | Digitized (10 keV, 8 ns); RecHit *x* is the cell center |
| SimHit cell center | `xR/yR/zR` | Same converter as RecHit position |
| First/last SimHit | `xFirstHit` / `xLastHit` (`xEntry`/`xExit`) | Min/max time, not silicon faces |
| Track *p*, θ, loc | `trk_*` | **IP perigee**, not a B0-plane fit |
| On-plane state | `trk_state_*`, `trk_x_on_plane` | `state_index` 0 = innermost |
| Sensor map | `trk_state_mapping_method` | 0 unresolved, 1 exact, 2 gated fallback |
| PDG | `trk_pdg` | From `edm4eic::Track` (CKF pion hypothesis, not 2212) |
| q/p = 0 | `p` is NaN, `momentum_resolved=0` | Not *p* = 0 |

Stage counters: `n_simhits`, `n_rawhits`, `n_rechits`, `n_measurements`,
`n_truth_seeds`, `n_stub_seeds`, `n_ts_unfiltered`, `n_ts_filtered`,
`n_ckf_unfiltered`, `n_ckf_filtered`, plus skip diagnostics
(`n_simhits_unresolved_cellid`, `n_missing_mc_relation`,
`n_sensor_map_*`, `n_pixel_snap_failed`).

`schema_version`, `geometry_name` (`$DETECTOR_CONFIG`), and
`detector_path` are written on every entry.
