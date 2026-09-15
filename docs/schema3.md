# B0Trackers schema 3

Schema 3 makes the B0 tracking diagnostics explicit enough to distinguish reconstruction failures from fit-quality problems.

## Selected-primary reconstruction stages

The event now separates geometric reachability from usable reconstructed measurements and the later tracking stages.

Truth reachability remains `n_stations_primary`: the number of distinct physical B0 stations crossed by the selected primary in simulation. It is an acceptance quantity, not yet a reconstructability definition.

Measurement-level reconstructability is described by:

- `min_measurement_stations_required`
- `n_measurements_selected_primary`
- `n_measurement_stations_selected_primary`
- `sel_primary_measurement_reconstructable`

A `B0TrackerMeasurement` is attributed to the selected primary when that particle is the largest summed truth-energy contributor across the measurement's constituent raw-hit cells. `sel_primary_measurement_reconstructable=1` requires attributed measurements in at least `B0Trackers:min_measurement_stations` distinct physical stations (default 3, matching the current B0 CKF station minimum). If RawHit truth associations are unavailable, these measurement-level quantities are `-1` rather than being reported as a reconstruction failure.

For production truth semantics, use this together with eic/EICrecon#2967: that upstream fix makes RawHit→SimHit associations contain only deposits that actually contributed to the digitized charge. Without it, older association collections can include same-cell subthreshold SimHits and therefore make a measurement's truth composition less faithful to the digitizer response.

The later selected-primary stages are:

- `sel_primary_has_seed`
- `sel_primary_has_unfiltered_track`
- `sel_primary_has_filtered_track`
- `sel_primary_has_truth_matched_track`

Stub seeds are attributed from the truth composition of their constituent hit cells. Unfiltered and filtered track stages use the corresponding EDM track-to-MC association collections. Together, the intended selected-primary funnel is therefore:

```text
generated primary
  -> truth reachable
  -> measurement reconstructable
  -> truth-attributed seed
  -> unfiltered CKF candidate
  -> ambiguity-surviving track
  -> truth-matched final track
```

Schema-2 `n_stub_seeds` / `n_ckf_*` remain useful event-level multiplicities, not primary efficiencies.

Seed truth diagnostics are available in `seed_assoc_mcIndex`, `seed_assoc_mcCollectionID`, and `seed_assoc_weight`.

## Stable identity

The legacy `*_index` fields remain for compatibility, but track summaries also expose stable EDM/PODIO identity:

- `*_object_index`
- `*_object_collectionID`
- `*_seed_index`
- `*_seed_collectionID`
- `*_identity_valid`

Track-to-trajectory matching is performed through the EDM `Track.trajectory` relation rather than collection position. Schema 3 deliberately has **no positional fallback** from `trajectory[i]` to `Track[i]` or `TrackParameters[i]`: if the explicit relation or trajectory-owned parameters are unavailable, identity/parameters remain unresolved rather than being guessed from parallel collection order.

CKFTracking's ACTS dynamic `seed` column stores the input seed **position** (`iseed`), not the seed PODIO ObjectID. B0Trackers therefore translates that position through the actual input `TrackSeed` collection before recording seed identity or resolving ACTS-state ownership. Seed ownership maps use the full `(collectionID,index)` key.

ACTS states additionally store their parent seed and, when a seed maps uniquely to one EDM track, the stable parent-track identity:

- `*_state_parent_seed_index`
- `*_state_parent_seed_collectionID`
- `*_state_parent_track_index`
- `*_state_parent_track_collectionID`
- `*_state_parent_identity_valid`

A seed may produce multiple CKF candidates. In that case parent-track identity is deliberately marked unresolved rather than guessing from collection order.

## State estimate semantics

`*_state_estimate_kind` records which estimate backs the legacy `*_state_loc*`, `*_state_theta`, etc. branches:

- `0`: none
- `1`: predicted
- `2`: filtered
- `3`: smoothed

The preference remains smoothed -> filtered -> predicted, but it is no longer implicit.

`*_state_time` is now **nanoseconds**, matching EDM4eic track-head time. This is the intentional schema-breaking unit correction that requires the schema bump from 2 to 3.

## Predicted innovation

For measurement states with a predicted ACTS state, schema 3 exports the pre-update prediction and the innovation used to judge CKF compatibility:

- predicted parameters: `*_state_pred_*`
- upper-triangular predicted 6x6 covariance: `*_state_pred_cov_ij`
- measurement dimension and projected bound indices: `*_state_meas_dim`, `*_state_proj_index{0,1}`
- measurement covariance: `*_state_meas_cov00`, `01`, `11`
- innovation components: `*_state_innov0`, `*_state_innov1`
- component-normalized innovations: `*_state_innov_pull0`, `*_state_innov_pull1`
- full innovation chi2: `*_state_innov_chi2`

For a measurement vector `m`, predicted bound state `x_pred`, measurement projection `H`, measurement covariance `V`, and predicted state covariance `C_pred`, the diagnostic is

```text
r = m - H x_pred
S = V + H C_pred H^T
chi2_innovation = r^T S^-1 r
```

This is the preferred CKF compatibility diagnostic. Smoothed residuals remain useful fit diagnostics, but should not be quoted as an unbiased intrinsic detector resolution.

## Covariances

Track-head covariance is exported as `*_cov_ij` for the 21 independent upper-triangular elements. Predicted state covariance is exported as `*_state_pred_cov_ij` in EDM-like units (mm, rad, 1/GeV, ns) for the same 21 upper-triangular elements.

## Mapping diagnostics

Each reconstruction chain has separate counters:

- `*_n_sensor_map_exact`
- `*_n_sensor_map_fallback`
- `*_n_sensor_map_failed_physics`
- `*_n_sensor_map_unmapped_nonphysics`

This prevents truth-seeded and stub-seeded tracks from being merged into one mapping number and stops passive/non-physics ACTS states from inflating the physics mapping-failure count.
