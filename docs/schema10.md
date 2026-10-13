# B0Trackers schema 10: SimHit-anchored local propagation

Schema 10 adds `tsloc_*` rows: the decisive local test between navigation
failure and mean-vs-stochastic trajectory divergence. Every row starts from
an actual selected-primary SimTrackerHit (truth position and truth momentum,
so all stochastic history up to the hit is preserved) and targets a mapped
ACTS B0 surface with the same propagator construction CKFTracking uses.
Opt-in via `B0Trackers:write_local_propagation=true` (default off);
`B0Trackers:local_epsilon_mm` sets the upstream offset (default 2 mm).

## Tests

- Test 0 (same-surface): start at
  `x_hit - epsilon * p_hat`, target the hit's own ACTS surface (matched by
  plane/module/sensor through the sensor map).
- Test 1 (next-station): start at `x_hit`, target every mapped surface of
  station `n+1`.
- Test 2 (no-crossing control): start exactly at `x_hit` targeting its own
  surface. Test 0 misses that test 2 reaches isolate the en-route
  volume-boundary/portal crossing as the killer.

## Branches (flat across the event, all parallel)

- `tsloc_test`, `tsloc_hitidx` (index into the SimHit `*_T` vectors),
  `tsloc_hitstation`, `tsloc_tgtstation`, `tsloc_geoid`, `tsloc_epsilon`.
- `tsloc_ixvalid` / `tsloc_ixpath`: bounded straight-ray
  `surface.intersect(..., BoundaryTolerance::None())` reachable (or
  on-surface) and its path length.
- `tsloc_reached`, `tsloc_loc0` / `tsloc_loc1`, `tsloc_gx` / `tsloc_gy` /
  `tsloc_gz`: `SurfaceReached` propagation outcome and arrival.
- `tsloc_dhit`: `|arrival - hit|` [mm]; `tsloc_edge`: truth hit distance to
  its sensor edge [mm] via `globalToLocal` + rectangle bounds (NaN when the
  surface is not rectangular or the projection fails).
- `tsloc_cell`: 0 reached, 1 ray-hits-but-propagation-misses
  (navigator/propagator-target bug), 2 ray-misses-and-propagation-misses
  (surface bounds/transform/mapping), 3 ray-misses-but-reaches
  (intersect/propagate inconsistency flag).
- `tsloc_err`: 0 on success, else `100 * category + value` with category
  1 propagator, 2 navigator, 3 stepper, 5 surface, 6 magnetic field,
  7 portal, 4 other. The dominant miss mode is 501 =
  `SurfaceError::GlobalPositionNotOnSurface` from
  `resolveLowestTrackingVolume`: the trajectory grazes a volume/portal edge
  within lookup tolerance while entering the sensor volume, so the crossing
  aborts instead of resolving into the neighboring volume.

## Diagnostic reading

Join with the schema-9 `tsprop_*` rows on event + surface for the full
four-way table: local-ok but IP→surface missed means long-range trajectory
divergence (mean vs stochastic path, possibly landing on a neighboring
finite sensor), not a basic surface-navigation bug. If test 0 comes out at
99-100%, the follow-up is which neighboring surface the mean trajectory
intersects and how far each truth hit sits from its sensor edge
(`tsloc_edge`), i.e. whether CKF navigation must be widened
B0-specifically to neighboring sensor instances.
