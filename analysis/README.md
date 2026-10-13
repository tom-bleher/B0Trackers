# B0 tracking validation tools

These scripts turn `B0Trackers/hits` into a small set of named, reproducible
metrics instead of requiring performance changes to be judged only from plots.

## Interpretable report

The report requires Python `numpy`, `awkward`, and `uproot`; `matplotlib` is
needed only for plots.

```bash
python3 analysis/b0_report.py b0trackers.root --output-dir report
```

The output contains `summary.json` plus plots for:

- event-level reconstruction flow: truth-reachable selected primary -> seed ->
  unfiltered CKF -> ambiguity-filtered CKF -> truth-matched track;
- relative momentum residual for the truth-matched track;
- q/p, theta, and phi pulls;
- seed survival through CKF and ambiguity resolution;
- ACTS-to-sensor mapping exact/fallback fractions among resolved states, plus
  raw failure counts;
- measurement-state local residual RMS by B0 station;
- geometry/link/pixel-snap diagnostic failure counts.

The default denominator is a selected primary crossing at least three truth B0
stations. This is intentionally called **truth reachability**, not measurement-
level reconstructability: a SimHit crossing does not guarantee a digitized,
usable measurement. It is best suited to controlled single-primary samples.

State residuals also need care: `*_state_resid_*` use the best available
state in the order smoothed -> filtered -> predicted, recorded per row in
`*_state_estimate_kind` (3 smoothed, 2 filtered, 1 predicted). Smoothed
residuals include the hit itself and are shrunk, so they are fit-quality
diagnostics, not an intrinsic AC-LGAD spatial resolution. The unbiased
quantities are the predicted innovations `*_state_innov*` and their pulls.

Definitions the report applies, all recorded in `summary.json`:

- **Truth match** (`selection.truth_match`): the track's dominant MC
  association is the selected primary with weight >= `--min-truth-weight`
  (default 0.5). The efficiency numerator, the funnel's last stage and the
  fake/duplicate counts share this one threshold, so no track is both
  efficient and fake.
- **Intervals** (`selection.efficiency_interval`): Wilson score, z = 1 (68.3%).
- **Phi** (`truth_matched.phi_frame`): azimuth around the status-4 beam
  proton. In the lab frame the 25 mrad crossing angle piles every B0 proton
  up near pi.
- **Out of range** (`truth_matched.out_of_range`): eligible events outside the
  bin edges or without a finite truth value, which the binned efficiencies
  would otherwise drop silently. The default momentum edges reach 300 GeV.
- **Residuals**: `*_by_face` splits each station into its back and front
  faces (local axes flip between them), and
  `estimate_kind_of_measurement_states` says how many residuals come from
  smoothed states. Use `predicted_innovation` for unbiased pulls.
- **Seed survival**: `all seeds` is dominated by combinatorial triplets;
  `selected_primary_seeds` restricts to seeds whose hits belong to the
  selected primary.
- Files from schema < 11 written on the per-face geometry (station = cellID
  layer) are refused.

Likewise, schema 2 increments `n_sensor_map_failed` for some non-physics states
that are not expected to map to a B0 sensor. The report preserves that raw count
for debugging, but regression gating uses the fallback fraction among resolved
sensor mappings instead of treating the raw failure fraction as a physics
quality metric.

## Regression comparison

Generate reports for a known-good baseline and a candidate, then compare their
JSON summaries:

```bash
python3 analysis/compare_b0_metrics.py \
  baseline/summary.json candidate/summary.json \
  --policy analysis/regression_policy.example.json
```

The comparison exits non-zero when a configured check fails, so it can be used
in CI or a local validation script. Supported checks include absolute efficiency
drops, absolute increases, fractional worsening, and absolute changes, plus
`max_drop_sigma`. That check takes the path of an efficiency object
(`value`/`low`/`high`) and fails when the drop exceeds `tolerance` combined
Wilson half-widths, so it scales with the sample size. The default and example
policies use it for efficiency, and `half_width_68` ((q84 - q16) / 2) rather
than `std` for the momentum resolution. A run from an uncommitted checkout
(`*_dirty` in the manifest) triggers a warning.

The example policy is deliberately conservative but is **not a physics-approved
acceptance specification**. Tune tolerances to the sample size, detector setup,
and expected statistical fluctuations before making it a merge gate.

## Runtime benchmark

Use the same input and reconstruction arguments to compare ordinary EICrecon
with EICrecon plus B0Trackers, while scanning thread counts:

```bash
python3 analysis/benchmark_b0_runtime.py sim.edm4hep.root \
  --threads 1,2,4,8 \
  --runs 3 \
  --events 1000 \
  --arg=-Pacts:MaterialMap=/path/to/validated/material-map.cbor
```

For every configuration the benchmark records wall/user/system time, peak RSS
when GNU `/usr/bin/time` is available, PODIO output size, B0Trackers histogram
size, and the B0Trackers wall-time overhead relative to the baseline. Results are
written to `b0-runtime/runs.csv` and `b0-runtime/summary.json`.

Run the benchmark on an otherwise quiet machine, use the same input file and
material map, and prefer several repetitions. Thread scaling is particularly
important for this plugin because event collection retrieval happens before the
fill mutex, while the current event analysis and `TTree::Fill()` share a
serialized critical section.

## Dependency-light regression tests

The numerical metric logic has no ROOT dependency. With the build-only cleanup
from PR #1, both C++ helper tests and Python metric tests can run without an
EICrecon installation:

```bash
cmake -S . -B build-helpers -DB0Trackers_BUILD_PLUGIN=OFF
cmake --build build-helpers
ctest --test-dir build-helpers --output-on-failure
```
