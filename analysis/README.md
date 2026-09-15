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
- ACTS-to-sensor mapping exact/fallback/failure fractions;
- measurement-state local residual RMS by B0 station;
- geometry/link/pixel-snap diagnostic failure counts.

The default denominator is a selected primary crossing at least three truth B0
stations. This is intentionally called **truth reachability**, not measurement-
level reconstructability: a SimHit crossing does not guarantee a digitized,
usable measurement. It is best suited to controlled single-primary samples.

Schema-2 state residuals also need care: the plugin currently stores the best
available state in the order smoothed -> filtered -> predicted, without a branch
recording which estimate was chosen. These residuals are therefore fit-quality
diagnostics, not an intrinsic AC-LGAD spatial-resolution measurement. A future
schema should export the estimate kind and predicted innovation/pull explicitly.

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
drops, absolute increases, fractional worsening, and absolute changes.

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
