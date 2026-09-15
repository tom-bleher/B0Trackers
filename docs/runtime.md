# Runtime and selective analysis modes

B0Trackers defaults to the same full diagnostic behavior as before, but expensive reconstruction products can now be disabled before `JEvent::Get()` requests them.

## Controls

```text
B0Trackers:enable_truth_seeded_chain=1
B0Trackers:enable_stub_seeded_chain=1
B0Trackers:write_track_states=1
```

The defaults are all enabled, preserving the full output.

### Stub-chain production study

For the usual B0 algorithm-performance study, the truth-seeded CKF chain is often only a reference. Skip it with:

```bash
eicrecon \
  -Pplugins=B0Trackers \
  -PB0Trackers:enable_truth_seeded_chain=0 \
  input.edm4hep.root
```

This prevents B0Trackers itself from requesting the truth-seeded seed/CKF products.

### Summary-only track study

If per-state ACTS diagnostics are not needed, keep the EDM track/trajectory summaries while skipping the large ACTS state containers:

```bash
eicrecon \
  -Pplugins=B0Trackers \
  -PB0Trackers:write_track_states=0 \
  input.edm4hep.root
```

This removes state traversal, surface mapping, and pixel snapping from the plugin analysis path and avoids requesting `*ActsTrackStates` / `*ActsTracks` solely for B0Trackers.

### Stub-only summary mode

For the lightest standard performance pass:

```bash
eicrecon \
  -Pplugins=B0Trackers \
  -PB0Trackers:enable_truth_seeded_chain=0 \
  -PB0Trackers:write_track_states=0 \
  input.edm4hep.root
```

The stub-seeded EDM tracks, associations, seed survival, hit/truth data, and event-level counters remain available.

## Reproducibility

Every tree entry records:

- `config_enable_truth_seeded_chain`
- `config_enable_stub_seeded_chain`
- `config_write_track_states`

A disabled chain intentionally reports its corresponding `has_*` branches as false; the config branches distinguish “disabled by this analysis mode” from “requested but no factory was registered.”

## Benchmarking modes

Use the paired randomized benchmark and pass selective controls only to the plugin run:

```bash
python3 analysis/benchmark_b0_runtime.py sim.edm4hep.root \
  --threads 1,2,4,8 \
  --runs 5 \
  --plugin-extra-arg=-PB0Trackers:enable_truth_seeded_chain=0 \
  --plugin-extra-arg=-PB0Trackers:write_track_states=0
```

Benchmark full and selective modes separately; do not compare files with different active analysis modes as if their output payloads were identical.

## Remaining serialization work

The current processor still uses member-backed ROOT branch buffers and therefore serializes the event-analysis/fill section with `m_fillMutex`. Moving CPU-heavy analysis outside that lock requires event-local buffers and verification that the DD4hep/ACTS geometry reads used by mapping and pixel snapping are safe concurrently. That refactor should be judged with the benchmark above rather than assumed to improve throughput.
